#' Causal inference with multiple treatments using VM for ATT effects
#'
#' The function \code{ce_estimate_vm_att} implements
#' VM to estimate ATT effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param reference_trt A numeric value indicating reference treatment group
#' for ATT effect.
#' @param caliper A numeric value denoting the caliper on the logit of
#' GPS within each cluster formed by K-means clustering.
#' The caliper is in standardized units.
#' For example, \code{caliper = 0.25} means that
#' all matches greater than 0.25 standard deviations of the
#' logit of GPS are dropped. The default value is 0.25.
#' @param n_cluster A numeric value denoting the number of clusters to
#' form using K means clustering on the logit of GPS.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function. The output also contains the number
#' of matched individuals.
#' @importFrom dplyr filter group_by summarise ungroup inner_join mutate select
#' @importFrom nnet multinom
#' @importFrom Matching Matchby
#' @references
#' Venables, W. N. & Ripley, B. D. (2002)
#' \emph{Modern Applied Statistics with S}.
#' Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Jasjeet S. Sekhon (2011).
#' Multivariate and Propensity Score Matching Software with A
#' utomated Balance Optimization: The Matching Package for R.
#' \emph{Journal of Statistical Software}, \strong{42}(7), 1-52
ce_estimate_vm_att <-
  function(y, x, w, reference_trt, caliper, n_cluster) {
    n_trt <- length(unique(w))
    if (n_trt > 3) {
      stop("We do not recommend using VM for more than 3 treatments")
    }
    # estimate generalized propensity scores using
    # multinomial logistic regression
    xwdata <- cbind(w, x)
    xwdata <- xwdata[order(xwdata$w), ]
    # Estimate the GPS from the full data using multinomial logistic regressio
    ps_fit <-
      nnet::multinom(as.factor(w) ~ ., data = xwdata, trace = FALSE)
    probs_logit1 <- data.frame(stats::fitted(ps_fit))
    colnames_probs_logit1 <- NULL
    for (i in 1:n_trt) {
      colnames_probs_logit1_once <- paste0("p", i)
      colnames_probs_logit1 <-
        c(colnames_probs_logit1, colnames_probs_logit1_once)
    }
    colnames(probs_logit1) <- colnames_probs_logit1
    xwdata <- cbind(xwdata, probs_logit1)

    # Determine eligibility
    min_max_ps <- NULL
    for (i in 1:n_trt) {
      xwdata_summarise_once <- xwdata %>%
        dplyr::group_by(w) %>%
        dplyr::summarise(min(eval(parse(text = paste0(
          "p", i
        )))), max(eval(parse(text = paste0(
          "p", i
        ))))) %>%
        dplyr::ungroup()
      names(xwdata_summarise_once)[c(2, 3)] <-
        c(paste0("min", i), paste0("max", i))
      if (i == 1) {
        min_max_ps <- xwdata_summarise_once
      } else {
        min_max_ps <-
          min_max_ps %>% dplyr::inner_join(xwdata_summarise_once, by = "w")
      }
    }

    # Discard units that fall outside the rectangular region (common
    # support region) defined by the maximum value of the smallest GPS
    # and the minimum value of the largest GPS in each treatment group.
    eligible <- TRUE
    for (i in 1:n_trt) {
      eligible_once <-
        xwdata[[paste0("p", i)]] >= max(min_max_ps[[paste0("min", i)]]) &
        xwdata[[paste0("p", i)]] <= min(min_max_ps[[paste0("max", i)]])
      eligible <- eligible & eligible_once
    }
    xwdata <- xwdata %>%
      dplyr::mutate(eligible = eligible)

    xwydata <- cbind(y, xwdata)
    xwydata <- dplyr::filter(xwydata, eligible)

    # Recalculate the GPS for the inferential units within
    # the common support region
    p1 <- NULL; p2 <- NULL; p3 <- NULL
    ps_fit_e <-
      nnet::multinom(w ~ ., data = xwydata %>%
                       dplyr::select(-p1, -p2, -p3, -eligible, -y),
                     trace = FALSE)
    probs_logit1_e <- stats::fitted(ps_fit_e)

    colnames_probs_logit1_e <- NULL
    for (i in 1:n_trt) {
      colnames_probs_logit1_e_once <- paste0("p", i)
      colnames_probs_logit1_e <-
        c(colnames_probs_logit1_e, colnames_probs_logit1_e_once)
    }
    colnames(probs_logit1_e) <- colnames_probs_logit1_e
    xwydata <- xwydata %>%
      dplyr::select(-"p1", -"p2", -"p3")
    xwydata <- cbind(xwydata, probs_logit1_e)

    for (i in 1:n_trt) {
      assign(paste0("n", i), sum(xwydata$w == i))
    }

    # k-means clustering based on the logit GPS
    for (i in 1:n_trt) {
      temp_once <-
        stats::kmeans(stats::qlogis(xwydata[[paste0("p", i)]]), n_cluster)
      xwydata <- xwydata %>%
        dplyr::mutate(temp_once$cluster)
      names(xwydata)[length(names(xwydata))] <- paste0("Quint", i)
    }
    treat <- NULL
    colnames(xwydata)[1:2] <- c("Y_obs", "treat")
    trt_indicator <- 1:n_trt
    trt_indicator_no_reference <-
      trt_indicator[trt_indicator != reference_trt]

    for (i in seq_len(length(trt_indicator_no_reference))) {
      assign(
        paste0("temp", reference_trt, trt_indicator_no_reference[i]),
        dplyr::filter(
          xwydata,
          treat %in% c(reference_trt, trt_indicator_no_reference[i])
        )
      )
    }
    assign(
      paste0(
        "temp",
        trt_indicator_no_reference[1],
        trt_indicator_no_reference[2]
      ),
      dplyr::filter(xwydata, treat %in% trt_indicator_no_reference)
    )

    # 1:1 matching based on the GPS
    # Those receiving w = 1 are matched to those receiving w = 2
    # using stats::qlogis(temp12$gps1)
    # within K-means strata of stats::qlogis(xwydata$gps3)
    # Those receiving w = 1 are matched to those receiving w = 3
    # using stats::qlogis(temp13$gps1)
    # within K-means strata of stats::qlogis(xwydata$gps2)
    for (i in seq_len(length(trt_indicator_no_reference))) {
      trt_indicator_left <-
        trt_indicator[!(trt_indicator %in% c(reference_trt,
                                             trt_indicator_no_reference[i]))]
      assign(
        paste0("match", reference_trt, trt_indicator_no_reference[i]),
        Matching::Matchby(
          Y = eval(parse(
            text = paste0("temp", reference_trt,
                          trt_indicator_no_reference[i])
          ))[["Y_obs"]],
          Tr = eval(parse(
            text = paste0("temp", reference_trt,
                          trt_indicator_no_reference[i])
          ))[["treat"]] == reference_trt,
          X = stats::qlogis(eval(parse(
            text = paste0("temp", reference_trt,
                          trt_indicator_no_reference[i])
          ))[[paste0("p", reference_trt)]]),
          by = eval(parse(
            text = paste0("temp", reference_trt,
                          trt_indicator_no_reference[i])
          ))[[paste0("Quint", trt_indicator_left)]],
          caliper = caliper,
          replace = T,
          estimand = "ATT",
          print.level = 0
        )
      )
    }

    # Extract the units receiving w = 1 who were matched to
    # units receiving w = 2 and w = 3 as well as their matches
    rownames(xwydata) <- seq_len(nrow(xwydata))
    xwydata$id <- seq_len(nrow(xwydata))
    eligible_matching <- TRUE
    for (i in seq_len(length(trt_indicator_no_reference))) {
      eligible_matching_once <-
        xwydata$id %in% eval(parse(
          text = paste0("match", reference_trt,
                        trt_indicator_no_reference[i])
        ))[["index.treated"]]
      eligible_matching <- eligible_matching & eligible_matching_once
    }
    if (sum(eligible_matching) == 0) {
      stop("No eligible matches could be found.")
    }
    xwydata$both_1 <- eligible_matching
    temp <- xwydata[xwydata$both_1 == "TRUE", ]


    for (i in seq_len(length(trt_indicator_no_reference))) {
      if (i == 1) {
        assign(
          paste0("m", reference_trt, trt_indicator_no_reference[i]),
          cbind(eval(parse(
            text = paste0("match", reference_trt,
                          trt_indicator_no_reference[i])
          ))[["index.treated"]], eval(parse(
            text = paste0("match", reference_trt,
                          trt_indicator_no_reference[i])
          ))[["index.control"]])
        )
      } else {
        trt_indicator_left <-
          trt_indicator[!(trt_indicator %in% c(reference_trt,
                                               trt_indicator_no_reference[i]))]
        assign(
          paste0("m", reference_trt, trt_indicator_no_reference[i]),
          cbind(
            eval(parse(
              text = paste0("match", reference_trt,
                            trt_indicator_no_reference[i])
            ))[["index.treated"]],
            eval(parse(
              text = paste0("match", reference_trt,
                            trt_indicator_no_reference[i])
            ))[["index.control"]] + sum(xwydata[["treat"]] ==
                                          trt_indicator_left)
          )
        )
      }
    }
    # Identify those who matched with w = 2 and w = 3
    for (i in seq_len(length(trt_indicator_no_reference))) {
      assign(paste0("m", reference_trt, trt_indicator_no_reference[i]),
             eval(parse(
               text = paste0("m", reference_trt,
                             trt_indicator_no_reference[i])
             ))[eval(parse(text = paste0(
               "m", reference_trt, trt_indicator_no_reference[i]
             )))[, 1] %in% rownames(temp), ])
    }
    triplets <- NULL
    for (i in seq_len(length(trt_indicator_no_reference))) {
      triplets_once <-
        eval(parse(text = paste0(
          "m", reference_trt, trt_indicator_no_reference[i]
        )))[order(eval(parse(
          text = paste0("m", reference_trt, trt_indicator_no_reference[i])
        ))[, 1]), ]
      triplets <- cbind(triplets, triplets_once)
    }
    triplets <- as.matrix(triplets[, c(1, 2, 4)])
    n_trip <- nrow(triplets)
    # Estimate the ATT effects based on the matched individuals:
    for (i in 1:3) {
      assign(paste0("Y", i, "_imp"),
             xwydata$Y_obs[triplets[, i]])
      assign(paste0("y", i, "_hat"),
             mean(eval(parse(
        text = paste0("Y", i, "_imp")
      ))))
    }

    for (i in seq_len(length(trt_indicator_no_reference))) {
      assign(
        paste0("RD", reference_trt,
               trt_indicator_no_reference[i], "_est"),
        eval(parse(text = paste0(
          "y", reference_trt, "_hat"
        ))) - eval(parse(
          text = paste0("y",
                        trt_indicator_no_reference[i], "_hat")
        ))
      )
      assign(
        paste0("RR", reference_trt,
               trt_indicator_no_reference[i], "_est"),
        eval(parse(text = paste0(
          "y", reference_trt, "_hat"
        ))) / eval(parse(
          text = paste0("y",
                        trt_indicator_no_reference[i], "_hat")
        ))
      )
      assign(
        paste0("OR", reference_trt,
               trt_indicator_no_reference[i], "_est"),
        (eval(parse(
          text = paste0("y", reference_trt, "_hat")
        )) / (1 - eval(
          parse(text = paste0("y", reference_trt, "_hat"))
        ))) / (eval(parse(
          text = paste0("y", trt_indicator_no_reference[i], "_hat")
        )) / (1 - eval(
          parse(text = paste0(
            "y", trt_indicator_no_reference[i], "_hat"
          ))
        )))
      )
    }
    result_list <- NULL
    for (j in seq_len(length(trt_indicator_no_reference))) {
      result_once <-
        rbind(eval(parse(
          text = paste0(
            "round(RD",
            reference_trt,
            trt_indicator_no_reference[j],
            "_est,2)"
          )
        )), eval(parse(
          text = paste0(
            "round(RR",
            reference_trt,
            trt_indicator_no_reference[j],
            "_est,2)"
          )
        )), eval(parse(
          text = paste0(
            "round(OR",
            reference_trt,
            trt_indicator_no_reference[j],
            "_est,2)"
          )
        )))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <-
        paste0("ATT", reference_trt, trt_indicator_no_reference[j])
      result_list <- c(result_list, result_once_list)
    }
    result_list <-
      c(
        result_list,
        number_matched = n_trip,
        method = "VM",
        estimand = "ATT"
      )
    class(result_list) <- "CIMTx_nonIPTW_once"
    return(result_list)
  }
