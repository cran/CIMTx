#' Flexible Monte Carlo sensitivity analysis for unmeasured confounding
#'
#' The function \code{sa} implements the flexible sensitivity analysis approach for unmeasured confounding with multiple treatments and a binary outcome.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param formula A \code{\link[stats]{formula}} object for the analysis. The default is to use all terms specified in \code{x}.
#' @param prior_c_function 1) A vector of characters indicating the prior distributions for the confounding functions. Each character contains the random number generation code from the standard probability \code{\link[stats:Distributions]{distributions}} in the \code{\link[stats:stats-package]{stats}} package. 2) A vector of characters including the grid specifications for the confounding functions. It should be used when users want to formulate the  confounding  functions as scalar values. 3) A matrix indicating the point mass prior for the confounding functions
#' @param M1 A numeric value indicating the number of draws of the GPS from the posterior predictive distribution
#' @param M2 A numeric value indicating the number of draws from the prior distributions of the confounding functions
#' @param nCores A numeric value indicating number of cores to use for parallel computing.
#' @param estimand A character string representing the type of causal estimand. Only \code{"ATT"} or \code{"ATE"} is allowed. When the \code{estimand = "ATT"}, users also need to specify the reference treatment group by setting the \code{reference_trt} argument.
#' @param reference_trt A numeric value indicating reference treatment group for ATT effect.
#' @param ... Other parameters that can be passed to BART functions
#'
#' @return A list of causal estimands including risk difference (RD) between different treatment groups.
#'
#' @export
#' @importFrom foreach %dopar%
#' @references
#'
#' Hu, L., Ji, J. (2021). CIMTx: An R package for causal inference with multiple treatments using observational data. arXiv:2110.10276
#' @examples
#' \donttest{
#'lp_w_all <-
#'  c(".4*x1 + .1*x2  - 1.1*x4 + 1.1*x5",    # w = 1
#'    ".2 * x1 + .2 * x2  - 1.2 * x4 - 1.3 * x5")  # w = 2
#'nlp_w_all <-
#'  c("-.5*x1*x4  - .1*x2*x5", # w = 1
#'    "-.3*x1*x4 + .2*x2*x5")# w = 2
#'lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - 1.1*x4 - 1.2*x5", 3)
#'nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#'X_all <- c(
#'  "rnorm(100, 0, 0.5)",# x1
#'  "rbeta(100, 2, .4)",   # x2
#'  "runif(100, 0, 0.5)",# x3
#'  "rweibull(100,1,2)",  # x4
#'  "rbinom(100, 1, .4)"# x5
#')
#'set.seed(1111)
#'data <- data_sim(
#'  sample_size = 100,
#'  n_trt = 3,
#'  X = X_all,
#'  lp_y = lp_y_all,
#'  nlp_y  = nlp_y_all,
#'  align = FALSE,
#'  lp_w = lp_w_all,
#'  nlp_w = nlp_w_all,
#'  tau = c(0.5,-0.5,0.5),
#'  delta = c(0.5,0.5),
#'  psi = 2
#')
#'c_grid <- c(
#'  "runif(-0.6, 0)",# c(1,2)
#'  "runif(0, 0.6)",# c(2,1)
#'  "runif(-0.6, 0)", # c(2,3)
#'  "seq(-0.6, 0, by = 0.3)", # c(1,3)
#'  "seq(0, 0.6, by = 0.3)", # c(3,1)
#'  "runif(0, 0.6)" # c(3,2)
#')
#'sensitivity_analysis_parallel_result <-
#'  sa(
#'    M1 = 1,
#'    x = data$covariates,
#'    y = data$y,
#'    w = data$w,
#'    prior_c_function = c_grid,
#'    nCores = 1,
#'    estimand = "ATE",
#'  )
#'  }
sa <- function(x, y, w, formula = NULL, prior_c_function, M1, M2 = NULL, nCores = 1, estimand, reference_trt,... ){
  if (!(estimand %in% c("ATE", "ATT"))) stop("Estimand only supported for \"ATT\" or \"ATE\"", call. = FALSE)
  if (estimand == "ATT" && !(reference_trt %in% unique(w))) stop(paste0("Please set the reference_trt from ", paste0(sort(unique(w)), collapse = ", "), "."), call. = FALSE)
  if (sum(c(length(w) == length(y), length(w) == nrow(x), length(y) == nrow(x))) != 3) stop(paste0("The length of y, the length of w and the nrow for x should be equal. Please double check the input."), call. = FALSE)
  if (!is.null(formula)){
    x <- as.data.frame(stats::model.matrix(object = formula, cbind(y,x)))
    x <- x[,!(names(x) == "(Intercept)")]
  }
  if (any(stringr::str_detect(prior_c_function, "seq")) == FALSE && is.numeric(prior_c_function) == FALSE){
    prior_c_function_all <- matrix(NA, ncol = length(prior_c_function), nrow = M2)
    for (i in 1:length(prior_c_function)){
      str_locate_parenthesis <- stringr::str_locate(prior_c_function[i],"\\(")
      # set.seed(seed)
      prior_c_function_all[,i] <- eval(parse(text = paste0(paste0(stringr::str_sub(prior_c_function[i], 1, str_locate_parenthesis[1]), M2, ",",stringr::str_sub(prior_c_function[i], str_locate_parenthesis[1]+1)))))
    }
    prior_c_function_used <- prior_c_function_all
  }
  if (any(stringr::str_detect(prior_c_function, "seq")) == TRUE){
    c_index_with_grid <- which(stringr::str_detect(prior_c_function, "seq"))
    c_index_without_grid <- which(!stringr::str_detect(prior_c_function, "seq"))
    n_c_with_grid <- length(prior_c_function[stringr::str_detect(prior_c_function, "seq")])
    grid_length <- length(eval(parse(text = prior_c_function[stringr::str_detect(prior_c_function, "seq")])))
    M2 <- grid_length^n_c_with_grid
    c_with_grid <- prior_c_function[c_index_with_grid]
    c_without_grid <- prior_c_function[c_index_without_grid]
    c_without_grid_all <- matrix(NA, ncol = length(c_without_grid), nrow = M2)
    for (i in 1:length(c_without_grid)){
      str_locate_parenthesis <- stringr::str_locate(c_without_grid[i],"\\(")
      # set.seed(seed)
      c_without_grid_all[,i] <- eval(parse(text = paste0(paste0(stringr::str_sub(c_without_grid[i], 1, str_locate_parenthesis[1]), M2, ",",stringr::str_sub(c_without_grid[i], str_locate_parenthesis[1]+1)))))
    }
    colnames(c_without_grid_all) <- c_index_without_grid
    c_with_grid_1 <- NULL
    for (i in 1:length(c_index_with_grid)){
      assign(paste0("c_with_grid_",i), eval(parse(text = c_with_grid[i])))
    }
    c_with_grid_all <- c_with_grid_1
    for (i in 1:(length(c_index_with_grid)-1)){
      c_with_grid_all <- tidyr::expand_grid(c_with_grid_all, eval(parse(text = paste0("c_with_grid_",(i+1)))))
    }
    colnames(c_with_grid_all) <- c_index_with_grid
    c_functions_grid_final <- cbind(as.data.frame(c_without_grid_all), c_with_grid_all)
    names(c_functions_grid_final) <- paste0("c", names(c_functions_grid_final) )
    c_functions_grid_final <- c_functions_grid_final %>%
      select(paste0("c",1:length(prior_c_function)))
    prior_c_function_used <- c_functions_grid_final
  }

  if (is.numeric(prior_c_function) == TRUE){
    prior_c_function_used <- t(apply(prior_c_function, 2, mean))
  }
  # change the type of y and w as the input parameter of bart function
  x = as.matrix(x)
  y = as.numeric(y)
  w = as.integer(w)
  n_trt <- length(unique(w))
  prior_c_function_used = as.matrix(prior_c_function_used)
  n_alpha = nrow(prior_c_function_used)
  # set.seed(seed)
  # fit the treatment assigment model, to use gap-sampling, we over sample n * 10 samples, and select a sample per 10 turns
  A_model <- BART::mbart2(x.train = x, as.integer(as.factor(w)), x.test = x, ndpost = M1 * 10, mc.cores = nCores)
  # assign the estimated assignment probability to each sample, the size is (n, #treatment, sample_size)
  gps = array(A_model$prob.test[seq(1, nrow(A_model$prob.test), 10),], dim = c(M1, length(unique(w)), length(w)))
  train_x = cbind(x, w)
  n_trt <- length(unique(w))
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  iterations <- n_alpha

  if (estimand == "ATE"){
    out <-
      foreach::foreach(
        i = 1:n_alpha,
        .combine = function(x, y) {
          result_list_final <- vector("list", length = (n_trt * (n_trt - 1)/2))
          counter <- 1
          for (k in 1:(n_trt-1)){
            for (m in (k + 1):n_trt){
              result_list_final[[counter]] <- rbind(x[[paste0("ATE_",k,m)]], y[[paste0("ATE_",k,m)]])
              names(result_list_final)[[counter]] <- paste0("ATE_",k,m)
              counter <- counter +1
            }
          }
          result_list_final
        }
      ) %dopar% {

        cat('Starting ', i, 'th job.\n', sep = '')
        for (k in 1:(n_trt-1)){
          for (m in (k + 1):n_trt){
            assign(paste0("ATE_",k,m), NULL)
          }
        }

        for (j in 1:M1) {
          # correct the binary outcome based on w, prior_c_function_used, gps
          train_y = ifelse(train_x[, "w"] == sort(unique(train_x[, "w"]))[1], y - (unlist(prior_c_function_used[i, 1]) * gps[j, 2, ] + unlist(prior_c_function_used[i, 4]) * gps[j, 3, ]),
                           ifelse(train_x[, "w"] == sort(unique(train_x[, "w"]))[2], y - (unlist(prior_c_function_used[i, 2]) * gps[j, 1, ] + unlist(prior_c_function_used[i, 3]) * gps[j, 3, ]),
                                  y - (unlist(prior_c_function_used[i, 5]) * gps[j, 1, ] + unlist(prior_c_function_used[i, 6]) * gps[j, 2, ])))

          # set.seed(seed)
          # fit the bart model to estimate causal effect
          bart_mod = BART::wbart(x.train = cbind(x, w),  y.train = train_y, printevery = 10000)
          n_trt <- length(unique(w))
          for (k in 1:n_trt){
            assign(paste0("predict_",k), BART::pwbart(cbind(x, w = k), bart_mod$treedraws))
          }

          for (k in 1:(n_trt-1)){
            for (m in (k + 1):n_trt){
              assign(paste0("ATE_",k,m), c(eval(parse(text = paste0("ATE_",k,m))), rowMeans(eval(parse(text = paste0("predict_",k))) - eval(parse(text = paste0("predict_",m))))))
            }
          }
        }
        result_list <- vector("list", length = (n_trt * (n_trt - 1)/2))
        counter <- 1
        for (k in 1:(n_trt-1)){
          for (m in (k + 1):n_trt){
            result_list[[counter]] <- eval(parse(text = paste0("ATE_",k,m)))
            names(result_list)[[counter]] <- paste0("ATE_",k,m)
            counter <- counter +1
          }
        }
        return(result_list)
      }
    parallel::stopCluster(cl)
    result_list_final <- vector("list", length = (n_trt * (n_trt - 1)/2))
    counter <- 1
    for (k in 1:(n_trt-1)){
      for (m in (k + 1):n_trt){
        result_list_final[[counter]] <- out[[counter]]
        names(result_list_final)[[counter]] <- paste0("ATE_",k,m)
        counter <- counter +1
      }
    }
    # result_list_final_2 <- result_list_final
    if (any(stringr::str_detect(prior_c_function, "seq")) == TRUE){
      result_final <- vector("list", length = (n_trt * (n_trt - 1)/2))
      counter <- 1
      for (k in 1:(n_trt-1)){
        for (m in (k + 1):n_trt){
          result_final[[counter]] <- apply(result_list_final[[paste0("ATE_",k,m)]],1,mean)
          names(result_final)[[counter]] <- paste0("ATE",k,m)
          counter <- counter +1
        }
      }
      result_final <- c(result_final, list(c_functions = prior_c_function_used, grid_index = c_index_with_grid))
      class(result_final) <- "CIMTx_sa_grid"
      return(result_final)
    }
    if (any(stringr::str_detect(prior_c_function, "seq")) == FALSE){
      result_final <- NULL
      counter <- 1
      for (k in 1:(n_trt-1)){
        for (m in (k + 1):n_trt){
          result_final <- c(result_final, list(result_list_final[[paste0("ATE_",k,m)]]))
          names(result_final)[[counter]] <- paste0("ATE_RD",k,m)
          counter <- counter +1
        }
      }
      # colnames(result_final) <- c("EST","SE","LOWER","UPPER")
      class(result_final) <- "CIMTx_ATE_sa"
      return(result_final)
    }
  }

  if (estimand == "ATT"){
    w_ind <- 1:n_trt
    w_ind_no_reference <- w_ind[w_ind!=reference_trt]
    out <-
      foreach::foreach(
        i = 1:n_alpha,
        .combine = function(x, y) {
          result_list_final <- vector("list", length = (n_trt - 1))
          counter <- 1
          for (k in 1:(n_trt-1)){
            result_list_final[[counter]] <- rbind(x[[paste0("ATT_",reference_trt,w_ind_no_reference[k])]], y[[paste0("ATT_",reference_trt,w_ind_no_reference[k])]])
            names(result_list_final)[[counter]] <- paste0("ATT_",reference_trt,w_ind_no_reference[k])
            counter <- counter +1
          }
          result_list_final
        }
      ) %dopar% {

        cat('Starting ', i, 'th job.\n', sep = '')
        for (k in 1:(n_trt-1)){
          assign(paste0("ATT_",reference_trt,w_ind_no_reference[k]), NULL)
        }

        for (j in 1:M1) {
          # correct the binary outcome based on w, prior_c_function, gps
          train_y = ifelse(train_x[, "w"] == sort(unique(train_x[, "w"]))[1], y - (unlist(prior_c_function[i, 1]) * gps[j, 2, ] + unlist(prior_c_function[i, 4]) * gps[j, 3, ]),
                           ifelse(train_x[, "w"] == sort(unique(train_x[, "w"]))[2], y - (unlist(prior_c_function[i, 2]) * gps[j, 1, ] + unlist(prior_c_function[i, 3]) * gps[j, 3, ]),
                                  y - (unlist(prior_c_function[i, 5]) * gps[j, 1, ] + unlist(prior_c_function[i, 6]) * gps[j, 2, ])))

          # set.seed(seed)
          # fit the bart model to estimate causal effect
          bart_mod = BART::wbart(x.train = cbind(x, w),  y.train = train_y, printevery = 10000,...)
          n_trt <- length(unique(w))
          for (k in 1:n_trt){
            assign(paste0("predict_",k), BART::pwbart(cbind(x[w == reference_trt,], w = k), bart_mod$treedraws))
          }

          for (k in 1:(n_trt-1)){
            assign(paste0("ATT_",reference_trt,w_ind_no_reference[k]), c(eval(parse(text = paste0("ATT_",reference_trt,w_ind_no_reference[k]))), rowMeans(eval(parse(text = paste0("predict_",reference_trt))) - eval(parse(text = paste0("predict_",w_ind_no_reference[k]))))))
          }
        }
        result_list <- vector("list", length = (n_trt * (n_trt - 1)/2))
        counter <- 1
        for (k in 1:(n_trt-1)){
          result_list[[counter]] <- eval(parse(text = paste0("ATT_",reference_trt,w_ind_no_reference[k])))
          names(result_list)[[counter]] <- paste0("ATT_",reference_trt,w_ind_no_reference[k])
          counter <- counter +1
        }
        return(result_list)
      }
    # close(pb)
    parallel::stopCluster(cl)
    result_list_final <- vector("list", length =  (n_trt - 1))
    counter <- 1
    for (k in 1:(n_trt-1)){
      result_list_final[[counter]] <- out[[counter]]
      names(result_list_final)[[counter]] <- paste0("ATT_",reference_trt,w_ind_no_reference[k])
      counter <- counter +1
    }
    if (any(stringr::str_detect(prior_c_function, "seq")) == TRUE){
      result_final <- vector("list", length = (n_trt - 1))
      counter <- 1
      for (k in 1:(n_trt-1)){
        result_final[[counter]] <- apply(result_list_final[[paste0("ATT_",reference_trt,w_ind_no_reference[k])]],1,mean)
        names(result_final)[[counter]] <- paste0("ATT",reference_trt,w_ind_no_reference[k])
        counter <- counter +1
      }
      result_final <- c(result_final, list(c_functions = prior_c_function_used, grid_index = c_index_with_grid))
      class(result_final) <- "CIMTx_sa_grid"
      return(result_final)
    }
    if (any(stringr::str_detect(prior_c_function, "seq")) == FALSE){
      result_final <- NULL
      counter <- 1
      for (k in 1:(n_trt-1)){
        result_final <- c(result_final, list(result_list_final[[paste0("ATT_",reference_trt,w_ind_no_reference[k])]]))

        names(result_final)[[counter]] <- paste0("ATT_RD",reference_trt,w_ind_no_reference[k])
        counter <- counter +1
      }
    }
    class(result_final) <- "CIMTx_ATT_sa"
    return(result_final)
  }
}
