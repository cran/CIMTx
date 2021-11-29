ce_estimate_bart_att <- function(y, x, w, discard = FALSE, ndpost = 1000,reference_trt,...) {

  n_trt <- length(unique(w))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]

  for (i in 1:n_trt){
    assign(paste0("n",i), sum(w==i))
  }

  xwdata = cbind(w,x)

  # Fit BART
  bart_mod = BART::pbart(x.train = xwdata, y.train = y, ndpost = ndpost, ...)

  # Predict potential outcomes for among those in reference_trt group

  assign(paste0("xp",reference_trt), xwdata[w==reference_trt,])
  for (i in trt_indicator[trt_indicator!=reference_trt]){
    assign(paste0("xp",i), xwdata[w==reference_trt,])
    assign(paste0("xp",i), as.data.frame(eval(parse(text = paste0("xp",i)))) %>%
             dplyr::mutate(w = i))
  }

  for (j in 1:(n_trt)){
      assign(paste0("bart_pred",reference_trt,j), BART::pwbart(eval(parse(text = paste0("xp",j))), bart_mod$treedraws, mu=mean(y)))
      assign(paste0("pred_prop",reference_trt,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference_trt,j)))))
    }

  if (discard == FALSE) {

    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference_trt,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference_trt,j)))))
    }
    n_discard_att <- 0
  } else if (discard == TRUE){
    for (j in 1:(n_trt)){
      assign(paste0("post.ind.sd",j), apply(eval(parse(text = paste0("pred_prop",reference_trt,j))), 2, stats::sd))
    }
    threshold <- max(eval(parse(text = paste0("post.ind.sd",reference_trt))))

    for (j in 1:length((trt_indicator_no_reference))) {
      assign(paste0("eligible",trt_indicator_no_reference[j]), (eval(parse(text = paste0("post.ind.sd",trt_indicator_no_reference[j]))) <= threshold))
    }
    eligible <- rep(TRUE, dim(eval(parse(text = paste0("pred_prop",reference_trt,1))))[2])
    for (j in 1:length((trt_indicator_no_reference))){
      eligible <- (eval(parse(text = paste0("eligible",trt_indicator_no_reference[j]))) & eligible)

    }
    n_discard_att <- sum(eligible == FALSE)
    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference_trt,j), (eval(parse(text = paste0("pred_prop",reference_trt,j))) %>%
               as.data.frame() %>%
               dplyr::select(which(eligible)) %>%
               as.matrix()))
    }

  }
  # Estimate causal effects
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est"), NULL)
    }

  for (i in 1:n_trt) {
    assign(paste0("y", i, "_pred"), NULL)
    assign(paste0("y", i, "_pred"), matrix(stats::rbinom(
      dim(eval(parse(text = (
        paste0("pred_prop",reference_trt, i)
      ))))[1] * dim(eval(parse(text = (
        paste0("pred_prop",reference_trt, i)
      ))))[2], 1, eval(parse(text = (
        paste0("pred_prop",reference_trt, i)
      )))
    ), nrow = ndpost))
  }



  for (j in 1:length((trt_indicator_no_reference))){
    assign(paste0("RD", reference_trt,trt_indicator_no_reference[j], "_est"), list(rowMeans(eval(parse(text = (
      paste0("y", reference_trt, "_pred")
    )))) - rowMeans(eval(parse(text = (
      paste0("y",trt_indicator_no_reference[j], "_pred")
    ))))))
    assign(paste0("RR", reference_trt,trt_indicator_no_reference[j], "_est"), list(rowMeans(eval(parse(text = (
      paste0("y", reference_trt, "_pred")
    )))) / rowMeans(eval(parse(text = (
      paste0("y",trt_indicator_no_reference[j], "_pred")
    ))))))
    assign(paste0("OR", reference_trt,trt_indicator_no_reference[j], "_est"), list(rowMeans((eval(parse(
      text = (paste0("y", reference_trt, "_pred"))
    ))) / (1 - rowMeans(eval(
      parse(text = (paste0(
        "y", reference_trt, "_pred"
      )))
    )))) / rowMeans((eval(parse(
      text = (paste0("y",trt_indicator_no_reference[j], "_pred"))
    ))) / (1 - rowMeans(eval(
      parse(text = (paste0(
        "y",trt_indicator_no_reference[j], "_pred"
      )))
    ))))))
  }

  # result <- NULL
  # for (j in 1:length((trt_indicator_no_reference))){
  #   assign(paste0("att",reference_trt,trt_indicator_no_reference[j]), posterior_summary(eval(parse(text =(paste0("RD",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference_trt,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference_trt,trt_indicator_no_reference[j], "_est"))))))
  #   assign(paste0("ATT",reference_trt,trt_indicator_no_reference[j]), list(round(eval(parse(text =(paste0("att",reference_trt,trt_indicator_no_reference[j])))), digits = 2)))
  #   assign(paste0("ATT",reference_trt,trt_indicator_no_reference[j]), stats::setNames(eval(parse(text =(paste0("ATT",reference_trt,trt_indicator_no_reference[j])))), paste0("ATT",reference_trt,trt_indicator_no_reference[j])))
  #   result <- c(result, (eval(parse(text =(paste0("ATT",reference_trt,trt_indicator_no_reference[j]))))))
  # }

  result <- NULL
  for (i in 1:(n_trt-1)){
    for (j in 1:length(trt_indicator_no_reference)){
      assign(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_RD",reference_trt, trt_indicator_no_reference[j])))
      assign(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_RR",reference_trt, trt_indicator_no_reference[j])))
      assign(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"), stats::setNames(eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est")))), paste0("ATT_OR",reference_trt, trt_indicator_no_reference[j])))
      result <- c(result, (eval(parse(text =(paste0("RD",reference_trt, trt_indicator_no_reference[j], "_est"))))), (eval(parse(text =(paste0("RR",reference_trt, trt_indicator_no_reference[j], "_est"))))), (eval(parse(text =(paste0("OR",reference_trt, trt_indicator_no_reference[j], "_est"))))))
    }
  }
  result <- c(result, list(n_discard = n_discard_att))
  class(result) <- "CIMTx_ATT_posterior"

  return(result)
}
