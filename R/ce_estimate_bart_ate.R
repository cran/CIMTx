
ce_estimate_bart_ate <- function(y, x, w, discard = FALSE,ndpost = 1000,... ) {

  n_trt <- length(unique(w))
  for (i in 1:n_trt){
      assign(paste0("n",i), sum(w==i))
    }
    xwdata = cbind(w,x)
    # Fit BART
    bart_mod = BART::pbart(x.train = xwdata, y.train = y, ndpost = ndpost,...)

    # Predict potential outcomes for w
    for (i in 1:n_trt){
      assign(paste0("xp",i), xwdata)
      for (j in 1:(n_trt)){
        assign(paste0("xp",j),  as.data.frame(eval(parse(text =paste0("xp",i)))) %>%
                 dplyr::mutate(w = j))
        assign(paste0("xp",i,j),  as.data.frame(eval(parse(text =paste0("xp",i)))) %>%
                 dplyr::filter(w == i) %>%
                 dplyr::mutate(w = j))
        assign(paste0("bart_pred",i,j), stats::predict(bart_mod, newdata = eval(parse(text = paste0("xp",i,j)))))
        assign(paste0("bart_pred",j), stats::predict(bart_mod, newdata = eval(parse(text = paste0("xp",j)))))
        assign(paste0("pred_prop",i,j), eval(parse(text = paste0("bart_pred",i,j)))[["prob.test"]])
        assign(paste0("pred_prop",j), eval(parse(text = paste0("bart_pred",j)))[["prob.test"]])
      }
    }
    if (discard == FALSE) {
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("bart_pred",i,j)))[["prob.test"]])
        }
      }
      discard_all <- 0
    } else if(discard == TRUE){
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("post.ind.sd",i,j), apply(eval(parse(text = paste0("pred_prop",i,j))), 2, stats::sd))
        }
        assign(paste0("threshold",i), max(eval(parse(text = paste0("post.ind.sd",i,i)))))
      }
   # Discarding rule
      for (i in 1:n_trt){
        n_trt_no_i <- unique(w)[unique(w)!=i]
        assign(paste0("eligible",i), TRUE)
        assign(paste0("criteria",i), TRUE)
        for (j in 1:length((n_trt_no_i))) {
          assign(paste0("criteria",i,n_trt_no_i[j]), eval(parse(text = paste0("post.ind.sd",i,n_trt_no_i[j]))) <= eval(parse(text = paste0("threshold",i))))
          assign(paste0("eligible",i), eval(parse(text = paste0("criteria",i,n_trt_no_i[j]))) & eval(parse(text = paste0("criteria",i))))
        }

        assign(paste0("n_",i, "_discard"), sum(eval(parse(text = paste0("eligible",i)))== FALSE))
      }
      discard_all <- NULL
      for (i in 1:n_trt){
        discard_all <- c(discard_all, eval(parse(text = paste0("n_",i, "_discard"))))
      }
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("pred_prop",i,j))) %>%
                   as.data.frame() %>%
                   dplyr::select(which(eval(parse(text = paste0("eligible",i))))) %>%
                   as.matrix())
        }
      }
    }

    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j, "_est"), NULL)
        assign(paste0("RR",i,j, "_est"), NULL)
        assign(paste0("OR",i,j, "_est"), NULL)
      }
    }

    for (i in 1:n_trt) {
      assign(paste0("y", i, "_pred"), NULL)
      assign(paste0("y", i, "_pred"), matrix(stats::rbinom(
        dim(eval(parse(text = (
          paste0("pred_prop", i)
        ))))[1] * dim(eval(parse(text = (
          paste0("pred_prop", i)
        ))))[2], 1, eval(parse(text = (
          paste0("pred_prop", i)
        )))
      ), nrow = ndpost))
    }
    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j, "_est"), NULL)
        assign(paste0("RR",i,j, "_est"), NULL)
        assign(paste0("OR",i,j, "_est"), NULL)
      }
    }

    for (i in 1:(n_trt - 1)) {
      for (j in (i + 1):n_trt) {
        assign(paste0("RD", i, j, "_est"), list(rowMeans(eval(parse(text = (
          paste0("y", i, "_pred")
        )))) - rowMeans(eval(parse(text = (
          paste0("y", j, "_pred")
        ))))))
        assign(paste0("RR", i, j, "_est"), list(rowMeans(eval(parse(text = (
          paste0("y", i, "_pred")
        )))) / rowMeans(eval(parse(text = (
          paste0("y", j, "_pred")
        ))))))
        assign(paste0("OR", i, j, "_est"), list(rowMeans((eval(parse(
          text = (paste0("y", i, "_pred"))
        ))) / (1 - rowMeans(eval(
          parse(text = (paste0(
            "y", i, "_pred"
          )))
        )))) / rowMeans((eval(parse(
          text = (paste0("y", j, "_pred"))
        ))) / (1 - rowMeans(eval(
          parse(text = (paste0(
            "y", j, "_pred"
          )))
        ))))))
      }
    }

    # result <- NULL
    # for (i in 1:(n_trt-1)){
    #   for (j in (i + 1):n_trt){
    #     assign(paste0("ate",i,j), posterior_summary(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est"))))))
    #     assign(paste0("ATE",i,j), list(round(eval(parse(text =(paste0("ate",i,j)))), digits = 2)))
    #     assign(paste0("ATE",i,j), stats::setNames(eval(parse(text =(paste0("ATE",i,j)))), paste0("ATE",i,j)))
    #     result <- c(result, (eval(parse(text =(paste0("ATE",i,j))))))
    #   }
    # }

    result <- NULL
    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        # assign(paste0("RD",i,j, "_est"), list(eval(parse(text =(paste0("y",i, "_pred")))) - eval(parse(text =(paste0("y",j, "_pred"))))))
        # assign(paste0("RR",i,j, "_est"), list(eval(parse(text =(paste0("y",i, "_pred")))) / eval(parse(text =(paste0("y",j, "_pred"))))))
        # assign(paste0("OR",i,j, "_est"), list((eval(parse(text =(paste0("y",i, "_pred")))) / (1 - eval(parse(text =(paste0("y",i, "_pred")))))) / (eval(parse(text =(paste0("y",j, "_pred")))) / (1 - eval(parse(text =(paste0("y",j, "_pred"))))))))
        assign(paste0("RD",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("RD",i,j, "_est")))), paste0("ATE_RD",i,j)))
        assign(paste0("RR",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("RR",i,j, "_est")))), paste0("ATE_RR",i,j)))
        assign(paste0("OR",i,j, "_est"), stats::setNames(eval(parse(text =(paste0("OR",i,j, "_est")))), paste0("ATE_OR",i,j)))
        result <- c(result, (eval(parse(text =(paste0("RD",i,j, "_est"))))), (eval(parse(text =(paste0("RR",i,j, "_est"))))), (eval(parse(text =(paste0("OR",i,j, "_est"))))))
      }
    }
    result <- c(result, list(n_discard = discard_all))
    class(result) <- "CIMTx_ATE_posterior"

return(result)
}
