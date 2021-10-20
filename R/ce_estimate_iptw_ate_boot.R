ce_estimate_iptw_ate_boot <- function(y, x, w, reference_trt, method, nboots,...) {
  trim_perc <- parent.frame()$trim_perc
  SL.library <- parent.frame()$SL.library

  n_trt <- length(unique(w))
  for (i in 1:n_trt){
    assign(paste0("iptw_multiTrt_ate_result_",i, "_all"), NULL)
  }
  number_matched_all <- NULL
  names_result <- NULL
  for (j in 1:nboots) {
    bootstrap_id <- sample(length(y), replace = T)
    y_boot <- y[bootstrap_id]
    w_boot <- w[bootstrap_id]
    x_boot <- x[bootstrap_id,]
    iptw_multiTrt_ate_result <- ce_estimate_iptw_ate(y = y_boot,
                                                  x = x_boot,
                                                  w = w_boot,
                                                  method = method,... )
    names_result <- names(iptw_multiTrt_ate_result)
    for (i in 1:n_trt){
      assign(paste0("iptw_multiTrt_ate_result_",i, "_all"), cbind(eval(parse(text = paste0("iptw_multiTrt_ate_result_",i, "_all"))), iptw_multiTrt_ate_result[[i]]))
    }
    print(paste0("Finish bootstrapping ", j))
  }

  for (i in 1:n_trt){
    assign(paste0("RD_mean_",i), mean(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[1,]))
    assign(paste0("RD_se_",i), stats::sd(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[1,]))
    assign(paste0("RD_lower_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[1,], probs=0.025, na.rm = T))
    assign(paste0("RD_upper_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[1,], probs=0.975, na.rm = T))
    assign(paste0("RR_mean_",i), mean(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[2,]))
    assign(paste0("RR_se_",i), stats::sd(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[2,]))
    assign(paste0("RR_lower_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[2,], probs=0.025, na.rm = T))
    assign(paste0("RR_upper_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[2,], probs=0.975, na.rm = T))
    assign(paste0("OR_mean_",i), mean(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[3,]))
    assign(paste0("OR_se_",i), stats::sd(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[3,]))
    assign(paste0("OR_lower_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[3,], probs=0.025, na.rm = T))
    assign(paste0("OR_upper_",i), stats::quantile(eval(parse(text = paste0("iptw_multiTrt_ate_result_", i, "_all")))[3,], probs=0.975, na.rm = T))
    # summarize results
    assign(paste0("RD_",i), c(eval(parse(text = paste0("RD_mean_",i))),
                              eval(parse(text = paste0("RD_se_",i))),
                              eval(parse(text = paste0("RD_lower_",i))),
                              eval(parse(text = paste0("RD_upper_",i)))
    ))
    assign(paste0("RR_",i), c(eval(parse(text = paste0("RR_mean_",i))),
                              eval(parse(text = paste0("RR_se_",i))),
                              eval(parse(text = paste0("RR_lower_",i))),
                              eval(parse(text = paste0("RR_upper_",i)))
    ))
    assign(paste0("OR_",i), c(eval(parse(text = paste0("OR_mean_",i))),
                              eval(parse(text = paste0("OR_se_",i))),
                              eval(parse(text = paste0("OR_lower_",i))),
                              eval(parse(text = paste0("OR_upper_",i)))
    ))
    assign(paste0("res_",i), rbind(eval(parse(text = paste0("RD_",i))),
                                   eval(parse(text = paste0("RR_",i))),
                                   eval(parse(text = paste0("OR_",i)))
    ))

  }
  result_list <- NULL
  for (i in 1:n_trt){
    result_list <- c(result_list, list(round(eval(parse(text = paste0("res_",i))),2)))
  }
  for (i in 1:n_trt){
    colnames(result_list[[i]]) <- c("EST","SE","LOWER","UPPER")
    rownames(result_list[[i]]) <- c("RD", "RR", "OR")
  }

  names(result_list) <- names_result[1:n_trt]

  return(result_list)
}
