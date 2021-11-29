ce_estimate_vm_att_boot <- function(y, x, w, reference_trt, caliper, n_cluster, nboots,verbose_boot) {

  n_trt <- length(unique(w))
  for (i in 1:(n_trt-1)){
    assign(paste0("vm_multiTrt_att_result_",i, "_all"), NULL)
  }
  number_matched_all <- NULL
  names_result <- NULL
  for (j in 1:nboots) {
    bootstrap_id <- sample(length(y), replace = T)
    y_boot <- y[bootstrap_id]
    trt_boot <- w[bootstrap_id]
    x_boot <- x[bootstrap_id,]
    vm_multiTrt_att_result <- ce_estimate_vm_att(y = y_boot,
                                              x = x_boot,
                                              w = trt_boot,
                                              reference_trt,
                                              caliper = caliper,
                                              n_cluster = n_cluster)
    names_result <- names(vm_multiTrt_att_result)
    for (i in 1:(length(names_result)-1)){
      assign(paste0("vm_multiTrt_att_result_",i, "_all"), cbind(eval(parse(text = paste0("vm_multiTrt_att_result_",i, "_all"))), vm_multiTrt_att_result[[i]]))
    }
    number_matched_all <-
      c(number_matched_all,
        vm_multiTrt_att_result$number_matched)
    if (verbose_boot == TRUE){
      print(paste0("Finish bootstrapping ", j))
    }
  }

  result <- NULL
  for (i in 1:(n_trt-1)){
    assign(paste0("RD_",i), list(as.double(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[1,])))
    assign(paste0("RR_",i), list(as.double(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[2,])))
    assign(paste0("OR_",i), list(as.double(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[3,])))

    assign(paste0("RD_",i), stats::setNames(eval(parse(text =(paste0("RD_",i)))), paste0("ATT_RD", stringr::str_sub(names_result[i], 4,5))))
    assign(paste0("RR_",i), stats::setNames(eval(parse(text =(paste0("RR_",i)))), paste0("ATT_RR", stringr::str_sub(names_result[i], 4,5))))
    assign(paste0("OR_",i), stats::setNames(eval(parse(text =(paste0("OR_",i)))), paste0("ATT_OR", stringr::str_sub(names_result[i], 4,5))))
    result <- c(result, (eval(parse(text =(paste0("RD_",i))))), (eval(parse(text =(paste0("RR_",i))))), (eval(parse(text =(paste0("OR_",i))))))
  }
  result <- c(result, list(number_matched_all = number_matched_all))
  class(result) <- "CIMTx_ATT_posterior"

  # for (i in 1:(n_trt-1)){
  #   assign(paste0("RD_mean_",i), mean(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[1,]))
  #   assign(paste0("RD_se_",i), stats::sd(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[1,]))
  #   assign(paste0("RD_lower_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[1,], probs=0.025, na.rm = T))
  #   assign(paste0("RD_upper_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[1,], probs=0.975, na.rm = T))
  #   assign(paste0("RR_mean_",i), mean(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[2,]))
  #   assign(paste0("RR_se_",i), stats::sd(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[2,]))
  #   assign(paste0("RR_lower_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[2,], probs=0.025, na.rm = T))
  #   assign(paste0("RR_upper_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[2,], probs=0.975, na.rm = T))
  #   assign(paste0("OR_mean_",i), mean(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[3,]))
  #   assign(paste0("OR_se_",i), stats::sd(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[3,]))
  #   assign(paste0("OR_lower_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[3,], probs=0.025, na.rm = T))
  #   assign(paste0("OR_upper_",i), stats::quantile(eval(parse(text = paste0("vm_multiTrt_att_result_", i, "_all")))[3,], probs=0.975, na.rm = T))
  #   # summarize results
  #   assign(paste0("RD_",i), c(eval(parse(text = paste0("RD_mean_",i))),
  #                             eval(parse(text = paste0("RD_se_",i))),
  #                             eval(parse(text = paste0("RD_lower_",i))),
  #                             eval(parse(text = paste0("RD_upper_",i)))
  #                             ))
  #   assign(paste0("RR_",i), c(eval(parse(text = paste0("RR_mean_",i))),
  #                             eval(parse(text = paste0("RR_se_",i))),
  #                             eval(parse(text = paste0("RR_lower_",i))),
  #                             eval(parse(text = paste0("RR_upper_",i)))
  #   ))
  #   assign(paste0("OR_",i), c(eval(parse(text = paste0("OR_mean_",i))),
  #                             eval(parse(text = paste0("OR_se_",i))),
  #                             eval(parse(text = paste0("OR_lower_",i))),
  #                             eval(parse(text = paste0("OR_upper_",i)))
  #   ))
  #   assign(paste0("res_",i), rbind(eval(parse(text = paste0("RD_",i))),
  #                                  eval(parse(text = paste0("RR_",i))),
  #                                  eval(parse(text = paste0("OR_",i)))
  #                                  ))
  # }
  #
  # result_list <- NULL
  # for (i in 1:(n_trt-1)){
  #   result_list <- c(result_list, list(round(eval(parse(text = paste0("res_",i))),2)))
  # }
  # for (i in 1:(n_trt-1)){
  #   colnames(result_list[[i]]) <- c("EST","SE","LOWER","UPPER")
  #   rownames(result_list[[i]]) <- c("RD", "RR", "OR")
  # }
  # result_list <- c(result_list, list(number_matched_all))
  # names(result_list) <- names_result

  return(result)
}
