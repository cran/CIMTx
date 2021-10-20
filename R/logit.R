#' Logit function
#'
#' @param x probability from 0 to 1
#'
#' @return Logit of x
#' @export
#'
#' @examples
#' logit(0.5)
logit <- function(x){
  log(x/(1-x))
}
