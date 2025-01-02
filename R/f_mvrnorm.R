#' Generate a multivate normal vector
#'
#' @param seed seed for the randomization
#' @param m mean of distribution
#' @param sd standard deviation of the distribution
#'
#' @return a multivate normal vector
#' 
#' @import MASS


f_mvrnorm<- function(seed, m, sd){
  set.seed(seed)
  return(MASS::mvrnorm(n=1, mu = m, Sigma = sd))
}
