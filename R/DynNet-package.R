#' A package for Causal Inference in a Latent Processes Network
#'
#' This package Contains functions to fit a dynamic model for multiple latent processes. 
#' The methodology is based on the combinaison 
#' of a linear mixed model for the latent processes trajectories and a system of 
#' differences equations for assessing their temporal influences. 
#' The estimation is done in the maximum likelihood framework.
#'
#' @importFrom graphics abline axis lines par plot points title
#' @importFrom stats as.formula model.matrix na.action na.omit pchisq pnorm printCoefmat quantile terms var
#' @importFrom Rcpp evalCpp 
#' @name DynNet-package
#' 
#' @keywords "Causality", "Dynamic model"," Latent processes"," multivariate longitudinal data"
#' @useDynLib DynNet, .registration = TRUE
"_PACKAGE"
