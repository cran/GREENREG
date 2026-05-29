# R/globals.R

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "ACF", "Ajustado", "Conteo", "Cook", "Distancia", "Etiqueta", "Freq", 
    "Fuente", "High", "ID", "IRR", "Imaginario", "Index", "Label", "Lag", 
    "Leverage", "Low", "Lower", "Media_Movil", "Observado", "PACF", "P_val", 
    "Predicho", "Real", "ResidStd", "Residuos", "SD_Movil", "Tiempo", "Tipo", 
    "Upper", "VIF", "Valor", "Variable", "X_val", "Y", "ajustados", 
    "conf.high", "conf.low", "density", "dnorm", "esp", "estimate", 
    "index", "leverage", "obs", "obs_x", "obs_y", "pearson_res", "pred", 
    "predicho", "prob", "prob_media", "res_sq", "resid", "resid_std", 
    "residuos", "term", "vjust_pos", "y"
  ))
}

#' @importFrom stats as.formula density dnorm printCoefmat reorder resid
#' @importFrom utils globalVariables
NULL