# Alias obsoletos (E-03)
#
# A partir de la version 0.1.2, las funciones exportadas se renombraron con el
# prefijo "greenreg_" para dejar de enmascarar funciones de paquetes base o
# recomendados (en particular, rlm() enmascaraba MASS::rlm() sin que el
# usuario recibiera ningun aviso: el orden de library() decidia en silencio
# que modelo se ajustaba). Estos alias se mantienen durante una version para
# no romper codigo existente, pero emiten un aviso con .Deprecated() y
# reenvian la llamada a la funcion nueva.

#' @rdname greenreg_rlm
#' @export
rlm <- function(formula, data, nivel_confianza = 0.95) {
  .Deprecated("greenreg_rlm", package = "GREENREG",
              msg = "rlm() esta obsoleta y ademas enmascara MASS::rlm(); usa greenreg_rlm().")
  greenreg_rlm(formula = formula, data = data, nivel_confianza = nivel_confianza)
}

#' @rdname greenreg_rls
#' @export
rls <- function(formula, data, nivel_confianza = 0.95) {
  .Deprecated("greenreg_rls", package = "GREENREG",
              msg = "rls() esta obsoleta; usa greenreg_rls().")
  greenreg_rls(formula = formula, data = data, nivel_confianza = nivel_confianza)
}

#' @rdname greenreg_logistico
#' @export
logistico <- function(formula, data, umbral = 0.5) {
  .Deprecated("greenreg_logistico", package = "GREENREG",
              msg = "logistico() esta obsoleta; usa greenreg_logistico().")
  greenreg_logistico(formula = formula, data = data, umbral = umbral)
}

#' @rdname greenreg_reg_poisson
#' @export
reg_poisson <- function(formula, data) {
  .Deprecated("greenreg_reg_poisson", package = "GREENREG",
              msg = "reg_poisson() esta obsoleta; usa greenreg_reg_poisson().")
  greenreg_reg_poisson(formula = formula, data = data)
}

#' @rdname greenreg_modelo_ar
#' @export
modelo_ar <- function(x, p = 1, nivel_confianza = 0.95) {
  .Deprecated("greenreg_modelo_ar", package = "GREENREG",
              msg = "modelo_ar() esta obsoleta; usa greenreg_modelo_ar().")
  greenreg_modelo_ar(x = x, p = p, nivel_confianza = nivel_confianza)
}

#' @rdname greenreg_modelo_ma
#' @export
modelo_ma <- function(x, q = 1, include_mean = TRUE) {
  .Deprecated("greenreg_modelo_ma", package = "GREENREG",
              msg = "modelo_ma() esta obsoleta; usa greenreg_modelo_ma().")
  greenreg_modelo_ma(x = x, q = q, include_mean = include_mean)
}

#' @rdname greenreg_modelo_arma
#' @export
modelo_arma <- function(x, p = 1, q = 1, include_mean = TRUE) {
  .Deprecated("greenreg_modelo_arma", package = "GREENREG",
              msg = "modelo_arma() esta obsoleta; usa greenreg_modelo_arma().")
  greenreg_modelo_arma(x = x, p = p, q = q, include_mean = include_mean)
}

#' @rdname greenreg_modelo_arima
#' @export
modelo_arima <- function(x, p = 1, d = 1, q = 1, include_mean = TRUE, include_drift = FALSE) {
  .Deprecated("greenreg_modelo_arima", package = "GREENREG",
              msg = "modelo_arima() esta obsoleta; usa greenreg_modelo_arima().")
  greenreg_modelo_arima(x = x, p = p, d = d, q = q,
                        include_mean = include_mean, include_drift = include_drift)
}

#' @rdname greenreg_analisis_datos
#' @export
analisis_datos <- function(data, variable_objetivo, verbose = TRUE) {
  .Deprecated("greenreg_analisis_datos", package = "GREENREG",
              msg = "analisis_datos() esta obsoleta; usa greenreg_analisis_datos().")
  greenreg_analisis_datos(data = data, variable_objetivo = variable_objetivo, verbose = verbose)
}
