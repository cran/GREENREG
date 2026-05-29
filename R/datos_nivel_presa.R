#' Niveles de agua en presas de Mexico
#'
#' Conjunto de datos que presenta los niveles de almacenamiento de agua
#' registrados en distintas presas y fechas, util para analisis de series de tiempo
#' y modelos ARIMA.
#'
#' @format Un data frame con 1825 observaciones y 2 variables:
#' \describe{
#'   \item{fecha}{Fecha del registro (AAAA-MM-DD)}
#'   \item{nivel_m}{Nivel de almacenamiento de la presa (metros)}
#' }
#'
#' @details Los datos fueron simulados a partir de patrones observados en reportes
#' de CONAGUA, manteniendo coherencia temporal y realismo. Se generaron con fines
#' didacticos para demostrar la aplicacion de modelos AR, MA y ARIMA en \pkg{GREENREG}.
#'
#' @examples
#' data("datos_nivel_presa")
#' ts_data <- datos_nivel_presa$nivel_m
#' modelo_ARIMA <- modelo_arima(ts_data, p = 1, d = 1, q = 1)
#' print(modelo_ARIMA)
#' plot(modelo_ARIMA)
#'
#' @source Elaboracion propia con fines didacticos basada en patrones de datos de CONAGUA.
"datos_nivel_presa"

