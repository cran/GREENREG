#' Anomalias de temperatura media diaria
#'
#' Conjunto de datos que muestra las anomalias (diferencias respecto a la media historica)
#' de la temperatura media diaria registradas a lo largo del tiempo. Este dataset puede
#' utilizarse para ilustrar analisis de series de tiempo y modelos autorregresivos o ARIMA.
#'
#' @format Un data frame con 730 observaciones y 2 variables:
#' \describe{
#'   \item{fecha}{Fecha del registro (AAAA-MM-DD)}
#'   \item{anomalia_c}{Anomalia de temperatura respecto al promedio historico (C)}
#' }
#'
#' @details Los datos fueron construidos de manera simulada a partir de patrones reales
#' observados en bases climaticas del INEGI y de la CONAGUA. La creacion propia de esta
#' base permite garantizar que los valores sean adecuados para fines didacticos y
#' compatibles con las funciones del paquete \pkg{GREENREG}.
#'
#' @examples
#' data("datos_anomalia_temperatura")
#' modelo_AR <- modelo_ar(datos_anomalia_temperatura$anomalia_c, p = 2)
#' modelo_AR
#' plot(modelo_AR)
#' 
#' @source Elaboracion propia con fines didacticos a partir de informacion publica de INEGI y CONAGUA.
"datos_anomalia_temperatura"

