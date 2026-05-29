#' Registros de tala ilegal en distintas regiones
#'
#' Conjunto de datos que contiene el numero de eventos de tala ilegal observados
#' junto con variables explicativas relacionadas con la densidad poblacional,
#' superficie forestal y esfuerzo de vigilancia. Ideal para ejemplificar la aplicacion
#' de modelos de regresion Poisson.
#'
#' @format Un data frame con 200 observaciones y 4 variables:
#' \describe{
#'   \item{numero_talas}{Numero de eventos de tala ilegal observados}
#'   \item{densidad_poblacion}{Densidad poblacional de la region (hab/km2)}
#'   \item{superficie_forestal_ha}{Superficie forestal en hectareas}
#'   \item{dias_vigilancia}{Numero de dias de vigilancia reportados}
#' }
#'
#' @details Los datos fueron generados de manera simulada con base en valores reportados
#' por CONAFOR y SEMARNAT. Esto permite realizar ejercicios reproducibles y seguros
#' sin exponer datos sensibles o confidenciales.
#'
#' @examples
#' data(datos_tala_ilegal_poisson)
#' modelo_poisson <- reg_poisson(numero_talas ~ densidad_poblacion + 
#'                               superficie_forestal_ha + dias_vigilancia, 
#'                               data = datos_tala_ilegal_poisson)
#'                               
#'
#' @source Elaboracion propia con fines didacticos tomando referencia de CONAFOR e INEGI.
"datos_tala_ilegal_poisson"

