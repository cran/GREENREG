#' Presencia de roya del cafe en diferentes localidades
#'
#' Base de datos simulada que presenta la presencia/ausencia de roya del cafe
#' junto con variables ambientales y de manejo agronomico relevantes.
#' Es util para ejercicios de regresion logistica.
#'
#' @format Un data frame con 300 observaciones y 4 variables:
#' \describe{
#'   \item{presencia_roya}{Presencia (1) o ausencia (0) de roya del cafe}
#'   \item{humedad_relativa_pct}{Humedad relativa promedio (%)}
#'   \item{altitud_msnm}{Altitud de la finca (m s.n.m.)}
#'   \item{manejo_sombra}{Porcentaje de cobertura de sombra en la plantacion}
#' }
#'
#' @details Los datos fueron generados de forma controlada utilizando parametros
#' similares a los observados por SIAP y estudios agronomicos del INIFAP.
#' Esta base fue creada para ser utilizada en funciones de regresion logistica del paquete \pkg{GREENREG}.
#'
#' @examples
#' data(datos_roya_cafe)
#' modelo_log <- logistico(presencia_roya ~ humedad_relativa_pct + 
#'                         altitud_msnm + manejo_sombra, 
#'                         data = datos_roya_cafe)
#'                         
#'
#' @source Elaboracion propia con fines didacticos basada en informacion de SIAP e INIFAP.
"datos_roya_cafe"
