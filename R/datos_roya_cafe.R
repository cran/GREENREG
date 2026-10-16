#' Presencia de roya del café en diferentes localidades
#'
#' Base de datos simulada que presenta la presencia/ausencia de roya del café
#' junto con variables ambientales y de manejo agronomico relevantes.
#' Es útil para ejercicios de regresión logística.
#'
#' @format Un data frame con 300 observaciones y 4 variables:
#' \describe{
#'   \item{presencia_roya}{Presencia (1) o ausencia (0) de roya del café}
#'   \item{humedad_relativa_pct}{Humedad relativa promedio (\%)}
#'   \item{altitud_msnm}{Altitud de la finca (m s.n.m.)}
#'   \item{manejo_sombra}{Porcentaje de cobertura de sombra en la plantacion}
#' }
#'
#' @details Los datos fueron generados de forma controlada utilizando parámetros
#' similares a los observados por SIAP y estudios agronomicos del INIFAP.
#' Esta base fue creada para ser utilizada en funciones de regresión logística del paquete \pkg{GREENREG}.
#'
#' @examples
#' data(datos_roya_cafe)
#' modelo_log <- greenreg_logistico(presencia_roya ~ humedad_relativa_pct +
#'                         altitud_msnm + manejo_sombra,
#'                         data = datos_roya_cafe)
#'
#'
#' @source Elaboracion propia con fines didacticos basada en información de SIAP e INIFAP.
"datos_roya_cafe"
