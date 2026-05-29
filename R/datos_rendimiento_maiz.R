#' Rendimiento de maiz bajo distintos niveles de fertilizacion
#'
#' Conjunto de datos que presenta el rendimiento del maiz en funcion
#' de diferentes variables climaticas y de manejo del cultivo.
#' Este dataset es util para ejemplificar modelos de regresion lineal simple y multiple.
#'
#' @format Un data frame con 250 observaciones y 4 variables:
#' \describe{
#'   \item{rendimiento_maiz_ton_ha}{Rendimiento del cultivo de maiz (ton/ha)}
#'   \item{precipitacion_mm}{Precipitacion registrada durante el ciclo agricola (mm)}
#'   \item{temperatura_c}{Temperatura media (C)}
#'   \item{fertilizante_kg_ha}{Cantidad de fertilizante aplicado (kg/ha)}
#' }
#'
#' @details Los datos fueron generados de forma controlada, tomando como referencia
#' distribuciones y rangos reales reportados por SIAP e INEGI. Esta base se diseno
#' especificamente para probar las funciones de regresion del paquete \pkg{GREENREG}.
#'
#' @examples
#' data(datos_rendimiento_maiz)
#' modelo_rlm <- rlm(rendimiento_maiz_ton_ha ~., data = datos_rendimiento_maiz)
#' modelo_rlm       
#' plot(modelo_rlm)
#'
#' @source Elaboracion propia con fines didacticos tomando referencia de datos del SIAP e INEGI.
"datos_rendimiento_maiz"
