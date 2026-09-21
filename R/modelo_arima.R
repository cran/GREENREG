#' Modelo Integrado de Media Movil Autorregresivo (ARIMA)
#'
#' Esta funcion ajusta el modelo mas versatil para series de tiempo univariadas:
#' el ARIMA(p, d, q). Combina la memoria de largo plazo (AR), la estabilizacion
#' mediante diferenciacion (I) para eliminar tendencias, y la respuesta a
#' choques aleatorios de corto plazo (MA).
#'
#' @details
#' El flujo interno de la funcion ejecuta una auditoria de cinco fases:
#' \enumerate{
#'   \item \strong{Diferenciacion (Fase I):} Si \eqn{d > 0}, el modelo resta los valores
#'         consecutivos para eliminar la tendencia y lograr que la serie sea
#'         "Estacionaria" (promedio constante).
#'   \item \strong{Estimacion ARMA:} Sobre la serie ya estabilizada, estima los
#'         parametros \eqn{\phi} (AR) y \eqn{\theta} (MA) mediante Maxima Verosimilitud.
#'   \item \strong{Doble Validacion de Estructura:} Analiza las raices del polinomio
#'         caracteristico para asegurar que la parte AR sea \strong{Estable} y la
#'         parte MA sea \strong{Invertible}.
#'   \item \strong{Auditoria de Residuos:} Verifica que el error final sea "Ruido Blanco"
#'         mediante las pruebas de Ljung-Box (independencia) y Shapiro-Wilk (normalidad).
#'   \item \strong{Diagnostico Dinamico:} Evalua si la media y la varianza se mantienen
#'         estables a lo largo del tiempo despues de la integracion.
#' }
#'
#' @section Guia de las 10 Graficas de Diagnostico:
#' Al ejecutar \code{plot(modelo)}, se despliega el diagnostico mas completo del paquete:
#' \itemize{
#'   \item \strong{1. Serie Original:} Visualiza la serie bruta. Fundamental para identificar la tendencia inicial que justifica el valor de \eqn{d}.
#'   \item \strong{2. Ajuste Integrado:} Superpone el modelo (purpura) sobre la realidad, incluyendo bandas de confianza del 95\%.
#'   \item \strong{3. Serie Diferenciada:} Muestra la serie despues de aplicar \eqn{d}. Si el modelo es correcto, esta serie debe verse plana y sin tendencia.
#'   \item \strong{4. Circulo Unitario:} Grafica de convergencia matematica. Los puntos rojos (AR) y azules (MA) deben estar dentro del circulo gris.
#'   \item \strong{5. Residuos en el Tiempo:} Busca patrones de error. Los puntos deben formar una nube aleatoria sin ciclos visibles.
#'   \item \strong{6. ACF de Residuos:} Prueba de independencia. Las barras deben estar dentro de los limites para confirmar el Ruido Blanco.
#'   \item \strong{7. PACF de Residuos:} Ayuda a identificar si el orden autorregresivo \eqn{p} es suficiente para explicar la inercia.
#'   \item \strong{8. Normalidad (Q-Q Plot):} Valida que los errores sigan una distribucion gaussiana para asegurar pronosticos confiables.
#'   \item \strong{9. Precision (Real vs Ajustado):} Evaluacion de la correlacion de Pearson entre los datos y la estimacion.
#'   \item \strong{10. Estacionariedad Final:} Muestra la media y varianza movil. Una franja purpura estable certifica un modelo robusto.
#' }
#'
#' @param x Para \code{modelo_arima()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arima_greenreg"} devuelto por \code{modelo_arima()}.
#' @param p Orden Autorregresivo (memoria del pasado).
#' @param d Orden de Integracion (numero de diferenciaciones).
#' @param q Orden de Media Movil (memoria de los errores).
#' @param include_mean Logico. Incluir constante? (Suele ser \code{FALSE} si \eqn{d > 0}).
#'
#' @return Un objeto de clase \code{"arima_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: El objeto \code{arima} original de \code{stats}.
#'   \item \code{ordenes}: Vector con los valores de p, d y q.
#'   \item \code{coeficientes_tabla}: Matriz con estimaciones, errores estandar y p-valores.
#'   \item \code{estatus}: Lista con el estado de estabilidad e invertibilidad.
#'   \item \code{diagnosticos}: Resultados de Ljung-Box, Shapiro-Wilk y AIC.
#'   \item \code{data_plot}: Dataframe listo para visualización con \code{ggplot2}.
#'   \item \code{data_diff}: Dataframe de la serie diferenciada (si \eqn{d > 0}).
#' }
#'
#' @examples
#' data("datos_nivel_presa")
#' ts_data <- datos_nivel_presa$nivel_m
#' modelo_ARIMA <- modelo_arima(ts_data, p = 1, d = 1, q = 1)
#' print(modelo_ARIMA)
#' plot(modelo_ARIMA)
#' 
#'
#' @seealso \code{\link[stats]{arima}}
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf filter sd
#' @import ggplot2
#' @export
modelo_arima <- function(x, p = 1, d = 1, q = 1, include_mean = TRUE) {

  # --- 1. Validacion y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numerico.")

  # Ajuste ARIMA(p, d, q)
  modelo_base <- stats::arima(x, order = c(p, d, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadisticas ---
  coefs <- modelo_base$coef

  # Si no hay coeficientes (ej: ARIMA(0,1,0)), manejamos el caso
  if (length(coefs) > 0) {
    se <- sqrt(diag(modelo_base$var.coef))
    z_val <- coefs / se
    p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
    tabla_coefs <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)
  } else {
    tabla_coefs <- matrix(NA, nrow=0, ncol=4)
  }

  # --- 3. Analisis de Raices (Estabilidad e Invertibilidad) ---

  # A) Raices AR (Estabilidad)
  ar_coefs <- coefs[grepl("^ar", names(coefs))]
  raices_ar_inv <- NULL
  es_estable <- TRUE
  if (length(ar_coefs) > 0) {
    raices_ar_inv <- 1 / polyroot(c(1, -ar_coefs))
    es_estable <- all(Mod(raices_ar_inv) < 1)
  }

  # B) Raices MA (Invertibilidad)
  ma_coefs <- coefs[grepl("^ma", names(coefs))]
  raices_ma_inv <- NULL
  es_invertible <- TRUE
  if (length(ma_coefs) > 0) {
    raices_ma_inv <- 1 / polyroot(c(1, ma_coefs))
    es_invertible <- all(Mod(raices_ma_inv) < 1)
  }

  # --- 4. Diagnosticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # Ljung-Box (Ruido Blanco)
  lb_test <- stats::Box.test(residuos, type = "Ljung-Box", lag = p + q + 5)

  # Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- stats::shapiro.test(residuos)
  } else {
    shapiro_test <- list(p.value = NA)
  }

  # --- 5. Datos para Graficas ---
  ajustados <- x - residuos # En ARIMA esto devuelve la serie original ajustada
  sigma <- sqrt(modelo_base$sigma2)
  # Intervalos 95%
  lower_ci <- ajustados - 1.96 * sigma
  upper_ci <- ajustados + 1.96 * sigma

  df_plot <- data.frame(
    Tiempo = 1:length(x),
    Observado = as.numeric(x),
    Ajustado = as.numeric(ajustados),
    Lower = as.numeric(lower_ci),
    Upper = as.numeric(upper_ci),
    Residuos = as.numeric(residuos)
  )

  # Datos de la serie diferenciada (para ver si d fue efectivo)
  df_diff <- NULL
  if (d > 0) {
    serie_diff <- diff(x, differences = d)
    df_diff <- data.frame(Tiempo = 1:length(serie_diff), Valor = as.numeric(serie_diff))
  }

  # Datos combinados para circulo unitario
  df_raices <- data.frame(Real = numeric(), Imaginario = numeric(), Tipo = character())
  if (!is.null(raices_ar_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ar_inv), Imaginario = Im(raices_ar_inv), Tipo = "AR (Estabilidad)"))
  }
  if (!is.null(raices_ma_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ma_inv), Imaginario = Im(raices_ma_inv), Tipo = "MA (Invertibilidad)"))
  }

  # --- 6. Objeto Final ---
  resultado <- list(
    modelo = modelo_base,
    ordenes = c(p=p, d=d, q=q),
    coeficientes_tabla = tabla_coefs,
    estatus = list(estable = es_estable, invertible = es_invertible),
    diagnosticos = list(lb = lb_test, shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_diff = df_diff,
    data_raices = df_raices
  )

  class(resultado) <- "arima_greenreg"
  return(resultado)
}

#' @rdname modelo_arima
#' @param x Para \code{modelo_arima()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arima_greenreg"} devuelto por \code{modelo_arima()}.
#' @param ... Argumentos adicionales.
#' @export
print.arima_greenreg <- function(x, ...) {
  
  p <- x$ordenes['p']
  d <- x$ordenes['d']
  q <- x$ordenes['q']
  
  cat("\n==========================================================\n")
  cat(sprintf("MODELO INTEGRADO ARIMA(%d, %d, %d) \n", p, d, q))
  cat("==========================================================\n\n")
  
  # --- 1. MODELO Y ECUACION ---
  cat("--- 1. MODELO Y ECUACION ---\n")
  
  # Componentes para la descripcion
  parte_ar <- if(p > 0) paste0("AR(", p, ")") else ""
  parte_i  <- if(d > 0) paste0("I(", d, ")") else ""
  parte_ma <- if(q > 0) paste0("MA(", q, ")") else ""
  
  # USAMOS Filter para limpiar los espacios vacios
  estructura <- Filter(function(z) nchar(z) > 0, c(parte_ar, parte_i, parte_ma))
  
  cat("Estructura:", paste(estructura, collapse = " + "), "\n")
  
  if(d > 0) {
    cat(sprintf("Ecuacion con d = %d: Y(t) = f(Pasado, Choques)\n", d))
    cat(" NOTA: El modelo trabaja sobre la serie diferenciada para eliminar la tendencia.\n")
  } else {
    cat("Ecuacion: Y(t) = f(Pasado, Choques)\n")
  }
  
  aic_val <- x$diagnosticos$aic
  aic_txt <- if(is.null(aic_val) || is.na(aic_val)) "NA" else round(aic_val, 2)
  cat("AIC (Criterio de Akaike):", aic_txt, "\n\n")
  
  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  
  if (!is.null(x$coeficientes_tabla) && nrow(x$coeficientes_tabla) > 0) {
    stats::printCoefmat(as.matrix(x$coeficientes_tabla), digits = 4, signif.stars = TRUE)
  } else {
    cat("  (Modelo sin coeficientes estimados: Caminata aleatoria pura)\n")
  }
  cat("\n GUIA: Los coeficientes explican la relacion entre los datos ya diferenciados.\n\n")
  
  # --- 3. IDENTIFICACION DEL MODELO (p, d, q) ---
  cat("--- 3. IDENTIFICACION DEL MODELO ---\n")
  cat(sprintf(" p (Auto-Regresivo): %d periodos de memoria propia.\n", p))
  cat(sprintf(" d (Integracion):    %d diferenciacion(es) para lograr estacionariedad.\n", d))
  cat(sprintf(" q (Media Movil):    %d choques pasados que afectan el presente.\n", q))
  cat(" NOTA: Un modelo ARIMA es exitoso si 'd' estabiliza la media y 'p,q' eliminan la autocorrelacion.\n\n")
  
  # --- 4. VERIFICACION DE SUPUESTOS ---
  cat("--- 4. VERIFICACION DE SUPUESTOS ---\n")
  
  # A) Estacionariedad y Estructura
  cat("A) Estacionariedad y Estabilidad:\n")
  if(d > 0) cat(sprintf("   Estacionariedad: Lograda mediante %d diferenciacion(es).\n", d))
  
  if(p > 0) {
    est <- x$estatus$estable
    status_ar <- if(!is.null(est) && !is.na(est) && est) "ESTABLE" else "INESTABLE o Indeterminado"
    cat(sprintf("  Parte AR: %s (Raices Rojas dentro del circulo).\n", status_ar))
  }
  
  if(q > 0) {
    inv <- x$estatus$invertible
    status_ma <- if(!is.null(inv) && !is.na(inv) && inv) "INVERTIBLE" else "NO INVERTIBLE o Indeterminado"
    cat(sprintf("  Parte MA: %s (Rombos Azules dentro del circulo).\n", status_ma))
  }
  
  # B) Residuos ~ Ruido Blanco (E)
  cat("\nB) Diagnostico de Residuos (E ~ Ruido Blanco):\n")
  
  # Media Cero
  m_res <- mean(x$data_plot$Residuos, na.rm = TRUE)
  m_res_txt <- if(is.na(m_res)) "NA" else sprintf("%.4f", m_res)
  cat(sprintf("  Media de residuos: %s (Esperado: 0)\n", m_res_txt))
  
  # No Autocorrelacion
  lb_p <- x$diagnosticos$lb$p.value
  lb_txt <- if(is.null(lb_p) || is.na(lb_p)) "NA" else sprintf("%.4f", lb_p)
  cat(sprintf("  Independencia (Ljung-Box p-valor): %s\n", lb_txt))
  
  if (!is.null(lb_p) && !is.na(lb_p) && lb_p > 0.05) {
    cat("     -> OKEY: Los residuos son independientes (Ruido Blanco).\n")
  } else {
    cat("     ->  ALERTA: El modelo no capturo toda la informacion disponible.\n")
  }
  
  # Normalidad y Varianza
  sh_p <- x$diagnosticos$shapiro$p.value
  sh_txt <- if(is.null(sh_p) || is.na(sh_p)) "N/D" else sprintf("%.4f", sh_p)
  cat(sprintf("  Normalidad (Shapiro p-valor): %s\n", sh_txt))
  
  # --- 5. RESUMEN DEL MODELO ---
  cat("\n--- 5. RESUMEN DEL MODELO ---\n")
  
  # Logica de validacion blindada contra NAs
  ar_ok <- if (p == 0) TRUE else (!is.null(x$estatus$estable) && !is.na(x$estatus$estable) && x$estatus$estable)
  ma_ok <- if (q == 0) TRUE else (!is.null(x$estatus$invertible) && !is.na(x$estatus$invertible) && x$estatus$invertible)
  lb_ok <- (!is.null(lb_p) && !is.na(lb_p) && lb_p > 0.05)
  
  if(ar_ok && ma_ok && lb_ok) {
    cat(" MODELO VALIDO: El ARIMA ha sido correctamente identificado y ajustado.\n")
  } else {
    cat(" REVISION NECESARIA: El modelo presenta problemas de estabilidad o residuos sucios.\n")
  }
  
  cat("\n Usa plot(modelo) para ver las 9 graficas, incluyendo la Serie Diferenciada.\n")
  cat("==========================================================\n")
  
  return(invisible(x))
}

#' @rdname modelo_arima 
#' @description Genera un panel visual para evaluar el ajuste del modelo.
#' @param x Para \code{modelo_arima()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arima_greenreg"} devuelto por \code{modelo_arima()}.
#' @param ... Argumentos adicionales.
#' @export
plot.arima_greenreg <- function(x, ...) {
  
  # --- 0. EXTRACCION DE PARAMETROS ---
  df <- x$data_plot
  p  <- x$ordenes['p']
  d  <- x$ordenes['d']
  q  <- x$ordenes['q']
  
  # TEMA
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, face = "italic",
                                           color = "#555555", margin = ggplot2::margin(t = 15),
                                           lineheight = 1.1),
      plot.margin = ggplot2::margin(10, 10, 20, 10),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
  
  # CORRECCION 1: Usar message()
  if (interactive()) {
    message("Iniciando diagnostico ARIMA... (10 graficas totales)")
  }
  
  # --- GRAFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---
  
  media_orig <- mean(df$Observado, na.rm = TRUE)
  
  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    # Linea de la serie original
    ggplot2::geom_line(color = "#34495E", linewidth = 0.7) +
    # Puntos sutiles
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    # Linea de la media (Referencia)
    ggplot2::geom_hline(yintercept = media_orig, linetype = "dashed",
                        color = "#8E44AD", linewidth = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualizacion de la serie temporal bruta (sin diferenciar)",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA PURPURA: Promedio historico general.\n",
                       " ANALISIS: Si la serie sube o baja constantemente (tendencia),\n",
                       "   el valor de d debe ser mayor a 0 para estabilizarla.\n",
                       " NOTA: El ARIMA usara d=", d, " para corregir lo que vemos aqui.")
    ) +
    mi_tema
  
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1: Inspeccion de la serie original. > ")}
  
  # --- GRAFICA 2: AJUSTE DEL MODELO ARIMA (Real vs Ajustado) ---
  
  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    # Sombra de Confianza (IC 95%) en Purpura
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#8E44AD", alpha = 0.15) +
    # Serie Observada (Realidad)
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), linewidth = 0.7) +
    # Serie Ajustada (Modelo ARIMA)
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo ARIMA"), linewidth = 0.8) +
    # Escala de colores
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo ARIMA" = "#8E44AD")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo ARIMA(", p, ",", d, ",", q, ")"),
      subtitle = "Serie Observada vs. Valores Ajustados (Fitted)",
      x = "Tiempo (Periodos)",
      y = "Valor",
      color = "Referencia:",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: La linea purpura debe seguir la trayectoria de la oscura.\n",
                       " SOMBRA: Rango donde el modelo espera que caigan los datos reales.\n",
                       " NOTA: El ARIMA integra la diferenciacion (d) para seguir tendencias\n",
                       "   y los componentes (p,q) para ajustar las oscilaciones.")
    ) +
    mi_tema
  
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2: Evaluacion del ajuste integrado. > ")}
  
  # --- GRAFICA 3: SERIE DIFERENCIADA (ESTACIONARIEDAD LOGRADA) ---
  
  if (!is.null(x$data_diff)) {
    # Usamos el dataframe de la serie diferenciada calculado en la funcion principal
    df_diff <- x$data_diff
    media_diff <- mean(df_diff$Valor, na.rm = TRUE)
    
    g3 <- ggplot2::ggplot(df_diff, ggplot2::aes(x = Tiempo, y = Valor)) +
      # Linea de la serie transformada (Diferenciada)
      ggplot2::geom_line(color = "#8E44AD", linewidth = 0.6) +
      # Linea de referencia en Cero
      ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
      # Linea de la media de la diferencia
      ggplot2::geom_hline(yintercept = media_diff, linetype = "dashed",
                          color = "#E74C3C", linewidth = 0.8) +
      ggplot2::labs(
        title = paste0("3. Serie Diferenciada (d = ", d, ")"),
        subtitle = "Transformacion para lograr Estacionariedad",
        x = "Tiempo (Periodos tras diferencia)",
        y = paste0("^", d, " Valor"),
        caption = paste0("INTERPRETACION:\n",
                         " BIEN: La serie debe oscilar aleatoriamente alrededor de la linea roja.\n",
                         " LINEA ROJA: Media de los cambios. Si es cercana a cero, no hay deriva.\n",
                         " NOTA: Esta es la serie que los componentes AR(", p, ") y MA(", q, ")\n",
                         "   estan intentando modelar realmente.")
      ) +
      mi_tema
    
    print(g3)
    
  } else {
    # Caso donde d = 0 (Modelo es un ARMA simple)
    # CORRECCION 1: message en lugar de cat
    if (interactive()) {
      message(" d = 0: No hay serie diferenciada para mostrar. La serie original ya es estacionaria.")
    }
  }
  if (interactive()) {readline(prompt = "Grafica 3: Verificacion del efecto de la diferenciacion. > ")}
  
  
  # --- GRAFICA 4: ESTRUCTURA MATEMATICA (CIRCULO UNITARIO) ---
  
  # CORRECCION 4: Validacion defensiva de data_raices
  if (!is.null(x$data_raices) && nrow(x$data_raices) > 0) {
    # 1. Definicion del circulo de referencia
    theta_circ <- seq(0, 2*pi, length.out = 200)
    df_circle_arima <- data.frame(x = cos(theta_circ), y = sin(theta_circ))
    
    # 2. Construccion del grafico dual de raices
    g4 <- ggplot2::ggplot() +
      # Circulo unitario (Limite de estabilidad e invertibilidad)
      ggplot2::geom_path(data = df_circle_arima, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", linewidth = 0.8) +
      # area sombreada interna (Zona de validez)
      ggplot2::geom_polygon(data = df_circle_arima, ggplot2::aes(x, y),
                            fill = "#8E44AD", alpha = 0.05) +
      # Ejes de coordenadas
      ggplot2::geom_vline(xintercept = 0, color = "gray85") +
      ggplot2::geom_hline(yintercept = 0, color = "gray85") +
      # Raices AR (Rojas) y MA (Azules)
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario, color = Tipo, shape = Tipo),
                          size = 4, stroke = 1.5) +
      # Colores y formas diferenciadas
      ggplot2::scale_color_manual(values = c("AR (Estabilidad)" = "#E74C3C",
                                             "MA (Invertibilidad)" = "#2980B9")) +
      ggplot2::scale_shape_manual(values = c("AR (Estabilidad)" = 16,
                                             "MA (Invertibilidad)" = 18)) +
      # Proporcion fija para que el circulo no parezca un ovalo
      ggplot2::coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
      ggplot2::labs(
        title = "4. Estructura Matematica (Circulo Unitario)",
        subtitle = paste0("Raices del Polinomio ARIMA(", p, ",", d, ",", q, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        color = "Componente:", shape = "Componente:",
        caption = paste0("INTERPRETACION:\n",
                         " PUNTOS ROJOS (AR): Estabilidad del sistema dinamico.\n",
                         " ROMBOS AZULES (MA): Invertibilidad (unicidad del modelo).\n",
                         " CRITICO: Si una raiz toca el borde gris, el modelo es casi\n",
                         "   no-estacionario (posiblemente d requiere ser mayor).")
      ) +
      mi_tema + ggplot2::theme(legend.position = "right")
    
    print(g4)
    
  } else {
    # CORRECCION 1: message en lugar de cat
    if (interactive()) {
      message(paste0(" ARIMA(", p, ",", d, ",", q, ") sin componentes AR/MA para graficar raices."))
    }
  }
  if (interactive()) {readline(prompt = "Grafica 4: Analisis de la convergencia del modelo. > ")}
  
  # --- GRAFICA 5: RESIDUOS EN EL TIEMPO (E vs t) ---
  
  # CORRECCION 2: stats::sd
  sd_res_arima <- stats::sd(df$Residuos, na.rm = TRUE)
  
  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    # Banda de confianza (+-2 SD) en Purpura tenue
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_arima, ymax = 2*sd_res_arima),
                         fill = "#8E44AD", alpha = 0.1) +
    # Linea de los residuos (errores)
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.5) +
    # Puntos para detectar shocks o outliers
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.4, size = 1) +
    # Linea central de referencia (Media Cero)
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Lineas criticas de control (+-2 SD)
    ggplot2::geom_hline(yintercept = c(-2*sd_res_arima, 2*sd_res_arima),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.7) +
    ggplot2::labs(
      title = "5. Residuos en el Tiempo (E vs t)",
      subtitle = paste0("Analisis de Independencia para ARIMA(", p, ",", d, ",", q, ")"),
      x = "Tiempo (Periodos)",
      y = "Error Residual",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos deben ser una nube aleatoria sin patron aparente.\n",
                       " ROJO: Puntos fuera de las lineas punteadas son 'shocks' imprevistos.\n",
                       " NOTA: Si ves que la nube se ensancha o se estrecha (Varianza no constante),\n",
                       "   el supuesto de Homocedasticidad se esta violando.")
    ) +
    mi_tema
  
  print(g5)
  if (interactive()) {readline(prompt = "Grafica 5: Inspeccion visual de los errores del modelo. > ")}
  
  # --- GRAFICA 6: AUTOCORRELACION DE RESIDUOS (ACF) ---
  
  # 1. Calculo de la ACF (Lags 1 a 20)
  acf_res_arima <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
  df_acf_arima <- data.frame(
    Lag = as.numeric(acf_res_arima$lag[-1]),
    ACF = as.numeric(acf_res_arima$acf[-1])
  )
  
  # Limite de significancia estadistica (95%)
  limite_arima <- 1.96 / sqrt(nrow(df))
  
  # Logica de deteccion: Resaltar fallos en Rojo
  df_acf_arima$Etiqueta <- ifelse(abs(df_acf_arima$ACF) > limite_arima, as.character(df_acf_arima$Lag), "")
  df_acf_arima$vjust_pos <- ifelse(df_acf_arima$ACF > 0, -0.5, 1.5)
  
  # 2. Construccion de la Grafica (Identidad Purpura ARIMA)
  g6 <- ggplot2::ggplot(df_acf_arima, ggplot2::aes(x = Lag, y = ACF)) +
    # Barras de correlacion (Gris oscuro para sobriedad tecnica)
    ggplot2::geom_col(fill = "#34495E", width = 0.5) +
    # Lineas de significancia (Purpura ARIMA)
    ggplot2::geom_hline(yintercept = c(limite_arima, -limite_arima),
                        linetype = "dashed", color = "#8E44AD", linewidth = 0.8) +
    # Linea base
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Alertas visuales (Numeros de Lag que fallan)
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelacion de Residuos (ACF)",
      subtitle = "Se ha extraido toda la dependencia temporal?",
      x = "Lag (Retraso)",
      y = "Coeficiente ACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Si todas las barras estan dentro de las lineas purpuras.\n",
                       " ROJO: El lag indicado aun tiene informacion no capturada por el modelo.\n",
                       " NOTA: El ARIMA(", p, ",", d, ",", q, ") busca que el error presente\n",
                       "   sea totalmente independiente de los errores pasados.")
    ) +
    mi_tema
  
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6: Validacion de independencia estadistica. > ")}
  
  # --- GRAFICA 7: AUTOCORRELACION PARCIAL DE RESIDUOS (PACF) ---
  
  # 1. Calculo de la PACF (Lags 1 a 20)
  pacf_res_arima <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)
  df_pacf_arima <- data.frame(
    Lag = as.numeric(pacf_res_arima$lag),
    PACF = as.numeric(pacf_res_arima$acf)
  )
  
  # Limite de significancia (95%)
  limite_arima <- 1.96 / sqrt(nrow(df))
  
  # Logica de deteccion: Alertas en Rojo
  df_pacf_arima$Etiqueta <- ifelse(abs(df_pacf_arima$PACF) > limite_arima, as.character(df_pacf_arima$Lag), "")
  df_pacf_arima$vjust_pos <- ifelse(df_pacf_arima$PACF > 0, -0.5, 1.5)
  
  # 2. Construccion de la Grafica
  g7 <- ggplot2::ggplot(df_pacf_arima, ggplot2::aes(x = Lag, y = PACF)) +
    # Barras de la PACF (Color Purpura ARIMA)
    ggplot2::geom_col(fill = "#8E44AD", width = 0.5) +
    # Lineas de significancia (Naranja para precaucion tecnica)
    ggplot2::geom_hline(yintercept = c(limite_arima, -limite_arima),
                        linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
    # Linea base
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Alertas visuales
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "7. Autocorrelacion Parcial de Residuos (PACF)",
      subtitle = "Identificacion de dependencias directas residuales",
      x = "Lag (Retraso)",
      y = "Coeficiente PACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Todas las barras dentro de las lineas naranjas.\n",
                       " NUMEROS ROJOS: Indican que el componente AR(", p, ") no fue suficiente.\n",
                       " NOTA: El PACF limpia el efecto de los lags intermedios para ver la\n",
                       "   conexion pura entre el presente y el pasado del error.")
    ) +
    mi_tema
  
  print(g7)
  if (interactive()) {readline(prompt = "Grafica 7: Diagnostico final de la estructura del error. > ")}
  
  # --- GRAFICA 8: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---
  
  g8 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
    # Puntos de cuantiles observados (Purpura ARIMA)
    ggplot2::stat_qq(color = "#8E44AD", alpha = 0.6, size = 2) +
    # Linea de referencia teorica (Distribucion Normal)
    ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "8. Normalidad de los Errores (Q-Q Plot)",
      subtitle = "Distribucion de Residuos vs. Cuantiles Teoricos Normales",
      x = "Cuantiles Teoricos (Normal)",
      y = "Cuantiles Observados (Residuos ARIMA)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos deben estar alineados sobre la diagonal punteada.\n",
                       " COLAS: Si los extremos se alejan, hay valores extremos no capturados.\n",
                       " NOTA: La normalidad permite que las inferencias estadisticas sean exactas.\n",
                       " CLAVE: Si hay mucha curvatura, el modelo ARIMA puede estar sesgado.")
    ) +
    mi_tema
  
  print(g8)
  if (interactive()) {readline(prompt = "Grafica 8: Verificacion de la distribucion de los errores. > ")}
  
  
  # --- GRAFICA 9: PRECISION DEL MODELO (REAL VS AJUSTADO) ---
  
  # Calculamos la correlacion de Pearson para el reporte de precision
  cor_arima <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")
  
  g9 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    # Puntos de dispersion (Color Purpura ARIMA)
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.5, size = 2) +
    # Linea de Identidad (Ideal: Real = Predicho)
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    # Linea de tendencia ajustada (Regresion de validacion)
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "9. Precision del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlacion de Pearson: ", round(cor_arima, 3)),
      x = "Valor Ajustado (Modelo ARIMA)",
      y = "Valor Observado (Realidad)",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA ROJA: El escenario ideal de prediccion perfecta.\n",
                       " BIEN: Los puntos purpuras deben seguir la trayectoria de la linea roja.\n",
                       " DISPERSION: Indica la incertidumbre del modelo ante choques aleatorios.\n",
                       " NOTA: El ARIMA integra la tendencia para mejorar esta alineacion.")
    ) +
    mi_tema
  
  print(g9)
  if (interactive()) {readline(prompt = "Grafica 9: Evaluacion de la capacidad de ajuste integrada. > ")}
  
  # --- GRAFICA 10: DIAGNOSTICO DE ESTACIONARIEDAD FINAL ---
  
  # 1. Definir ventana (10% de los datos o minimo 5)
  ventana_arima <- max(5, round(nrow(df) * 0.1))
  
  # 2. Calculos sobre la Serie Original (para ver la estabilidad del ajuste)
  df_est_arima <- df
  df_est_arima$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_arima, ventana_arima), sides = 2))
  
  # CORRECCION 2: stats::sd
  df_est_arima$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_arima):min(nrow(df), i+ventana_arima)
    stats::sd(df$Observado[idx], na.rm = TRUE)
  })
  
  # 3. FILTRO (Eliminamos NAs de los extremos para un reporte limpio)
  df_final_arima <- df_est_arima[!is.na(df_est_arima$Media_Movil), ]
  
  # 4. CONSTRUCCION DE LA GRAFICA (Estilo Purpura ARIMA)
  g10 <- ggplot2::ggplot(df_final_arima, ggplot2::aes(x = Tiempo)) +
    # Serie original completa de fondo (Atenuada)
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray85", linewidth = 0.5) +
    
    # Franja de Varianza Movil (Purpura ARIMA)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#8E44AD", alpha = 0.15) +
    
    # Linea de la Media Movil (Nivel local)
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#8E44AD", linewidth = 1.2) +
    
    # Linea de la Media Global (Referencia fija en Rojo)
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
    
    ggplot2::labs(
      title = "10. Diagnostico de Estacionariedad Final",
      subtitle = paste0("Analisis Dinamico de Media y Varianza (Ventana k = ", ventana_arima, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA PURPURA: Media local. En un ARIMA con d=", d, ", debe seguir la tendencia.\n",
                       " LINEA ROJA: Promedio global de la serie original.\n",
                       " BIEN: Si la franja purpura tiene un ancho constante (Varianza estable).\n",
                       " NOTA: Si el ancho de la franja crece, el modelo ARIMA requiere transformacion Log.")
    ) +
    mi_tema
  
  print(g10)
  
  # CORRECCION 1: Usar message() en lugar de cat()
  if (interactive()) {
    message("Diagnostico ARIMA finalizado con exito (10/10 graficas).")
  }
  
  # CORRECCION 3: Retornar el objeto original (x)
  return(invisible(x))
}