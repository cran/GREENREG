#' Modelo Autorregresivo (AR)
#'
#' Esta funcion ajusta un modelo de series de tiempo donde el valor actual \eqn{(Y_t)}
#' se explica mediante una combinacion lineal de sus propios valores pasados
#' \eqn{(Y_{t-1}, Y_{t-2}, \dots)}. Es la herramienta estandar para analizar la
#' "persistencia" o "memoria" de un fenomeno a traves del tiempo.
#'
#' @details
#' El flujo interno de la funcion se divide en cinco etapas diagnosticas:
#' \enumerate{
#'   \item **Estimacion (MLE):** Utiliza Maxima Verosimilitud para encontrar los
#'         coeficientes \eqn{\phi} (Phi) que mejor describen la dependencia temporal.
#'   \item **Analisis de Estabilidad:** Calcula las raices inversas del polinomio
#'         caracteristico. Un modelo AR solo es valido si es "Estacionario", es decir,
#'         si las raices se encuentran dentro del circulo unitario.
#'   \item **Inferencia de Parametros:** Calcula errores estandar y p-valores para
#'         determinar cuantos rezagos (lags) son realmente significativos.
#'   \item **Auditoria de Residuos:** Evalua si el error restante es "Ruido Blanco"
#'         mediante la prueba de Ljung-Box (independencia) y Shapiro-Wilk (normalidad).
#'   \item **Traduccion Didactica:** Genera un reporte que explica si el modelo
#'         "explota" (inestable) o si ha capturado exitosamente la memoria de la serie.
#' }
#'
#' @section Guia de las 9 Graficas de Diagnostico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una auditoria temporal completa:
#' \describe{
#'   \item{\strong{1. Inspeccion Visual}}{Muestra la serie bruta con su media historica. Ayuda a detectar tendencias que el modelo AR no puede manejar por si solo.}
#'   \item{\strong{2. Ajuste Temporal}}{Superpone el modelo sobre la realidad. La sombra roja indica la incertidumbre (Intervalo de Confianza).}
#'   \item{\strong{3. Residuos en el Tiempo}}{Un "electrocardiograma" del error. No debe presentar ondas, ciclos ni tendencias.}
#'   \item{\strong{4. Circulo Unitario}}{Grafica matematica critica. Si los rombos rojos salen del circulo, el modelo es inestable y no sirve para pronosticar.}
#'   \item{\strong{5. ACF de Residuos}}{Busca correlaciones remanentes. Las barras deben ser pequenas para confirmar que el error es ruido aleatorio.}
#'   \item{\strong{6. PACF de Residuos}}{Ayuda a confirmar si el orden \code{p} seleccionado fue suficiente para limpiar la dependencia directa.}
#'   \item{\strong{7. Normalidad (Q-Q Plot)}}{Valida que los choques aleatorios sigan una distribucion normal, requisito para que los intervalos de confianza sean exactos.}
#'   \item{\strong{8. Precision (Real vs Ajustado)}}{Muestra la correlacion directa entre la prediccion y el dato observado.}
#'   \item{\strong{9. Estacionariedad Movil}}{Analiza si la media y la varianza cambian a lo largo de la serie. Una linea verde horizontal confirma un modelo robusto.}
#' }
#'
#' @param x Vector numerico o serie de tiempo (\code{ts}). No debe contener valores \code{NA}.
#' @param p Orden del modelo (numero de rezagos a considerar). Por defecto es 1.
#' @param nivel_confianza Nivel de significancia para las bandas de error (default 0.95).
#'
#' @return Un objeto de clase \code{"ar_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{arima} original.
#'   \item \code{orden_p}: El orden del modelo ajustado.
#'   \item \code{coeficientes}: Matriz con estimaciones, errores estandar y p-valores.
#'   \item \code{estabilidad}: Lista con las raaces inversas y diagnostico de estacionariedad.
#'   \item \code{diagnosticos}: Resultados de las pruebas de Ljung-Box, Shapiro-Wilk y AIC.
#'   \item \code{notas}: Reporte interpretado sobre estabilidad y ruido blanco.
#'   \item \code{data_plot}: Dataframe listo para visualización avanzada.
#' }
#' 
#' @examples
#'
#' data("datos_anomalia_temperatura")
#' modelo_AR <- modelo_ar(datos_anomalia_temperatura$anomalia_c, p = 2)
#' modelo_AR
#' plot(modelo_AR)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf filter sd
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
modelo_ar <- function(x, p = 1, nivel_confianza = 0.95) {

  # --- 1. Validacion y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numerico o una serie de tiempo (ts).")

  # Ajuste usando ARIMA con q=0, d=0 (es decir, un AR puro)
  modelo_base <- stats::arima(x, order = c(p, 0, 0))

  # --- 2. Coeficientes e Interpretacion ---
  coefs <- coef(modelo_base)
  # Error estandar aproximado (matriz de varianza-covarianza)
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- data.frame(
    Estimacion = coefs,
    Std.Error = se,
    z.value = z_val,
    P.valor = p_val
  )

  # --- 3. Analisis de Estabilidad (Raices) ---
  # El modelo AR es estable si las raices inversas del polinomio estan dentro del circulo unitario
  # Extraemos solo los coeficientes AR (sin intercepto)
  ar_coefs <- coefs[grepl("^ar", names(coefs))]

  if (length(ar_coefs) > 0) {
    poly_roots <- polyroot(c(1, -ar_coefs))
    inverse_roots <- 1 / poly_roots
    modulos <- Mod(inverse_roots)
    es_estable <- all(modulos < 1)
  } else {
    inverse_roots <- NULL
    es_estable <- TRUE # AR(0) es ruido blanco, siempre estable
  }

  # --- 4. Diagnosticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # A) Test Ljung-Box (Independencia / Ruido Blanco)
  lb_test <- stats::Box.test(residuos, type = "Ljung-Box", lag = p + 5)
  es_ruido_blanco <- lb_test$p.value > 0.05

  # B) Test Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- stats::shapiro.test(residuos)
    es_normal <- shapiro_test$p.value > 0.05
  } else {
    shapiro_test <- list(p.value = NA)
    es_normal <- NA
  }

  # --- 5. Notas ---
  nota_estabilidad <- if (es_estable) {
    " CONCLUSION: El modelo es ESTABLE (Estacionario). Puede usarse para predecir."
  } else {
    " CONCLUSION: El modelo es INESTABLE (No Estacionario). No es confiable a largo plazo."
  }

  nota_lb <- if (es_ruido_blanco) {
    " CONCLUSION: Ruido Blanco (OK). No queda informacion util en los errores."
  } else {
    " CONCLUSION: Existe Autocorrelacion. El modelo no capturo toda la dinamica (Faltan lags)."
  }

  nota_norm <- if (!is.na(es_normal) && es_normal) {
    " CONCLUSION: Los errores son Normales (OK). Los intervalos de confianza son fiables."
  } else {
    " CONCLUSION: Los errores NO son normales. Precaucion con los intervalos."
  }

  # Interpretacion del primer lag
  phi1 <- if(length(ar_coefs) > 0) ar_coefs[1] else 0
  nota_phi <- paste0("Persistencia: El coeficiente AR1 es ", round(phi1, 2), ". ",
                     "Esto indica que tanto del valor de ayer (t-1) se transfiere a hoy.")

  # --- 6. Preparar Datos para Graficas ---
  # Valores ajustados: Real - Residuo
  ajustados <- x - residuos

  # Calculo de Intervalos de Confianza (Shadows)
  sigma <- sqrt(modelo_base$sigma2)
  z_score <- stats::qnorm(1 - (1 - nivel_confianza) / 2)
  lower_ci <- ajustados - z_score * sigma
  upper_ci <- ajustados + z_score * sigma

  df_plot <- data.frame(
    Tiempo = 1:length(x),
    Observado = as.numeric(x),
    Ajustado = as.numeric(ajustados),
    Lower = as.numeric(lower_ci),
    Upper = as.numeric(upper_ci),
    Residuos = as.numeric(residuos)
  )

  # Datos para circulo unitario (si hay raices)
  df_raices <- NULL
  if (!is.null(inverse_roots)) {
    df_raices <- data.frame(
      Real = Re(inverse_roots),
      Imaginario = Im(inverse_roots)
    )
  }

  resultado <- list(
    x_original = x,
    modelo = modelo_base,
    orden_p = p,
    coeficientes = tabla_coefs,
    estabilidad = list(es_estable = es_estable, raices_inv = inverse_roots),
    diagnosticos = list(LjungBox = lb_test, Shapiro = shapiro_test, AIC = modelo_base$aic),
    notas = list(est = nota_estabilidad, lb = nota_lb, norm = nota_norm, phi = nota_phi),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "ar_greenreg"
  return(resultado)
}

#' @rdname modelo_ar
#' @param x Objeto de clase.
#' @param ... Argumentos adicionales.
#' @export
print.ar_greenreg <- function(x, ...) {
  
  cat("\n==========================================================\n")
  cat("     MODELO AUTORREGRESIVO AR(", x$orden_p, ") \n", sep = "")
  cat("==========================================================\n\n")
  
  # --- 1. MODELO Y ECUACION ---
  cat("--- 1. MODELO Y ECUACION ---\n")
  phi <- round(x$modelo$coef[grepl("^ar", names(x$modelo$coef))], 3)
  c_int <- round(x$modelo$coef["intercept"], 3)
  
  # Construccion visual de la ecuacion de diferencias
  eq <- paste0("Y(t) = ", if(is.na(c_int)) 0 else c_int)
  if(length(phi) > 0) {
    for(i in 1:length(phi)) {
      eq <- paste0(eq, " + (", phi[i], " * Y(t-", i, "))")
    }
  }
  eq <- paste0(eq, " + E(t)")
  
  cat("Ecuacion de Diferencias Estimada:\n", eq, "\n")
  
  aic_val <- x$diagnosticos$AIC
  aic_txt <- if(is.null(aic_val) || is.na(aic_val)) "NA" else round(aic_val, 2)
  cat("AIC (Criterio de Akaike):", aic_txt, "\n")
  cat(" NOTA: El intercepto representa el valor medio de la serie si es estacionaria.\n\n")
  
  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS (Parametros Phi) ---\n")
  stats::printCoefmat(as.matrix(x$coeficientes), digits = 4, signif.stars = TRUE)
  cat("\n", x$notas$phi, "\n")
  cat(" GUIA: Si |Phi| es cercano a 1, la serie tiene mucha memoria (persistencia).\n\n")
  
  # --- 3. IDENTIFICACION DEL MODELO (ORDEN p) ---
  cat("--- 3. IDENTIFICACION DEL MODELO (Orden p) ---\n")
  cat(" Orden seleccionado: p =", x$orden_p, "\n")
  cat(" Criterio ACF: Las autocorrelaciones decaen exponencialmente o de forma sinoidal.\n")
  cat(" Criterio PACF: Los coeficientes de autocorrelacion parcial se cortan despues del lag", x$orden_p, ".\n")
  cat(" NOTA: Si ves barras significativas en PACF despues de p, podrias necesitar un orden mayor.\n\n")
  
  # --- 4. VERIFICACION DE SUPUESTOS ---
  cat("--- 4. VERIFICACION DE SUPUESTOS ---\n")
  
  # A) Estacionariedad (Estabilidad)
  cat("A) Estacionariedad (Raices Unitarias):\n")
  if (!is.null(x$estabilidad$raices_inv)) {
    modulos <- Mod(x$estabilidad$raices_inv)
    cat("   - Modulos de raices inversas:", paste(round(modulos, 3), collapse = ", "), "\n")
  }
  cat("   -", x$notas$est, "\n")
  
  # B) Errores ~ Ruido Blanco
  cat("\nB) Diagnostico de Errores (E ~ RB):\n")
  
  # Media Cero (Protegido)
  mean_res <- mean(x$data_plot$Residuos, na.rm = TRUE)
  mean_txt <- if(is.na(mean_res)) "NA" else sprintf("%.4f", mean_res)
  cat(sprintf("  Media de los residuos: %s (Esperado: 0)\n", mean_txt))
  
  # No Autocorrelacion (Ljung-Box Protegido)
  lb_pval <- x$diagnosticos$LjungBox$p.value
  lb_txt <- if(is.null(lb_pval) || is.na(lb_pval)) "NA" else sprintf("%.4f", lb_pval)
  cat(sprintf("  No Autocorrelacion (Ljung-Box p-valor): %s\n", lb_txt))
  cat("     ->", x$notas$lb, "\n")
  
  # Varianza Constante y Normalidad (Protegido)
  sh_pval <- x$diagnosticos$Shapiro$p.value
  sh_txt <- if(is.null(sh_pval) || is.na(sh_pval)) "NA" else sprintf("%.4f", sh_pval)
  cat(sprintf("  Normalidad (Shapiro p-valor): %s\n", sh_txt))
  cat("     ->", x$notas$norm, "\n\n")
  
  # --- 5. RESUMEN DEL MODELO ---
  cat("--- 5. RESUMEN DEL MODELO ---\n")
  
  # Validacion logica segura para evitar colapsos por NA
  es_valido <- FALSE
  if (!is.null(x$estabilidad$es_estable) && !is.null(lb_pval) && !is.na(lb_pval)) {
    if (x$estabilidad$es_estable && lb_pval > 0.05) {
      es_valido <- TRUE
    }
  }
  
  if(es_valido) {
    cat("RESULTADO FINAL: El modelo AR(", x$orden_p, ") es VALIDO.\n", sep="")
    cat("   Los residuos son ruido blanco y la serie es estacionaria.\n")
  } else {
    cat(" ADVERTENCIA: El modelo presenta deficiencias tecnicas.\n")
    cat("   Revisa si es necesario diferenciar la serie (d=1) o cambiar el orden p.\n")
  }
  
  cat("\n Usa plot(modelo) para visualizar la serie, ACF y raices unitarias.\n")
  cat("==========================================================\n")
  
  return(invisible(x))
}
#' @rdname modelo_ar
#' @description Genera una secuencia de 9 graficas de diagnostico para validar el modelo AR.
#' @param x Objeto de clase ar_greenreg.
#' @param ... Argumentos adicionales.
#' @export
plot.ar_greenreg <- function(x, ...) {
  
  df <- x$data_plot
  
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
    message("Generando graficas de diagnostico AR... (Presiona [Enter] para avanzar)\n")
  }
  
  # --- GRAFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---
  
  # CORRECCION 2: Usar stats::sd
  media_serie <- mean(df$Observado, na.rm = TRUE)
  sd_serie <- stats::sd(df$Observado, na.rm = TRUE)
  
  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = media_serie - sd_serie, ymax = media_serie + sd_serie,
                      fill = "#34495E", alpha = 0.1) +
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.7) +
    ggplot2::geom_point(color = "#2C3E50", alpha = 0.3, size = 1) +
    ggplot2::geom_hline(yintercept = media_serie, linetype = "dashed",
                        color = "#E74C3C", linewidth = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Analisis de Trayectoria y Estacionariedad Visual",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA ROJA: Es el valor promedio historico de tu serie.\n",
                       " ESTACIONARIA: Si la serie cruza constantemente la linea roja y no tiene tendencia.\n",
                       " TENDENCIA: Si la serie se aleja de la linea roja de forma sostenida (requiere d=1).\n",
                       " NOTA: El area gris representa la zona de fluctuacion normal (+-1 Desviacion Estandar).")
    ) +
    mi_tema
  
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1: Inspeccion visual de la serie bruta > ")}
  
  
  # --- GRAFICA 2: AJUSTE DEL MODELO ---
  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#E74C3C", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo AR"), linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo AR" = "#E74C3C")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo AR(", x$orden_p, ")"),
      subtitle = "Serie Observada vs. Valores Ajustados en el Tiempo",
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      color = "Leyenda:",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: La linea roja (Modelo) debe seguir la trayectoria de la oscura (Realidad).\n",
                       " SOMBRA ROJA: Representa el intervalo de confianza al 95%.\n",
                       " NOTA: Si el modelo esta muy 'liso' y los datos son muy volatiles, falta orden p.")
    ) +
    mi_tema
  
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2: Comparacion temporal del ajuste > ")}
  
  
  # --- GRAFICA 3: RESIDUOS EN EL TIEMPO (E vs t) ---
  
  # CORRECCION 2: Usar stats::sd
  sd_res <- stats::sd(df$Residuos, na.rm = TRUE)
  
  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res, ymax = 2*sd_res),
                         fill = "#3498DB", alpha = 0.1) +
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.5) +
    ggplot2::geom_point(color = "#3498DB", alpha = 0.4, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(-2*sd_res, 2*sd_res),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.6) +
    ggplot2::labs(
      title = "3. Residuos en el Tiempo",
      subtitle = "Buscando Aleatoriedad (E ~ Ruido Blanco)",
      x = "Tiempo (Periodos)",
      y = "Error (Residuo)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos oscilan al azar alrededor de cero (sin patrones).\n",
                       " BANDAS ROJAS: El 95% de los errores deberian estar dentro de estos limites.\n",
                       " MAL: Si ves una 'U', una tendencia o rachas largas arriba/abajo del cero.\n",
                       " CLAVE: Si hay patrones, el modelo AR(", x$orden_p, ") es insuficiente.")
    ) +
    mi_tema
  
  print(g3)
  if (interactive()) {readline(prompt = "Grafica 3: Analisis de aleatoriedad de los errores > ")}
  
  
  # --- GRAFICA 4: ESTABILIDAD (CIRCULO UNITARIO) ---
  if (!is.null(x$data_raices)) {
    theta <- seq(0, 2*pi, length.out = 200)
    df_circle <- data.frame(x = cos(theta), y = sin(theta))
    
    g4 <- ggplot2::ggplot() +
      ggplot2::geom_path(data = df_circle, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, color = "#BDC3C7") +
      ggplot2::geom_hline(yintercept = 0, color = "#BDC3C7") +
      ggplot2::geom_polygon(data = df_circle, ggplot2::aes(x, y),
                            fill = "#2ECC71", alpha = 0.05) +
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario),
                          color = "#E74C3C", size = 4, shape = 18) +
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::labs(
        title = "4. Estabilidad (Circulo Unitario)",
        subtitle = paste0("Raices Inversas del Polinomio AR(", x$orden_p, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        caption = paste0("INTERPRETACION:\n",
                         " ESTABLE: Todas las raices (rombos rojos) estan DENTRO del circulo.\n",
                         " ALERTA: Si una raiz toca el borde, la serie tiene una raiz unitaria (necesita d=1).\n",
                         " INESTABLE: Si una raiz sale del circulo, el modelo no converge y fallara al predecir.")
      ) +
      mi_tema
    
    print(g4)
    
  } else {
    # CORRECCION 1: Usar message() en lugar de cat()
    message(" El modelo AR(0) no tiene raices; es ruido blanco puro y siempre es estable.")
  }
  if (interactive()) {readline(prompt = "Grafica 4: Analisis de raices > ")}
  
  
  # --- GRAFICA 5: AUTOCORRELACION DE RESIDUOS (ACF) ---
  acf_res <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
  
  df_acf <- data.frame(
    Lag = as.numeric(acf_res$lag[-1]),
    ACF = as.numeric(acf_res$acf[-1])
  )
  
  n_obs <- nrow(df)
  limite <- 1.96 / sqrt(n_obs)
  
  df_acf$Etiqueta <- ifelse(abs(df_acf$ACF) > limite, as.character(df_acf$Lag), "")
  df_acf$vjust_pos <- ifelse(df_acf$ACF > 0, -0.5, 1.5)
  
  g5 <- ggplot2::ggplot(df_acf, ggplot2::aes(x = Lag, y = ACF)) +
    ggplot2::geom_col(fill = "#34495E", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite, -limite),
                        linetype = "dashed", color = "#2980B9", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "5. Autocorrelacion de Residuos (ACF)",
      subtitle = "Quedo informacion sin capturar en los errores?",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente de Autocorrelacion",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN (RUIDO BLANCO): Todas las barras estan dentro de las lineas azules.\n",
                       " NUMEROS ROJOS: Indican retrasos donde el error sigue un patron.\n",
                       " MAL: Si el Lag 1 es alto, el modelo necesita un orden 'p' mayor.\n",
                       " NOTA: El Lag 0 se omite por ser la correlacion del dato consigo mismo (1.0).")
    ) +
    mi_tema
  
  print(g5)
  if (interactive()) {readline(prompt = "Grafica 5: Evaluacion de la independencia de los errores > ")}
  
  
  # --- GRAFICA 6: AUTOCORRELACION PARCIAL (PACF) DE RESIDUOS ---
  pacf_res <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)
  
  df_pacf <- data.frame(
    Lag = as.numeric(pacf_res$lag),
    PACF = as.numeric(pacf_res$acf)
  )
  
  n_obs <- nrow(df)
  limite <- 1.96 / sqrt(n_obs)
  
  df_pacf$Etiqueta <- ifelse(abs(df_pacf$PACF) > limite, as.character(df_pacf$Lag), "")
  df_pacf$vjust_pos <- ifelse(df_pacf$PACF > 0, -0.5, 1.5)
  
  g6 <- ggplot2::ggplot(df_pacf, ggplot2::aes(x = Lag, y = PACF)) +
    ggplot2::geom_col(fill = "#2980B9", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite, -limite),
                        linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelacion Parcial (PACF)",
      subtitle = "Identificacion de lags residuales directos",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente PACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Ninguna barra sobresale de las lineas naranjas.\n",
                       " NUMEROS ROJOS: Indican que ese lag especifico aun tiene informacion util.\n",
                       " ACCION: Si el Lag 1 o 2 es rojo, aumenta el orden 'p' en tu funcion modelo_ar().\n",
                       " NOTA: La PACF ayuda a confirmar si el proceso es realmente un AR.")
    ) +
    mi_tema
  
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6: Verificacion de la estructura del orden p > ")}
  
  
  # --- GRAFICA 7: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---
  df_limpio <- df[!is.na(df$Residuos), ]
  
  g7 <- ggplot2::ggplot(df_limpio, ggplot2::aes(sample = Residuos)) +
    ggplot2::stat_qq(color = "#8E44AD", alpha = 0.6, size = 2) +
    ggplot2::stat_qq_line(color = "#E74C3C", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "7. Normalidad (Q-Q Plot)",
      subtitle = "Comparacion de Cuantiles de Residuos vs. Normal Teorica",
      x = "Cuantiles Teoricos (Normal)",
      y = "Cuantiles Observados (Residuos)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos morados deben seguir fielmente la linea roja diagonal.\n",
                       " COLAS: Si los extremos se alejan de la linea, hay valores atipicos (outliers).\n",
                       " NOTA: Si los puntos forman una 'S', la distribucion es asimetrica.\n",
                       " CLAVE: La normalidad asegura que tus pronosticos sean estadisticamente confiables.")
    ) +
    mi_tema
  
  print(g7)
  if (interactive()) {readline(prompt = "Grafica 7: Verificacion de la distribucion gaussiana > ")}
  
  
  # --- GRAFICA 8: PRECISION DEL MODELO (REAL VS AJUSTADO) ---
  cor_val <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")
  
  g8 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    ggplot2::geom_point(color = "#2980B9", alpha = 0.5, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    ggplot2::geom_smooth(method = "lm", color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "8. Precision del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlacion de Pearson: ", round(cor_val, 3)),
      x = "Valor Ajustado (Modelo AR)",
      y = "Valor Observado (Realidad)",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA ROJA: Diagonal de perfeccion (Realidad = Modelo).\n",
                       " BIEN: Los puntos deben estar cerca y paralelos a la linea roja.\n",
                       " DISPERSION: Si los puntos se alejan mucho, el modelo tiene poca capacidad predictiva.\n",
                       " NOTA: El modelo AR suele ser conservador y 'suavizar' los extremos.")
    ) +
    mi_tema
  
  print(g8)
  if (interactive()) {readline(prompt = "Grafica 8: Comparacion directa de valores reales y estimados > ")}
  
  
  # --- GRAFICA 9: DIAGNOSTICO DE ESTACIONARIEDAD (MOVIL) ---
  ventana_ar <- max(5, round(nrow(df) * 0.1))
  
  df_est_ar <- df
  df_est_ar$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_ar, ventana_ar), sides = 2))
  
  # CORRECCION 2: Usar stats::sd
  df_est_ar$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_ar):min(nrow(df), i+ventana_ar)
    stats::sd(df$Observado[idx], na.rm = TRUE)
  })
  
  df_final_ar <- df_est_ar[!is.na(df_est_ar$Media_Movil), ]
  
  g9 <- ggplot2::ggplot(df_final_ar, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray80", linewidth = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#27AE60", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#27AE60", linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
    ggplot2::labs(
      title = "9. Diagnostico de Estacionariedad",
      subtitle = paste0("Analisis de Media y Varianza Movil (Ventana k = ", ventana_ar, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA VERDE: Media local movil. Debe ser horizontal y estable.\n",
                       " LINEA ROJA: Media global (Referencia constante).\n",
                       " BIEN: Si la franja verde es estable y no tiene tendencia.\n",
                       " NOTA: Se han filtrado los extremos para un reporte libre de avisos de error.")
    ) +
    mi_tema
  
  print(g9)
  
  # CORRECCION 1: Usar message()
  if (interactive()) {
    message("Secuencia completa de 9 graficas para Series de Tiempo (AR) finalizada.\n")
  }
  
  # CORRECCION 3: Retornar x
  return(invisible(x))
}