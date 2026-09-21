#' Modelo de Media Movil (MA)
#'
#' Esta funcion ajusta un modelo de series de tiempo donde el valor actual \eqn{(Y_t)}
#' se explica mediante una combinacion lineal de los errores (shocks) pasados
#' (\eqn{\epsilon_{t-1}}, \eqn{\epsilon_{t-2}}, ...). Es ideal para modelar fenomenos con
#' efectos de memoria corta o impactos transitorios.
#'
#' @details
#' El flujo interno de la funcion se divide en cinco etapas fundamentales:
#' \enumerate{
#'   \item **Estimacion (MLE):** Utiliza Maxima Verosimilitud para estimar los
#'         parametros \eqn{\theta} (Theta). A diferencia del AR, el MA requiere
#'         algoritmos iterativos ya que los errores no son directamente observables.
#'   \item **Analisis de Invertibilidad:** Calcula las raices inversas del polinomio
#'         caracteristico. Un modelo MA debe ser "Invertible" para que el error
#'         pueda expresarse como una funcion convergente de los datos observados.
#'   \item **Identificacion Estructural:** Evalua la ACF (Autocorrelacion). En un
#'         proceso MA(q), la ACF debe "cortarse" abruptamente despues del lag q.
#'   \item **Auditoria de Residuos:** Verifica que el error remanente sea Ruido Blanco
#'         mediante la prueba de Ljung-Box y la distribucion Gaussiana (Shapiro-Wilk).
#'   \item **Traduccion Didactica:** Genera un reporte que explica la persistencia
#'         de los choques y valida si el modelo es unico y matematicamente solido.
#' }
#'
#' @section Guia de las 9 Graficas de Diagnostico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una auditoria de "shocks" temporales:
#' \describe{
#'   \item{\strong{1. Inspeccion Visual}}{Muestra la serie bruta con su media historica. El modelo MA supone que la serie siempre es atraida hacia este nivel central.}
#'   \item{\strong{2. Ajuste Temporal}}{Superpone el modelo (azul) sobre la realidad. Permite ver como el modelo "suaviza" los picos de la serie.}
#'   \item{\strong{3. Precision (Real vs Ajustado)}}{Grafica de dispersion para evaluar la correlacion. Los puntos deben seguir la diagonal roja.}
#'   \item{\strong{4. Circulo Unitario}}{Prueba de Invertibilidad. Los rombos azules deben estar dentro del circulo para asegurar que el modelo sea unico.}
#'   \item{\strong{5. Residuos en el Tiempo}}{Evaluacion de aleatoriedad. No debe haber patrones; los errores deben ser impredecibles.}
#'   \item{\strong{6. ACF de Residuos}}{Prueba critica de Ruido Blanco. Si hay barras largas (numeros rojos), el orden \code{q} es insuficiente.}
#'   \item{\strong{7. PACF de Residuos}}{Complemento de la ACF para verificar que no existan dependencias directas remanentes.}
#'   \item{\strong{8. Normalidad (Q-Q Plot)}}{Valida que los choques sigan una campana de Gauss, asegurando la validez de los intervalos de confianza.}
#'   \item{\strong{9. Estacionariedad Movil}}{Confirma que la media y varianza no se desplazan. Vital para asegurar que el modelo MA sea aplicable.}
#' }
#'
#' @param x Para \code{modelo_ma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"ma_greenreg"} devuelto por \code{modelo_ma()}.
#' @param q Orden del modelo (numero de choques pasados a considerar). Por defecto es 1.
#' @param include_mean Logico. Incluir el promedio constante de la serie? Por defecto \code{TRUE}.
#'
#' @return Un objeto de clase \code{"ma_greenreg"}, que consiste en una lista con:
#' \itemize{
#'   \item \code{modelo}: El objeto \code{arima} original (orden 0,0,q).
#'   \item \code{orden_q}: El orden del modelo ajustado.
#'   \item \code{coeficientes_tabla}: Matriz con estimaciones, errores estandar y p-valores.
#'   \item \code{invertibilidad}: Lista con las raices inversas y diagnostico de invertibilidad.
#'   \item \code{diagnosticos}: Resultados de las pruebas de Ljung-Box, Shapiro-Wilk y AIC.
#'   \item \code{data_plot}: Dataframe listo para visualizacion con \code{ggplot2}.
#' }
#' @examples
#' data("datos_anomalia_temperatura")
#' modelo_MA <- modelo_ma(datos_anomalia_temperatura$anomalia_c, q = 1, include_mean = TRUE)
#' modelo_MA
#' plot(modelo_MA)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf logLik filter sd
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
modelo_ma <- function(x, q = 1, include_mean = TRUE) {

  # --- 1. Validacion y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numerico.")

  # Ajuste usando ARIMA(0,0,q) que equivale a MA(q)
  modelo_base <- stats::arima(x, order = c(0, 0, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadisticas ---
  coefs <- modelo_base$coef
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)

  # --- 3. Analisis de Invertibilidad (Raices Unitarias) ---
  # La invertibilidad en MA es equivalente a la estabilidad en AR
  # Extraemos solo los coeficientes MA (ignoramos intercepto)
  ma_coefs <- coefs[grepl("^ma", names(coefs))]

  if (length(ma_coefs) > 0) {
    # Polinomio caracteristico: 1 + theta1*z + theta2*z^2 ...
    poly_roots <- polyroot(c(1, ma_coefs))
    inverse_roots <- 1 / poly_roots
    modulos <- Mod(inverse_roots)
    # Regla: Invertible si TODOS los modulos son < 1
    es_invertible <- all(modulos < 1)
  } else {
    inverse_roots <- NULL
    es_invertible <- TRUE
  }

  # --- 4. Diagnosticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # Ljung-Box (Ruido Blanco)
  lb_test <- stats::Box.test(residuos, type = "Ljung-Box", lag = q + 5)

  # Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- stats::shapiro.test(residuos)
  } else {
    shapiro_test <- list(p.value = NA)
  }

  # --- 5. Datos para Graficas ---
  ajustados <- x - residuos
  sigma <- sqrt(modelo_base$sigma2)
  # Intervalos de confianza 95%
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

  # Datos Circulo Unitario
  df_raices <- NULL
  if (!is.null(inverse_roots)) {
    df_raices <- data.frame(Real = Re(inverse_roots), Imaginario = Im(inverse_roots))
  }

  # --- 6. Objeto Final ---
  resultado <- list(
    x_original = x,
    modelo = modelo_base,
    orden_q = q,
    coeficientes_tabla = tabla_coefs,
    invertibilidad = list(es_invertible = es_invertible, raices = inverse_roots),
    diagnosticos = list(lb = lb_test, shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "ma_greenreg"
  return(resultado)
}

#' @rdname modelo_ma
#' @param x Para \code{modelo_ma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"ma_greenreg"} devuelto por \code{modelo_ma()}.
#' @param ... Argumentos adicionales.
#' @export
print.ma_greenreg <- function(x, ...) {
  
  cat("\n==========================================================\n")
  cat("      MODELO DE MEDIA MOVIL MA(", x$orden_q, ") \n", sep = "")
  cat("==========================================================\n\n")
  
  # --- 1. MODELO Y ECUACION ---
  cat("--- 1. MODELO Y ECUACION ---\n")
  theta <- round(x$modelo$coef[grepl("^ma", names(x$modelo$coef))], 3)
  c_int <- round(x$modelo$coef["intercept"], 3)
  
  # Proteccion contra interceptos NA
  mu <- if(is.na(c_int)) 0 else c_int
  
  # Construccion visual de la ecuacion de choques
  eq <- paste0("Y(t) = ", mu, " + E(t)")
  if(length(theta) > 0) {
    for(i in 1:length(theta)) {
      eq <- paste0(eq, " + (", theta[i], " * E(t-", i, "))")
    }
  }
  
  cat("Ecuacion de Medias Moviles Estimada:\n", eq, "\n")
  
  # Proteccion NA para AIC
  aic_val <- x$diagnosticos$aic
  aic_txt <- if(is.null(aic_val) || is.na(aic_val)) "NA" else round(aic_val, 2)
  cat("AIC (Criterio de Akaike):", aic_txt, "\n")
  cat(" NOTA: El modelo MA(", x$orden_q, ") describe como los choques (E) del pasado afectan al presente.\n\n", sep = "")
  
  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS (Parametros Theta) ---\n")
  # Obligatorio usar as.matrix para evitar desajustes
  stats::printCoefmat(as.matrix(x$coeficientes_tabla), digits = 4, signif.stars = TRUE)
  cat("\n GUIA: Los coeficientes Theta miden la persistencia de los choques externos.\n\n")
  
  # --- 3. IDENTIFICACION DEL MODELO (ORDEN q) ---
  cat("--- 3. IDENTIFICACION DEL MODELO (Orden q) ---\n")
  cat(" Orden seleccionado: q =", x$orden_q, "\n")
  cat(" Criterio ACF: La autocorrelacion debe cortarse abruptamente despues del lag", x$orden_q, ".\n")
  cat(" Criterio PACF: Los coeficientes parciales deben decaer exponencialmente hacia cero.\n")
  cat(" NOTA: El modelo MA es el 'espejo' del AR. Aqui la ACF es la que manda para elegir q.\n\n")
  
  # --- 4. VERIFICACION DE SUPUESTOS ---
  cat("--- 4. VERIFICACION DE SUPUESTOS ---\n")
  
  # A) Estacionariedad e Invertibilidad
  cat("A) Estacionariedad e Invertibilidad:\n")
  cat("  Estacionariedad: Un modelo MA puro es SIEMPRE estacionario por definicion.\n")
  
  if (!is.null(x$invertibilidad$raices)) {
    modulos <- Mod(x$invertibilidad$raices)
    cat("  Modulos de raices inversas (Invertibilidad):", paste(round(modulos, 3), collapse = ", "), "\n")
    
    es_inv <- x$invertibilidad$es_invertible
    if(!is.null(es_inv) && !is.na(es_inv) && es_inv) {
      cat("     -> OKEY: El modelo es Invertible (Raices < 1). Los errores son estimables.\n")
    } else {
      cat("     ->  ALERTA: Modelo NO Invertible. Los choques pasados tienen un peso infinito.\n")
    }
  }
  
  # B) Errores ~ Ruido Blanco
  cat("\nB) Diagnostico de Errores (E ~ RB):\n")
  
  # Media Cero (Protegido)
  mean_res <- mean(x$data_plot$Residuos, na.rm = TRUE)
  mean_txt <- if(is.na(mean_res)) "NA" else sprintf("%.4f", mean_res)
  cat(sprintf("  Media de los residuos: %s (Esperado: 0)\n", mean_txt))
  
  # No Autocorrelacion Ljung-Box (Protegido)
  lb_p <- x$diagnosticos$lb$p.value
  lb_txt <- if(is.null(lb_p) || is.na(lb_p)) "NA" else sprintf("%.4f", lb_p)
  cat(sprintf("  No Autocorrelacion (Ljung-Box p-valor): %s\n", lb_txt))
  
  if (!is.null(lb_p) && !is.na(lb_p) && lb_p > 0.05) {
    cat("     -> OKEY: Los residuos son ruido blanco. Toda la senal fue capturada.\n")
  } else {
    cat("     ->  ALERTA: Queda estructura en los residuos. Considera aumentar el orden q.\n")
  }
  
  # Varianza y Normalidad (Protegido)
  sh_p <- x$diagnosticos$shapiro$p.value
  sh_txt <- if(is.null(sh_p) || is.na(sh_p)) "NA" else sprintf("%.4f", sh_p)
  cat(sprintf("  Normalidad (Shapiro p-valor): %s\n", sh_txt))
  cat("     ->", if(!is.null(sh_p) && !is.na(sh_p) && sh_p > 0.05) "OKEY: Errores normales.\n" else " Errores no normales o indeterminados.\n")
  
  # --- 5. RESUMEN DEL MODELO ---
  cat("\n--- 5. RESUMEN DEL MODELO ---\n")
  
  # Validacion logica segura para evitar colapso por TRUE/FALSE faltante
  es_valido <- FALSE
  es_inv <- x$invertibilidad$es_invertible
  
  if (!is.null(es_inv) && !is.na(es_inv) && !is.null(lb_p) && !is.na(lb_p)) {
    if (es_inv && lb_p > 0.05) {
      es_valido <- TRUE
    }
  }
  
  if(es_valido) {
    cat("RESULTADO FINAL: El modelo MA(", x$orden_q, ") es VALIDO y UNICO.\n", sep="")
  } else {
    cat(" RECOMENDACION: Revisa el orden q o la invertibilidad para asegurar un ajuste correcto.\n")
  }
  
  cat("\n Usa plot(modelo) para visualizar la serie y la ACF de los residuos.\n")
  cat("==========================================================\n")
  
  return(invisible(x))
}

#' @rdname modelo_ma
#' @description Genera una secuencia de 9 graficas de diagnostico para validar el modelo MA.
#' @param x Para \code{modelo_ma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"ma_greenreg"} devuelto por \code{modelo_ma()}.
#' @param ... Argumentos adicionales.
#' @export
plot.ma_greenreg <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instala ggplot2.")
  
  df <- x$data_plot
  
  # TEMA ESTANDAR
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
  
  # CORRECCION 1: message en lugar de cat
  if (interactive()) {
    message("Generando graficas de diagnostico MA... (Presiona [Enter] para avanzar)")
  }
  
  # --- GRAFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---
  
  # Calculamos la media para la linea de referencia
  media_serie <- mean(df$Observado, na.rm = TRUE)
  
  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    # Linea de la serie original
    ggplot2::geom_line(color = "#34495E", linewidth = 0.7) +
    # Puntos sutiles
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    # Linea de la media historica (Referencia fundamental en MA)
    ggplot2::geom_hline(yintercept = media_serie, linetype = "dashed",
                        color = "#2980B9", linewidth = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualizacion de la serie temporal bruta",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA AZUL: Es el valor promedio historico.\n",
                       " MEDIA MOVIL: El modelo MA supone que la serie siempre regresa a esta media.\n",
                       " NOTA: Si observas una tendencia marcada (subida o bajada constante),\n",
                       "   un modelo MA puro podria no ser suficiente sin diferenciar primero.")
    ) +
    mi_tema
  
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1: Inspeccion de la serie pura > ")}
  
  # --- GRAFICA 2: AJUSTE DEL MODELO (Real vs Ajustado) ---
  
  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    # Intervalo de confianza (Sombra azul)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#2980B9", alpha = 0.15) +
    # Serie Observada (Realidad)
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), linewidth = 0.7) +
    # Serie Ajustada (Modelo MA)
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo MA"), linewidth = 0.8) +
    # Configuracion de colores
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo MA" = "#2980B9")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo MA(", x$orden_q, ")"),
      subtitle = "Serie Observada vs. Valores Ajustados (Fitted Values)",
      x = "Tiempo (Periodos)",
      y = "Valor",
      color = "Leyenda:",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: La linea azul debe capturar la tendencia central de la oscura.\n",
                       " SOMBRA AZUL: Margen de error esperado (IC 95%).\n",
                       " NOTA: El modelo MA suele ser mas 'suave' que la serie original,\n",
                       "   ya que filtra el ruido aleatorio para quedarse con el promedio movil.")
    ) +
    mi_tema
  
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2: Comparacion del modelo contra los datos reales > ")}
  
  # --- GRAFICA 3: PRECISION DEL MODELO (REAL VS AJUSTADO) ---
  
  # Calculamos la correlacion para el reporte visual
  cor_ma <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")
  
  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    # Puntos de dispersion (Color azul MA)
    ggplot2::geom_point(color = "#2980B9", alpha = 0.5, size = 2) +
    # Linea de Identidad (El ideal de perfeccion)
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    # Linea de tendencia real con formula explicita para evitar avisos
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "3. Precision del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlacion de Pearson: ", round(cor_ma, 3)),
      x = "Valor Ajustado (Modelo MA)",
      y = "Valor Observado (Realidad)",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA ROJA: Diagonal de 45 (Donde el modelo es perfecto).\n",
                       " BIEN: Los puntos deben agruparse cerca de la linea roja.\n",
                       " NOTA: Si la nube de puntos es muy ancha, el modelo MA tiene poca precision.\n",
                       " CLAVE: Un valor de correlacion cercano a 1 indica un ajuste excelente.")
    ) +
    mi_tema
  
  print(g3)
  if (interactive()) {readline(prompt = "Grafica 3: Analisis de correlacion y precision. > ")}
  
  # --- GRAFICA 4: INVERTIBILIDAD (CIRCULO UNITARIO) ---
  
  if (!is.null(x$data_raices)) {
    # 1. Crear el circulo de referencia
    theta <- seq(0, 2*pi, length.out = 200)
    df_circle <- data.frame(x = cos(theta), y = sin(theta))
    
    # 2. Construccion de la Grafica con Estilo MA (Azul)
    g4 <- ggplot2::ggplot() +
      # Circulo unitario punteado
      ggplot2::geom_path(data = df_circle, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", linewidth = 0.8) +
      # Ejes cartesianos
      ggplot2::geom_vline(xintercept = 0, color = "#BDC3C7") +
      ggplot2::geom_hline(yintercept = 0, color = "#BDC3C7") +
      # area de Invertibilidad (Verde tenue)
      ggplot2::geom_polygon(data = df_circle, ggplot2::aes(x, y),
                            fill = "#27AE60", alpha = 0.05) +
      # Raices inversas (Puntos azules para MA)
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario),
                          color = "#2980B9", size = 4, shape = 18) +
      # Ajuste de proporciones
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::labs(
        title = "4. Invertibilidad (Circulo Unitario)",
        subtitle = paste0("Raices Inversas del Polinomio MA(", x$orden_q, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        caption = paste0("INTERPRETACION:\n",
                         " INVERTIBLE: Todas las raices (rombos azules) estan DENTRO del circulo.\n",
                         " ALERTA: Si una raiz toca o sale del circulo, el modelo no es unico.\n",
                         " NOTA: La invertibilidad asegura que los errores del modelo sean estimables\n",
                         "   a partir de los datos observados de forma convergente.")
      ) +
      mi_tema
    
    print(g4)
    
  } else {
    # CORRECCION 1: Usar message()
    message("El modelo MA(0) no tiene raices; es ruido blanco y siempre es invertible.")
  }
  if (interactive()) {readline(prompt = "Grafica 4: Analisis de raices del proceso de media movil > ")}
  
  # --- GRAFICA 5: RESIDUOS EN EL TIEMPO (E vs t) ---
  
  # CORRECCION 2: stats::sd()
  sd_res_ma <- stats::sd(df$Residuos, na.rm = TRUE)
  
  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    # Banda de confianza (+-2 desviaciones estandar) en color Morado/Azul
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_ma, ymax = 2*sd_res_ma),
                         fill = "#8E44AD", alpha = 0.1) +
    # Linea de los residuos
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.5) +
    # Puntos para identificar errores individuales
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.4, size = 1) +
    # Linea central (Media Cero)
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Umbrales criticos
    ggplot2::geom_hline(yintercept = c(-2*sd_res_ma, 2*sd_res_ma),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.6) +
    ggplot2::labs(
      title = "5. Residuos en el Tiempo",
      subtitle = "Analisis de Aleatoriedad (E ~ Ruido Blanco)",
      x = "Tiempo (Periodos)",
      y = "Error (Residuo)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos deben oscilar al azar alrededor de la linea negra.\n",
                       " BANDAS ROJAS: El 95% de los errores deben estar dentro de estos limites.\n",
                       " MAL: Si ves ciclos o tendencias, el orden 'q' es insuficiente.\n",
                       " CLAVE: Si la variabilidad aumenta con el tiempo, el modelo es inestable.")
    ) +
    mi_tema
  
  print(g5)
  if (interactive()) {readline(prompt = "Grafica 5: Verificacion de la independencia temporal de los errores. > ")}
  
  # --- GRAFICA 6: AUTOCORRELACION DE RESIDUOS (ACF) ---
  
  # 1. Calculo de la ACF de los residuos (maximo 20 lags)
  acf_res_ma <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
  
  # Creamos el dataframe omitiendo el Lag 0 (que siempre es 1)
  df_acf_ma <- data.frame(
    Lag = as.numeric(acf_res_ma$lag[-1]),
    ACF = as.numeric(acf_res_ma$acf[-1])
  )
  
  # Limite de significancia estadistica (95%)
  n_obs_ma <- nrow(df)
  limite_ma <- 1.96 / sqrt(n_obs_ma)
  
  # Logica de etiquetas rojas para alertar sobre fallos en lags especificos
  df_acf_ma$Etiqueta <- ifelse(abs(df_acf_ma$ACF) > limite_ma, as.character(df_acf_ma$Lag), "")
  df_acf_ma$vjust_pos <- ifelse(df_acf_ma$ACF > 0, -0.5, 1.5)
  
  # 2. Construccion de la Grafica
  g6 <- ggplot2::ggplot(df_acf_ma, ggplot2::aes(x = Lag, y = ACF)) +
    # Barras de correlacion (Color oscuro para contraste)
    ggplot2::geom_col(fill = "#34495E", width = 0.4) +
    # Lineas de significancia (Azules para MA)
    ggplot2::geom_hline(yintercept = c(limite_ma, -limite_ma),
                        linetype = "dashed", color = "#2980B9", linewidth = 0.8) +
    # Linea base en cero
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Etiquetas de error (Numeros Rojos)
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    # Escala del eje X
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelacion de Residuos (ACF)",
      subtitle = "Quedo memoria de choques sin capturar?",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente ACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Si todas las barras estan dentro de las lineas azules.\n",
                       " NUMEROS ROJOS: Indican que el error en ese lag NO es aleatorio.\n",
                       " MAL: Si el Lag 1 es rojo, tu modelo MA(", x$orden_q, ") es insuficiente.\n",
                       " CLAVE: Los residuos deben ser 'Ruido Blanco' (puro azar).")
    ) +
    mi_tema
  
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6: Diagnostico de independencia de los errores. > ")}
  
  # --- GRAFICA 7: AUTOCORRELACION PARCIAL (PACF) DE RESIDUOS ---
  
  # 1. Calculo de la PACF de los residuos
  pacf_res_ma <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)
  
  df_pacf_ma <- data.frame(
    Lag = as.numeric(pacf_res_ma$lag),
    PACF = as.numeric(pacf_res_ma$acf)
  )
  
  # Limite de significancia (Mismo que en ACF)
  n_obs_ma <- nrow(df)
  limite_ma <- 1.96 / sqrt(n_obs_ma)
  
  # Logica de etiquetas rojas para fallas
  df_pacf_ma$Etiqueta <- ifelse(abs(df_pacf_ma$PACF) > limite_ma, as.character(df_pacf_ma$Lag), "")
  df_pacf_ma$vjust_pos <- ifelse(df_pacf_ma$PACF > 0, -0.5, 1.5)
  
  # 2. Construccion de la Grafica con Estilo MA (Azul/Naranja)
  g7 <- ggplot2::ggplot(df_pacf_ma, ggplot2::aes(x = Lag, y = PACF)) +
    # Barras de la PACF (Color Azul MA)
    ggplot2::geom_col(fill = "#2980B9", width = 0.4) +
    # Lineas de significancia (Cambiamos a Naranja para diferenciar de ACF)
    ggplot2::geom_hline(yintercept = c(limite_ma, -limite_ma),
                        linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
    # Linea base
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Etiquetas de error (Chivatos en Rojo)
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    # Escala del eje X
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "7. Autocorrelacion Parcial (PACF)",
      subtitle = "Identificacion de dependencias directas en residuos",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente PACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Ninguna barra sobresale de las lineas naranjas.\n",
                       " NUMEROS ROJOS: Indican que ese retraso especifico aun tiene informacion.\n",
                       " NOTA: En un modelo MA bien ajustado, tanto ACF como PACF de residuos \n",
                       "   deberian verse como Ruido Blanco (barras cortas).")
    ) +
    mi_tema
  
  print(g7)
  if (interactive()) {readline(prompt = "Grafica 7: Verificacion de la limpieza estructural del error > ")}
  
  # --- GRAFICA 8: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---
  
  g8 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
    # Puntos de los cuantiles observados (Color Naranja para MA)
    ggplot2::stat_qq(color = "#E67E22", alpha = 0.6, size = 2) +
    # Linea de referencia teorica (Normalidad perfecta)
    ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "8. Normalidad de los Errores (Q-Q Plot)",
      subtitle = "Comparacion de Cuantiles de Residuos vs. Normal Teorica",
      x = "Cuantiles Teoricos (Distribucion Normal)",
      y = "Cuantiles Observados (Residuos MA)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos deben seguir la linea diagonal punteada.\n",
                       " EXTREMOS: Si los puntos se alejan en las puntas, hay valores atipicos.\n",
                       " NOTA: La normalidad valida que el error tipico (Sigma) sea constante.\n",
                       " CLAVE: Si los puntos forman una curva, la serie tiene asimetria.")
    ) +
    mi_tema
  
  print(g8)
  if (interactive()) {readline(prompt = "Grafica 8: Evaluacion de la distribucion gaussiana del error > ")}
  
  # --- GRAFICA 9: DIAGNOSTICO DE ESTACIONARIEDAD (MOVIL) ---
  
  # 1. Definimos la ventana (10% de los datos o minimo 5)
  ventana_ma <- max(5, round(nrow(df) * 0.1))
  
  df_est_ma <- df
  
  # 2. CALCULOS (Es vital forzar as.numeric para ggplot2)
  df_est_ma$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_ma, ventana_ma), sides = 2))
  
  # CORRECCION 2: stats::sd()
  df_est_ma$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_ma):min(nrow(df), i+ventana_ma)
    stats::sd(df$Observado[idx], na.rm = TRUE)
  })
  
  # 3. FILTRADO (Eliminamos los NAs antes de entrar a ggplot)
  df_final_ma <- df_est_ma[!is.na(df_est_ma$Media_Movil), ]
  
  # 4. CONSTRUCCION DE LA GRAFICA
  g9 <- ggplot2::ggplot(df_final_ma, ggplot2::aes(x = Tiempo)) +
    # Serie original completa de fondo (usando df original para no perder datos visuales)
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray80", linewidth = 0.5) +
    
    # Franja de Varianza Movil (Azul MA)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#2980B9", alpha = 0.15) +
    
    # Linea de la Media Movil (Nivel local)
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#2980B9", linewidth = 1.2) +
    
    # Linea de la Media Global (Referencia fija en Rojo)
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
    
    ggplot2::labs(
      title = "9. Diagnostico de Estacionariedad",
      subtitle = paste0("Analisis de Media y Varianza Movil (Ventana k = ", ventana_ma, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA AZUL: Media local movil. Debe ser horizontal y estable.\n",
                       " LINEA ROJA: Media global (El objetivo del modelo MA).\n",
                       " BIEN: Si la franja azul mantiene su ancho y nivel constante.\n",
                       " NOTA: Se han filtrado los extremos para evitar avisos de valores nulos.")
    ) +
    mi_tema
  
  print(g9)
  
  # CORRECCION 1: message en lugar de cat
  if (interactive()) {
    message("Grafica 9: Verificacion final de estabilidad estadistica completada.")
  }
  
  # CORRECCION 4: Retorno del objeto x original
  return(invisible(x))
}