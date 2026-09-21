#' Modelo Autorregresivo de Media Movil (ARMA)
#'
#' Esta funcion ajusta un modelo hibrido ARMA(p, q) que combina la dependencia
#' de los valores pasados (AutoRegresivo) con la influencia de los errores
#' pasados (Media Movil). Es la herramienta estandar para modelar series de
#' tiempo estacionarias con estructuras de dependencia complejas.
#'
#' @details
#' El flujo interno de la funcion realiza una auditoria estructural en cinco pasos:
#' \enumerate{
#'   \item **Estimacion Conjunta:** Utiliza Maxima Verosimilitud (MLE) para encontrar
#'         simultaneamente los coeficientes \eqn{\phi} (AR) y \eqn{\theta} (MA).
#'   \item **Estabilidad (Parte AR):** Verifica que el componente autorregresivo no sea
#'         explosivo. Las raices inversas del polinomio AR deben estar dentro del circulo.
#'   \item **Invertibilidad (Parte MA):** Verifica que el componente de media movil sea
#'         univoco. Las raices inversas del polinomio MA deben estar dentro del circulo.
#'   \item **Diagnostico de Residuos:** Evalua si el error resultante es Ruido Blanco
#'         mediante el test de Ljung-Box, asegurando que no quede senal sin capturar.
#'   \item **Criterio de Informacion:** Proporciona el AIC para facilitar la comparacion
#'         entre diferentes combinaciones de (p, q), buscando el equilibrio entre
#'         complejidad y precision.
#' }
#'
#' @section Guia de las 9 Graficas de Diagnostico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validacion mixta:
#' \describe{
#'   \item{\strong{1. Inspeccion de Serie}}{Muestra la serie bruta con su media. Fundamental para confirmar visualmente la estacionariedad.}
#'   \item{\strong{2. Ajuste Mixto}}{Superpone el modelo (verde) sobre la realidad. Permite evaluar como interactuan la inercia y los shocks.}
#'   \item{\strong{3. Residuos en el Tiempo}}{Evaluacion de la aleatoriedad del error. Los puntos deben formar una nube sin patrones.}
#'   \item{\strong{4. Circulo Unitario Dual}}{Grafica critica exclusiva de ARMA. Muestra simultaneamente la Estabilidad (Puntos Rojos) e Invertibilidad (Rombos Azules).}
#'   \item{\strong{5. ACF de Residuos}}{Busca autocorrelaciones remanentes. Las barras deben ser minimas para confirmar que el modelo es optimo.}
#'   \item{\strong{6. PACF de Residuos}}{Ayuda a identificar si el orden autorregresivo (p) es suficiente.}
#'   \item{\strong{7. Normalidad (Q-Q Plot)}}{Valida que los errores sigan una campana de Gauss, requisito para la validez de los intervalos de confianza.}
#'   \item{\strong{8. Precision (Real vs Ajustado)}}{Muestra la correlacion directa. Un ARMA exitoso suele tener una alineacion estrecha sobre la diagonal.}
#'   \item{\strong{9. Estacionariedad Movil}}{Prueba final de estabilidad. Una franja verde horizontal confirma que el modelo no se desplaza en el tiempo.}
#' }
#'
#' @param x Para \code{modelo_arma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{modelo_arma()}.
#' @param p Orden Autorregresivo (pasados de la variable).
#' @param q Orden de Media Movil (pasados del error).
#' @param include_mean Logico. Incluir constante/intercepto? Por defecto \code{TRUE}.
#'
#' @return Un objeto de clase \code{"arma_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{arima} original (orden p, 0, q).
#'   \item \code{ordenes}: Vector con los valores de p y q.
#'   \item \code{coeficientes_tabla}: Matriz con estimaciones, errores estandar y p-valores.
#'   \item \code{estatus}: Lista con el estado de estabilidad e invertibilidad.
#'   \item \code{diagnosticos}: Resultados de Ljung-Box, Shapiro-Wilk y AIC.
#'   \item \code{data_plot}: Dataframe listo para visualizacion con \code{ggplot2}.
#' }
#'
#' @examples
#' data("datos_anomalia_temperatura")
#' ts_data <- datos_anomalia_temperatura$anomalia_c
#' modelo <- modelo_arma(ts_data, p = 1, q = 1)
#' modelo
#' plot(modelo)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf filter sd
#' @import ggplot2
#' @export
modelo_arma <- function(x, p = 1, q = 1, include_mean = TRUE) {

  # --- 1. Validacion y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numerico.")

  # Ajuste ARIMA(p, 0, q)
  modelo_base <- stats::arima(x, order = c(p, 0, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadisticas ---
  coefs <- modelo_base$coef
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)

  # --- 3. Analisis de Raices (Doble Verificacion) ---

  # A) Estabilidad (Parte AR) - Raices en ROJO
  ar_coefs <- coefs[grepl("^ar", names(coefs))]
  raices_ar_inv <- NULL
  es_estable <- TRUE

  if (length(ar_coefs) > 0) {
    poly_roots_ar <- polyroot(c(1, -ar_coefs))
    raices_ar_inv <- 1 / poly_roots_ar
    es_estable <- all(Mod(raices_ar_inv) < 1)
  }

  # B) Invertibilidad (Parte MA) - Raices en AZUL
  ma_coefs <- coefs[grepl("^ma", names(coefs))]
  raices_ma_inv <- NULL
  es_invertible <- TRUE

  if (length(ma_coefs) > 0) {
    poly_roots_ma <- polyroot(c(1, ma_coefs))
    raices_ma_inv <- 1 / poly_roots_ma
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
  ajustados <- x - residuos
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
    ordenes = c(p = p, q = q),
    coeficientes_tabla = tabla_coefs,
    estatus = list(estable = es_estable, invertible = es_invertible),
    diagnosticos = list(lb = lb_test, shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "arma_greenreg"
  return(resultado)
}

#' @rdname modelo_arma
#' @param x Para \code{modelo_arma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{modelo_arma()}.
#' @param ... Argumentos adicionales.
#' @export
print.arma_greenreg <- function(x, ...) {
  
  # Uso de sep="" para que los numeros p y q queden pegados a las comas
  cat("\n==========================================================\n")
  cat(" Modelo Mixto ARMA(", x$ordenes['p'], ",", x$ordenes['q'], ") \n", sep="")
  cat("==========================================================\n\n")
  
  cat("--- 1. COEFICIENTES ESTIMADOS ---\n")
  stats::printCoefmat(as.matrix(x$coeficientes_tabla), digits = 4, signif.stars = TRUE)
  
  # Proteccion contra NA en AIC
  aic_val <- x$diagnosticos$aic
  aic_txt <- if(is.null(aic_val) || is.na(aic_val)) "NA" else round(aic_val, 2)
  cat("\nAIC:", aic_txt, "(Menor es mejor)\n\n")
  
  cat("--- 2. VERIFICACION ESTRUCTURAL (Raices Unitarias) ---\n")
  
  # Check AR (Protegido contra NA)
  if (x$ordenes['p'] > 0) {
    est <- x$estatus$estable
    if (!is.null(est) && !is.na(est) && est) {
      cat(" Parte AR: ESTABLE (Raices dentro del circulo).\n")
    } else { 
      cat(" Parte AR: INESTABLE (Explosivo o Indeterminado).\n")
    }
  } else { 
    cat(" Parte AR: No aplica (p=0).\n") 
  }
  
  # Check MA (Protegido contra NA)
  if (x$ordenes['q'] > 0) {
    inv <- x$estatus$invertible
    if (!is.null(inv) && !is.na(inv) && inv) {
      cat(" Parte MA: INVERTIBLE (Modelo unico).\n")
    } else { 
      cat(" Parte MA: NO INVERTIBLE (Problemas de identificacion).\n")
    }
  } else { 
    cat(" Parte MA: No aplica (q=0).\n") 
  }
  cat("\n")
  
  cat("--- 3. DIAGNOSTICO DE RESIDUOS ---\n")
  
  # Proteccion contra NA en Ljung-Box y sprintf
  lb_p <- x$diagnosticos$lb$p.value
  lb_txt <- if(is.null(lb_p) || is.na(lb_p)) "NA" else sprintf("%.4f", lb_p)
  
  cat(sprintf(" Independencia (Ljung-Box): p-valor = %s\n", lb_txt))
  
  # Evaluacion logica blindada
  if (!is.null(lb_p) && !is.na(lb_p) && lb_p > 0.05) {
    cat(" Ruido Blanco: El modelo capturo la dinamica temporal.\n")
  } else {
    cat(" Autocorrelacion: El modelo ARMA actual es insuficiente o indeterminado.\n")
  }
  
  cat("\n Usa plot(modelo) para graficas de diagnostico del modelo ARMA.\n")
  cat("==========================================================\n")
  
  # Retorno obligatorio S3
  return(invisible(x))
}

#' @rdname modelo_arma
#' @description Genera una secuencia de 9 graficas de diagnostico para validar el modelo ARMA.
#' @param x Para \code{modelo_arma()}: vector numerico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{modelo_arma()}.
#' @param ... Argumentos adicionales.
#' @export
plot.arma_greenreg <- function(x, ...) {
  
  df <- x$data_plot
  p  <- x$ordenes['p'] # Extraemos p del objeto x
  q  <- x$ordenes['q'] # Extraemos q del objeto x
  
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
  
  # CORRECCION: Usar message
  if (interactive()) {
    message("Generando graficas ARMA... (Presiona [Enter] para avanzar)")
  }
  
  # --- GRAFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---
  
  media_global <- mean(df$Observado, na.rm = TRUE)
  
  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    ggplot2::geom_line(color = "#34495E", linewidth = 0.7) +
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    ggplot2::geom_hline(yintercept = media_global, linetype = "dashed",
                        color = "#16A085", linewidth = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualizacion de la serie temporal antes del ajuste mixto",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA VERDE: Promedio historico de la serie.\n",
                       " OBJETIVO: El modelo ARMA asume que la serie fluctua alrededor\n",
                       "   de este nivel constante (Estacionariedad en media).\n",
                       " NOTA: Si ves una 'escalera' o tendencia, el modelo ARMA(p,q)\n",
                       "   podria requerir una diferencia previa (Modelo ARIMA).")
    ) +
    mi_tema
  
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1: Inspeccion de los datos brutos > ")}
  
  
  # --- GRAFICA 2: AJUSTE DEL MODELO ARMA (Real vs Ajustado) ---
  
  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#16A085", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo ARMA"), linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo ARMA" = "#16A085")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo ARMA(", p, ",", q, ")"),
      subtitle = "Comparacion de valores reales vs. estimacion mixta",
      x = "Tiempo (Periodos)",
      y = "Valor",
      color = "Referencia:",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: La linea verde debe 'abrazar' los movimientos de la oscura.\n",
                       " SOMBRA: Zona de confianza. El 95% de los datos deben caer aqui.\n",
                       " NOTA: El ARMA es mas potente que el AR o MA por separado,\n",
                       "   logrando un equilibrio entre tendencia y choques aleatorios.")
    ) +
    mi_tema
  
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2: Evaluacion visual del ajuste mixto > ")}
  
  
  # --- GRAFICA 3: RESIDUOS (ERRORES) VS TIEMPO ---
  
  # CORRECCION: stats::sd
  sd_res_arma <- stats::sd(df$Residuos, na.rm = TRUE)
  
  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_arma, ymax = 2*sd_res_arma),
                         fill = "#16A085", alpha = 0.1) +
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.5) +
    ggplot2::geom_point(color = "#16A085", alpha = 0.4, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(-2*sd_res_arma, 2*sd_res_arma),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.7) +
    ggplot2::labs(
      title = "3. Analisis de Residuos (E vs t)",
      subtitle = "Verificacion de Aleatoriedad y Media Cero",
      x = "Tiempo (Periodos)",
      y = "Residuo (Error)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los errores deben ser una 'nube' aleatoria centrada en cero.\n",
                       " ROJO: Puntos fuera de estas lineas son shocks no explicados (Outliers).\n",
                       " MAL: Si ves ciclos, 'olas' o tendencias, el modelo no es suficiente.\n",
                       " NOTA: El ARMA busca que el error de hoy no dependa del error de ayer.")
    ) +
    mi_tema
  
  print(g3)
  if (interactive()) {readline(prompt = "Grafica 3: Inspeccion de la calidad del error residual > ")}
  
  
  # --- GRAFICA 4: ESTRUCTURA (CIRCULO UNITARIO AR Y MA) ---
  
  # CORRECCION defensiva: validar is.null
  if (!is.null(x$data_raices) && nrow(x$data_raices) > 0) {
    theta_circ <- seq(0, 2*pi, length.out = 200)
    df_circle_arma <- data.frame(x = cos(theta_circ), y = sin(theta_circ))
    
    g4 <- ggplot2::ggplot() +
      ggplot2::geom_path(data = df_circle_arma, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_polygon(data = df_circle_arma, ggplot2::aes(x, y),
                            fill = "#BDC3C7", alpha = 0.1) +
      ggplot2::geom_vline(xintercept = 0, color = "gray80") +
      ggplot2::geom_hline(yintercept = 0, color = "gray80") +
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario, color = Tipo, shape = Tipo),
                          size = 4, stroke = 1.2) +
      ggplot2::scale_color_manual(values = c("AR (Estabilidad)" = "#E74C3C",
                                             "MA (Invertibilidad)" = "#2980B9")) +
      ggplot2::scale_shape_manual(values = c("AR (Estabilidad)" = 16,
                                             "MA (Invertibilidad)" = 18)) +
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::labs(
        title = "4. Estructura de Raices (Circulo Unitario)",
        subtitle = paste0("Diagnostico Dual: ARMA(", p, ",", q, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        color = "Componente:", shape = "Componente:",
        caption = paste0("INTERPRETACION:\n",
                         " PUNTOS ROJOS (AR): Si estan dentro, el modelo es ESTABLE.\n",
                         " ROMBOS AZULES (MA): Si estan dentro, el modelo es INVERTIBLE.\n",
                         " CRITICO: Si algun punto toca o sale del circulo gris,\n",
                         "   los resultados estadisticos NO son confiables.")
      ) +
      mi_tema + ggplot2::theme(legend.position = "right")
    
    print(g4)
    
  } else {
    # CORRECCION: Usar message
    message("El modelo ARMA(0,0) no tiene raices; es ruido blanco.")
  }
  if (interactive()) {readline(prompt = "Grafica 4: Verificacion de estabilidad e invertibilidad mixta > ")}
  
  
  # --- GRAFICA 5: AUTOCORRELACION DE RESIDUOS (ACF) ---
  
  acf_res_arma <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
  df_acf_arma <- data.frame(
    Lag = as.numeric(acf_res_arma$lag[-1]),
    ACF = as.numeric(acf_res_arma$acf[-1])
  )
  
  limite_arma <- 1.96 / sqrt(nrow(df))
  
  df_acf_arma$Etiqueta <- ifelse(abs(df_acf_arma$ACF) > limite_arma, as.character(df_acf_arma$Lag), "")
  df_acf_arma$vjust_pos <- ifelse(df_acf_arma$ACF > 0, -0.5, 1.5)
  
  g5 <- ggplot2::ggplot(df_acf_arma, ggplot2::aes(x = Lag, y = ACF)) +
    ggplot2::geom_col(fill = "#34495E", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                        linetype = "dashed", color = "#16A085", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "5. Autocorrelacion de Residuos (ACF)",
      subtitle = "Quedo estructura temporal en el error?",
      x = "Lag (Retraso)",
      y = "Coeficiente ACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Si todas las barras estan 'atrapadas' en las lineas verdes.\n",
                       " NUMEROS ROJOS: Indican que lags NO son ruido blanco.\n",
                       " NOTA: El ARMA busca 'limpiar' la ACF. Si el Lag 1 o 2 es rojo,\n",
                       "   el modelo mixto aun no es optimo.")
    ) +
    mi_tema
  
  print(g5)
  if (interactive()) {readline(prompt = "Grafica 5: Diagnostico de independencia de los errores mixtos > ")}
  
  
  # --- GRAFICA 6: AUTOCORRELACION PARCIAL (PACF) DE RESIDUOS ---
  
  pacf_res_arma <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)
  
  df_pacf_arma <- data.frame(
    Lag = as.numeric(pacf_res_arma$lag),
    PACF = as.numeric(pacf_res_arma$acf)
  )
  
  limite_arma <- 1.96 / sqrt(nrow(df))
  
  df_pacf_arma$Etiqueta <- ifelse(abs(df_pacf_arma$PACF) > limite_arma, as.character(df_pacf_arma$Lag), "")
  df_pacf_arma$vjust_pos <- ifelse(df_pacf_arma$PACF > 0, -0.5, 1.5)
  
  g6 <- ggplot2::ggplot(df_pacf_arma, ggplot2::aes(x = Lag, y = PACF)) +
    ggplot2::geom_col(fill = "#16A085", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                        linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelacion Parcial (PACF)",
      subtitle = "Identificacion de dependencias directas en errores mixtos",
      x = "Lag (Retraso)",
      y = "Coeficiente PACF",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Ninguna barra sobresale de las lineas naranjas.\n",
                       " NUMEROS ROJOS: Indican que ese lag especifico aun tiene 'senal'.\n",
                       " NOTA: Si el Lag 1 es significativo, el componente AR(", p, ") \n",
                       "   no esta capturando bien la inercia inmediata.")
    ) +
    mi_tema
  
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6: Auditoria final de la estructura del error > ")}
  
  
  # --- GRAFICA 7: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---
  
  g7 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
    ggplot2::stat_qq(color = "#16A085", alpha = 0.6, size = 2) +
    ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "7. Normalidad de los Errores (Q-Q Plot)",
      subtitle = "Comparacion de Cuantiles de Residuos Mixtos vs. Normal Teorica",
      x = "Cuantiles Teoricos (Normal)",
      y = "Cuantiles Observados (Residuos ARMA)",
      caption = paste0("INTERPRETACION:\n",
                       " BIEN: Los puntos deben seguir la linea diagonal punteada.\n",
                       " COLAS: Si los extremos se curvan, hay eventos extremos (outliers).\n",
                       " NOTA: La normalidad asegura que el error tipico sea consistente.\n",
                       " CLAVE: Si los puntos se alejan mucho, los p-valores del modelo son dudosos.")
    ) +
    mi_tema
  
  print(g7)
  if (interactive()) {readline(prompt = "Grafica 7: Verificacion de la distribucion gaussiana del error > ")}
  
  
  # --- GRAFICA 8: PRECISION DEL MODELO (REAL VS AJUSTADO) ---
  
  cor_arma <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")
  
  g8 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    ggplot2::geom_point(color = "#16A085", alpha = 0.5, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "8. Precision del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlacion de Pearson: ", round(cor_arma, 3)),
      x = "Valor Ajustado (Modelo ARMA)",
      y = "Valor Observado (Realidad)",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA ROJA: El ideal donde el modelo predice exactamente la realidad.\n",
                       " BIEN: Los puntos verdes deben seguir la trayectoria de la linea roja.\n",
                       " DISPERSION: Si los puntos estan muy lejos de la linea, el modelo es pobre.\n",
                       " CLAVE: Un ARMA exitoso debe tener una correlacion superior a 0.7.")
    ) +
    mi_tema
  
  print(g8)
  if (interactive()) {readline(prompt = "Grafica 8: Analisis de correlacion y capacidad predictiva > ")}
  
  
  # --- GRAFICA 9: DIAGNOSTICO DE ESTACIONARIEDAD (MOVIL) ---
  
  ventana_arma <- max(5, round(nrow(df) * 0.1))
  
  df_est_arma <- df
  
  df_est_arma$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_arma, ventana_arma), sides = 2))
  
  # CORRECCION: stats::sd
  df_est_arma$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_arma):min(nrow(df), i+ventana_arma)
    stats::sd(df$Observado[idx], na.rm = TRUE)
  })
  
  df_final_grafica <- df_est_arma[!is.na(df_est_arma$Media_Movil), ]
  
  g9 <- ggplot2::ggplot(df_final_grafica, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray85", linewidth = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#16A085", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#16A085", linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
    ggplot2::labs(
      title = "9. Diagnostico de Estacionariedad",
      subtitle = paste0("Analisis de Media y Varianza Movil (Ventana k = ", ventana_arma, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("INTERPRETACION:\n",
                       " LINEA VERDE: Media local. Debe ser horizontal y estable.\n",
                       " LINEA ROJA: Promedio global. La base del modelo.\n",
                       " BIEN: Si la franja verde mantiene un ancho constante.\n",
                       " NOTA: Se han omitido los extremos para evitar errores de calculo movil.")
    ) +
    mi_tema
  
  print(g9)
  
  # CORRECCION: Usar message y devolver el objeto x
  if (interactive()) {
    message("Recuerda que un buen modelo ARMA es el equilibrio entre parsimonia y precision.")
    message("Diagnostico ARMA completado.")
  }
  
  return(invisible(x))
}