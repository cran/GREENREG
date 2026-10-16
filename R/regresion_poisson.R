#' Regresión de Poisson (Modelos de Conteo)
#'
#' Esta función ajusta un Modelo Lineal Generalizado (GLM) de la familia Poisson con
#' función de enlace logaritmica \eqn{\ln(\lambda)}. Esta disenada para modelar variables
#' de respuesta que representan conteos enteros no negativos (0, 1, 2, 3...) y tasas
#' de incidencia en estudios agronomicos, forestales y sociales.
#'
#' @details
#' El flujo interno de la función se divide en seis etapas criticas:
#' \enumerate{
#'   \item **Ajuste Log-Lineal:** Estima la relación entre los predictores y el logaritmo
#'         de la tasa de eventos (\eqn{\lambda}). Utiliza Maxima Verosimilitud (MLE) para el ajuste.
#'   \item **Interpretación de Tasas (IRR):** Los coeficientes originales (\eqn{\beta}) se
#'         exponencian automaticamente \eqn{e^{\beta}} para obtener los **Ratios de Tasas de Incidencia**.
#'         Esto permite concluir en que porcentaje aumenta o disminuye la frecuencia de eventos
#'         por cada unidad de cambio en X.
#'   \item **Diagnóstico de Equidispersion:** Calcula el Ratio de Dispersión basado en el
#'         estadístico de Pearson \eqn{(X^2 / gl)}. En una distribución Poisson perfecta, este
#'         valor debe ser 1. Valores > 1.5 indican **Sobredispersión**, sugiriendo que los
#'         datos varian más de lo que el modelo permite.
#'   \item **Inferencia Robusta:** Genera intervalos de confianza para los IRR, permitiendo
#'         identificar factores de riesgo (IRR > 1) y factores protectores (IRR < 1).
#'   \item **Evaluación de Residuos:** Calcula residuos de devianza y Pearson para identificar
#'         patrones de error sistematicos o falta de ajuste en conteos especificos (como el cero).
#'   \item **Traduccion Didactica:** Genera un reporte que advierte si el modelo es válido
#'         o si se debe migrar a un modelo Binomial Negativa debido a la dispersión.
#' }
#'
#' @section Guia de las 8 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validación para conteos:
#' \describe{
#'   \item{\strong{1. Distribución de Conteos (Real vs Modelo)}}{Compara las frecuencias observadas contra las que el modelo "esperaba" predecir. Crucial para detectar el exceso de ceros.}
#'   \item{\strong{2. Precisión (Predicho vs Real)}}{Scatterplot entre el Lambda estimado y el conteo real. Los puntos deben seguir la tendencia de la línea roja diagonal.}
#'   \item{\strong{3. Importancia de Variables (IRR)}}{Visualiza el impacto relativo de cada predictor. Los puntos verdes son significativos; si el intervalo cruza el 1, la variable es ruido.}
#'   \item{\strong{4. Verificación de Sobredispersión}}{Gráfica la relación Media-Varianza. Si la línea verde (datos) se eleva sobre la negra (teoria), el supuesto Poisson se rompe.}
#'   \item{\strong{5. Distancia de Cook}}{Identifica casos atípicos en el conteo que tienen un peso excesivo en la estimación de la tasa promedio.}
#'   \item{\strong{6. Residuos de Devianza}}{Analiza la calidad del ajuste en escala logaritmica. Busca una nube uniforme sin formas de "U" o de embudo.}
#'   \item{\strong{7. Apalancamiento (Leverage)}}{Detecta observaciones con valores de X inusuales que controlan la pendiente del modelo.}
#'   \item{\strong{8. Densidad de Ajuste Global}}{Superpone las curvas de probabilidad del modelo y la realidad para evaluar el solapamiento total.}
#' }
#'
#' @param formula Una fórmula de R (ej. \code{numero_frutos ~ riego + fertilizante}).
#' @param data Un \code{data.frame} con las observaciones (el conteo debe ser entero).
#'
#' @return Un objeto de clase \code{"poisson_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{glm} original de la familia Poisson.
#'   \item \code{irr_table}: Tabla con los Ratios de Tasas de Incidencia, intervalos y p-valores.
#'   \item \code{dispersión}: El ratio de dispersión calculado (Pearson / GL).
#'   \item \code{metricas}: Lista con AIC, Log-Verosimilitud y Devianzas.
#'   \item \code{predicciones}: Dataframe con valores observados, lambda estimado y diagnósticos (leverage/residuos).
#' }
#'
#' @examples
#' data(datos_tala_ilegal_poisson)
#' modelo_poisson <- greenreg_reg_poisson(numero_talas ~ densidad_poblacion +
#'                               superficie_forestal_ha + dias_vigilancia,
#'                               data = datos_tala_ilegal_poisson)
#' plot(modelo_poisson)
#' @importFrom stats glm poisson coef confint formula na.omit predict fitted residuals df.residual pnorm qnorm logLik pchisq cooks.distance hatvalues rstandard dpois quantile
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
greenreg_reg_poisson <- function(formula, data) {
  # --- 1. Ajuste del Modelo ---
  if (!inherits(formula, "formula")) stop("El argumento 'f\u00f3rmula' debe ser una f\u00f3rmula (ej. y ~ x1 + x2).")
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")
  var_y_check <- all.vars(formula)[1]
  y_check <- data[[var_y_check]]
  if (!is.numeric(y_check) || any(y_check < 0, na.rm = TRUE) || any(y_check != floor(y_check), na.rm = TRUE)) {
    stop(sprintf(
      "La variable respuesta '%s' debe ser un conteo: num\u00e9rica, entera y no negativa (0, 1, 2, ...).",
      var_y_check))
  }
  modelo_base <- stats::glm(formula, data = data, family = stats::poisson(link = "log"), na.action = stats::na.omit)

  # --- 2. Verificacion de Equidispersion (Ratio de Pearson) ---
  # Si el ratio es > 1.5 o 2, hay sobredispersion
  pearson_chisq <- sum(stats::residuals(modelo_base, type = "pearson")^2)
  ratio_dispersion <- pearson_chisq / modelo_base$df.residual

  # --- 3. Coeficientes e IRR ---
  resumen <- summary(modelo_base)
  coefs_log <- stats::coef(modelo_base)
  irr <- exp(coefs_log)
  ci_irr <- exp(stats::confint(modelo_base))
  tabla_irr <- cbind(IRR = irr, `2.5 %` = ci_irr[,1], `97.5 %` = ci_irr[,2], `P-valor` = resumen$coefficients[,4])

  df_preds <- data.frame(
    obs = modelo_base$y,
    pred = stats::fitted(modelo_base),
    leverage = stats::hatvalues(modelo_base),
    resid_std = stats::rstandard(modelo_base)
  )
  # --- 4. Empaquetado ---
  resultado <- list(
    formula = formula,
    modelo = modelo_base,
    predicciones = df_preds,
    coeficientes = resumen$coefficients,
    irr_table = tabla_irr,
    dispersion = ratio_dispersion,
    metricas = list(
      aic = modelo_base$aic,
      loglik = as.numeric(stats::logLik(modelo_base)),
      null_dev = modelo_base$null.deviance,
      res_dev = modelo_base$deviance
    )
  )
  class(resultado) <- "poisson_greenreg"
  return(resultado)
}

#' @rdname greenreg_reg_poisson
#' @param x Objeto de clase.
#' @param ... Argumentos adicionales.
#' @export
print.poisson_greenreg <- function(x, ...) {

  cat("\n==============================================================\n")
  cat("            REGRESI\u00d3N POISSON (MODELO DE CONTEO) \n")
  cat("==============================================================\n\n")

  # --- 1. MODELO Y ECUACION ---
  cat("--- 1. MODELO Y ECUACI\u00d3N ---\n")
  cat("F\u00f3rmula:", deparse(x$formula), "\n")

  # Espacio de nombres stats::
  b <- signif(stats::coef(x$modelo), 4)
  vars <- names(b)
  ecuacion <- paste0("log(mu) = ", b[1])
  if(length(b) > 1) {
    for(i in 2:length(b)) ecuacion <- paste0(ecuacion, " + (", b[i], " * ", vars[i], ")")
  }
  cat("Ecuaci\u00f3n estimada:\n", ecuacion, "\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS (IRR) ---\n")
  cat("Interpretaci\u00f3n en Tasas de Incidencia (Incidence Rate Ratios):\n")
  print(round(x$irr_table, 4))
  cat("\n NOTA: Los IRR indican cuanto se multiplica la tasa de eventos.\n")
  cat("   - IRR > 1: Efecto multiplicativo (Riesgo/Aumento).\n")
  cat("   - IRR < 1: Efecto de reducci\u00f3n (Protecci\u00f3n/Disminuci\u00f3n).\n\n")

  # --- 3. EVALUACION DEL MODELO ---
  cat("--- 3. EVALUACI\u00d3N DEL MODELO (Prueba de Raz\u00f3n de Verosimilitud) ---\n")

  chi_cuad <- x$metricas$null_dev - x$metricas$res_dev
  df_chi <- x$modelo$df.null - x$modelo$df.residual

  # Blindaje logico contra valores nulos
  if (!is.null(chi_cuad) && !is.na(chi_cuad) && !is.null(df_chi) && !is.na(df_chi) && df_chi > 0) {
    p_val_lrt <- 1 - stats::pchisq(chi_cuad, df_chi)
    p_val_txt <- if(is.na(p_val_lrt)) "NA" else format.pval(p_val_lrt)
    concl_txt <- if(!is.na(p_val_lrt) && p_val_lrt < 0.05) "El modelo explica significativamente los datos" else "No se rechaza que el modelo sea equivalente al modelo nulo (sin evidencia suficiente de mejora)"
  } else {
    chi_cuad <- NA
    df_chi <- NA
    p_val_txt <- "NA"
    concl_txt <- "No disponible o indeterminado"
  }

  chi_txt <- if(is.na(chi_cuad)) "NA" else sprintf("%.3f", chi_cuad)
  df_txt <- if(is.na(df_chi)) "NA" else as.character(df_chi)

  cat(sprintf(" Chi-cuadrado Global: %s | GL: %s\n", chi_txt, df_txt))
  cat(" P-valor (LRT):", p_val_txt, "\n")
  cat(" Conclusi\u00f3n:", concl_txt, "\n\n")

  # --- 4. BONDAD DE AJUSTE ---
  cat("--- 4. BONDAD DE AJUSTE ---\n")

  pseudo_r2 <- 1 - (x$metricas$res_dev / x$metricas$null_dev)

  # Proteccion estricta para sprintf
  loglik_txt <- if(is.null(x$metricas$loglik) || is.na(x$metricas$loglik)) "NA" else sprintf("%10.2f", x$metricas$loglik)
  resdev_txt <- if(is.null(x$metricas$res_dev) || is.na(x$metricas$res_dev)) "NA" else sprintf("%10.2f", x$metricas$res_dev)
  pseudo_txt <- if(is.null(pseudo_r2) || is.na(pseudo_r2)) "NA" else sprintf("%7.2f%%", pseudo_r2 * 100)

  cat(sprintf(" Log-Verosimilitud: %s\n", loglik_txt))
  cat(sprintf(" Devianza Residual: %s\n", resdev_txt))
  cat(sprintf(" Pseudo R2 (McFadden): %s\n", pseudo_txt))
  cat("  (Representa la reducci\u00f3n de la incertidumbre/error respecto al modelo nulo).\n\n")

  # --- 5. SIGNIFICANCIA DE COEFICIENTES ---
  cat("--- 5. SIGNIFICANCIA DE COEFICIENTES (Prueba Wald) ---\n")

  # Blindaje extra convirtiendo a matrix
  if (!is.null(x$coeficientes) && nrow(x$coeficientes) > 0) {
    stats::printCoefmat(as.matrix(x$coeficientes), digits = 4, signif.stars = TRUE)
  } else {
    cat(" No se pudieron estimar coeficientes validos.\n")
  }
  cat("\n")

  # --- 6. VERIFICACION DE SUPUESTOS ---
  cat("--- 6. VERIFICACI\u00d3N DE SUPUESTOS ---\n")

  # A) Equidispersion
  disp <- x$dispersion
  disp_txt <- if(is.null(disp) || is.na(disp)) "NA" else sprintf("%.2f", disp)

  cat(sprintf(" Equidispersion (Ratio): %s\n", disp_txt))

  if(!is.null(disp) && !is.na(disp)) {
    if(disp > 1.5) {
      cat(" ALERTA: Hay SOBREDISPERSI\u00d3N (Ratio > 1.5). Los p-valores podr\u00edan estar subestimados.\n")
    } else {
      cat(" OKEY: Se cumple razonablemente el supuesto de equidispersion (Media = Varianza).\n")
    }
  }

  # B) Otros supuestos
  cat(" Independencia: Asumida por el dise\u00f1o de muestreo.\n")
  cat(" Linealidad: Se asume que log(Media) es lineal respecto a las variables X.\n\n")

  # --- 7. RESUMEN DEL MODELO ---
  cat("--- 7. RESUMEN DEL MODELO ---\n")

  # Validacion estructural de la tabla antes de hacer el filtro
  if (!is.null(x$irr_table) && "P-valor" %in% colnames(x$irr_table)) {
    p_vals <- x$irr_table[, "P-valor"]
    sig_vars <- rownames(x$irr_table)[!is.na(p_vals) & p_vals < 0.05]

    if(length(sig_vars) > 0) {
      cat(" Variables significativas identificadas:", paste(sig_vars, collapse = ", "), "\n")
    } else {
      cat(" No se encontraron variables con impacto estad\u00edstico significativo.\n")
    }
  }

  if(!is.null(disp) && !is.na(disp) && disp > 2) {
    cat("\n RECOMENDACI\u00d3N: Debido a la alta sobredispersi\u00f3n, se sugiere \n probar un modelo Binomial Negativa o usar errores est\u00e1ndar robustos.\n")
  }

  cat("\n==============================================================\n")

  return(invisible(x))
}


#' @rdname greenreg_reg_poisson
#' @description Genera una secuencia de 8 gráficas de diagnóstico para validar el modelo Poisson.
#' @param x Objeto de clase.
#' @param ... Argumentos adicionales.
#' @param pausar Lógico. Si \code{TRUE}, pausa con \code{readline()} entre cada
#'   gráfica (por defecto, solo en sesión interactiva).
#' @export
plot.poisson_greenreg <- function(x, ..., pausar = interactive()) {
  # --- 1. Extraccion y Validacion de Datos ---
  modelo <- x$modelo
  if (is.null(x$predicciones)) stop("Error interno: No se encontraron predicciones.")
  df_preds <- x$predicciones

  # Tema visual consistente
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10))
    )

  # Usar message() en lugar de cat()
  if (pausar) {
    message("Generando 8 gr\u00e1ficas de diagn\u00f3stico... (Presiona [Enter] para avanzar)")
  }

  # --- GRAFICA 1: CONTEOS OBSERVADOS VS PREDICHOS ---

  # Determinamos el limite del eje X (evitando valores extremos que arruinen la visualizacion)
  max_val_obs <- max(df_preds$obs, na.rm = TRUE)
  max_count <- round(stats::quantile(df_preds$obs, 0.99, na.rm = TRUE))
  if (is.na(max_count) || max_count < 1) max_count <- max_val_obs

  # A) Frecuencia Real (Observada)
  obs_limpios <- df_preds$obs[!is.na(df_preds$obs)]
  freq_real <- table(factor(obs_limpios, levels = 0:max_count))
  df_counts <- data.frame(Conteo = 0:max_count, Freq = as.numeric(freq_real), Tipo = "Realidad")

  # B) Frecuencia Esperada (Modelo Poisson)
  # Calculamos la probabilidad de cada conteo k dado el lambda (pred) de cada fila
  probs_matrix <- sapply(0:max_count, function(k) {
    sum(stats::dpois(k, lambda = df_preds$pred), na.rm = TRUE)
  })

  df_expected <- data.frame(Conteo = 0:max_count, Freq = probs_matrix, Tipo = "Modelo (Esperado)")

  # Combinamos ambos para ggplot
  df_comparacion <- rbind(df_counts, df_expected)

  g1 <- ggplot2::ggplot(df_comparacion, ggplot2::aes(x = factor(Conteo), y = Freq, fill = Tipo)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Realidad" = "#34495E", "Modelo (Esperado)" = "#E67E22")) +
    ggplot2::labs(
      title = "1. Distribuci\u00f3n de Conteos (Real vs Modelo)",
      subtitle = "Predice el modelo la cantidad correcta de eventos?",
      x = "N\u00famero de Eventos (Conteo)",
      y = "Frecuencia Total",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Las barras naranjas (Modelo) tienen casi la misma altura que las grises (Realidad).\n",
                       " ALERTA: Si la barra gris del '0' es gigante y la naranja peque\u00f1a, hay exceso de ceros.\n",
                       " DISPERSI\u00d3N: Si las barras reales son m\u00e1s anchas que las del modelo, hay sobredispersi\u00f3n.")
    ) +
    mi_tema

  print(g1)
  if (pausar) {readline(prompt = "Gr\u00e1fica 1/8 > ")}


  # --- GRAFICA 2: VALORES PREDICHOS VS VALORES REALES ---

  # Calculamos la correlacion para mostrarla en el subtitulo
  cor_val <- stats::cor(df_preds$obs, df_preds$pred, method = "pearson")

  g2 <- ggplot2::ggplot(df_preds, ggplot2::aes(x = pred, y = obs)) +
    # Puntos con ligera transparencia para ver densidad (overplotting)
    ggplot2::geom_point(color = "#2E86C1", alpha = 0.5, size = 2) +
    # Linea de Identidad (Lo que seria un modelo perfecto)
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    # Suavizado greenreg_logistico/tendencia para ver el comportamiento real
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "poisson"),
                         color = "#F1C40F", se = TRUE, fill = "#F1C40F", alpha = 0.2) +
    ggplot2::labs(
      title = "2. Precisi\u00f3n del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlaci\u00f3n de Pearson: ", round(cor_val, 3)),
      x = "Valor Predicho (Lambda estimado)",
      y = "Valor Real (Conteo observado)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " L\u00cdNEA ROJA: Es la diagonal de perfeccion (donde Predicho = Real).\n",
                       " L\u00cdNEA AMARILLA: Es la tendencia real de tus datos.\n",
                       " BIEN: Si los puntos y la l\u00ednea amarilla siguen de cerca a la roja.\n",
                       " ALERTA: Si los puntos se alejan mucho hacia arriba conforme aumenta X, hay subestimaci\u00f3n.")
    ) +
    mi_tema

  print(g2)
  if (pausar) { readline(prompt = "Gr\u00e1fica 2/8 > ")}


  # --- GRAFICA 3: IMPORTANCIA DE VARIABLES (IRR) ---

  # 1. Preparar los datos desde la tabla de IRR generada en la funcion principal
  df_irr <- as.data.frame(x$irr_table)
  df_irr$Variable <- rownames(df_irr)

  # 2. Limpieza: Excluimos el Intercepto para no distorsionar la escala visual
  df_irr <- df_irr[df_irr$Variable != "(Intercept)", ]

  # Renombrar columnas para facilitar el manejo en ggplot
  colnames(df_irr) <- c("IRR", "Low", "High", "P_val", "Variable")

  # 3. Construccion de la Grafica
  g3 <- ggplot2::ggplot(df_irr, ggplot2::aes(x = reorder(Variable, IRR), y = IRR)) +
    # Linea de referencia en 1 (indica que no hay efecto)
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Barras de error (Intervalos de Confianza al 95%)
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Low, ymax = High), width = 0.2, color = "#7F8C8D") +
    # Puntos del IRR colorizados por significancia
    ggplot2::geom_point(ggplot2::aes(color = P_val < 0.05), size = 4) +
    # Etiquetas con el valor del IRR sobre cada punto
    ggplot2::geom_text(ggplot2::aes(label = round(IRR, 2)), vjust = -1.2, size = 3.5, fontface = "bold") +
    # Invertimos ejes para que los nombres de variables sean legibles (horizontal)
    ggplot2::coord_flip() +
    # Colores semanticos: Verde para significativo, Rojo para no significativo
    ggplot2::scale_color_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#E74C3C"),
                                labels = c("No Significativo", "Significativo (p < 0.05)")) +
    ggplot2::labs(
      title = "3. Importancia de Variables (IRR)",
      subtitle = "Impacto relativo en la tasa de eventos (Incidence Rate Ratio)",
      x = "Variables Predictoras",
      y = "IRR (Razon de Tasas de Incidencia)",
      color = "Estado del Coeficiente:",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " VERDE: La variable tiene un impacto real (p < 0.05).\n",
                       " L\u00cdNEA EN 1: Si el intervalo cruza esta l\u00ednea, la variable no influye.\n",
                       " DERECHA (>1): La variable INCREMENTA la frecuencia de eventos.\n",
                       " IZQUIERDA (<1): La variable REDUCE la frecuencia de eventos.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "top")

  print(g3)
  if (pausar) {readline(prompt = "Gr\u00e1fica 3/8 > ")}


  # --- GRAFICA 4: VERIFICACION DE SOBREDISPERSION ---

  # Preparar datos: La varianza empirica se aproxima con el cuadrado de los residuos
  df_disp <- data.frame(
    predicho = df_preds$pred,
    res_sq = (df_preds$obs - df_preds$pred)^2
  )

  g4 <- ggplot2::ggplot(df_disp, ggplot2::aes(x = predicho, y = res_sq)) +
    # Puntos de dispersion
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.4, size = 2) +
    # Linea de Equidispersion (Varianza = Media)
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    # Tendencia de la varianza real
    ggplot2::geom_smooth(method = "gam", color = "#2ECC71", fill = "#2ECC71", alpha = 0.2) +
    ggplot2::labs(
      title = "4. Verificaci\u00f3n de Sobredispersi\u00f3n",
      subtitle = "Se cumple el supuesto de Media = Varianza?",
      x = "Media Predicha ",
      y = "Varianza Observada (Residuos al cuadrado)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " L\u00cdNEA NEGRA: Representa la Equidispersion teorica (Poisson perfecta).\n",
                       " L\u00cdNEA VERDE: Es la relaci\u00f3n real entre media y varianza en tus datos.\n",
                       " BIEN: Si la l\u00ednea verde sigue de cerca a la negra punteada.\n",
                       " ALERTA: Si la l\u00ednea verde sube mucho m\u00e1s rapido que la negra, hay SOBREDISPERSI\u00d3N.")
    ) +
    mi_tema

  print(g4)
  if (pausar) { readline(prompt = "Gr\u00e1fica 4/8 > ")}


  # --- GRAFICA 5: DATOS INFLUYENTES (DISTANCIA DE COOK) ---

  # 1. Calculo de la Distancia de Cook
  cooks_d <- stats::cooks.distance(x$modelo)
  n_obs <- length(cooks_d)

  # Umbral critico comun: 4/n
  umbral_cook <- 4 / n_obs

  df_cook <- data.frame(
    ID = 1:n_obs,
    Distancia = cooks_d
  )

  # Etiquetamos solo los puntos que superan el umbral
  df_cook$Label <- ifelse(df_cook$Distancia > umbral_cook, as.character(df_cook$ID), "")

  g5 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = ID, y = Distancia)) +
    # Lineas verticales
    ggplot2::geom_segment(ggplot2::aes(xend = ID, yend = 0), color = "darkgrey", alpha = 0.5) +
    # Puntos destacados
    ggplot2::geom_point(ggplot2::aes(color = Distancia > umbral_cook), size = 2) +
    # Linea de umbral critico
    ggplot2::geom_hline(yintercept = umbral_cook, linetype = "dashed", color = "#E74C3C", linewidth = 0.8)

  # Check defensivo de ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g5 <- g5 + ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold")
  } else {
    g5 <- g5 + ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -1, size = 3.5, fontface = "bold")
  }

  g5 <- g5 +
    # Colores: Rojo para los que superan el umbral
    ggplot2::scale_color_manual(values = c("FALSE" = "#3498DB", "TRUE" = "#E74C3C")) +
    ggplot2::labs(
      title = "5. Datos Influyentes (Distancia de Cook)",
      subtitle = paste0("Umbral cr\u00edtico (4/n): ", round(umbral_cook, 4)),
      x = "\u00cdndice de la Observaci\u00f3n (Fila)",
      y = "Distancia de Cook",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " PUNTOS ROJOS: Son observaciones que tienen un impacto excesivo en el modelo.\n",
                       " UMBRAL: Si un punto est\u00e1 muy por encima de la l\u00ednea roja, revisa sus datos.\n",
                       " ACCI\u00d3N: Podr\u00edan ser errores de dedo o casos excepcionales que distorsionan los IRR.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g5)
  if (pausar) {readline(prompt = "Gr\u00e1fica 5/8 > ")}


  # --- GRAFICA 6: ANALISIS DE RESIDUOS (DEVIANZA VS AJUSTADOS) ---

  # Extraemos residuos de devianza (los estandar para GLM Poisson)
  df_res <- data.frame(
    ajustados = log(df_preds$pred), # Usamos el log(lambda) para linearizar la escala
    residuos = stats::residuals(x$modelo, type = "deviance")
  )

  g6 <- ggplot2::ggplot(df_res, ggplot2::aes(x = ajustados, y = residuos)) +
    # Puntos de los residuos
    ggplot2::geom_point(color = "#16A085", alpha = 0.5, size = 2) +
    # Linea central en 0
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Lineas de umbral (+-2)
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    # Tendencia suavizada (Loess)
    ggplot2::geom_smooth(method = "loess", color = "#E67E22", fill = "#E67E22", alpha = 0.2, formula = y ~ x) +
    ggplot2::labs(
      title = "6. An\u00e1lisis de Residuos de Devianza",
      subtitle = "Residuos vs. Valores Ajustados (escala log)",
      x = "Valores Ajustados (Log-Lambda)",
      y = "Residuos de Devianza",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Los puntos se reparten al azar arriba y abajo del cero (nube uniforme).\n",
                       " CURVA: Si la l\u00ednea naranja hace una 'U', te falta una variable cuadratica o interacci\u00f3n.\n",
                       " EMBUDO: Si la dispersi\u00f3n aumenta hacia la derecha, confirma la Sobredispersi\u00f3n.\n",
                       " OUTLIERS: Puntos fuera de las l\u00edneas rojas (+-2) est\u00e1n muy mal predichos.")
    ) +
    mi_tema

  print(g6)
  if (pausar) {readline(prompt = "Gr\u00e1fica 6/8 > ")}


  # --- GRAFICA 7: APALANCAMIENTO (LEVERAGE) ---

  # 1. Preparar datos de Leverage (h)
  lev <- x$predicciones$leverage
  n_obs <- length(lev)
  p_vars <- length(stats::coef(x$modelo)) # Numero de parametros (k)

  # Umbral comun: 2k/n
  umbral_lev <- (2 * p_vars) / n_obs

  df_lev <- data.frame(
    ID = 1:n_obs,
    Leverage = lev
  )

  # Etiquetar solo puntos con alto apalancamiento
  df_lev$Label <- ifelse(df_lev$Leverage > umbral_lev, as.character(df_lev$ID), "")

  g7 <- ggplot2::ggplot(df_lev, ggplot2::aes(x = ID, y = Leverage)) +
    # Lineas verticales para cada punto
    ggplot2::geom_segment(ggplot2::aes(xend = ID, yend = 0), color = "lightgrey", alpha = 0.6) +
    # Puntos con tamano basado en su propio valor de leverage
    ggplot2::geom_point(ggplot2::aes(color = Leverage > umbral_lev, size = Leverage)) +
    # Linea de umbral critico (2k/n)
    ggplot2::geom_hline(yintercept = umbral_lev, linetype = "dashed", color = "#2980B9", linewidth = 0.8)

  # Check defensivo de ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g7 <- g7 + ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold")
  } else {
    g7 <- g7 + ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -1, size = 3.5, fontface = "bold")
  }

  g7 <- g7 +
    # Escala de colores (Azul para normal, Naranja para alto leverage)
    ggplot2::scale_color_manual(values = c("FALSE" = "#BDC3C7", "TRUE" = "#E67E22")) +
    ggplot2::labs(
      title = "7. An\u00e1lisis de Apalancamiento (Leverage)",
      subtitle = paste0("Umbral de dise\u00f1o (2k/n): ", round(umbral_lev, 4)),
      x = "\u00cdndice de la Observaci\u00f3n (Fila)",
      y = "Valor de Leverage (h)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " PUNTOS NARANJA: Son observaciones con valores de X muy poco comunes.\n",
                       " UMBRAL AZUL: Marca el limite donde un punto empieza a tener mucho peso.\n",
                       " NOTA: Un alto leverage no es malo por si solo, pero si ese punto \n",
                       "   tambien tiene un residuo alto (Gr\u00e1fica 6), est\u00e1 da\u00f1ando tu modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g7)
  if (pausar) {readline(prompt = "Gr\u00e1fica 7/8 > ")}


  # --- GRAFICA 8: AJUSTE DEL MODELO (DENSIDAD REAL VS PREDICHA) ---

  # Creamos un dataframe largo para facilitar el mapeo de colores en ggplot
  df_ajuste <- data.frame(
    Valor = c(df_preds$obs, df_preds$pred),
    Fuente = rep(c("Observado (Real)", "Ajustado (Modelo)"), each = nrow(df_preds))
  )

  g8 <- ggplot2::ggplot(df_ajuste, ggplot2::aes(x = Valor, fill = Fuente)) +
    # Usamos geom_density con alfa para ver la superposicion
    ggplot2::geom_density(alpha = 0.4, color = "white") +
    # Anadimos lineas de media para cada grupo
    ggplot2::geom_vline(xintercept = mean(df_preds$obs), color = "#2C3E50", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = mean(df_preds$pred), color = "#D35400", linetype = "dashed", linewidth = 1) +
    ggplot2::scale_fill_manual(values = c("Observado (Real)" = "#34495E", "Ajustado (Modelo)" = "#E67E22")) +
    ggplot2::labs(
      title = "8. Ajuste Global del Modelo",
      subtitle = "Comparaci\u00f3n de Densidades: Realidad vs. Predicci\u00f3n",
      x = "N\u00famero de Eventos",
      y = "Densidad (Probabilidad suavizada)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " AJUSTE \u00d3PTIMO: Si ambas 'monta\u00f1as' coinciden en posici\u00f3n y ancho.\n",
                       " DESPLAZAMIENTO: Si las l\u00edneas punteadas (medias) est\u00e1n lejos, el modelo es sesgado.\n",
                       " FORMA: Si la curva naranja es mucho m\u00e1s delgada que la gris, confirma la Sobredispersi\u00f3n.\n",
                       " NOTA: La l\u00ednea naranja representa la 'esperanza' del modelo Poisson.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "top")

  print(g8)

  # E-04: se devuelve la lista de graficas ggplot en vez de invisible(x)
  if (pausar) {
    message("Diagn\u00f3stico Poisson completado.")
  }
  return(invisible(list(g1 = g1, g2 = g2, g3 = g3, g4 = g4, g5 = g5,
                        g6 = g6, g7 = g7, g8 = g8)))
}