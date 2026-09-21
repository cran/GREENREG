#' Regresion Logistica Binaria (MLG Logit)
#'
#' Esta funcion ajusta un Modelo Lineal Generalizado (GLM) de la familia binomial para
#' predecir la probabilidad de un evento dicotomico (Si/No, 0/1). A diferencia de la
#' regresion lineal, el motor de calculo estima los coeficientes mediante Maxima
#' Verosimilitud (MLE) y ofrece un diagnostico especializado en clasificacion.
#'
#' @details
#' El flujo interno de la funcion se divide en cinco etapas fundamentales:
#' \enumerate{
#'   \item **Ajuste Logit:** Transforma la variable respuesta mediante la funcion de enlace
#'         \eqn{logit(p) = ln(p / (1-p))}, permitiendo que la probabilidad se mantenga entre 0 y 1.
#'   \item **Conversion de Impacto:** Los coeficientes estimados se transforman automaticamente
#'         a **Odds Ratios (OR)** mediante \eqn{e^\beta}. Esto permite concluir cuanto
#'         multiplica cada variable la posibilidad de que ocurra el evento.
#'   \item **Evaluacion de Bondad (McFadden):** Calcula el Pseudo-\eqn{R^2} comparando la
#'         verosimilitud del modelo ajustado contra un modelo nulo (sin predictores).
#'   \item **Metricas de Clasificacion:** Genera la Matriz de Confusion basada en el \code{umbral}
#'         especificado y calcula la Precision Global (Accuracy).
#'   \item **Diagnostico de Influencia:** Evalua si existen observaciones que deforman
#'         las probabilidades mediante el calculo de la Distancia de Cook y el Apalancamiento (Leverage).
#' }
#'
#' @section Guia de las 9 Graficas de Diagnostico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validacion predictiva:
#' \describe{
#'   \item{\strong{1. Capacidad de Separacion}}{Muestra si el modelo realmente distingue entre grupos. Buscamos "montanas" de densidad que no se traslapen.}
#'   \item{\strong{2. Matriz de Confusion Visual}}{Mapa de calor que resume aciertos (diagonal) y errores (falsos positivos/negativos).}
#'   \item{\strong{3. Odds Ratios (Forest Plot)}}{Visualiza la importancia de las variables. Los puntos a la derecha del 1 son factores de riesgo; a la izquierda, protectores.}
#'   \item{\strong{4. Curva de Calibracion}}{Compara la probabilidad que "dice" el modelo contra la frecuencia real observada. Debe seguir la diagonal.}
#'   \item{\strong{5. Curva ROC (AUC)}}{Mide la potencia global de diagnostico. Un AUC de 1.0 es una clasificacion perfecta; 0.5 es puro azar.}
#'   \item{\strong{6. Distancia de Cook}}{Detecta casos atipicos en la clasificacion que tienen un peso excesivo en los coeficientes.}
#'   \item{\strong{7. Trade-off de Umbrales}}{Permite ver como cambia la Sensibilidad y Especificidad al mover el punto de corte.}
#'   \item{\strong{8. Residuos de Pearson}}{Busca patrones sistematicos en el error. Idealmente, la linea roja debe mantenerse cerca del cero.}
#'   \item{\strong{9. Leverage vs Residuos}}{Identifica observaciones con combinaciones de variables X muy raras que "jalonean" las probabilidades.}
#' }
#'
#' @param formula Una formula de R (ej. \code{presencia_plaga ~ humedad + altitud}).
#' @param data Un \code{data.frame} con las observaciones.
#' @param umbral Punto de corte para clasificar como "Evento (1)" (por defecto 0.5).
#'
#' @return Un objeto de clase \code{"logistico_greenreg"}, que consiste en una lista con los siguientes componentes:
#' \itemize{
#'   \item \code{formula}: La formula de R aplicada al modelo.
#'   \item \code{modelo}: El objeto de clase \code{glm} original con el ajuste binomial.
#'   \item \code{variables}: Una lista que contiene el nombre de la variable respuesta (\code{y}).
#'   \item \code{coeficientes}: Matriz con las estimaciones de los coeficientes, errores estandar, estadisticos z y p-valores.
#'   \item \code{odds_ratios}: Tabla con los Odds Ratios (OR) calculados, sus intervalos de confianza y p-valores asociados.
#'   \item \code{metricas}: Una lista que incluye el Pseudo-R2 de McFadden, el Criterio de Informacion de Akaike (AIC) y la precision global (\emph{Accuracy}).
#'   \item \code{confusion}: Matriz de confusion que cruza los valores reales frente a los predichos segun el umbral especificado.
#'   \item \code{predicciones}: Un \code{data.frame} con los valores observados, probabilidades predichas, clasificacion final y el apalancamiento (\emph{leverage}).
#'   \item \code{diag_leverage}: Lista con informacion diagnostica sobre el apalancamiento y deteccion de puntos extremos.
#'   \item \code{notas}: Lista de mensajes interpretativos sobre la bondad de ajuste y la precision.
#'   \item \code{supuestos}: Un vector de caracteres con las validaciones tecnicas (semaforos) sobre significancia, multicolinealidad e influencia.
#' }
#'
#'
#' @examples
#' data(datos_roya_cafe)
#' modelo_log <- logistico(presencia_roya ~ humedad_relativa_pct + 
#'                         altitud_msnm + manejo_sombra, 
#'                         data = datos_roya_cafe)
#'                         
#' @importFrom stats glm binomial fitted coef confint formula na.omit predict model.frame logLik aggregate pchisq cooks.distance hatvalues rstandard
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
# Asegurate de usar esta version actualizada de tu funcion principal
logistico <- function(formula, data, umbral = 0.5) {

  # --- 1. Preparacion y Ajuste ---
  # Dejamos que glm maneje el punto (.) y los NAs automaticamente
  modelo_base <- stats::glm(formula, data = data, family = stats::binomial(link = "logit"), na.action = stats::na.omit)
  resumen <- summary(modelo_base)

  # Extraemos los datos EXACTOS que uso el modelo (sin NAs)
  data_model <- modelo_base$model

  # Identificar la variable respuesta (siempre es la 1ra columna del model frame)
  var_y_name <- names(data_model)[1]
  y_real <- modelo_base$y # vector 0/1 real usado en el ajuste

  # --- 2. Interpretacion de Coeficientes (Odds Ratios) ---
  coefs_log <- coef(modelo_base)
  coefs_or <- exp(coefs_log)
  cis_or <- exp(stats::confint(modelo_base)) # Intervalos de confianza para OR

  tabla_or <- cbind(OR = coefs_or, `2.5 %` = cis_or[,1], `97.5 %` = cis_or[,2], `P-valor` = resumen$coefficients[,4])

  # --- 3. Metricas de Bondad de Ajuste ---
  # Pseudo R2 de McFadden: 1 - (LogLik_Modelo / LogLik_Nulo)
  ll_modelo <- as.numeric(logLik(modelo_base))

  # Ajustamos modelo nulo sobre LOS MISMOS DATOS exactos (data_model)
  modelo_nulo <- stats::glm(as.formula(paste(var_y_name, "~ 1")), data = data_model, family = binomial)
  ll_nulo <- as.numeric(logLik(modelo_nulo))

  pseudo_r2 <- 1 - (ll_modelo / ll_nulo)

  # Matriz de Confusion y Accuracy
  probs_pred <- fitted(modelo_base)
  clase_pred <- ifelse(probs_pred > umbral, 1, 0)

  tabla_confusion <- table(Real = y_real, Predicho = clase_pred)
  accuracy <- sum(diag(tabla_confusion)) / sum(tabla_confusion)

  # --- 4. Notas Didacticas ---
  nota_r2 <- paste0("Pseudo R2 (McFadden) = ", round(pseudo_r2, 3),
                    ". (Valores entre 0.2 y 0.4 se consideran excelentes en logistica).")

  nota_acc <- paste0("Precision Global (Accuracy): ", round(accuracy * 100, 1),
                     "%. El modelo clasifico correctamente este porcentaje de casos.")

  # Mensaje sobre Odds Ratios significativos
  significativos <- rownames(tabla_or)[tabla_or[, "P-valor"] < 0.05 & rownames(tabla_or) != "(Intercept)"]
  if (length(significativos) > 0) {
    ejemplo_var <- significativos[1]
    or_val <- tabla_or[ejemplo_var, "OR"]
    direccion <- ifelse(or_val > 1, "AUMENTA", "DISMINUYE")
    nota_interp <- paste0("Ejemplo de interpretacion: La variable '", ejemplo_var,
                          "' es significativa. Su OR es ", round(or_val, 2),
                          ", lo que significa que ", direccion, " la probabilidad del evento.")
  } else {
    nota_interp <- "Ninguna variable parece tener un efecto estadisticamente significativo."
  }
  # --- 4.5 Semaforo de Validacion (Diagnostico Tecnico) ---
  # A) Significancia del Modelo (Likelihood Ratio Test)
  p_modelo <- 1 - pchisq(modelo_base$null.deviance - modelo_base$deviance,
                         df = modelo_base$df.null - modelo_base$df.residual)

  nota_signif <- if(p_modelo < 0.05) "OKEY: El modelo es significativamente mejor que el nulo." else " ALERTA: El modelo no aporta informacion relevante (p > 0.05)."

  # B) Multicolinealidad (VIF)
  vif_vals <- if(requireNamespace("car", quietly = TRUE)) car::vif(modelo_base) else NA
  nota_vif <- if(any(vif_vals > 10, na.rm = TRUE)) " ALERTA: Multicolinealidad alta (VIF > 10). Revisa las variables redundantes." else "OK: No hay multicolinealidad severa."

  # C) Influencia (Puntos Cook > Umbral)
  cooksd <- stats::cooks.distance(modelo_base)
  umbral_cook <- 4 / length(cooksd)
  hay_influencia <- any(cooksd > umbral_cook)
  nota_influencia <- if(hay_influencia) " ALERTA: Se detectaron casos influyentes (Distancia de Cook). Revisa los casos atipicos." else "OKEY: No hay puntos influyentes criticos."

  # --- NUEVO: D) Leverage (Apalancamiento) ---
  leverage <- stats::hatvalues(modelo_base) # Calculamos h
  num_vars <- length(coefs_log) # numero de parametros (incluyendo intercepto)
  n_obs <- length(y_real)
  umbral_leverage <- 2 * num_vars / n_obs # Umbral comun: 2*k/n
  hay_leverage_alto <- any(leverage > umbral_leverage)
  nota_leverage <- if(hay_leverage_alto) " ALERTA: Hay observaciones con Leverage alto (puntos atipicos en X). Usa plot() para ver el grafico." else "OKEY: No hay puntos con apalancamiento extremo."

  # --- 5. Construccion del Objeto Final ---
  # Consolidamos TODO en una sola lista
  resultado <- list(
    formula = formula, # <--- Guardamos la formula (necesario para el print seccion 1)
    modelo = modelo_base,
    variables = list(y = var_y_name),
    coeficientes = resumen$coefficients,
    odds_ratios = tabla_or,
    metricas = list(pseudo_r2 = pseudo_r2, aic = modelo_base$aic, accuracy = accuracy),
    confusion = tabla_confusion,
    predicciones = data.frame(
      obs = y_real,
      prob = probs_pred,
      clasif = clase_pred,
      leverage = leverage # <--- Guardamos leverage para el plot
    ),
    diag_leverage = list(umbral = umbral_leverage, hay_alto = hay_leverage_alto), # Info diagnostica
    notas = list(
      r2 = nota_r2,
      acc = nota_acc,
      interp = nota_interp
    ),
    # Aqui unificamos el semaforo que creaste arriba, incluyendo el nuevo de leverage
    supuestos = c(nota_signif, nota_vif, nota_influencia, nota_leverage)
  )

  class(resultado) <- "logistico_greenreg"
  return(resultado)
}



#' @rdname logistico
#' @param x Objeto de clase logistico_greenreg.
#' @param ... Argumentos adicionales.
#' @export
print.logistico_greenreg <- function(x, ...) {
  
  cat("\n==========================================================\n")
  cat("               REGRESION LOGISTICA \n")
  cat("==========================================================\n\n")
  
  # --- 1. MODELO Y ECUACION ---
  cat("--- 1. MODELO Y ECUACION ---\n")
  cat("Formula aplicada:", deparse(x$formula), "\n")
  
  # CORRECCION CRAN: Uso explicito de stats::
  b <- round(stats::coef(x$modelo), 3)
  vars <- names(b)
  ecuacion <- paste0("logit(p) = ", b[1])
  if(length(b) > 1) {
    for(i in 2:length(b)) {
      ecuacion <- paste0(ecuacion, " + (", b[i], " * ", vars[i], ")")
    }
  }
  cat("Ecuacion estimada:\n", ecuacion, "\n\n")
  
  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  cat("Estimaciones en escala Logit (log-odds):\n")
  stats::printCoefmat(x$coeficientes, digits = 4, signif.stars = TRUE)
  cat("\nInterpretacion en Odds Ratios (OR):\n")
  print(round(x$odds_ratios, 4))
  cat(" GUIA: OR > 1 Aumenta prob. | OR < 1 Disminuye prob.\n\n")
  
  # --- 3. EVALUACION DEL MODELO ---
  cat("--- 3. EVALUACION DEL MODELO ---\n")
  
  # CORRECCION CRAN: Uso explicito de stats::
  chi_cuad <- x$modelo$null.deviance - x$modelo$deviance
  df_chi <- x$modelo$df.null - x$modelo$df.residual
  
  if (!is.na(chi_cuad) && !is.na(df_chi) && df_chi > 0) {
    p_val_lrt <- 1 - stats::pchisq(chi_cuad, df_chi)
    p_val_txt <- if(is.na(p_val_lrt)) "NA" else sprintf("%.4f", p_val_lrt)
    resultado_txt <- if(!is.na(p_val_lrt) && p_val_lrt < 0.05) "Modelo Globalmente Significativo" else "Modelo No Significativo"
  } else {
    chi_cuad <- NA
    df_chi <- NA
    p_val_txt <- "NA"
    resultado_txt <- "No disponible"
  }
  
  cat("Prueba de Razon de Verosimilitud (LRT):\n")
  cat(sprintf(" Chi-cuadrado: %s | GL: %s | p-valor: %s\n", 
              if(is.na(chi_cuad)) "NA" else round(chi_cuad, 3),
              if(is.na(df_chi)) "NA" else df_chi, 
              p_val_txt))
  cat(" Resultado:", resultado_txt, "\n\n")
  
  # --- 4. BONDAD DE AJUSTE ---
  cat("--- 4. BONDAD DE AJUSTE ---\n")
  
  # CORRECCION CRAN: Uso explicito de stats::
  loglik_val <- try(stats::logLik(x$modelo), silent = TRUE)
  loglik_txt <- if(inherits(loglik_val, "try-error") || is.na(loglik_val)) "NA" else round(as.numeric(loglik_val), 3)
  
  cat(" Log-Verosimilitud (LogLik):", loglik_txt, "\n")
  cat(" AIC (Criterio de Akaike): ", round(x$metricas$aic, 3), "\n")
  cat("*", x$notas$r2, "\n")
  cat(" Prueba de Wald: Revisar asteriscos en la seccion 2 (Significancia individual).\n\n")
  
  # --- 5. MATRIZ DE CONFUSION ---
  cat("--- 5. MATRIZ DE CONFUSION (Real vs Predicho) ---\n")
  print(x$confusion)
  cat("\n INTERPRETACION:\n")
  cat("Buscamos una diagonal principal (Real y Predicho coinciden) con valores altos.\n")
  cat("Los valores fuera de la diagonal representan Falsos Positivos y Falsos Negativos.\n")
  cat("*", x$notas$acc, "\n\n")
  
  # --- 6. VERIFICACION DE SUPUESTOS ---
  cat("--- 6. VERIFICACION DE SUPUESTOS ---\n")
  cat(" Linealidad en el Logit: Supuesta (Usar plot(x) para verificar graficos de residuos).\n")
  cat(" Independencia: Observaciones independientes asumidas por diseno.\n")
  
  # Mostrar las alertas de forma limpia y segura
  if (!is.null(x$supuestos)) {
    for (i in 1:length(x$supuestos)) {
      cat(" - ", x$supuestos[i], "\n")
    }
  }
  cat("\n")
  
  # --- 7. RESUMEN DEL MODELO ---
  cat("--- 7. RESUMEN DEL MODELO ---\n")
  cat("Conclusion rapida:", x$notas$interp, "\n")
  
  hay_alerta <- FALSE
  if (!is.null(x$supuestos)) {
    hay_alerta <- any(grepl("ALERTA", x$supuestos))
  }
  
  if(hay_alerta) {
    cat("\n ALERTA: El modelo presenta advertencias tecnicas (Multicolinealidad, Puntos influyentes o Leverage).\n")
    cat(" Se recomienda usar plot(nombre_del_objeto) para ver los graficos de diagnostico.\n")
  }
  
  cat("\n==========================================================\n")
  
  return(invisible(x))
}
#' @rdname logistico
#' @description Genera una secuencia de 9 graficas para validar el modelo de regresion logistica.
#' @param x Objeto de clase logistico_greenreg.
#' @param ... Argumentos adicionales.
#' @export
plot.logistico_greenreg <- function(x, ...) {
  
  df_preds <- x$predicciones 
  modelo   <- x$modelo
  umbral   <- if(!is.null(x$umbral)) x$umbral else 0.5
  
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
  
  # CORRECCION 1: Usar message() en lugar de cat()
  if (interactive()) {
    message("Generando graficas de diagnostico... (Presione [Enter] para avanzar)\n")
  }
  
  # --- GRAFICA 1: SEPARACION DE CLASES (DENSIDAD) ---
  g1 <- ggplot2::ggplot(df_preds, ggplot2::aes(x = prob, fill = factor(obs))) +
    ggplot2::geom_density(alpha = 0.6, color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "#3498DB", "1" = "#E74C3C"), name = "Clase Real") +
    ggplot2::labs(
      title = "1. Capacidad de Separacion (Discriminacion)",
      subtitle = "Distribucion de las probabilidades predichas para cada grupo",
      x = "Probabilidad Predicha por el Modelo (0 a 1)",
      y = "Densidad",
      caption = "INTERPRETACION:\n BIEN: Dos montanas separadas. La azul a la izquierda (cerca de 0) y la roja a la derecha (cerca de 1).\n MAL: Si las montanas estan una encima de la otra (solapadas), el modelo esta confundido."
    ) + mi_tema
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1/9 > ")}
  
  # --- GRAFICA 2: MATRIZ DE CONFUSION VISUAL ---
  pred_clase <- ifelse(df_preds$prob > umbral, 1, 0)
  tabla <- table(Real = df_preds$obs, Predicho = pred_clase)
  df_matriz <- as.data.frame(tabla)
  
  g2 <- ggplot2::ggplot(df_matriz, ggplot2::aes(x = Predicho, y = Real, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 8, color = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#95A5A6", high = "#2C3E50") +
    ggplot2::labs(
      title = paste("2. Matriz de Confusion (Umbral:", umbral, ")"),
      subtitle = "Conteo de Aciertos y Errores",
      caption = "INTERPRETACION:\n DIAGONAL (0-0 y 1-1): Son los Aciertos correctos.\n CRUZADA (0-1 y 1-0): Son los Errores (Falsos Positivos y Falsos Negativos).\nObjetivo: Que los cuadros oscuros esten en la diagonal."
    ) + mi_tema
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2/9 > ")}
  
  # --- GRAFICA 3: ODDS RATIOS (FOREST PLOT MEJORADO) ---
  if (requireNamespace("broom", quietly = TRUE)) {
    coefs <- broom::tidy(modelo, conf.int = TRUE, exponentiate = TRUE)
    coefs <- coefs[coefs$term != "(Intercept)", ]
    
    coefs$Tipo <- "No Significativo"
    coefs$Tipo[coefs$p.value < 0.05 & coefs$estimate > 1] <- "Riesgo (Aumenta Prob)"
    coefs$Tipo[coefs$p.value < 0.05 & coefs$estimate < 1] <- "Protector (Baja Prob)"
    
    g3 <- ggplot2::ggplot(coefs, ggplot2::aes(x = estimate, y = reorder(term, estimate), color = Tipo)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = conf.low, xmax = conf.high), size = 0.8) +
      ggplot2::scale_color_manual(values = c("Riesgo (Aumenta Prob)" = "#E74C3C",
                                             "Protector (Baja Prob)" = "#27AE60",
                                             "No Significativo" = "#95A5A6")) +
      ggplot2::labs(
        title = "3. Importancia de Variables (Odds Ratios)",
        subtitle = "Impacto de cada variable en la probabilidad del evento",
        x = "Odds Ratio (Escala Logaritmica)", y = "Variables",
        caption = "INTERPRETACION:\n  Derecha del 1: AUMENTA la probabilidad del evento (Riesgo).\n Izquierda del 1: DISMINUYE la probabilidad (Proteccion).\n GRIS (Cruza el 1): La variable NO sirve, estadisticamente no afecta."
      ) + mi_tema
    print(g3)
  } else {
    # CORRECCION 1: Usar message()
    message("Instala el paquete 'broom' para ver la grafica de Odds Ratios: install.packages('broom')")
  }
  if (interactive()) {readline(prompt = "Grafica 3/9 > ")}
  
  # --- GRAFICA 4: CURVA DE CALIBRACION ---
  df_preds$bin <- cut(df_preds$prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  
  calibracion <- aggregate(obs ~ bin, data = df_preds, mean, na.rm=TRUE)
  prob_medias <- aggregate(prob ~ bin, data = df_preds, mean, na.rm=TRUE)
  calibracion$prob_media <- prob_medias$prob
  
  g4 <- ggplot2::ggplot(calibracion, ggplot2::aes(x = prob_media, y = obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_line(color = "#8E44AD", linewidth = 1) +
    ggplot2::geom_point(color = "#8E44AD", size = 3) +
    ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(
      title = "4. Calibracion del Modelo",
      subtitle = "Dice la verdad el modelo sobre sus probabilidades?",
      x = "Probabilidad que dijo el modelo",
      y = "Porcentaje Real de Eventos",
      caption = "INTERPRETACION:\n BIEN: Los puntos morados siguen la linea gris punteada.\n MAL: Si la curva esta muy por debajo de la linea, el modelo es 'exagerado'."
    ) + mi_tema
  print(g4)
  if (interactive()) {readline(prompt = "Grafica 4/9 > ")}
  
  # --- GRAFICA 5: CURVA ROC ---
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(df_preds$obs, df_preds$prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    roc_df <- data.frame(esp = roc_obj$specificities, sens = roc_obj$sensitivities)
    
    g5 <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - esp, y = sens)) +
      ggplot2::geom_line(color = "#D35400", linewidth = 1.2) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      ggplot2::labs(
        title = paste0("5. Curva ROC (AUC = ", round(auc_val, 2), ")"),
        subtitle = "Capacidad global de diagnostico",
        x = "Tasa de Falsos Positivos", y = "Sensibilidad (Tasa de Aciertos)",
        caption = "INTERPRETACION:\n BIEN: La curva naranja se acerca a la esquina superior izquierda.\n MAL: Si la curva esta cerca de la diagonal gris, el modelo es puro azar."
      ) + mi_tema
    print(g5)
  } else {
    # CORRECCION 1: Usar message()
    message("Instala el paquete 'pROC' para ver la Curva ROC: install.packages('pROC')")
  }
  if (interactive()) {readline(prompt = "Grafica 5/9 > ")}
  
  # --- GRAFICA 6: VALORES INFLUYENTES (COOK) ---
  cooksd <- stats::cooks.distance(modelo)
  n_datos <- length(cooksd)
  umbral_cook <- 4 / n_datos
  df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)
  df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")
  
  g6 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
    ggplot2::geom_bar(stat = "identity", fill = "#34495E", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "6. Datos Influyentes (Distancia de Cook)",
      subtitle = "Deteccion de casos que distorsionan el modelo",
      x = "Numero de Caso", y = "Influencia",
      caption = "INTERPRETACION:\n ALERTA: Las barras altas que cruzan la linea roja son casos atipicos.\nRevisalos en tu base de datos."
    ) + mi_tema
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6/9 > ")}
  
  # --- GRAFICA 7: BARRIDO DE UMBRALES (Trade-off) ---
  thresholds <- seq(0, 1, length = 100)
  sens <- sapply(thresholds, function(t) mean(df_preds$prob[df_preds$obs == 1] >= t, na.rm=TRUE))
  spec <- sapply(thresholds, function(t) mean(df_preds$prob[df_preds$obs == 0] < t, na.rm=TRUE))
  df_sweep <- data.frame(thresholds, sens, spec)
  
  g7 <- ggplot2::ggplot(df_sweep) +
    ggplot2::geom_line(ggplot2::aes(x = thresholds, y = sens, color = "Sensibilidad"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(x = thresholds, y = spec, color = "Especificidad"), linewidth = 1) +
    ggplot2::geom_vline(xintercept = umbral, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("Sensibilidad" = "#E74C3C", "Especificidad" = "#2980B9")) +
    ggplot2::labs(
      title = "7. Trade-off: Umbral de Clasificacion",
      subtitle = "Como el umbral afecta los errores",
      x = "Umbral", y = "Tasa", color = "Metrica",
      caption = "INTERPRETACION:\nEl punto de cruce es el umbral 'balanceado'.\nSi necesitas detectar todos los casos, mueve el umbral a la izquierda (mas sensibilidad)."
    ) + mi_tema
  print(g7)
  
  if (interactive()) {readline(prompt = "Grafica 7/9 > ")}
  
  # --- GRAFICA 8: RESIDUOS DE PEARSON ---
  df_preds$pearson_res <- (df_preds$obs - df_preds$prob) / sqrt(df_preds$prob * (1 - df_preds$prob))
  
  g8 <- ggplot2::ggplot(df_preds, ggplot2::aes(x = prob, y = pearson_res)) +
    ggplot2::geom_point(alpha = 0.5, color = "#2C3E50") +
    ggplot2::geom_smooth(method = "loess", color = "#E74C3C", se = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "8. Analisis de Residuos (Pearson)",
      subtitle = "Existen patrones ocultos en el error?",
      x = "Probabilidad Predicha", y = "Residuo de Pearson",
      caption = "INTERPRETACION:\nSi la linea roja se aleja del 0, el modelo tiene problemas estructurales."
    ) + mi_tema
  print(g8)
  if (interactive()) {readline(prompt = "Grafica 8/9 > ")}
  
  # --- GRAFICA 9: RESIDUOS VS APALANCAMIENTO ---
  lev <- stats::hatvalues(x$modelo)
  res_std <- stats::rstandard(x$modelo) 
  n_datos <- length(lev)
  
  # CORRECCION 2: stats::coef
  p_vars <- length(stats::coef(x$modelo)) 
  umbral_lev <- (2 * p_vars) / n_datos
  
  df_diag <- data.frame(
    index = 1:n_datos,
    leverage = lev,
    resid_std = res_std
  )
  
  df_diag$Label <- ifelse(df_diag$leverage > umbral_lev | abs(df_diag$resid_std) > 2,
                          as.character(df_diag$index), "")
  
  g9 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = leverage, y = resid_std)) +
    ggplot2::geom_point(ggplot2::aes(size = leverage), color = "#C0392B", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5)
  
  # CORRECCION 4: Proteccion defensiva para ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g9 <- g9 + ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold")
  } else {
    g9 <- g9 + ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.5, vjust = -0.5)
  }
  
  g9 <- g9 +
    ggplot2::labs(
      title = "9. Diagnostico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificando puntos que 'jalonean' el ajuste logistico",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("INTERPRETACION:\n",
                       " EJE X: Puntos a la derecha de la linea azul tienen valores de X muy atipicos.\n",
                       " EJE Y: Puntos fuera de las lineas rojas punteadas son Outliers (errores en clasificacion).\n",
                       " CLAVE: Puntos arriba/abajo a la derecha estan deformando las probabilidades del modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")
  
  print(g9)
  
  # CORRECCION 1: Usar message()
  if (interactive()) {
    message("Secuencia de graficas completada. Se ha generado la grafica 9 de influencia.\n")
  }
  
  # CORRECCION 3: Retornar el objeto original de forma invisible
  return(invisible(x))
}