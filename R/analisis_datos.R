#' Asistente Interactivo de Análisis de Datos
#'
#'
#' Esta es la función principal y punto de entrada del paquete. Actua como un consultor
#' automatizado que guia al usuario desde la carga de datos hasta la generacion del código final,
#' seleccionando la metodologia adecuada segun la naturaleza de la variable objetivo.
#'
#' @section Como funciona este Asistente paso a paso?:
#' \enumerate{
#'   \item **Escaneo y Validación:**
#'         Analiza la base de datos para asegurar que la variable objetivo existe.
#'         Ademas, realiza un test interno de autocorrelación para verificar si tiene sentido aplicar modelos de Series de Tiempo (evitando análisis espurios).
#'   \item **Menu Interactivo:**
#'         Despliega opciones en la consola basandose en los datos detectados, preguntando al usuario que enfoque desea (Regresión Clasica, Modelos Generalizados o Series de Tiempo).
#'   \item **Generacion de Estrategias (El "Motor"):**
#'         \itemize{
#'            \item \emph{Para Regresión y GLM:} Crea dos versiones del código.
#'            \enumerate{
#'               \item **Modelo Full:** Incluye todas las variables disponibles.
#'               \item **Modelo Optimizado:** Ejecuta internamente un algoritmo "Stepwise" basado en AIC para limpiar el modelo y quedarse solo con los predictores relevantes.
#'            }
#'            \item \emph{Para Series de Tiempo:} Ejecuta una "Busqueda en Rejilla" (Grid Search) probando combinaciones de parámetros (p, d, q). Identifica automaticamente si la serie necesita diferenciación (ARIMA), si es puramente autoregresiva (AR), de media móvil (MA) o mixta (ARMA), y genera el código para la función especifica.
#'         }
#' }
#'
#' @section Guia de la Gráfica Exploratoria Generada:
#' Dependiendo de la opción elegida, la función genera una visualización preliminar diferente para entender los datos antes de modelar:
#' \describe{
#'   \item{\strong{Escenario 1: Regresión Lineal (Histograma + Gauss)}}{
#'     Muestra un histograma azul de tus datos superpuesto con una curva roja (Campana de Gauss teorica).
#'     \emph{Objetivo:} Evaluar visualmente si la variable respuesta se parece a una distribución Normal (supuesto clave de la regresión lineal).
#'   }
#'   \item{\strong{Escenario 2: GLM (Barras o Histograma)}}{
#'     Detecta el tipo de variable:
#'     \itemize{
#'       \item *Si es Binaria (0/1):* Gráfico de barras naranja. Ayuda a ver si las clases están desbalanceadas (ej. muchos "0" y pocos "1").
#'       \item *Si es Conteo:* Histograma verde. Muestra la frecuencia de eventos (típico de Poisson).
#'     }
#'   }
#'   \item{\strong{Escenario 3: Serie de Tiempo (Línea de Tendencia)}}{
#'     Gráfica la evolución de la variable a lo largo del índice temporal.
#'     \emph{Objetivo:} Identificar tendencias (subidas/bajadas a largo plazo) o ciclos repetitivos antes de ajustar el modelo.
#'   }
#' }
#'
#' @param data Data frame o base de datos que contiene las variables.
#' @param variable_objetivo Cadena de texto con el nombre exacto de la columna que deseas analizar (ej. "Ventas").
#' @param verbose Lógico. Si es TRUE, despliega los menus. Por defecto es TRUE.
#'
#' @return Devuelve de forma invisible un objeto de clase \code{ggplot2}.
#' La función se utiliza principalmente por su efecto secundario de imprimir
#' código sugerido en la consola y generar visualizaciones.
#'
#' @examples
#'   data(datos_rendimiento_maiz)
#'   greenreg_analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")
#'
#' @importFrom stats lm glm step cor shapiro.test arima gaussian binomial poisson acf
#' @import ggplot2
#' @export
greenreg_analisis_datos <- function(data, variable_objetivo, verbose = TRUE) {

  # --- 1. Validaciones e Informacion de la Base ---

  # Verificamos si la columna objetivo realmente existe en el data frame proporcionado.
  if (!(variable_objetivo %in% names(data))) stop("Error: La variable objetivo no existe en la base de datos.")
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("El argumento 'data' debe ser un data.frame con al menos una fila.")
  }
  if (ncol(data) < 2) {
    stop("'data' necesita al menos una variable predictora adem\u00e1s de 'variable_objetivo'.")
  }

  # Extraccion de la variable respuesta (y) y lista de variables predictoras (vars_x).
  y <- data[[variable_objetivo]]
  vars_x <- setdiff(names(data), variable_objetivo)
  n <- nrow(data)

  # Obtenemos el nombre del objeto 'data' tal como lo escribio el usuario (ej. "mi_base")
  # para usarlo textualmente en la generacion automatica del codigo sugerido.
  nombre_data <- deparse(substitute(data))
  # --- 2A. Deteccion Automatica y Recomendacion ---

  # --- CORRECCIoN: Definimos p_externas aqui ---
  p_externas <- length(vars_x)
  # ---------------------------------------------

  # --- 2. Deteccion Automatica y Menu Interactivo ---
  col_fecha <- NULL

  # 1. Buscar si ya existe columna fecha nativa, o si la variable objetivo (u
  # otra columna) ya es una serie de tiempo (is.ts), sin exigir columna de fecha.
  cols_nativas <- names(data)[sapply(data, inherits, what = c("Date", "POSIXct", "POSIXt"))]
  es_ts_nativo <- stats::is.ts(data[[variable_objetivo]]) || any(sapply(data, stats::is.ts))

  if (length(cols_nativas) > 0) {
    col_fecha <- cols_nativas[1]
  } else if (!es_ts_nativo) {
    # 2. Intentar convertir columnas de texto a fecha
    formatos_prueba <- c("%d/%m/%Y", "%Y-%m-%d", "%m/%d/%Y", "%d-%m-%Y", "%Y/%m/%d")
    cols_texto <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]

    for (col in cols_texto) {
      vals <- as.character(data[[col]])
      exito <- FALSE
      for (fmt in formatos_prueba) {
        fechas_test <- try(as.Date(vals, format = fmt), silent = TRUE)
        # Criterio: Si no da error y tiene pocos valores nulos (<20%), es fecha valida
        if (!inherits(fechas_test, "try-error") && sum(is.na(fechas_test)) < (0.2 * length(vals))) {
          data[[col]] <- fechas_test # Auto-conversion!
          col_fecha <- col
          exito <- TRUE
          message(sprintf(" Auto-correcci\u00f3n: Columna '%s' detectada como FECHA (formato '%s').\n", col, fmt))
          break
        }
      }
      if (exito) break
    }
  }

  tiene_fecha <- !is.null(col_fecha) || es_ts_nativo
  hay_autocorr <- FALSE

  # 3. Verificar Autocorrelacion (solo si hay fecha o es serie ts)
  if (tiene_fecha) {
    if (!is.null(col_fecha)) {
      data <- data[order(data[[col_fecha]]), ] # Ordenar cronologicamente
    }
    y_ord <- data[[variable_objetivo]]

    # Reasignamos 'y' con la serie ya ordenada cronologicamente.
    # Sin esto, el grid search de ARIMA/AR/MA/ARMA (Opcion 3) y el codigo
    # sugerido al usuario terminaban usando la serie en su orden ORIGINAL
    # (potencialmente desordenada), invalidando la deteccion de autocorrelacion
    # y el ajuste de los modelos de series de tiempo.
    y <- y_ord

    # Ljung-Box global en vez de "algun rezago individual supera 2/sqrt(n)":
    # con ~20 rezagos ese criterio da falsos positivos frecuentes en ruido blanco.
    lb_global <- tryCatch(
      stats::Box.test(y_ord, type = "Ljung-Box", lag = min(10, round(length(y_ord) / 5))),
      error = function(e) list(p.value = NA)
    )
    hay_autocorr <- !is.na(lb_global$p.value) && lb_global$p.value < 0.05
  }
  es_ts <- tiene_fecha && hay_autocorr

  # --- 2B. Logica de Recomendacion ---
  y_temp <- data[[variable_objetivo]]

  es_binaria <- all(y_temp %in% c(0, 1), na.rm=TRUE) || (is.factor(y_temp) && length(levels(y_temp))==2)
  es_conteo <- .greenreg_es_conteo(y_temp, es_binaria)

  opcion_sugerida <- "1"
  razon_sugerencia <- "Variable continua (Regresi\u00f3n Lineal)."

  if (es_ts) {
    opcion_sugerida <- "3"
    razon_sugerencia <- if (!is.null(col_fecha)) {
      sprintf("Serie de Tiempo (Fecha detectada: '%s' + Autocorrelaci\u00f3n).", col_fecha)
    } else {
      "Serie de Tiempo (clase ts detectada + Autocorrelaci\u00f3n)."
    }
  } else if (es_binaria) {
    opcion_sugerida <- "2"
    razon_sugerencia <- "Variable Binaria (Log\u00edstica)."
  } else if (es_conteo) {
    opcion_sugerida <- "2"
    razon_sugerencia <- "Variable de Conteo (Poisson)."
  }

  # --- 2C. Menu Visual ---
  # Se controla unicamente con 'verbose': un usuario que ejecute la funcion
  # con verbose = TRUE desde un script (Rscript, .Rmd, etc.) debe ver la
  # salida igual que en una sesion interactiva. 'interactive()' solo se usa
  # mas abajo para decidir si se puede pedir input con readline().
  if (verbose) {
    cat("\nAsistente de An\u00e1lisis de Datos\n")
    cat("===========================================\n")
    cat(sprintf("   Base de datos: %s\n", nombre_data))
    cat(sprintf("   Objetivo: '%s'\n", variable_objetivo))
    cat(sprintf("   Observaciones: %d | Predictoras: %d\n\n", n, p_externas))

    cat("RECOMENDACI\u00d3N AUTOM\u00c1TICA:\n")
    cat(sprintf("   -> Sugerencia: Opci\u00f3n %s\n", opcion_sugerida))
    cat(sprintf("   -> Motivo: %s\n\n", razon_sugerencia))

    cat("Opciones Disponibles:\n")

    # Marcas visuales para guiar al usuario
    mark1 <- if(opcion_sugerida == "1") "  * (Recomendado)" else ""
    mark2 <- if(opcion_sugerida == "2") "  * (Recomendado)" else ""
    mark3 <- if(opcion_sugerida == "3") "  * (Recomendado)" else ""

    cat(sprintf("   1. Regresi\u00f3n Lineal%s\n", mark1))
    cat(sprintf("   2. Modelos GLM (Log\u00edstica/Poisson)%s\n", mark2))
    cat(sprintf("   3. Series de Tiempo%s\n", mark3))
  }

  # 1. Manejo de la seleccion
  if (interactive() && verbose) {
    seleccion <- readline(prompt = "Selecciona una opci\u00f3n (1-3) [Enter para auto]: ")
  } else {
    # En sesiones no interactivas (como los tests de CRAN),
    # forzamos la sugerencia sin preguntar.
    seleccion <- ""
  }

  # 2. Aplicar la logica de seleccion
  if (seleccion == "") {
    seleccion <- opcion_sugerida

    # Usamos message() en lugar de cat() para que sea suprimible.
    # Se muestra siempre que verbose = TRUE (no solo en sesiones interactivas).
    if (verbose) {
      message(sprintf("-> Seleccion autom\u00e1tica aplicada: %s", seleccion))
    }
  }

  # --- 3. Procesamiento segun Seleccion ---
  # Inicializacion de variables para almacenar el codigo generado y explicaciones.
  codigo_full <- ""
  codigo_opt <- ""
  explicacion <- ""
  p_graph <- NULL

  # ==============================================================================
  # OPCIoN 1: REGRESIoN LINEAL
  # Estrategia: Ajustar minimos cuadrados ordinarios (RLS/RLM).
  # ==============================================================================
  if (seleccion == "1") {
    message("\n Analizando para Regresi\u00f3n Lineal...\n")

    # A) Codigo Modelo Completo
    # Genera el string del codigo para un modelo que incluye TODAS las variables (formula = y ~ .).
    codigo_full <- paste0(
      sprintf("modelo_rlm_full <- greenreg_rlm(%s ~ ., data = %s)\n", variable_objetivo, nombre_data),
      "modelo_rlm_full\n",
      "plot(modelo_rlm_full)"
    )

    # B) Codigo Modelo Optimizado (Stepwise AIC)
    # Ejecuta internamente un 'stepwise selection' para encontrar que variables
    # realmente aportan valor al modelo y descartar las redundantes.
    formula_opt <- "."
    vars_seleccionadas <- "Todas"

    if (length(vars_x) > 0) {
      # Ajuste temporal para evaluacion (no se muestra al usuario).
      mod_temp <- stats::lm(as.formula(paste(variable_objetivo, "~ .")), data = data)
      # Seleccion de variables basada en AIC (trace=0 silencia la salida tecnica).
      mod_step <- stats::step(mod_temp, trace = 0)
      vars_opt <- names(coef(mod_step))[-1] # Extraemos nombres de coeficientes (excluyendo intercepto).

      if (length(vars_opt) > 0) {
        formula_opt <- paste(vars_opt, collapse = " + ")
        vars_seleccionadas <- paste(vars_opt, collapse = ", ")
      }
      if (length(vars_opt) == length(vars_x)) {
        explicacion <- "El proceso stepwise (AIC) conserv\u00f3 todas las variables; no se elimin\u00f3 ninguna."
      } else if (length(vars_opt) > 0) {
        explicacion <- sprintf("El proceso stepwise (AIC) elimin\u00f3 %d variable(s); quedaron: %s.",
                                length(vars_x) - length(vars_opt), vars_seleccionadas)
      }
    }

    # Genera el string del codigo para el modelo optimizado.
    codigo_opt <- paste0(
      sprintf("modelo_rlm_opt <- greenreg_rlm(%s ~ %s, data = %s)\n", variable_objetivo, formula_opt, nombre_data),
      "modelo_rlm_opt\n",
      "plot(modelo_rlm_opt)"
    )

    # C) Grafica Exploratoria con NOTAS
    # Crea un histograma de la variable objetivo y superpone una curva normal teorica
    # para ensenar al usuario sobre el supuesto de normalidad visualmente.
    media_y <- mean(y, na.rm=TRUE); sd_y <- sd(y, na.rm=TRUE)
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable_objetivo]])) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, fill = "#3498DB", color = "white", alpha = 0.7) +
      ggplot2::stat_function(fun = dnorm, args = list(mean = media_y, sd = sd_y), color = "red", linewidth = 1.2) +
      ggplot2::labs(
        title = paste("Distribuci\u00f3n de", variable_objetivo),
        subtitle = "Histograma (Datos) vs Campana de Gauss (Te\u00f3rico)",
        x = variable_objetivo, y = "Densidad",
        caption = "NOTA:\n La l\u00ednea ROJA representa la normalidad perfecta.\n Si las barras azules se alejan mucho de la forma de campana roja,\n  los supuestos de normalidad podr\u00edan no cumplirse (considera transformar variables)."
      ) + ggplot2::theme_minimal()
  }

  # ==============================================================================
  # OPCIoN 2: GLM (Logistica / Poisson)
  # Estrategia: Detectar automaticamente si es binaria (0/1) o conteo (enteros positivos).
  # ==============================================================================
  else if (seleccion == "2") {
    message("\n Analizando para GLM...\n")

    # Deteccion heuristica del tipo de variable respuesta
    es_binaria <- all(y %in% c(0, 1), na.rm=TRUE) || (is.factor(y) && length(levels(y))==2)
    es_conteo <- .greenreg_es_conteo(y, es_binaria)

    # Configuracion de parametros segun el tipo detectado
    if (es_binaria) {
      fun_glm <- "greenreg_logistico"
      nombre_mod <- "modelo_log"
      familia <- "binomial"

      if (verbose) {
        cat("   -> Detectada Variable Binaria. Sugerencia: Regresi\u00f3n Log\u00edstica.\n")
      }

      nota_grafica <- paste("NOTA:\n Revisa el balance entre las barras (0 y 1).\n",
                            "Si una barra es diminuta comparada con la otra (Desbalance),\n",
                            "el modelo tendra dificultades para predecir la clase minoritaria.")

    } else if (es_conteo) {
      fun_glm <- "greenreg_reg_poisson"
      nombre_mod <- "modelo_pois"
      familia <- "poisson"

      if (verbose) {
        cat("   -> Detectada Variable de Conteo. Sugerencia: Regresi\u00f3n Poisson.\n")
      }

      nota_grafica <- paste("NOTA:\n En Poisson, es comun ver barras altas a la izquierda (valores bajos)\n",
                            "que descienden hacia la derecha. Esto es normal y esperado.")

    } else {
      # Fallback a gaussiana si no es clara la deteccion
      fun_glm <- "glm_gen"
      nombre_mod <- "modelo_glm"
      familia <- "gaussian"

      if (verbose) {
        cat("Tipo ambiguo. Usando configuraci\u00f3n gen\u00e9rica.\n")
      }

      nota_grafica <- "NOTA: Distribuci\u00f3n de la variable respuesta."
    }
    # A) Codigo Modelo Completo
    codigo_full <- paste0(
      sprintf("%s_full <- %s(%s ~ ., data = %s)\n", nombre_mod, fun_glm, variable_objetivo, nombre_data),
      sprintf("%s_full\n", nombre_mod),
      sprintf("plot(%s_full)", nombre_mod)
    )

    # B) Codigo Modelo Optimizado
    # Intenta realizar una seleccion de variables stepwise para GLM.
    formula_opt <- "."
    if (length(vars_x) > 0) {
      try({
        mod_temp <- stats::glm(as.formula(paste(variable_objetivo, "~.")), data = data, family = familia)
        mod_step <- stats::step(mod_temp, trace = 0)
        vars_opt <- names(coef(mod_step))[-1]
        if (length(vars_opt) > 0) {
          formula_opt <- paste(vars_opt, collapse = " + ")
        }
        if (length(vars_opt) == length(vars_x)) {
          explicacion <- "El proceso stepwise (AIC) conserv\u00f3 todas las variables; no se elimin\u00f3 ninguna."
        } else if (length(vars_opt) > 0) {
          explicacion <- sprintf("El proceso stepwise (AIC) elimin\u00f3 %d variable(s); quedaron: %s.",
                                  length(vars_x) - length(vars_opt), paste(vars_opt, collapse = ", "))
        }
      }, silent = TRUE)
    }

    codigo_opt <- paste0(
      sprintf("%s_opt <- %s(%s ~ %s, data = %s)\n", nombre_mod, fun_glm, variable_objetivo, formula_opt, nombre_data),
      sprintf("%s_opt\n", nombre_mod),
      sprintf("plot(%s_opt)", nombre_mod)
    )

    # C) Grafica
    # Grafico de barras para binaria o histograma para conteo/continuo
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable_objetivo]])) +
      (if(es_binaria) ggplot2::geom_bar(fill = "orange", alpha=0.8) else ggplot2::geom_histogram(fill="#2ECC71", bins=20, color="white")) +
      ggplot2::labs(
        title = "Distribuci\u00f3n de la Variable Respuesta",
        caption = nota_grafica
      ) + ggplot2::theme_minimal()
  }

  # ==============================================================================
  # OPCIoN 3: SERIES DE TIEMPO (Sin Alertas)
  # Estrategia: Grid Search (fuerza bruta controlada) para hallar (p, d, q) optimos por AIC.
  # ==============================================================================
  else if (seleccion == "3") {
    message(" An\u00e1lisis de Serie de Tiempo (AR / MA / ARMA / ARIMA)...")

    if (p_externas > 0) {
      message(" Los modelos AR/MA/ARMA/ARIMA son univariados.")
    }


    mejor_aic <- Inf
    mejor_orden <- c(0, 0, 0) # Orden inicial (p, d, q)

    # Grid Search con silenciador de warnings (suppressWarnings)
    # Se iteran combinaciones de p (AR), d (Diferencias) y q (MA).
    for (d_test in 0:1) {
      for (p_test in 0:3) {
        for (q_test in 0:3) {
          if (p_test==0 && q_test==0 && d_test==0) next # Saltar modelo nulo
          try({
            # Ajuste de prueba. suppressWarnings() es clave aqui para evitar que el usuario
            # se asuste con alertas de "posible no convergencia" durante la busqueda automatica.
            mod <- suppressWarnings(stats::arima(y, order = c(p_test, d_test, q_test), method="CSS-ML"))

            # Guardamos el modelo si tiene el menor AIC encontrado hasta ahora
            if (!is.null(mod$aic) && mod$aic < mejor_aic) {
              mejor_aic <- mod$aic
              mejor_orden <- c(p_test, d_test, q_test)
            }
          }, silent = TRUE)
        }
      }
    }

    p_opt <- mejor_orden[1]; d_opt <- mejor_orden[2]; q_opt <- mejor_orden[3]

    # Seleccion de Funcion de la paqueteria segun el orden optimo encontrado
    if (d_opt > 0) {
      # Si d > 0, requiere diferenciacion -> ARIMA
      funcion_chapi <- "greenreg_modelo_arima"; nombre_modelo <- "modelo_ARIMAchapi"
      args <- sprintf("p = %d, d = %d, q = %d", p_opt, d_opt, q_opt)
    } else {
      # Si d = 0, es estacionaria (teoricamente). Verificamos componentes AR y MA.
      if (p_opt > 0 && q_opt == 0) {
        funcion_chapi <- "greenreg_modelo_ar"; nombre_modelo <- "modelo_AR"; args <- sprintf("p = %d", p_opt)
      } else if (p_opt == 0 && q_opt > 0) {
        funcion_chapi <- "greenreg_modelo_ma"; nombre_modelo <- "modelo_MA"; args <- sprintf("q = %d, include_mean = TRUE", q_opt)
      } else {
        funcion_chapi <- "greenreg_modelo_arma"; nombre_modelo <- "modelo_ARMA"; args <- sprintf("p = %d, q = %d, include_mean = TRUE", p_opt, q_opt)
      }
    }

    explicacion <- sprintf("AIC \u00f3ptimo encontrado con orden (%d, %d, %d).", p_opt, d_opt, q_opt)

    # Codigo unico: En series de tiempo univariadas no seleccionamos variables,
    # sino que seleccionamos el orden del modelo, por eso solo hay un bloque de codigo.
    codigo_salida <- paste0(
      sprintf("ts_data <- %s$%s\n", nombre_data, variable_objetivo),
      sprintf("%s <- %s(ts_data, %s)\n", nombre_modelo, funcion_chapi, args),
      sprintf("%s\n", nombre_modelo),
      sprintf("plot(%s)", nombre_modelo)
    )

    # Grafica de Serie de Tiempo con Notas
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = 1:n, y = .data[[variable_objetivo]])) +
      ggplot2::geom_line(color = "#8E44AD") +
      ggplot2::labs(
        title = paste("Evoluci\u00f3n Temporal:", variable_objetivo),
        x = "Tiempo / \u00edndice", y = "Valor",
        caption = "NOTA:\n Observa la gr\u00e1fica: Los datos suben o bajan con el tiempo (Tendencia)?\n Se repiten patrones regularmente (Ciclos)?\n El modelo sugerido abajo ha sido seleccionado para manejar estas caracteristicas."
      ) + ggplot2::theme_minimal()

    codigo_full <- codigo_salida # Para unificar la salida final
    codigo_opt <- "" # No aplica segunda opcion en TS univariado
  } else {
    stop("Opci\u00f3n invalida.")
  }

  # --- 4. Salida Final ---
  if (verbose) {
    cat("\n AN\u00c1LISIS COMPLETADO.\n")

    if (explicacion != "") {
      cat(paste("NOTA:", explicacion, "\n"))
    }
    cat("--------------------------------------------------\n")

    # Imprimir Opcion 1 (Codigo sugerido para Series de Tiempo)
    if (seleccion == "3") {
      cat(" C\u00f3digo Sugerido (Mejor Ajuste AIC):\n")
      cat("--------------------------------------------------\n")
      cat(codigo_full, "\n")
      cat("--------------------------------------------------\n")
    } else {
      # Opciones para Regresion y GLM
      cat(" OPCI\u00d3N A: Modelo con TODAS las variables\n")
      cat("--------------------------------------------------\n")
      cat(codigo_full, "\n")
      cat("--------------------------------------------------\n\n")

      cat(" OPCI\u00d3N B: Modelo OPTIMIZADO (Variables Seleccionadas)\n")
      cat("   (Elimina variables que aportan ruido o son redundantes)\n")
      cat("--------------------------------------------------\n")
      cat(codigo_opt, "\n")
      cat("--------------------------------------------------\n")
    }
    if (!is.null(p_graph)) print(p_graph)
  }

  return(invisible(p_graph))
}

# Heuristica interna para detectar variables de conteo (Poisson). No
# exportada (por eso usa comentarios normales, no roxygen): unifica el
# criterio que antes estaba duplicado (y desalineado) en dos puntos de
# greenreg_analisis_datos() [hallazgo C-11]. Exige enteros no negativos, que
# no sea binaria, y ademas pocos valores distintos relativos a n o una media
# baja; sin esto ultimo cualquier variable continua medida en enteros (edad,
# altura en cm, ingresos) se clasificaria como conteo.
.greenreg_es_conteo <- function(y, es_binaria) {
  if (!is.numeric(y)) return(FALSE)
  y_ok <- y[!is.na(y)]
  if (length(y_ok) == 0) return(FALSE)
  es_entero_no_neg <- all(y_ok == floor(y_ok)) && all(y_ok >= 0)
  if (!es_entero_no_neg || isTRUE(es_binaria)) return(FALSE)
  prop_distintos <- length(unique(y_ok)) / length(y_ok)
  (mean(y_ok) < 30) || (prop_distintos <= 0.3)
}
