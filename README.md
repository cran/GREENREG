# GREENREG

<!-- badges: start -->
<!-- badges: end -->

**GREENREG** ofrece un conjunto de funciones accesibles y automatizadas para
aplicar modelos estadisticos comunes en ciencias agropecuarias y ambientales:
Regresion Lineal Simple (RLS), Regresion Lineal Multiple (RLM), Modelos
Lineales Generalizados (logistica y Poisson), y modelos de series de tiempo
(AR, MA, ARMA, ARIMA). Ademas de ajustar el modelo, cada funcion valida
automaticamente los supuestos estadisticos, traduce los resultados a un
veredicto interpretable, y genera graficas de diagnostico listas para usarse.

Desarrollado como parte de una tesis de licenciatura en la Universidad
Autonoma Chapingo.

## Instalacion

Version publicada en CRAN:

```r
install.packages("GREENREG")
```

Version en desarrollo (con las correcciones mas recientes):

```r
# install.packages("remotes")
remotes::install_local("ruta/al/paquete/GREENREG")
```

## Ejemplo rapido

```r
library(GREENREG)

data(datos_rendimiento_maiz)

modelo <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
print(modelo)   # Coeficientes, VIF, supuestos y veredicto interpretado
plot(modelo)    # Secuencia de 10 graficas de diagnostico
```

## Funciones incluidas

| Funcion         | Modelo                                             |
|-----------------|-----------------------------------------------------|
| `greenreg_rls()`         | Regresion Lineal Simple                              |
| `greenreg_rlm()`         | Regresion Lineal Multiple                             |
| `greenreg_logistico()`   | Regresion Logistica (binomial)                        |
| `greenreg_reg_poisson()` | Regresion de Poisson                                   |
| `greenreg_modelo_ar()`   | Modelo Autorregresivo (AR)                             |
| `greenreg_modelo_ma()`   | Modelo de Media Movil (MA)                             |
| `greenreg_modelo_arma()` | Modelo Autorregresivo de Media Movil (ARMA)            |
| `greenreg_modelo_arima()`| Modelo Autorregresivo Integrado de Media Movil (ARIMA) |
| `greenreg_analisis_datos()` | Sugiere que modelo usar segun el tipo de variable   |

## Mas informacion

Consulta la viñeta `vignette("introduccion-greenreg")` para un recorrido
completo, y `NEWS.md` para el historial de cambios.

## Reportar problemas

Este paquete se mantiene activamente; si encuentras un comportamiento
inesperado, por favor documenta el caso (datos, codigo y salida obtenida) y
contacta al mantenedor listado en `DESCRIPTION`.
