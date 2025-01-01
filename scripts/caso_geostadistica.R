#Librerias:
library(spBayes)
library(MBA)
library(geoR)
library(fields)
library(sp)
library(classInt)
library(lattice)
library(MASS)
library(coda)
library(spdep)
remove(list=ls())
ls()
library(INLA)
inla.setOption(scale.model.default=FALSE)

# Cargamos los archivos “LDNSuicides.shp” (polígono espacial) y “LondonSuicides.RData” (datos). 
# y Luego unimos el polígono espacial con el data.frame que contiene las variables.

load("London Suicides/LondonSuicides.RData")
ls()
london.gen <- st_read("London Suicides/LDNSuicides.shp")
# Ver los nombres de las columnas de los datos
names(london.gen)
# Ver las primeras filas de los datos
head(london.gen)

# Acceder a la geometría de los polígonos
london_geometry <- st_geometry(london.gen)

# Mostrar las geometrías
print(london_geometry)

# Graficar la geometría
plot(london_geometry)

if (!dir.exists("Graficas_Espacial")) {
  dir.create("Graficas_Espacial")
}

temp <- poly2nb(london_geometry) #Esto calcula las vecindades espaciales basadas en las geometrías.
nb2INLA("LDN.graph", temp) #Genera el archivo de adyacencia requerido para el análisis con INLA.

LDN.adj <- paste(getwd(),"/LDN.graph",sep="") #Contiene la ruta completa del archivo
#de adyacencia que describe la conexión espacial entre las áreas.

H <- inla.read.graph(filename="LDN.graph")
H

# inla.read.graph(): Esta función de la librería INLA se utiliza para leer un archivo 
# de adyacencia previamente creado (en este caso, LDN.graph).
# El archivo contiene la información de vecindad espacial entre las áreas (polígonos).
# Describe:
# Número de nodos (áreas geográficas).
# Vecinos: Indica qué áreas están conectadas (adyacentes) entre sí.
#$nnbs
# [1] 4 4 6 4 5 3 3 4 7 5 6 4 4 6 4 3 5 4 3 4 5 7 6 4 5 5 5 3 5 6 6 5
# Sitio 1 esta conectado con 4 sitios, el 2 tambien con 4 sitios


## Se une el objeto espacial con los datos
names <- sort(london.gen$NAME) 
data.suicides <- data.frame(NAME=names, y=y, E=E, x1=x1, x2=x2) # datos de los suicidios
Nareas <- length(data.suicides[,1])

# Reordenamos los datos
#data.boroughs <- attr(london.gen, "data") # En sf no es necesario usar attr
data.boroughs <- london.gen  # london.gen es ya un `data.frame` con geometría en `sf`
order <- match(data.boroughs$NAME,data.suicides$NAME)
data.suicides <- data.suicides[order,]
data.suicides$ID <- seq(1,Nareas)

# Incluimos en el dataframe de london.gen la variable ID
london.gen <- merge(data.boroughs,data.suicides,by="NAME")
names(london.gen)


#Especifico para sf
#install.packages("ggplot2")
library(ggplot2)

# Creamos el gráfico de distribución de suicidios
ggplot(data = london.gen) +
  geom_sf(aes(fill = y)) +                  # Mapea la columna 'y' como relleno de color
  labs(title = "Suicidios en Londres \n entre 1989 y 1993") +
  theme_minimal()                           # Estilo limpio de ggplot2



# Aplicamos los test de Moran y Geary para la variable SMR, con el objetivo de analizar si se justifica
# la utilización de un modelo espacial para estos datos (utilizaremos la matriz W estandarizada por fila).

# Calcular el SMR
london.gen$SMR <- london.gen$y / london.gen$E
#y: Número de suicidios observados.
#E: Número de suicidios esperados.
# SMR>1: más eventos observados de lo esperado
# SMR<1: menos eventos observados de lo esperado


# Creamos la matriz de pesos espaciales basada en la matriz de adyacencia
library(spdep)

nb <- poly2nb(london.gen)  # Crear la lista de vecinos
listw <- nb2listw(nb, style = "W")  # Matriz de pesos espaciales estandarizada por filas ("nb")

# Test de Moran
moran.test(london.gen$SMR, listw)
#Test de Moran:
#0.257
#p-valor: 
#0.0039
#Interpretación: Altamente significativa, lo que indica una autocorrelación espacial positiva.

# Test de Geary
geary.test(london.gen$SMR, listw)
#Test de Geary:
#0.726
#p-valor: 
#0.00599
#Interpretación: Altamente significativa, confirmando también autocorrelación positiva.



# Con el objetivo de modelar de la mejor manera el SMR, propondremos tres modelos: M1 (lineal), M2 (Modelo SAR),
# M3 (modelo CAR) y eligiremos el modelo cuyo resultado AIC sea menor:

#(M1) : y = β0 + β1x1 + β2x2 + β3x1x2 + η, η ∼ N(0, σ2);
#(M2) : y = β0 + β1x1 + β2x2 + β3x1x2 + η, η ∼ SAR(ρ);
#(M3) : y = β0 + β1x1 + β2x2 + β3x1x2 + η, η ∼ CAR(ρ)

# Para los modelos SAR y CAR utilizaremos la libreria "spatialreg".

# Modelo 1:
# Ajuste del modelo M1
M1 <- lm(SMR ~ x1 + x2 + I(x1*x2), data = london.gen)
summary(M1)

# Modelo 2:
library(spatialreg)

# Ajuste del modelo M2 (SAR)
M2 <- errorsarlm(SMR ~ x1 + x2 +I(x1*x2), data = london.gen, listw = listw)
summary(M2)


# Modelo 3:
# Ajuste del modelo M3 (CAR)

listw_binary <- nb2listw(nb, style = "B")  # Matriz binaria
M3 <- spautolm(SMR ~ x1 + x2 + I(x1*x2), data = london.gen, listw = listw_binary, family = "CAR")
summary(M3)


AIC(M1, M2, M3)

# Resultado de los modelos:
# Teniendo en cuenta el AIC se elegiría el modelo de regresión lineal (M1), sin embargo,
# existe el antecedente de los test de Moran y geary el cual nos señala que estos datos
# tienen una autocorrelación espacial positiva, por lo que nos inclinamos a utilizar un modelo
# que capture mejor esa dependencia espacial, SAR o CAR. Dentro de estos dos modelos el 
# que presenta un mejor AIC es el modelo CAR (M3).

# Observación:
# Para el modelo seleccionado las variables significativas son x1 y x2, pero no asi la
# interacción entre estas variables x1*x2, por lo que se procede a eliminar del modelo:

# Nuevo modelo 2 sin la multiplicación entre las variables
M3_2 <- spautolm(SMR ~ x1 + x2, data = london.gen, listw = listw_binary, family = "CAR")
summary(M3_2)


# Calculamos los valores ajustados y residuos con el objetivo de comparar el SMR real vs SMR estimado:
london.gen$SMR_estimado_M3 <- fitted(M3_2)  # Valores ajustados del modelo
london.gen$residuos_M3 <- residuals(M3_2)   # Residuos del modelo


# Graficamos el SMR real y el SMR estimado
# Mapa del SMR real
g1 <-ggplot(data = london.gen) +
  geom_sf(aes(fill = SMR)) +
  scale_fill_viridis_c() +
  labs(title = "SMR Real (CAR)", fill = "SMR") +
  theme_minimal()

# Mapa del SMR estimado
g2 <-ggplot(data = london.gen) +
  geom_sf(aes(fill = SMR_estimado_M3)) +
  scale_fill_viridis_c() +
  labs(title = "SMR Estimado (CAR)", fill = "SMR Estimado") +
  theme_minimal()

# Mapa de los residuos
g3 <-ggplot(data = london.gen) +
  geom_sf(aes(fill = residuos_M3)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Residuos del Modelo CAR", fill = "Residuos") +
  theme_minimal()


# Visualizamos los graficos SMR Real vs Estimado:
library(gridExtra)
library(grid)

grid.arrange(
  g1, g2, 
  nrow = 2, 
  heights = c(1, 1), 
  top = textGrob("Comparación SMR Real vs Estimado", gp = gpar(fontsize = 20)))


# Estos graficos nos permiten nos permiten comparar el smr real vs el estimado y visualizar que tan bien el
# modelo captura los patrones espaciales.
# Las discrepancias entre el SMR real y estimado pueden indicar áreas donde el modelo no se ajusta bien.
# Si el modelo logra replicar bien los valores reales en términos espaciales:
# Esto valida la elección de las covariables y la estructura espacial del modelo (como el uso de un modelo CAR).
# Si no lo hace, podría ser necesario ajustar el modelo, incluir nuevas variables o reconsiderar la matriz de pesos espaciales.


# Con el smr estimado podemos Evaluar el impacto de factores explicativos:
# Identificar factores socioeconómicos que contribuyen al riesgo en cada área,
# (Indice de privación social(x1), Indice de fragmentación social(x2)).

# Informar decisiones futuras:
# Si las condiciones cambian (por ejemplo, aumenta la privación social), 
# el modelo puede predecir cómo esto afectará el SMR, ayudando a planificar
# intervenciones proactivas.



# Como análisis final, aplicaremos los test de Moran y geary a los residuos del modelo final y evaluaremos
# si existe dependencia espacial:

# Test de Moran para residuos del modelo CAR
moran_car <- moran.test(london.gen$residuos_M3, listw_binary)
moran_car
#Test de Moran:
# −0.122 (cercano a 0, indicando falta de autocorrelación).
# valor p: 0.8045


# Test de Geary para residuos del modelo CAR
geary_car <- geary.test(london.gen$residuos_M3, listw_binary)
geary_car
# Test de Geary:
# 1.018 (cercano a 1, lo que indica independencia espacial).
# valor p:  0.5589

#No hay evidencia de autocorrelación espacial en los residuos del modelo CAR.
#Esto significa que el modelo CAR ha capturado correctamente la dependencia espacial 
#en los datos. Los residuos restantes son independientes y no presentan patrones espaciales
#significativos.

# Mapa de los residuos
g3
# Objetivo de visualizar el mapa de los residuos: 
# Localizar regiones donde el modelo subestima o sobreestima el SMR.

# Residuos positivos (en rojo): El modelo subestima el SMR real en esas áreas.
# Residuos negativos (en azul): El modelo sobreestima el SMR real.
# Esto puede indicar características locales o factores no incluidos en el 
# modelo que afectan la variable de interés.

# Verificamos la normalidad de los residuos
# Aunque los modelos espaciales como CAR o SAR tienen en cuenta la autocorrelación
# espacial, es igualmente importante verificar la normalidad de los residuos, 
# ya que una distribución no normal podría indicar problemas en el ajuste o 
# en la especificación del modelo.

# Test de normalidad
# H0: residuos siguen una dist normal
# H1: residuos no siguen una dist normal
shapiro.test(london.gen$residuos)
# Como no estamos en la region de rechazo del valor p (5%) No se rechaza H0





#En un modelo espacial bien ajustado:
#Los datos originales pueden tener autocorrelación espacial (capturada por el modelo).
#Los residuos no deben tener autocorrelación espacial significativa, ya que esto indicaría 
#que el modelo no ha capturado completamente la dependencia espacial.






