# Caso práctico de Geoestadística en R

Este proyecto aborda un caso práctico de geoestadística utilizando la librería `geoR` en R. El objetivo es analizar y modelar el número de suicidios en 32 municipios de
Londres (excluyendo City London) en el periodo de 1989-1993 para hombres y mujeres combinados. Las
variables que se registran son:
- Número de suicidios observados en el periodo estudiado (y),
- Número de casos de suicidio esperados (E),
- Indice de privación social (x1),
- Indice de fragmentación social (x2), que refleja la falta de conexiones sociales y de sentido de
comunidad.

Además, de las variables del sector geográfico, con el objetivo de explicar "y"", a través del SMR = y/E.

## Contenido del Proyecto
- `data/`: Contiene los datos necesarios para el análisis.
- `scripts/`: Scripts en R que reproducen cada etapa del caso.
- `results/`: Resultados generados (gráficos, mapas, tablas).

## Librerías: `geoR`, `spBayes`, `MBA`, `fields`, `sp`, `classInt`, `lattice`, `MASS`, `coda`, `spdep`.

## Scripts
El análisis completo se encuentra en el archivo `caso_geoestadistica.R`, ubicado en la carpeta `scripts/`. El script incluye:
1. Carga de librerías y datos: Importación de archivos como LDNSuicides.shp y LondonSuicides.RData, así como la     unión del polígono espacial con el conjunto de datos.

2. Análisis espacial inicial: Aplicación de los test de Moran y Geary para la variable SMR y evaluación de la       necesidad de un modelo espacial.

3. Ajuste de modelos espaciales: Estimación de los modelos M1 (lineal), M2 (SAR) y M3 (CAR), comparación mediante   el criterio AIC y selección del modelo óptimo.

4. Visualización y análisis de resultados:
    - Gráfico del SMR estimado vs. real.
    - Análisis de los residuos del modelo final y aplicación de los test de Moran y Geary.
    - Generación de gráficos de los residuos.


# Datos
Esta carpeta contiene los siguientes archivos:
- `LDNSuicides.shp`: Archivo de polígonos espaciales de los municipios.
- `LondonSuicides.RData`: Datos del estudio.


Para ejecutarlo:
```bash
Rscript scripts/caso_geoestadistica.R

