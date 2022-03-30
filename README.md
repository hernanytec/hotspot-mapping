# Grupo de Estudos de Aprendizado de Máquina (GEAM)
### Hotspot Mapping: identificando locais de risco para a ocorrência de crimes

## Sobre o problema
 - ### **Contexto**
    > #### Imagine que você foi encarregado para decidir quais áreas da região de Manhattan devem receber mais patrulhamento policial. 
   > **Como você resolveria esse problema?**

 - ### **Informações**
   - Conjunto de dados com registros de crimes de toda a cidade de NY
   - A força policial disponível é capaz de cobrir apenas 5% da região
 
 - ### **Objetivo**
  - Identificar áreas de risco na região de **Manhattan**
  - Determinar o nome das ruas da cidade que deverão receber reforço no patrulhamento policial

# Conjunto de dados

> [New York City Crimes](
https://www.kaggle.com/datasets/adamschroeder/crimes-new-york-city?resource=download)

 - Crimes reportados em todos os 5 bairros da cidade de NY (2014-2015)
 - Subconjunto dos dados disponíveis em [NYPD Complaint Data Historic, 2006 - 2019](https://data.cityofnewyork.us/Public-Safety/NYPD-Complaint-Data-Historic/qgea-i56i)

# O que foi visto com esse exercício?

  - [x] Operações com dados geoespaciais usando GeoPandas
  - [x] Tipos de Sistemas de projeção de coordenadas 
  - [x] Visualização dos dados usando contextily e matplotlib
  - [x] Geocoding e obtenção de dados do openstreetmap usando OSMNX
  - [x] Aplicação do KDE com grade de células para a geração de Hotspots

## Referências
 - [DOC] [Geopandas](https://geopandas.org/en/stable/docs.html)
 - [DOC] [OSMNX](https://osmnx.readthedocs.io/en/stable/)
 - [PAPER] [The Utility of Hotspot Mapping for Predicting Spatial Patterns of Crime](https://link.springer.com/article/10.1057/palgrave.sj.8350066)
 - [SITE] [Kernel Density Estimation](https://mathisonian.github.io/kde/)
