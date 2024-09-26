# Metelmann2019Model
The present code is a unofficial R version of the model by Metelmann et al 2019 (https://doi.org/10.1098/rsif.2018.0761). It is a mechanistic model describing the population dynamics of Ae. albopictus of different sites using ordinary differential equations. The Ae. albopictus populations are divided into five classes (eggs, diapausing eggs, juveniles, immature adults, adults). Specific environmental (human density) and meteorological (daily maximal, mean and minimum temperature, daily rain, photoperiod) drivers affect their survival, fertility and development rates. 


The author of the code rejects all responsibility about any discrepancies from the model developed in the article.

The model needs environmental and metereological  variables to run, formatted in a dataframe as the one specified in "00_Create_W_df.R". In absence of data, this script needs to be run first to generate synthetic series. The model itself is "01_Model.R", to be run afterwards. This will autonomously call the script "Integration_functions.R", which contains the ODE systmes (normal and log-transformed).
