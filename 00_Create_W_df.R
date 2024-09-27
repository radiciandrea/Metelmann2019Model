# README

# code to setup W_tot_df with synthetic data
# the dataframe which contains, for each simulated day and region, 
# population density (/km2), temperatures (°C), rain (mm)

library(dplyr)

rm(list = ls())

# spatal
regions = 1:3
pop = c(10, 200, 3000) # rural, intermediate, city
lat = c(40, 45, 50)
lon = c(4, 2, 0)

# temporal
years = c(2001, 2002, 2003)
DOY = 1:365
DOS = DOY

#size
n_r = length(regions)
n_d = 365

# weather is pseudorandomly generated
T_min = 10-10*cos(2*pi*DOS/365)
T_max = 18-8*cos(2*pi*DOS/365)
T_av = (T_min+T_max)*0.5
P = 2*sqrt((183 - abs(DOY - 182)) *(DOS %% 4 > DOS %% 23))

for (year in years){
  
  W_tot_df = data.frame(region = rep(regions, each = n_d),
                        lat = rep(lat, each = n_d),
                        lon = rep(lon, each = n_d),
                        pop = rep(pop, each = n_d),
                        date = as.Date(DOY, origin = paste0(year-1, "-12-31")),
                        DOS = DOS,
                        DOY = DOY,
                        T_m = 2*runif(n_r*n_d)+rep(T_min, n_r),
                        T_M = 2*runif(n_r*n_d)+rep(T_min, n_r),
                        T_av = 2*runif(n_r*n_d)+rep(T_min, n_r),
                        P = runif(n_r*n_d)*rep(P, n_r))
  
  DOS = DOY+max(DOS)
  
  save(W_tot_df, file = paste0("W_tot_", year, ".RData"))
}