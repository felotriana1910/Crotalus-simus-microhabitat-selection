### Libraries

library(ggplot2)
library(dplyr)
library(corrplot)
library(readxl)
library(adehabitatHR)
library(adehabitatHS)
library(sp)
library(sf)
library(lme4)
library(broom.mixed)
library(popbio)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(Hmisc)
library(MuMIn)
library(ggeffects)
library(corrplot)
library(glmmTMB)
library(DHARMa)
library(cowplot)
library(tidyr)

### We call the database identified by the DOI: 10.17632/s94ccbvzyj.1

# URL of the database
data_url <- "https://data.mendeley.com/public-files/datasets/s94ccbvzyj/files/429ceb5d-32ed-47a1-89c5-1d6b103245d7/file_downloaded"

# We download the file to a temporary location and read it
Csimus_microhabitat_data <- tempfile(fileext = ".csv")
download.file(data_url, destfile = Csimus_microhabitat_data, mode = "wb")

Habitat <- read.csv2(Csimus_microhabitat_data,
                     stringsAsFactors = FALSE,
                     fileEncoding = "latin1")

### Data pre-processing

# We convert the coordinates from degrees-minutes-seconds format to decimal degrees
dms_a_decimal <- function(coord) {
  # Remove spaces and separate degrees, minutes, and seconds
  partes <- regmatches(coord, regexec("([0-9]+)°([0-9]+)'([0-9]+[,.]?[0-9]*)", coord))[[1]]
  
  grados <- as.numeric(partes[2])
  minutos <- as.numeric(partes[3])
  segundos <- as.numeric(gsub(",", ".", partes[4]))  # Replace comma with period if necessary
  
  # Convert to decimal
  decimal <- grados + minutos / 60 + segundos / 3600
  return(decimal)
}

# Remove the negative sign before converting
Habitat$lat_decimal <- sapply(Habitat$Latitude, dms_a_decimal)
Habitat$lon_decimal <- -sapply(gsub("^-", "", Habitat$Longitude), dms_a_decimal)  # Preserve original negative sign

# Pseudo-absences are already in decimal degrees, so they are copied directly
# (replacing the decimal comma with a period)
Habitat$lat_decimal[Habitat$ID2 == "PS"] <- gsub(",", ".", Habitat$Latitude[Habitat$ID2 == "PS"])
Habitat$lon_decimal[Habitat$ID2 == "PS"] <- gsub(",", ".", Habitat$Longitude[Habitat$ID2 == "PS"])

# We set the variables as we need them
Habitat$lon_decimal <- as.numeric(Habitat$lon_decimal)
Habitat$lat_decimal <- as.numeric(Habitat$lat_decimal)
Habitat$Diel_period <- as.factor(Habitat$Diel_period)
Habitat$Substrate[Habitat$Substrate == ""] <- NA
Habitat$Substrate <- as.factor(Habitat$Substrate)

# We filter out the individuals with 23 or fewer data points (threshold 
# we decided on in order to analyze them)

Habitat_filtrada <- Habitat %>% filter(!ID %in% c("R2", "R6", "R7"))

## Pearson correlation between environmental variables

# We select only the variables of interest
variables_amb <- Habitat_filtrada[, c("Dist_large_tree", "Dist_small_tree", 
                                      "Dist_fallen_log", "Canopy_cover", 
                                      "Str_understory", "Temperature", 
                                      "Humidity", "Dist_water_source")]

# Correlation matrix
matriz_cor <- cor(variables_amb, method = "pearson", use = "complete.obs")
print(round(matriz_cor, 2)) # View matrix

# Visualization
corrplot(matriz_cor, 
         method = "color",        
         type = "upper",          
         addCoef.col = "black",   
         number.cex = 0.7,        
         tl.col = "black",        
         tl.srt = 45,             
         diag = FALSE)            


# Matrix with "r" y "p-values"
matriz_completa <- rcorr(as.matrix(variables_amb), type = "pearson")

# View r coefficients
print(round(matriz_completa$r, 2))

# View p-values
print(round(matriz_completa$P, 3))

# We set the data bases for analyses
# We obtain the home ranges (100% MCP) of each individual in order to compare the 
# locations selected by snakes only with the pseudo-absences within that individual's 
# home range, and thus truly have a measure of the habitat available to each individual

# We separate presences and absences
Presence <- Habitat_filtrada %>% filter(ID2 != "PS")
pseudoausencia <- Habitat_filtrada %>% filter(ID2 == "PS")

# Convert pseudo-absences to a spatial object
coordinates(pseudoausencia) <- ~lon_decimal+lat_decimal
proj4string(pseudoausencia) <- CRS("+proj=longlat +datum=WGS84")

DFfinal1=data.frame(matrix(ncol=(ncol(Habitat)+1),nrow=0)) # New database
colnames(DFfinal1)=c(colnames(Habitat),"IDF") # Add new variable for random effect in the model
for (i in unique(Presence$ID)){
  df <- Presence[Presence$ID==i,]
  coordinates(df) <- ~lon_decimal+lat_decimal
  proj4string(df) <- CRS("+proj=longlat +datum=WGS84")  # Geographic coordinate system
  df$IDF=i
  # 🔹 3. Calculate MCP 100%
  mcp <- mcp(df, percent = 100)
  
  esta_dentro <- !is.na(over(pseudoausencia, mcp)) # if the point is inside
  puntos_en_mcp <- pseudoausencia[esta_dentro[,1], ]
  puntos_en_mcp$IDF=i
  # 🔹 4. Plot
  
  plot(mcp, col = "lightblue", main = "MCP100 - Área mínima convexa")
  points(puntos_en_mcp, col = "red", pch = 20)
  #print(length(puntos_en_mcp))
  DFfinalPre1=rbind(df,puntos_en_mcp)
  DFfinal1=rbind(DFfinal1,DFfinalPre1@data)
}

### We rescale variables for analyses

DFfinal1_scaled <- DFfinal1 %>%
  dplyr::mutate(
    Dist_large_tree      = scale(Dist_large_tree),
    Dist_small_tree      = scale(Dist_small_tree),
    Dist_fallen_log      = scale(Dist_fallen_log),
    Dist_water_source  = scale(Dist_water_source),
    Canopy_cover         = scale(Canopy_cover),
    Str_understory       = scale(Str_understory),
    Temperature          = scale(Temperature),
    Humidity             = scale(Humidity)
  )


### Compositional analysis to examine substrate selection, based on Aebischer et al. (1993)

## We first create a dataset of the used and available substrates

#Create the data bases
pres_df <- DFfinal1_scaled %>% filter(Presence == 1)
aus_df <- DFfinal1_scaled %>% filter(Presence == 0)

# Matrix of USE per individual
used <- pres_df %>%
  group_by(IDF, Substrate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(IDF) %>%
  mutate(pct = n / sum(n) * 100) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = Substrate, values_from = pct, values_fill = 0) %>%
  as.data.frame()

rownames(used) <- used$IDF
used <- used[, -1]

# Matrix of availability per individual
avail <- aus_df %>%
  group_by(IDF, Substrate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(IDF) %>%
  mutate(pct = n / sum(n) * 100) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = Substrate, values_from = pct, values_fill = 0) %>%
  as.data.frame()

rownames(avail) <- avail$IDF
avail <- avail[, -1]


# We make sure the columns are the same
sustratos <- union(colnames(used), colnames(avail))
used[setdiff(sustratos, colnames(used))] <- 0
avail[setdiff(sustratos, colnames(avail))] <- 0
used <- used[, sustratos]
avail <- avail[, sustratos]

# Remove the "NA" column as a string
used <- used[, colnames(used) != "NA"]
avail <- avail[, colnames(avail) != "NA"]
colnames(used)
colnames(avail)

# We exclude only Rock because it is a missing available habitat type for most snakes
used_sub2 <- used[, colnames(used) != "Rock"]
avail_sub2 <- avail[, colnames(avail) != "Rock"]

# Check
colSums(avail_sub2 > 0)
colSums(used_sub2 > 0)

# We run the compositional analysis
resultado_comp <- compana(used_sub2, avail_sub2, 
                           test = "randomisation", 
                           nrep = 999,
                           rnv = 0.01)
print(resultado_comp)


## Models to assess microhabitat selection: generalized linear mixed model (glmm) 

# We add a variable for Translocated vs Natives
# Modelo total viendo la interacción del día/noche y translocada/residente por variable

DFfinal1_scaled$ConditionSnk="TRASLOCADAS"
DFfinal1_scaled$ConditionSnk[DFfinal1_scaled$IDF %in% c("N1", "N2", "N3", "N4", "N5")]="NATIVAS"

# Identify complete rows
DFfinal1_scaled_completo <- DFfinal1_scaled %>%
  filter(complete.cases(Dist_large_tree, Dist_small_tree, Dist_fallen_log,
                        Dist_water_source, Canopy_cover, Str_understory,
                        Temperature, Humidity)) 

## glmm with package glmmTMB

# Full model

DFfinal1_scaled_completo$ConditionSnk <- factor(DFfinal1_scaled_completo$ConditionSnk)
DFfinal1_scaled_completo$Diel_period <- factor(DFfinal1_scaled_completo$Diel_period)
 
modeloT_inter_completo <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + 
                                          Dist_fallen_log + Dist_water_source + Canopy_cover + 
                                          Temperature + Str_understory + Humidity +
                                          (Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                                             Dist_water_source + Canopy_cover + Temperature + 
                                             Str_understory + Humidity):Diel_period +
                                          (Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                                             Dist_water_source + Canopy_cover + Temperature + 
                                             Str_understory + Humidity):ConditionSnk +
                                          (1 | IDF),
                                        data = DFfinal1_scaled_completo,
                                        family = binomial)

summary(modeloT_inter_completo)


#### Step 1: We remove non-significant variables from the model

# Model without Dist_large_tree or Canopy_cover, which are the ones that are not significant at all
model1 <- glmmTMB(Presence ~ Dist_small_tree + Dist_fallen_log + 
                    Dist_water_source + Temperature + 
                    Str_understory + Humidity +
                    (Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                       Temperature + Str_understory + Humidity):Diel_period +
                    (Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                       Temperature + Str_understory + Humidity):ConditionSnk +
                    (1 | IDF),
                  data = DFfinal1_scaled_completo,
                  family = binomial)


summary(model1)

# Model without large tree
model2 <- glmmTMB(Presence ~ Dist_small_tree + Dist_fallen_log + 
                    Dist_water_source + Temperature + Canopy_cover +
                    Str_understory + Humidity +
                    (Dist_small_tree + Dist_fallen_log + Dist_water_source + Canopy_cover +
                       Temperature + Str_understory + Humidity):Diel_period +
                    (Dist_small_tree + Dist_fallen_log + Dist_water_source + Canopy_cover +
                       Temperature + Str_understory + Humidity):ConditionSnk +
                    (1 | IDF),
                  data = DFfinal1_scaled_completo,
                  family = binomial)

# Modelo without canopy cover
model3 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                    Dist_water_source + Temperature + 
                    Str_understory + Humidity +
                    (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                       Temperature + Str_understory + Humidity):Diel_period +
                    (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                       Temperature + Str_understory + Humidity):ConditionSnk +
                    (1 | IDF),
                  data = DFfinal1_scaled_completo,
                  family = binomial)


# Model removing non-significant interactions with Diel_period or Translocation
model4 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                    Dist_water_source + Canopy_cover + Temperature + 
                    Str_understory + Humidity +
                    Temperature:Diel_period +
                    Str_understory:Diel_period +
                    Humidity:Diel_period +
                    Dist_small_tree:ConditionSnk +
                    Dist_fallen_log:ConditionSnk +
                    Str_understory:ConditionSnk +
                    Humidity:ConditionSnk +
                    (1 | IDF),
                  data = DFfinal1_scaled_completo,
                  family = binomial)

summary(model4)


# PHASE 1: Remove one interaction with Night at a time

f1_sinDCT_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_small_tree + Dist_fallen_log + Dist_water_source + Canopy_cover + 
                            Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDST_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_large_tree + Dist_fallen_log + Dist_water_source + Canopy_cover + 
                            Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDT_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_water_source + Canopy_cover + 
                           Temperature + Str_understory + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDA_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Canopy_cover + 
                           Temperature + Str_understory + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinCD_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Temperature + Str_understory + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinT_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                       Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                       (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Str_understory + Humidity):Diel_period +
                       (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                       (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinES_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinH_N <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                       Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                       (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Temperature + Str_understory):Diel_period +
                       (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                       (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

# PHASE 1: Remove one interaction with Translocated at a time

f1_sinDCT_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                          Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                          (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                             Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                          (Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                             Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                          (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDST_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                          Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                          (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                             Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                          (Dist_large_tree + Dist_fallen_log + Dist_water_source + 
                             Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                          (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDT_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinDA_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                            Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinCD_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Temperature + Str_understory + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinT_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Str_understory + Humidity):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinES_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                         (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                            Canopy_cover + Temperature + Humidity):ConditionSnk +
                         (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f1_sinH_TR <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                        Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory + Humidity):Diel_period +
                        (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                           Canopy_cover + Temperature + Str_understory):ConditionSnk +
                        (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

# PHASE 2: Remove an entire variable (main effect + both interactions)

f2_sinDCT <- glmmTMB(Presence ~ Dist_small_tree + Dist_fallen_log + 
                       Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                       (Dist_small_tree + Dist_fallen_log + Dist_water_source + Canopy_cover + 
                          Temperature + Str_understory + Humidity):Diel_period +
                       (Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                       (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinDST <- glmmTMB(Presence ~ Dist_large_tree + Dist_fallen_log + 
                       Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                       (Dist_large_tree + Dist_fallen_log + Dist_water_source + Canopy_cover + 
                          Temperature + Str_understory + Humidity):Diel_period +
                       (Dist_large_tree + Dist_fallen_log + Dist_water_source + 
                          Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                       (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinDT <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + 
                      Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                      (Dist_large_tree + Dist_small_tree + Dist_water_source + Canopy_cover + 
                         Temperature + Str_understory + Humidity):Diel_period +
                      (Dist_large_tree + Dist_small_tree + Dist_water_source + 
                         Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                      (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinDA <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                      Canopy_cover + Temperature + Str_understory + Humidity +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Canopy_cover + 
                         Temperature + Str_understory + Humidity):Diel_period +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                         Canopy_cover + Temperature + Str_understory + Humidity):ConditionSnk +
                      (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinCD <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                      Dist_water_source + Temperature + Str_understory + Humidity +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                         Temperature + Str_understory + Humidity):Diel_period +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                         Temperature + Str_understory + Humidity):ConditionSnk +
                      (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinT <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                     Dist_water_source + Canopy_cover + Str_understory + Humidity +
                     (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                        Canopy_cover + Str_understory + Humidity):Diel_period +
                     (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                        Canopy_cover + Str_understory + Humidity):ConditionSnk +
                     (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinES <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                      Dist_water_source + Canopy_cover + Temperature + Humidity +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                         Canopy_cover + Temperature + Humidity):Diel_period +
                      (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                         Canopy_cover + Temperature + Humidity):ConditionSnk +
                      (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f2_sinH <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                     Dist_water_source + Canopy_cover + Temperature + Str_understory +
                     (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                        Canopy_cover + Temperature + Str_understory):Diel_period +
                     (Dist_large_tree + Dist_small_tree + Dist_fallen_log + Dist_water_source + 
                        Canopy_cover + Temperature + Str_understory):ConditionSnk +
                     (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)


# Phase 3: Combinations of two additional interactions starting from model4
# model4 already has: T:N, Str:N, H:N, DST:TR, DFL:TR, Str:TR, H:TR

f3_1 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + Humidity:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  Dist_water_source:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_2 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + Humidity:Diel_period +
                  Canopy_cover:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_3 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + Humidity:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  Canopy_cover:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_4 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + Humidity:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  Temperature:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_5 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + 
                  Humidity:Diel_period + Dist_small_tree:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_6 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + 
                  Humidity:Diel_period + Dist_fallen_log:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_7 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + 
                  Humidity:Diel_period + Dist_water_source:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

f3_8 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + 
                  Humidity:Diel_period + Dist_large_tree:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)


f4_1 <- glmmTMB(Presence ~ Dist_large_tree + Dist_small_tree + Dist_fallen_log + 
                  Dist_water_source + Canopy_cover + Temperature + Str_understory + Humidity +
                  Temperature:Diel_period + Str_understory:Diel_period + 
                  Humidity:Diel_period +
                  Dist_small_tree:ConditionSnk + Dist_fallen_log:ConditionSnk + 
                  Str_understory:ConditionSnk + Humidity:ConditionSnk +
                  Temperature:ConditionSnk + Canopy_cover:ConditionSnk +
                  (1 | IDF), data = DFfinal1_scaled_completo, family = binomial)

## We make the AIC table

modelos_lista <- list(
  "Completo" = modeloT_inter_completo,
  "model1" = model1,
  "model2" = model2,
  "model3" = model3,
  "model4" = model4,
  "f1_sinDCT_N" = f1_sinDCT_N,
  "f1_sinDST_N" = f1_sinDST_N,
  "f1_sinDT_N" = f1_sinDT_N,
  "f1_sinDA_N" = f1_sinDA_N,
  "f1_sinCD_N" = f1_sinCD_N,
  "f1_sinT_N" = f1_sinT_N,
  "f1_sinES_N" = f1_sinES_N,
  "f1_sinH_N" = f1_sinH_N,
  "f1_sinDCT_TR" = f1_sinDCT_TR,
  "f1_sinDST_TR" = f1_sinDST_TR,
  "f1_sinDT_TR" = f1_sinDT_TR,
  "f1_sinDA_TR" = f1_sinDA_TR,
  "f1_sinCD_TR" = f1_sinCD_TR,
  "f1_sinT_TR" = f1_sinT_TR,
  "f1_sinES_TR" = f1_sinES_TR,
  "f1_sinH_TR" = f1_sinH_TR,
  "f2_sinDCT" = f2_sinDCT,
  "f2_sinDST" = f2_sinDST,
  "f2_sinDT" = f2_sinDT,
  "f2_sinDA" = f2_sinDA,
  "f2_sinCD" = f2_sinCD,
  "f2_sinT" = f2_sinT,
  "f2_sinES" = f2_sinES,
  "f2_sinH" = f2_sinH,
  "f3_1" = f3_1,
  "f3_2" = f3_2,
  "f3_3" = f3_3,
  "f3_4" = f3_4,
  "f3_5" = f3_5,
  "f3_6" = f3_6,
  "f3_7" = f3_7,
  "f3_8" = f3_8,
  "f4_1" = f4_1
)

# AIC values
aicc_vals <- sapply(modelos_lista, AICc)
delta_aicc <- aicc_vals - min(aicc_vals)
w_aicc <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Table with AIC
tabla_modelos <- data.frame(
  AICc = round(aicc_vals, 1),
  Delta_AICc = round(delta_aicc, 1),
  wAIC = round(w_aicc, 4)
)

tabla_modelos <- tabla_modelos[order(tabla_modelos$AICc), ]
print(tabla_modelos)

# The model we are going to select is model 4. It is within Delta_AICc < 2 and 
# is the one with the fewest parameters.

### We are going to check the overdispersion of the model residuals with DHARMa.

# We calculate the residuals
simulationOutput <- simulateResiduals(fittedModel = model4, plot = FALSE)
residuals(simulationOutput)

# Plot scaled residuals
plot(simulationOutput)

# Test for over/underdispersion
testDispersion(simulationOutput)

# Test for sero-inflation
testZeroInflation(simulationOutput)


### Plots for publication 

### Substrate

# Put observations in a database
substrate_table <- table(
  group = c(rep("Presence", nrow(Presence)), rep("Pseudoabsence", nrow(pseudoausencia))),
  substrate = c(Presence$Substrate, pseudoausencia$Substrate)
)

substrate_table <- as.data.frame(substrate_table)
substrate_table$substrate <- factor(substrate_table$substrate, levels = c("Rock", "Vegetation", "Small_stones",
                                                                          "Log", "Bare_ground", "Leaf_litter"))

# Plot
ggplot(substrate_table, aes(x = factor(group, labels = c("Habitat used", "Available habitat")), 
                         y = Freq,
                         fill = substrate)) +
  geom_bar(position = "fill", stat = "identity", color = "black") + 
  scale_fill_manual(values = c("Rock"         = "#E7298A",
                               "Vegetation"   = "#66A61E",
                               "Small_stones" = "#E6AB02",
                               "Log"          = "#7570B3",
                               "Bare_ground"  = "#D95F02",
                               "Leaf_litter"  = "#1B9E77"),
                    labels = c("Rock"         = "Rock",
                               "Vegetation"   = "Vegetation",
                               "Small_stones" = "Small stones",
                               "Log"          = "Log",
                               "Bare_ground"  = "Bare ground",
                               "Leaf_litter"  = "Leaf litter")) +
  labs(y = "Proportion", fill = "Substrate") + 
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank())


### Numerical variables

## Distance large tree

# We save the original mean and standard deviation used to scale Dist_large_tree
mean_DCT <- attr(DFfinal1_scaled$Dist_large_tree, "scaled:center")
sd_DCT <- attr(DFfinal1_scaled$Dist_large_tree, "scaled:scale")
pred_DCT <- ggpredict(model4, 
                      terms = "Dist_large_tree [all]",
                      bias_correction = TRUE)

# Unscale the predictions
pred_DCT$x_real <- pred_DCT$x * sd_DCT + mean_DCT

# Unscale the points
DFfinal1_scaled_completo$DCT_real <- DFfinal1_scaled_completo$Dist_large_tree * sd_DCT + mean_DCT

# Plot
graf_DCT <- ggplot(pred_DCT, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(size = 1) +
  geom_jitter(data = DFfinal1_scaled_completo, 
              aes(x = DCT_real, y = Presence),
              height = 0.02, alpha = 0.3, size = 1) +
  labs(x = "Distance to the nearest large tree (m)",
       y = "Probability of presence") +
  theme_bw(base_size = 14)

## Distance small tree

# Same process as above
mean_DST <- attr(DFfinal1_scaled$Dist_small_tree, "scaled:center")
sd_DST <- attr(DFfinal1_scaled$Dist_small_tree, "scaled:scale")

pred_DST <- ggpredict(model4, 
                      terms = c("Dist_small_tree [all]", "ConditionSnk"),
                      bias_correction = TRUE)
pred_DST$x_real <- pred_DST$x * sd_DST + mean_DST
DFfinal1_scaled_completo$DST_real <- DFfinal1_scaled_completo$Dist_small_tree * sd_DST + mean_DST

# Plot
graf_DST <- ggplot() +
  geom_ribbon(data = pred_DST, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high,
                  group = group, fill = group), alpha = 0.2) +
  geom_line(data = pred_DST, 
            aes(x = x_real, y = predicted, color = group), 
            linewidth = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = DST_real, y = Presence, color = ConditionSnk),
              height = 0.02, alpha = 0.3, size = 1) +
  scale_color_manual(values = c("NATIVAS" = "#1b9e77", "TRASLOCADAS" = "#7570b3")) +
  scale_fill_manual(values = c("NATIVAS" = "#1b9e77", "TRASLOCADAS" = "#7570b3")) +
  labs(x = "Distance to the nearest small tree (m)",
       y = "Probability of presence") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

## Distance fallen log

mean_DFL <- attr(DFfinal1_scaled$Dist_fallen_log, "scaled:center")
sd_DFL <- attr(DFfinal1_scaled$Dist_fallen_log, "scaled:scale")
pred_DFL <- ggpredict(model4, 
                      terms = c("Dist_fallen_log [all]", "ConditionSnk"),
                      bias_correction = TRUE)
pred_DFL$x_real <- pred_DFL$x * sd_DFL + mean_DFL
DFfinal1_scaled_completo$DFL_real <- DFfinal1_scaled_completo$Dist_fallen_log * sd_DFL + mean_DFL

# Plot
graf_DFL <- ggplot() +
  geom_ribbon(data = pred_DFL, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high,
                  group = group, fill = group), alpha = 0.2) +
  geom_line(data = pred_DFL, 
            aes(x = x_real, y = predicted, color = group), 
            linewidth = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = DFL_real, y = Presence, color = ConditionSnk),
              height = 0.02, alpha = 0.3, size = 1) +
  scale_color_manual(values = c("NATIVAS" = "#1b9e77", "TRASLOCADAS" = "#7570b3")) +
  scale_fill_manual(values = c("NATIVAS" = "#1b9e77", "TRASLOCADAS" = "#7570b3")) +
  labs(x = "Distance to the nearest fallen log (m)",
       y = "Probability of presence") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

## Distance to the water

mean_DWS <- attr(DFfinal1_scaled$Dist_water_source, "scaled:center")
sd_DWS <- attr(DFfinal1_scaled$Dist_water_source, "scaled:scale")
pred_DWS <- ggpredict(model4, 
                      terms = "Dist_water_source [all]",
                      bias_correction = TRUE)
pred_DWS$x_real <- pred_DWS$x * sd_DWS + mean_DWS
DFfinal1_scaled_completo$DWS_real <- DFfinal1_scaled_completo$Dist_water_source * sd_DWS + mean_DWS

# Plot
graf_DWS <- ggplot() +
  geom_ribbon(data = pred_DWS, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey70") +
  geom_line(data = pred_DWS, 
            aes(x = x_real, y = predicted), 
            color = "black", size = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = DWS_real, y = Presence),
              height = 0.02, width = 1.8, alpha = 0.3, size = 1, color = "black") +
  labs(x = "Distance to the nearest water source (m)",
       y = "Probability of presence") +
  theme_bw(base_size = 14)  

## Canopy cover

mean_CC <- attr(DFfinal1_scaled$Canopy_cover, "scaled:center")
sd_CC <- attr(DFfinal1_scaled$Canopy_cover, "scaled:scale")
pred_CC <- ggpredict(model4, 
                     terms = "Canopy_cover [all]",
                     bias_correction = TRUE)
pred_CC$x_real <- pred_CC$x * sd_CC + mean_CC
DFfinal1_scaled_completo$CC_real <- DFfinal1_scaled_completo$Canopy_cover * sd_CC + mean_CC

# Plot
graf_CC <- ggplot() +
  geom_ribbon(data = pred_CC, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey70") +
  geom_line(data = pred_CC, 
            aes(x = x_real, y = predicted), 
            color = "black", size = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = CC_real, y = Presence),
              height = 0.02, alpha = 0.3, size = 1, color = "black") +
  labs(x = "Canopy cover (%)",
       y = "Probability of presence") +
  theme_bw(base_size = 14)

## Structure understory

mean_SU <- attr(DFfinal1_scaled$Str_understory, "scaled:center")
sd_SU <- attr(DFfinal1_scaled$Str_understory, "scaled:scale")
pred_SU <- ggpredict(model4, 
                     terms = c("Str_understory [all]", "Diel_period", "ConditionSnk"),
                     bias_correction = TRUE)
pred_SU$x_real <- pred_SU$x * sd_SU + mean_SU
DFfinal1_scaled_completo$SU_real <- DFfinal1_scaled_completo$Str_understory * sd_SU + mean_SU

# Plot
graf_SU <- ggplot() +
  geom_ribbon(data = pred_SU, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high,
                  group = interaction(group, facet), fill = group), alpha = 0.2) +
  geom_line(data = pred_SU, 
            aes(x = x_real, y = predicted, color = group, linetype = facet), 
            size = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = SU_real, y = Presence, color = Diel_period),
              height = 0.02, alpha = 0.3, size = 1) +
  scale_linetype_manual(values = c("NATIVAS" = "solid", "TRASLOCADAS" = "dashed")) +
  scale_color_manual(values = c("Day" = "#d95f02", "Night" = "#0072B2")) +
  scale_fill_manual(values = c("Day" = "#d95f02", "Night" = "#0072B2")) +
  labs(x = "Number of structures in the understory",
       y = "Probability of presence",
       color = "Diel period", fill = "Diel period", linetype = "Status") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

## Temperature

mean_T <- attr(DFfinal1_scaled$Temperature, "scaled:center")
sd_T <- attr(DFfinal1_scaled$Temperature, "scaled:scale")
pred_T <- ggpredict(model4, 
                    terms = c("Temperature [all]", "Diel_period"),
                    bias_correction = TRUE)
pred_T$x_real <- pred_T$x * sd_T + mean_T
DFfinal1_scaled_completo$T_real <- DFfinal1_scaled_completo$Temperature * sd_T + mean_T

# Plot
graf_T <- ggplot() +
  geom_ribbon(data = pred_T, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high,
                  group = group, fill = group), alpha = 0.2) +
  geom_line(data = pred_T, 
            aes(x = x_real, y = predicted, color = group), 
            linewidth = 1.2) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = T_real, y = Presence, color = Diel_period),
              height = 0.02, alpha = 0.3, size = 3) +
  scale_color_manual(values = c("Day" = "#D55E00", "Night" = "#0072B2")) +
  scale_fill_manual(values = c("Day" = "#D55E00", "Night" = "#0072B2")) +
  labs(x = "Temperature (°C)",
       y = "Probability of presence") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# Humidity

mean_H <- attr(DFfinal1_scaled$Humidity, "scaled:center")
sd_H <- attr(DFfinal1_scaled$Humidity, "scaled:scale")
pred_H <- ggpredict(model4, 
                    terms = c("Humidity [all]", "Diel_period", "ConditionSnk"),
                    bias_correction = TRUE)
pred_H$x_real <- pred_H$x * sd_H + mean_H
DFfinal1_scaled_completo$H_real <- DFfinal1_scaled_completo$Humidity * sd_H + mean_H

# Plot
graf_H <- ggplot() +
  geom_ribbon(data = pred_H, 
              aes(x = x_real, y = predicted, ymin = conf.low, ymax = conf.high,
                  group = interaction(group, facet), fill = group), alpha = 0.2) +
  geom_line(data = pred_H, 
            aes(x = x_real, y = predicted, color = group, linetype = facet), 
            size = 1) +
  geom_jitter(data = DFfinal1_scaled_completo,
              aes(x = H_real, y = Presence, color = Diel_period),
              height = 0.02, alpha = 0.3, size = 1) +
  scale_linetype_manual(values = c("NATIVAS" = "solid", "TRASLOCADAS" = "dashed")) +
  scale_color_manual(values = c("Day" = "#d95f02", "Night" = "#0072B2")) +
  scale_fill_manual(values = c("Day" = "#d95f02", "Night" = "#0072B2")) +
  labs(x = "Humidity (%)",
       y = "Probability of presence",
       color = "Diel period", fill = "Diel period", linetype = "Status") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

## Final Plot

# Leyend Panel 1: Diel period + Origin status
leyenda1_plot <- ggplot(pred_H, aes(x = x_real, y = predicted)) +
  geom_line(aes(color = group, linetype = facet), linewidth = 1) +
  scale_color_manual(values = c("Day" = "#d95f02", "Night" = "#0072B2"),
                     name = "Diel period",
                     labels = c("Day" = "Daytime", "Night" = "Nighttime")) +
  scale_linetype_manual(values = c("NATIVAS" = "solid", "TRASLOCADAS" = "dashed"),
                        name = "Origin status",
                        labels = c("NATIVAS" = "Resident", "TRASLOCADAS" = "Translocated")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", legend.box = "vertical")

leyenda1 <- get_legend(leyenda1_plot)

# Leyend Panel 2: All individuals
leyenda2_plot <- ggplot() +
  geom_line(data = data.frame(x = 1, y = 1), aes(x = x, y = y, linetype = "All individuals"),
            color = "black", linewidth = 1) +
  scale_linetype_manual(values = c("All individuals" = "solid"), name = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

leyenda2 <- get_legend(leyenda2_plot)

# Leyend Panel 3: Resident vs Translocated
leyenda3_plot <- ggplot(pred_DST, aes(x = x_real, y = predicted)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("NATIVAS" = "#1b9e77", "TRASLOCADAS" = "#7570b3"),
                     name = "Origin status",
                     labels = c("NATIVAS" = "Resident", "TRASLOCADAS" = "Translocated")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

leyenda3 <- get_legend(leyenda3_plot)

# Panel 1
panel1 <- plot_grid(graf_T, graf_H, graf_SU, ncol = 3)
panel1_con_leyenda <- plot_grid(panel1, leyenda1, ncol = 2, rel_widths = c(1, 0.15))

# Panel 2
panel2 <- plot_grid(graf_DWS, graf_DCT, graf_CC, ncol = 3)
panel2_con_leyenda <- plot_grid(panel2, leyenda2, ncol = 2, rel_widths = c(1, 0.15))

# Panel 3
panel3_graficos <- plot_grid(NULL, graf_DST, graf_DFL, NULL, 
                             ncol = 4, 
                             rel_widths = c(0.5, 1, 1, 0.5))

panel3_con_leyenda <- plot_grid(panel3_graficos, leyenda3, 
                                ncol = 2, rel_widths = c(1, 0.15))

### Final plot
plot_grid(panel1_con_leyenda, panel2_con_leyenda, panel3_con_leyenda, 
          ncol = 1)

