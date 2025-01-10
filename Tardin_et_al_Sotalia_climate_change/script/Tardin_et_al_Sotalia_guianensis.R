if(!require(rgdal)) install.packages("rgdal", dependencies = T)
if(!require(raster)) install.packages("raster", dependencies = T)
if(!require(spThin)) install.packages("spThin", dependencies = T)
if(!require(usdm)) install.packages("usdm", dependencies = T)
if(!require(biomod2)) install.packages("biomod2", dependencies = T)
if(!require(dplyr)) install.packages("dplyr", dependencies = T)
if(!require(ggplot2)) install.packages("ggplot2", dependencies = T)
if(!require(ncdf4)) install.packages("ncdf4", dependencies = T)

Species <- "Sotalia_guianensis"

#### Getting data ####

# Range data
RangeData <- readOGR(paste("Input/Range/", "Sotalia_50m_1degree.shp", sep = ""))
plot(RangeData)

#### Environmental variables - Present ####

distcoast <- raster("Input/Layers/Present/Distcoast.tif")
cur_dir_s <- raster("Input/Layers/Present/Cur_dir.nc")
cur_vel_s <- raster("Input/Layers/Present/Cur_vel.nc")
iron_s <- raster("Input/Layers/Present/Iron.nc")
nitrate_s <- raster("Input/Layers/Present/Nitrate.nc")
o2_s <- raster("Input/Layers/Present/o2.nc")
PAR_s <- raster("Input/Layers/Present/PAR.nc")
ph_s <- raster("Input/Layers/Present/ph.nc")
phosphate_s <- raster("Input/Layers/Present/Phosphate.nc")
primprod_s <- raster("Input/Layers/Present/PrimProd.nc")
Slope_s <- raster("Input/Layers/Present/Slope.nc")
sss_s <- raster("Input/Layers/Present/SSS.nc")
sst_s <- raster("Input/Layers/Present/SST.nc")

# Cropping layers
distcoast <- mask(crop(distcoast, RangeData), RangeData)
cur_dir_s <- mask(crop(cur_dir_s, RangeData), RangeData)
cur_vel_s <- mask(crop(cur_vel_s, RangeData), RangeData)
iron_s <- mask(crop(iron_s, RangeData), RangeData)
nitrate_s <- mask(crop(nitrate_s, RangeData), RangeData)
o2_s <- mask(crop(o2_s, RangeData), RangeData)
PAR_s <- mask(crop(PAR_s, RangeData), RangeData)
ph_s <- mask(crop(ph_s, RangeData), RangeData)
phosphate_s <- mask(crop(phosphate_s, RangeData), RangeData)
primprod_s <- mask(crop(primprod_s, RangeData), RangeData)
Slope_s <- mask(crop(Slope_s, RangeData), RangeData)
sss_s <- mask(crop(sss_s, RangeData), RangeData)
sst_s <- mask(crop(sst_s, RangeData), RangeData)

# Resampling layers
distcoast <- resample(distcoast, cur_dir_s, method = "ngb")

# Checking the layers' extension
compareRaster(distcoast,cur_dir_s, cur_vel_s, iron_s, nitrate_s, o2_s, PAR_s, ph_s,phosphate_s,primprod_s,Slope_s,sst_s)

# Merging layers
biostack <- stack(distcoast,cur_dir_s,cur_vel_s,iron_s,nitrate_s,o2_s,ph_s,phosphate_s,primprod_s,Slope_s,sss_s,sst_s)

# Plotting layers
pal <- colorRampPalette(c("yellow", "blue"))
plot(biostack, col = pal(20))

dir.create("Output")
dir.create(paste("Output/", Species, sep = ""))
setwd(paste("Output/", Species, sep = ""))

# Correlation tests

pairs(biostack) # Correlation
biostack_df <- as.data.frame(biostack)
vif(biostack_df) # Multicollinearity
vifstep(biostack_df)
vifcor(biostack_df, th = 0.7)

# Merging layers without collinearity problem
biostack <- stack(distcoast,cur_dir_s,cur_vel_s,primprod_s,Slope_s,sss_s,sst_s)

# Renaming biostack
names(biostack)[1] <- "Distance to coast"
names(biostack)[2] <- "Current direction"
names(biostack)[3] <- "Current speed"
names(biostack)[4] <- "Primary productivity"
names(biostack)[5] <- "Slope"
names(biostack)[6] <- "Salinity"
names(biostack)[7] <- "SST"

png("1. Biostack.png", width = 3200, height = 2400, res = 300); par(oma = c(1,1,1,1)); plot(biostack, col = pal(20)); dev.off()

# Presence data

PresenceData <- read.csv("D:/OneDrive/R/Sotalia_final_climatechange/Output/Thinned/sotcurrent.csv" , sep=";")
head(PresenceData)

# Cropping to the study area

coordinates(PresenceData) <- c("decimallon", "decimallat")
PresenceData <- SpatialPointsDataFrame(coords = PresenceData, data = data.frame(extract(sss_s, PresenceData)))
PresenceData <- subset(PresenceData, extract(sss_s, PresenceData) >= 0)
PresenceData <- as.data.frame(PresenceData)
PresenceData["Species"] <- c(paste(Species))

# Plotting presence data

plot(RangeData)
points(PresenceData$decimallon, PresenceData$decimallat, col="red", cex=1)

# Thinning presence data

dir.create("Thinned")

PresenceDataThinned <- thin( loc.data = PresenceData,
                             lat.col = "decimallat",
                             long.col = "decimallon",
                             spec.col = "Species",
                             thin.par = 15, # 1st approach: 50 / 2nd approach: 30
                             reps = 5,
                             locs.thinned.list.return = TRUE,
                             write.files = TRUE,
                             max.files = 5,
                             out.dir = "Thinned",
                             out.base = "Presence_thinned",
                             write.log.file = TRUE,
                             log.file = "Thinned/Presence_thinned_log_file.txt" )

PresenceData <- read.csv("Thinned/Presence_thinned_thin1.csv" , sep=";")
PresenceData["Occ"] <- c("1")
nrow(PresenceData)

plot(RangeData)
points(PresenceData$decimallon, PresenceData$decimallat, col="red", cex=1)

#### Prepare data & parameters ####

# Select the name of the studied species

myRespName <- "Sotalia_guianensis"

# Get corresponding presence/absence data

myResp <- as.numeric(PresenceData$Occ)
str(myResp)

myRespXY <- PresenceData[which(myResp==1),c("decimallon", "decimallat")]
colnames(myRespXY)

# Pseudo-absences extraction

bm.Species <- BIOMOD_FormatingData(resp.var = myResp,
                                   expl.var = biostack,
                                   resp.xy = myRespXY,
                                   resp.name = myRespName,
                                   PA.nb.rep = 3,
                                   PA.nb.absences = (nrow(PresenceData)*3), # 1st approach: 100 / 2nd approach: (nrow(PresenceData)*3)
                                   PA.strategy = "disk",
                                   PA.dist.min = 30000)

# Summary
bm.Species
sum(bm.Species@PA.table[["PA1"]] == TRUE)-nrow(PresenceData) # n of PAs

# Plotting pseudo-absences
plot(bm.Species)
png("2. Pseudo-absences.png", width = 2400, height = 2400, res = 300); plot(bm.Species); dev.off()

#### Parameterize modeling options ####

# Preparing for Maxent

# Saving explanatory variables in .ascii for MaxEnt algorithm

dir.create("maxent_bg")
maxent.background.dat.dir <- "maxent_bg"

for(var_ in names(biostack)){
  cat("/n> saving", paste0(var_, ".asc"))
  writeRaster(subset(biostack, var_), 
              filename = file.path(maxent.background.dat.dir, paste0(var_, ".asc")),
              overwrite = TRUE)
}

# Defining the path for maxent.jar file 

path.to.maxent.jar <- file.path("D:/OneDrive/Documents/R/win-library/4.1/dismo/java", "maxent.jar")

# Defining Models Options

bm.opt <- BIOMOD_ModelingOptions(GAM = list(k = 4), # 1st and 2nd approach = please remove this argument
                                 MAXENT = list(path_to_maxent.jar = path.to.maxent.jar,
                                               background_data_dir = maxent.background.dat.dir,
                                               maximumbackground = 10000))

save.image(file = (paste("../../", Species, ".RData", sep = "")))
getwd()
#### Running single models ####

{sotaliamodel <- BIOMOD_Modeling(bm.Species,
                                 modeling.id = paste("model_", Species, sep=""),
                                 models = c("MAXENT","GLM", "GAM", "RF", "GBM"),
                                 bm.options = bm.opt,
                                 CV.nb.rep = 10,
                                 CV.strategy = "kfold", # 1st and 2nd approach = please remove this argument
                                 CV.k = 3, # 1st and 2nd approach = please remove this argument
                                 CV.perc = 0.7,
                                 do.full.models = F,
                                 prevalence = 0.5,
                                 var.import = 3,
                                 metric.eval = c("TSS", "ROC"))
beepr::beep(8)}

# Get evaluation scores & variables importance

ModelEvalsot <- get_evaluations(sotaliamodel)
write.csv(ModelEvalsot, "3. ModelEval_single.csv")

VarImportsot <- get_variables_importance(sotaliamodel)
write.csv(VarImportsot, "4. VarImport_single.csv")

# Represent evaluation scores & variables importance

Graph_EvalSinglesot <- bm_PlotEvalMean(bm.out = sotaliamodel, dataset = "validation")
png("5. EvalSingle.png", width = 2400, height = 2400, res = 300); Graph_EvalSinglesot; dev.off()

Graph_BPEvalSinglesot <- bm_PlotEvalBoxplot(bm.out = sotaliamodel, dataset = "validation", group.by = c('algo', 'algo'))
png("6. BoxplotEvalSingle.png", width = 2400, height = 2400, res = 300); Graph_BPEvalSinglesot; dev.off()

ModelEval_TSS <- ModelEvalsot %>% filter(metric.eval == "TSS")
newalgo <- factor(ModelEval_TSS$algo, levels = c("GLM", "GAM", "RF", "GBM", "MAXENT"))

NewGraph_EvalSinglesot <-  ggplot(ModelEval_TSS, aes(x = newalgo, y = validation, fill = algo)) +
  stat_boxplot(geom = "errorbar") + 
  geom_hline(yintercept = 0.7, linetype='longdash', col='red') +
  labs(title = paste(Species), x = "Algorithms", y = "TSS")+
  guides(fill = "none", color = "none") +
  theme_classic() +
  geom_boxplot() +
  ylim(0.15, 1) +
  theme(legend.title = element_blank(),
        title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_EvalSinglesot)

png("6. NewGraph_EvalSingle.png", width = 1800, height = 1800, res = 300); NewGraph_EvalSinglesot; dev.off()

#### Running ensemble model ####

{sot_EM <- BIOMOD_EnsembleModeling(sotaliamodel,
                                            models.chosen = "all",
                                            em.by = "all",
                                            metric.select = c("TSS","ROC"),
                                            em.algo = c("EMmean", "EMwmean", "EMca"),
                                            metric.select.thresh = c(0.7, 0.7),
                                            metric.eval = c("TSS","ROC"),
                                            var.import = 3)
beepr::beep(8)}

# Get evaluation scores & variables importance

ModelEval_EM <- get_evaluations(sot_EMtestepa_run)
write.csv(ModelEval_EM,paste0("7. ModelEval_EMPA+run.csv"))

VarImport_EM <- get_variables_importance(sot_EM)
write.csv(VarImport_EM, "8. VarImport_EM.csv")

# Represent evaluation scores & variables importance

Graph_BoxplotVarAlgo <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('expl.var', 'algo', 'algo'))
png("9. BoxplotVarAlgo.png", width = 2400, height = 2400, res = 300); Graph_BoxplotVarAlgo; dev.off()

Graph_BoxplotAlgoVar <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('algo', 'expl.var', 'algo'))
png("10. BoxplotAlgoVar.png", width = 2400, height = 2400, res = 300); Graph_BoxplotAlgoVar; dev.off()

NewGraph_BoxplotVarImp <-  Graph_BoxplotVarAlgo$tab %>%
  filter(algo == "EMwmean") %>%
  ggplot() +
  geom_boxplot(aes(x = expl.var, y = var.imp, fill = factor(expl.var), color = factor(expl.var))) + scale_fill_manual(values = c("white", "white", "white", "white", "white", "white","white")) +
  scale_color_manual(values = c("red", "blue", "darkgreen", "orange", "cyan3", "brown4" ,"black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1) )+ 
  labs(title = paste(Species), x = "Explanatory variables", y = "Variable importance") +
  guides(fill = "none", color = "none") +
  ylim(0, 1) +
  theme(title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_BoxplotVarImp)

png("11. NewBoxplot_wmean.png", width = 1800, height = 1800, res = 300); NewGraph_BoxplotVarImp; dev.off()

# Represent response curves

Graph_ResponseCurves_wmean <- bm_PlotResponseCurves(bm.out = sot_EM,
                                                    models.chosen = get_built_models(sot_EM)[2],
                                                    fixed.var = 'median')
png("12. ResponseCurves_wmean.png", width = 2400, height = 2400, res = 300); Graph_ResponseCurves_wmean; dev.off()

ggdat <- Graph_ResponseCurves_wmean$tab
ggdat$algo <- ifelse(grepl("_GLM", ggdat$pred.name, ignore.case = T), "GLM",
                     ifelse(grepl("_GAM", ggdat$pred.name, ignore.case = T ), "GAM",
                            ifelse(grepl("_GBM", ggdat$pred.name, ignore.case = T ), "GBM",
                                   ifelse(grepl("_RF", ggdat$pred.name, ignore.case = T ), "RF",
                                          ifelse(grepl("_MAXENT", ggdat$pred.name, ignore.case = T ), "MAXENT", "NA")))))

ggdat$algo <- as.factor(ggdat$algo)

ggdat$id <- as.factor(ggdat$id)

NewResponseCurve <- ggplot(ggdat, aes(x = expl.val, y = pred.val)) +
  geom_smooth(aes(group = algo, color = algo), fill = "indianred1", se = TRUE, method = "gam", alpha = 0.01) +
  geom_smooth(color = "red", fill = "indianred1", se = TRUE, method = "gam") +
  facet_wrap(~ expl.name, scales = "free_x") +
  labs(title = paste(Species), x = "Explanatory variables value", y = "Suitability value") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "",
        title = element_text(size = 16, color = "gray33"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "gray33"),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines"))

plot(NewResponseCurve)

png("13. NewResponseCurves_wmean.png", width = 2400, height = 1800, res = 300); NewResponseCurve; dev.off()

#### Projecting single models ####

{Proj_sot <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                        proj.name = "Current",
                                        new.env = biostack,
                                        models.chosen = "all",
                                        metric.binary = "all",
                                        metric.filter = "all",
                                        build.clamping.mask = F,
                                        keep.in.memory = F)
beepr::beep(8)}

Proj_sot

#### Projecting ensemble model ####

{Proj_sot_EM <- BIOMOD_EnsembleForecasting(bm.em = sot_EM,
                                                    bm.proj = Proj_sot,
                                                    proj.name = "Current_EM",
                                                    models.chosen = "all",
                                                    metric.binary = "all",
                                                    metric.filter = "all",
                                                    output.format = ".img")
beepr::beep(8)}

Proj_sot_EM

#### Future projection - SSP1-2.6 2060-2070 ####

cur_dir_s26 <- raster("Input/Layers/SSP26_2070/Cur_dir.nc")
cur_vel_s26 <- raster("Input/Layers/SSP26_2070/Cur_vel.nc")
primprod_s26 <- raster("Input/Layers/SSP26_2070/PrimProd.nc")
sss_s26 <- raster("Input/Layers/SSP26_2070/SSS.nc")
sst_s26 <- raster("Input/Layers/SSP26_2070/SST.nc")

# Cropping layers
cur_dir_s26 <- mask(crop(cur_dir_s26, RangeData), RangeData)
cur_vel_s26 <- mask(crop(cur_vel_s26, RangeData), RangeData)
primprod_s26 <- mask(crop(primprod_s26, RangeData), RangeData)
sss_s26 <- mask(crop(sss_s26, RangeData), RangeData)
sst_s26 <- mask(crop(sst_s26, RangeData), RangeData)

# Merging layers
biostack_futuro <- stack(distcoast, cur_dir_s26, cur_vel_s26, primprod_s26, Slope_s, sss_s26, sst_s26)
pal <- colorRampPalette(c("yellow", "blue"))
plot(biostack_futuro, col = pal(20))

# Renaming biostack future
names(biostack_futuro)[1] <- "Distance to coast"
names(biostack_futuro)[2] <- "Current direction"
names(biostack_futuro)[3] <- "Current speed"
names(biostack_futuro)[4] <- "Primary productivity"
names(biostack_futuro)[5] <- "Slope"
names(biostack_futuro)[6] <- "Salinity"
names(biostack_futuro)[7] <- "SST"

#### Projecting individual models in geographic space in 2070 SSP1-2.6 - Optimistic scenario ####

project_sot26 <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                   proj.name = 'Mitigated_SSP26_2070',
                                   new.env = biostack_futuro,
                                   models.chosen = get_built_models(sotaliamodel),
                                   metric.binary = c("TSS"),
                                   metric.filter = 'all',
                                   build.clamping.mask = F,keep.in.memory = F, overwrite = T)
project_sot26

#### Projecting ensemble model in geographic space in 2070 SSP1-2.6 - Optimistic scenario ####

project_ensemble_sot26 <- BIOMOD_EnsembleForecasting(
  bm.em=sot_EM,
  bm.proj = project_sot26,
  proj.name = "Mitigated_EM_SSP26_2070",
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = c("ROC", "TSS"),
  metric.filter = c("ROC", "TSS"),
  compress = TRUE, output.format = ".img",
  nb.cpu = 1)

#### Comparing range sizes ####
# Load current and future binary projections

SotCurrentProj26 <- get_predictions(Proj_sot_EM, metric.binary = "TSS")
SotFutureProj26 <- get_predictions(project_ensemble_sot26, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize26 <- BIOMOD_RangeSize(proj.current = SotCurrentProj26, 
                                      proj.future = SotFutureProj26)

myBiomodRangeSize26$Compt.By.Models
plot(myBiomodRangeSize26$Diff.By.Pixel)


# Represent main results 
gg26 = bm_PlotRangeSize(bm.range = myBiomodRangeSize26, 
                      do.count = TRUE,
                      do.perc = TRUE,
                      do.maps = TRUE,
                      do.mean = TRUE,
                      do.plot = TRUE,
                      row.names = c("Species", "Dataset", "Run", "Algo"))
png("Rangechange_SSP1-26.png", width = 2400, height = 1800, res = 300); gg26; dev.off()

#### Future projection - SSP2-4.5 2060-2070 ####

project_sot45 <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                   proj.name = 'Middle_SSP45_2070',
                                   new.env = biostack_futuro,
                                   models.chosen = get_built_models(sotaliamodel),
                                   metric.binary = c("TSS"),
                                   metric.filter = 'all',
                                   build.clamping.mask = F,keep.in.memory = F, overwrite = T)
project_sot45

cur_dir_s45 <- raster("Input/Layers/SSP45_2070/Cur_dir.nc")
cur_vel_s45 <- raster("Input/Layers/SSP45_2070/Cur_vel.nc")
primprod_s45 <- raster("Input/Layers/SSP45_2070/PrimProd.nc")
sss_s45 <- raster("Input/Layers/SSP45_2070/SSS.nc")
sst_s45 <- raster("Input/Layers/SSP45_2070/SST.nc")

# Cropping layers
cur_dir_s45 <- mask(crop(cur_dir_s45, RangeData), RangeData)
cur_vel_s45 <- mask(crop(cur_vel_s45, RangeData), RangeData)
primprod_s45 <- mask(crop(primprod_s45, RangeData), RangeData)
sss_s45 <- mask(crop(sss_s45, RangeData), RangeData)
sst_s45 <- mask(crop(sst_s45, RangeData), RangeData)

# Merging layers
biostack_futuro <- stack(distcoast, cur_dir_s45, cur_vel_s45, primprod_s45, Slope_s, sss_s45, sst_s45)
pal <- colorRampPalette(c("yellow", "blue"))
plot(biostack_futuro, col = pal(20))

# Renaming biostack future
names(biostack_futuro)[1] <- "Distance to coast"
names(biostack_futuro)[2] <- "Current direction"
names(biostack_futuro)[3] <- "Current speed"
names(biostack_futuro)[4] <- "Primary productivity"
names(biostack_futuro)[5] <- "Slope"
names(biostack_futuro)[6] <- "Salinity"
names(biostack_futuro)[7] <- "SST"

#### Projecting individual models in geographic space in 2070 SSP2-4.5 - Intermediate scenario ####

project_sot45 <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                   proj.name = 'Mitigated_SSP45_2070',
                                   new.env = biostack_futuro,
                                   models.chosen = get_built_models(sotaliamodel),
                                   metric.binary = c("TSS"),
                                   metric.filter = 'all',
                                   build.clamping.mask = F,keep.in.memory = F, overwrite = T)
project_sot45

#### Projecting ensemble model in geographic space in 2070 SSP2-4.5 - Intermediate scenario ####

project_ensemble_sot45 <- BIOMOD_EnsembleForecasting(
  bm.em=sot_EM,
  bm.proj = project_sot45,
  proj.name = "Middle_EM_SSP45_2070",
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = c("ROC", "TSS"),
  metric.filter = c("ROC", "TSS"),
  compress = TRUE, output.format = ".img",
  nb.cpu = 1)

#### Comparing range sizes ####
# Load current and future binary projections
SotCurrentProj45 <- get_predictions(Proj_sot_EM, metric.binary = "TSS")
SotFutureProj45 <- get_predictions(project_ensemble_sot45, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize45 <- BIOMOD_RangeSize(proj.current = SotCurrentProj45, 
                                        proj.future = SotFutureProj45)

myBiomodRangeSize45$Compt.By.Models
plot(myBiomodRangeSize45$Diff.By.Pixel)

# Represent main results 
gg45 = bm_PlotRangeSize(bm.range = myBiomodRangeSize45, 
                        do.count = TRUE,
                        do.perc = TRUE,
                        do.maps = TRUE,
                        do.mean = TRUE,
                        do.plot = TRUE,
                        row.names = c("Species", "Dataset", "Run", "Algo"))

#### Future projection - SSP5-8.5 2060-2070 ####

cur_dir_s85 <- raster("Input/Layers/SSP85_2070/Cur_dir.nc")
cur_vel_s85 <- raster("Input/Layers/SSP85_2070/Cur_vel.nc")
primprod_s85 <- raster("Input/Layers/SSP85_2070/PrimProd.nc")
sss_s85 <- raster("Input/Layers/SSP85_2070/SSS.nc")
sst_s85 <- raster("Input/Layers/SSP85_2070/SST.nc")

# Cropping layers
cur_dir_s85 <- mask(crop(cur_dir_s85, RangeData), RangeData)
cur_vel_s85 <- mask(crop(cur_vel_s85, RangeData), RangeData)
primprod_s85 <- mask(crop(primprod_s85, RangeData), RangeData)
sss_s85 <- mask(crop(sss_s85, RangeData), RangeData)
sst_s85 <- mask(crop(sst_s85, RangeData), RangeData)

# Merging layers
biostack_futuro85 <- stack(distcoast, cur_dir_s85, cur_vel_s85,primprod_s85,  Slope_s, sss_s85, sst_s85)
pal <- colorRampPalette(c("yellow", "blue"))
plot(biostack_futuro85, col = pal(20))

# Renaming biostack future
names(biostack_futuro85)[1] <- "Distance to coast"
names(biostack_futuro85)[2] <- "Current direction"
names(biostack_futuro85)[3] <- "Current speed"
names(biostack_futuro85)[4] <- "Primary productivity"
names(biostack_futuro85)[5] <- "Slope"
names(biostack_futuro85)[6] <- "Salinity"
names(biostack_futuro85)[7] <- "SST"

plot(biostack_futuro85, col = pal(20))

#### Projecting individual models in geographic space in 2070 SSP5-8.5 - Pessimistic scenario ####

project_sot85 <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                   proj.name = 'Mitigated_SSP85_2070',
                                   new.env = biostack_futuro85,
                                   models.chosen = get_built_models(sotaliamodel),
                                   metric.binary = c("TSS"),
                                   metric.filter = 'all',
                                   build.clamping.mask = F,keep.in.memory = F, overwrite = T)
project_sot85

#### Projecting ensemble model in geographic space in 2070 SSP5-8.5 - Pessimistic scenario ####

project_ensemble_sot85 <- BIOMOD_EnsembleForecasting(
  bm.em=sot_EM,
  bm.proj = project_sot85,
  proj.name = "BAU_EM_SSP85_2070",
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = c("ROC", "TSS"),
  metric.filter = c("ROC", "TSS"),
  compress = TRUE, output.format = ".img",
  nb.cpu = 1)

#### Comparing range sizes ####
# Load current and future binary projections
SotCurrentProj85 <- get_predictions(Proj_sot_EM, metric.binary = "TSS")
SotFutureProj85 <- get_predictions(project_ensemble_sot85, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize85 <- BIOMOD_RangeSize(proj.current = SotCurrentProj85, 
                                        proj.future = SotFutureProj85)

myBiomodRangeSize85$Compt.By.Models
plot(myBiomodRangeSize85$Diff.By.Pixel)

# Represent main results 
gg85 = bm_PlotRangeSize(bm.range = myBiomodRangeSize85, 
                        do.count = TRUE,
                        do.perc = TRUE,
                        do.maps = TRUE,
                        do.mean = TRUE,
                        do.plot = TRUE,
                        row.names = c("Species", "Dataset", "Run", "Algo"))

save.image(file = (paste("../../", Species, ".RData", sep = "")))
