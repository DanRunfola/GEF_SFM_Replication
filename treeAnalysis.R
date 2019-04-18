library(rpart)
library(doBy)
library(Rcpp)
library(MatchIt)
library(rpart.plot)
library(rgdal)
library(stargazer)
set.seed("1269")
source("causal_tree_functions.R")
sourceCpp("split.cpp")

gefTrtSatData <- read.csv("data/merge_global_treatment.csv")
allGefSfmIds <- read.csv("data/all_sfm_sccf_projects.csv")
round2Dta <- read.csv("data/CombinedGEF.csv", header=FALSE)
round1Dta <- read.csv("data/GlobalEnvironmentFacility_GeocodedResearchRelease_Level1_v1.1/data/level_1a.csv", header=TRUE, stringsAsFactors=FALSE)
round1Anc <- read.csv("data/GlobalEnvironmentFacility_GeocodedResearchRelease_Level1_v1.1/data/projects_ancillary.csv", header=TRUE, stringsAsFactors=FALSE)
round1<- intersect(unlist(round1Anc["gef_id"][[1]]), unlist(allGefSfmIds))
allIDs <- append(unlist(round2Dta["V1"][[1]]), round1)

######################### Load Ancillary Year Data for Treatment Cases
#NOTE: Making assumption that if there is no year value, year = 2013 (as a reasonably conservative assumption)
ancillary_data <- read.csv("./data/combo_ancillary.csv", stringsAsFactors=FALSE, header=TRUE)
ancillary_data <- ancillary_data[unique(ancillary_data$GEF_ID, fromLast=TRUE),]
gefDta <- merge(gefTrtSatData, ancillary_data, by.x="project_id", by.y="GEF_ID", all.x=TRUE)
gefDta[is.na(gefDta$Date_start_CALC),]["Date_start_year"] <- 2013
gefDta[gefDta$Date_start_year == "1900",]["Date_start_year"] <- 2013


######################### Load Control Data and assign random years
control_dta <- read.csv("./data/merge_global_control.csv")
control_dta["Date_start_year"] <- sample(1999:2013, nrow(control_dta), replace=TRUE)
control_dta["Treatment"] <- 0
control_dta$ID <- seq.int(nrow(control_dta))

######################### Construct pre-trend and post-trend for each outcome variable - NTL
gefDta$ID <- seq.int(nrow(gefDta))
gefDta["Treatment"] <- 1
d <- data.frame('pre_year' = gefDta$Date_start_year - 1, hansenpost = '2014', ltdr_postyear = '2015', 'pre_year_B' = gefDta$Date_start_year - 2,ID = gefDta$ID)
d$NTL_prefix = "gefDta$v4composites."
d$NTL_post = ".mean"
NTL_preB_cols <- c('NTL_prefix', 'pre_year_B', 'NTL_post')
NTL_pre_cols <- c('NTL_prefix', 'pre_year', 'NTL_post')

d$LTDR_prefix = "gefDta$ltdr_avhrr_yearly_ndvi."
d$LTDR_post = ".mean"
LTDR_pre_cols <- c('LTDR_prefix', 'pre_year', 'LTDR_post')
LTDR_post_cols <- c('LTDR_prefix', 'ltdr_postyear', 'LTDR_post')

d$hansen_prefix = "gefDta$lossyear.na.categorical_"
d$hansen_post = "/gefDta$lossyear.na.categorical_count"
hansen_pre_cols <- c('hansen_prefix', 'pre_year', 'hansen_post')
hansen_post_cols <- c('hansen_prefix', 'hansenpost', 'hansen_post')
hansen_pre_absolute_cols <- c('hansen_prefix', 'pre_year')

d$pre_NTL_vars = do.call(paste, c(d[ , NTL_pre_cols], list(sep = "")))
d$pre_NTLB_vars = do.call(paste, c(d[ , NTL_preB_cols], list(sep = "")))

d$pre_LTDR_vars = do.call(paste, c(d[ , LTDR_pre_cols], list(sep = "")))
d$post_LTDR_vars = do.call(paste, c(d[ , LTDR_post_cols], list(sep = "")))

d$pre_hansen_vars = do.call(paste, c(d[ , hansen_pre_cols], list(sep = "")))
d$post_hansen_vars = do.call(paste, c(d[ , hansen_post_cols], list(sep = ""))) 
d$pre_hansen_absolute_vars = do.call(paste, c(d[ , hansen_pre_absolute_cols], list(sep = "")))

d$preNTLA <- eval(parse(text=d$pre_NTL_vars))
d$preNTLB <- eval(parse(text=d$pre_NTLB_vars))
d$preNTL = d$preNTLA - d$preNTLB
d$postNTL <- gefDta$viirs_ntl_yearly.2015.mean - gefDta$viirs_ntl_yearly.2014.mean

d$preLTDR <- eval(parse(text=d$pre_LTDR_vars))
d$postLTDR <- eval(parse(text=d$post_LTDR_vars))

d$prehansen <- eval(parse(text=d$pre_hansen_vars))
d$posthansen <- eval(parse(text=d$post_hansen_vars))
d$prehansen_absolute <- eval(parse(text=d$pre_hansen_absolute_vars))

gefDta_m <- merge(gefDta, d, by="ID")

#Remove NA values
gefDta_m <- gefDta_m[!is.na(gefDta_m$preNTL),]
gefDta_m <- gefDta_m[!is.na(gefDta_m$postNTL),]
gefDta_m <- gefDta_m[!is.na(gefDta_m$preLTDR),]
gefDta_m <- gefDta_m[!is.na(gefDta_m$postLTDR),]
gefDta_m <- gefDta_m[!is.na(gefDta_m$prehansen),]
gefDta_m <- gefDta_m[!is.na(gefDta_m$posthansen),]


#Control Cases
d <- data.frame('pre_year' = control_dta$Date_start_year - 1, hansenpost = '2014', ltdr_postyear = '2015', 'pre_year_B' = control_dta$Date_start_year - 2,ID = control_dta$ID)
d$NTL_prefix = "control_dta$v4composites."
d$NTL_post = ".mean"
NTL_preB_cols <- c('NTL_prefix', 'pre_year_B', 'NTL_post')
NTL_pre_cols <- c('NTL_prefix', 'pre_year', 'NTL_post')

d$LTDR_prefix = "control_dta$ltdr_avhrr_yearly_ndvi."
d$LTDR_post = ".mean"
LTDR_pre_cols <- c('LTDR_prefix', 'pre_year', 'LTDR_post')
LTDR_post_cols <- c('LTDR_prefix', 'ltdr_postyear', 'LTDR_post')

d$hansen_prefix = "control_dta$lossyear.na.categorical_"
d$hansen_post = "/control_dta$lossyear.na.categorical_count"
hansen_pre_cols <- c('hansen_prefix', 'pre_year', 'hansen_post')
hansen_post_cols <- c('hansen_prefix', 'hansenpost', 'hansen_post')
hansen_pre_absolute_cols <- c('hansen_prefix', 'pre_year')

d$pre_NTL_vars = do.call(paste, c(d[ , NTL_pre_cols], list(sep = "")))
d$pre_NTLB_vars = do.call(paste, c(d[ , NTL_preB_cols], list(sep = "")))

d$pre_LTDR_vars = do.call(paste, c(d[ , LTDR_pre_cols], list(sep = "")))
d$post_LTDR_vars = do.call(paste, c(d[ , LTDR_post_cols], list(sep = "")))

d$pre_hansen_vars = do.call(paste, c(d[ , hansen_pre_cols], list(sep = "")))
d$post_hansen_vars = do.call(paste, c(d[ , hansen_post_cols], list(sep = ""))) 
d$pre_hansen_absolute_vars = do.call(paste, c(d[ , hansen_pre_absolute_cols], list(sep = "")))

d$preNTLA <- eval(parse(text=d$pre_NTL_vars))
d$preNTLB <- eval(parse(text=d$pre_NTLB_vars))
d$preNTL = d$preNTLA - d$preNTLB
d$postNTL <- control_dta$viirs_ntl_yearly.2015.mean - control_dta$viirs_ntl_yearly.2014.mean

d$preLTDR <- eval(parse(text=d$pre_LTDR_vars))
d$postLTDR <- eval(parse(text=d$post_LTDR_vars))

d$prehansen <- eval(parse(text=d$pre_hansen_vars))
d$posthansen <- eval(parse(text=d$post_hansen_vars))
d$prehansen_absolute <- eval(parse(text=d$pre_hansen_absolute_vars))
ctrl_back <- control_dta
control_dta <- merge(control_dta, d, by="ID")

#Remove NA values
control_dta <- control_dta[!is.na(control_dta$preNTL),]
control_dta <- control_dta[!is.na(control_dta$postNTL),]
control_dta <- control_dta[!is.na(control_dta$preLTDR),]
control_dta <- control_dta[!is.na(control_dta$postLTDR),]
control_dta <- control_dta[!is.na(control_dta$prehansen),]
control_dta <- control_dta[!is.na(control_dta$posthansen),]

#Add empty location ID to facilitate merging later
control_dta["project_location_id"] <- 1
#Merge treatment and control into common dataframe
names_common <- intersect(names(control_dta), names(gefDta_m))

tDta <- rbind(control_dta[,names_common], gefDta_m[,names_common])

tDta$years_implemented <- 2015 - tDta$Date_start_year
tDta <- tDta[tDta$years_implemented > 0,]

#Choose pscore variables
pVars <- c("prehansen_absolute", "prehansen",
           "preLTDR", "preNTL", 
           "distance_to_drugdata_201708.na.mean",
           "diamond_distance_201708.na.mean",
           "dist_to_onshore_petroleum_v12.na.mean",
           "distance_to_coast.na.mean",
           "distance_to_gold_v12.na.mean",
           "distance_to_gemdata_201708.na.mean",
           "dist_to_groads.na.mean",
           "dist_to_water.na.mean",
           "dist_to_gadm28_borders.na.mean",
           "wdpa_iucn.na.categorical_mix",
           "gdp_grid.na.mean",
           "udel_precip_v4_01_yearly_max.2000.mean",
           "udel_air_temp_v4_01_yearly_max.2000.mean" ,
           "udel_air_temp_v4_01_yearly_mean.2000.mean" ,
           "udel_precip_v4_01_yearly_min.2000.mean" ,
           "modis_lst_day_yearly_mean.2001.mean"  ,
           "srtm_slope_500m.na.mean"  ,
           "srtm_elevation_500m.na.mean"   ,
           "air_pollution_o3.2000.mean",
           "gpw_v3_density.2000.mean" ,
           "gpw_v3_count.2000.sum",
           "accessibility_map.na.mean",
           "Treatment"
           )

#Choose analysis variables
aVars <- c("years_implemented", "Date_start_year", 
           "prehansen_absolute", "prehansen",
           "preLTDR", "preNTL", 
           "distance_to_drugdata_201708.na.mean",
           "diamond_distance_201708.na.mean",
           "dist_to_onshore_petroleum_v12.na.mean",
           "distance_to_coast.na.mean",
           "distance_to_gold_v12.na.mean",
           "distance_to_gemdata_201708.na.mean",
           "dist_to_groads.na.mean",
           "dist_to_water.na.mean",
           "dist_to_gadm28_borders.na.mean",
           "wdpa_iucn.na.categorical_mix",
           "gdp_grid.na.mean",
           "udel_precip_v4_01_yearly_max.2000.mean",
           "udel_air_temp_v4_01_yearly_max.2000.mean" ,
           "udel_air_temp_v4_01_yearly_mean.2000.mean" ,
           "udel_precip_v4_01_yearly_min.2000.mean" ,
           "modis_lst_day_yearly_mean.2001.mean"  ,
           "srtm_slope_500m.na.mean"  ,
           "srtm_elevation_500m.na.mean"   ,
           "air_pollution_o3.2000.mean",
           "gpw_v3_density.2000.mean" ,
           "gpw_v3_count.2000.sum",
           "accessibility_map.na.mean", "Treatment")

pscore.Calc <- matchit(Treatment ~ prehansen_absolute + prehansen +
                       preLTDR + preNTL +
                       distance_to_drugdata_201708.na.mean +
                       diamond_distance_201708.na.mean+
                       dist_to_onshore_petroleum_v12.na.mean+
                       distance_to_coast.na.mean+
                       distance_to_gold_v12.na.mean+
                       distance_to_gemdata_201708.na.mean+
                       dist_to_groads.na.mean+
                       dist_to_water.na.mean+
                       dist_to_gadm28_borders.na.mean+
                       wdpa_iucn.na.categorical_mix+
                       gdp_grid.na.mean+
                       udel_precip_v4_01_yearly_max.2000.mean+
                       udel_air_temp_v4_01_yearly_max.2000.mean+
                       udel_air_temp_v4_01_yearly_mean.2000.mean+
                       udel_precip_v4_01_yearly_min.2000.mean+
                       modis_lst_day_yearly_mean.2001.mean+
                       srtm_slope_500m.na.mean+
                       srtm_elevation_500m.na.mean+
                       air_pollution_o3.2000.mean+
                       gpw_v3_density.2000.mean+
                       gpw_v3_count.2000.sum+
                       accessibility_map.na.mean,
                       data= na.omit(tDta[pVars], cols=pVars),
                       method="nearest", distance="logit")

covariate_labs <- c("Pre-Implementation Deforestation Rate (Absolute)", 
                    "Pre-Implementation Deforestation Rate (Relative)",
           "Pre-Implementation NDVI", "Pre-Implementation NTL", 
           "Distance to Drug Cultivation Site",
           "Distance to Diamond Resources",
           "Distance to Petroleum Resources",
           "Distance to Coast",
           "Distance to Gold Deposit",
           "Distance to Any Gem Deposit",
           "Distance to Roads",
           "Distance to Water",
           "Distance to Country Border",
           "WDPA Inclusion",
           "GDP Estimate",
           "Pre-Implementation Maximum Yearly Precipitation",
           "Pre-Implementation Maximum Yearly Air Temperature" ,
           "Pre-Implementation Mean Yearly Air Temperature",
           "Pre-Implementation Minimum Yearly Precipitation",
           "Pre-Implementation Land Surface Temperature",
           "Slope",
           "Elevation",
           "Ozone Concentration (ug/m3)",
           "Pre-Implementation Population Density" ,
           "Pre-Implementation Absolute Population",
           "Travel Time to Urban Areas"
)
stargazer(pscore.Calc$model, style="qje", 
          covariate.labels=covariate_labs,
          digits.extra=10,
          font.size = "footnotesize",
          no.space = TRUE,
          out="propensity.html")

matched.dta <- match.data(pscore.Calc)


matched.dta.ntl <- merge(matched.dta, tDta[c("postNTL")], by="row.names")


  ###########################################
  ###########################################
  ###########################################
  ###########################################
  ###########################################
  #Tree - Hansen
  ###########################################
  ###########################################
  ###########################################
  ###########################################
  ###########################################
  matched.dtaB <- matched.dta[as.numeric(matched.dta$distance) < 0.99,]
  matched.dtaB <- matched.dtaB[as.numeric(matched.dta$distance) > 0.01,]
  tDta$hansenOutcome <- tDta$posthansen - tDta$prehansen
  matched.dta.hansen <- merge(matched.dtaB, tDta[c("hansenOutcome")], by="row.names")
  
  transOutcome <- list(rep(0,nrow(matched.dta.hansen)))
  
  for(i in 1:nrow(matched.dta.hansen))
  {
    if(matched.dta.hansen$Treatment[i] == 1)
    {
      #Treated
      transOutcome[i] = matched.dta.hansen$hansenOutcome[i] * 
        (1 / matched.dta.hansen$distance[i])
    }
    else
    {
      #Untreated
      transOutcome[i] = -1 * (matched.dta.hansen$hansenOutcome[i] * 
                                ((1-0) / (1 - matched.dta.hansen$distance[i])))
    }
  }
  matched.dta.hansen$transOutcome <- unlist(transOutcome)
  
  alist <- list(eval=ctev, split=ctsplit, init=ctinit)
  dbb = matched.dta.hansen
  k = 10
  n = dim(dbb)[1]
  crxvdata = dbb
  crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
  list = 1:k
  m.split = 200
  
  errset = list()
  
  
  
  for (i in 1:k){
    errset[[i]] = list()
    trainingset <- subset(crxvdata, id %in% list[-i])
    sub.fit = rpart(cbind(trainingset$hansenOutcome,trainingset$Treatment,
                          trainingset$distance,trainingset$transOutcome) ~
                      prehansen_absolute + prehansen +
                      preLTDR + preNTL +
                      distance_to_coast.na.mean+
                      dist_to_groads.na.mean+
                      dist_to_water.na.mean+
                      dist_to_gadm28_borders.na.mean+
                      wdpa_iucn.na.categorical_mix+
                      udel_precip_v4_01_yearly_max.2000.mean+
                      udel_air_temp_v4_01_yearly_max.2000.mean+
                      udel_air_temp_v4_01_yearly_mean.2000.mean+
                      udel_precip_v4_01_yearly_min.2000.mean+
                      modis_lst_day_yearly_mean.2001.mean+
                      srtm_slope_500m.na.mean+
                      srtm_elevation_500m.na.mean+
                      air_pollution_o3.2000.mean+
                      gpw_v3_density.2000.mean+
                      gpw_v3_count.2000.sum+
                      accessibility_map.na.mean,
                    trainingset,
                    control = rpart.control(cp = 0,minsplit = m.split),
                    method=alist)
    sub.fit.dm = data.matrix(sub.fit$frame)
    index = as.numeric(rownames(sub.fit$frame))
    removed_nodes = 0
    #fit1$frame$var = as.numeric(fit1$frame$var)
    removed_nodes = cross_validate(sub.fit.dm, index,removed_nodes)
    removed_nodes = removed_nodes[-1]
    #errset[i] = rep(0,length(removed_nodes))
    for(l in 1:length(removed_nodes)){
      error = 0
      sub.fit.pred = snip.rpart(sub.fit, removed_nodes[1:l])
      
      #Subset Fit
      testset <- subset(crxvdata, id %in% c(i))
      pt = predict(sub.fit.pred,testset,type = "matrix")
      y = data.frame(pt)
      val = data.matrix(y)
      idx = as.numeric(rownames(testset))
      dbidx = as.numeric(rownames(dbb))
      
      for(pid in 1:(dim(y)[1])){
        id = match(idx[pid],dbidx)
        error = error + (dbb$transOutcome[id] - val[pid])^2
      }
      
      if(error == 0){
        errset[[i]][l] = 1000000
      }
      else{
        errset[[i]][l] = error/k
      }
    }
  }
  
  #Identify the average error to depth ratio across all cross-validations
  avg.index <- vector()
  for(e in 1:length(errset))
  {
    avg.index[e] <- which.min(errset[[e]])
  }
  
  #---------------
  #Build Final Tree
  #---------------
  fit1 = rpart(cbind(crxvdata$hansenOutcome,crxvdata$Treatment,
                     crxvdata$distance,crxvdata$transOutcome) ~
                 prehansen_absolute + prehansen +
                 preLTDR + preNTL +
                 distance_to_coast.na.mean+
                 dist_to_groads.na.mean+
                 dist_to_water.na.mean+
                 dist_to_gadm28_borders.na.mean+
                 wdpa_iucn.na.categorical_mix+
                 udel_precip_v4_01_yearly_max.2000.mean+
                 udel_air_temp_v4_01_yearly_max.2000.mean+
                 udel_air_temp_v4_01_yearly_mean.2000.mean+
                 udel_precip_v4_01_yearly_min.2000.mean+
                 modis_lst_day_yearly_mean.2001.mean+
                 srtm_slope_500m.na.mean+
                 srtm_elevation_500m.na.mean+
                 air_pollution_o3.2000.mean+
                 gpw_v3_density.2000.mean+
                 gpw_v3_count.2000.sum+
                 accessibility_map.na.mean,
               crxvdata,
               control = rpart.control(cp = 0,minsplit = m.split),
               method=alist)
  fit = data.matrix(fit1$frame)
  index = as.numeric(rownames(fit1$frame))
  
  
  removed_nodes = 0
  removed_nodes = cross_validate(fit, index,removed_nodes)
  removed_nodes = removed_nodes[-1]
  pruned_nodes = removed_nodes[1:round(mean(avg.index))]
  final.tree <- snip.rpart(fit1, pruned_nodes)
  
  #Prep for output
  print.tree <- final.tree
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen"] <- "Base Defor. (Rel)"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_count.2000.sum"] <- "Absolute Population"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_precip_v4_01_yearly_max.2000.mean"] <- "Max Precip"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="air_pollution_o3.2000.mean"] <- "Ozone Concentration"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="preLTDR"] <- "Vegetation Density"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_density.2000.mean"] <- "Population Density"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_elevation_500m.na.mean"] <- "Slope"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen_absolute"] <- "Base Defor. (Abs)"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="distance_to_coast.na.mean"] <- "Dist. to Coast"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_precip_v4_01_yearly_min.2000.mean"] <- "Mean Precip"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_slope_500m.na.mean"] <- "Slope"
  levels(print.tree $frame$var)[levels(print.tree $frame$var)=="dist_to_groads.na.mean"] <- "Dist. to Roads"
  
  
  
  
  png("./GEF_LandCover.png", width = 1280, height = 720, res=300)
  rpart.plot(print.tree , cex=0.3, extra=1, branch=1, type=4, tweak=0.9, clip.right.labs=FALSE,
             box.col=c("palegreen3", "pink")[findInterval(print.tree $frame$yval, v = c(-1,0))],
             faclen=0,
             varlen=0
  )
  dev.off()


tDta$pred_hansen <- predict(final.tree, tDta)


###########################################
###########################################
###########################################
###########################################
###########################################
#Tree - LTDR
###########################################
###########################################
###########################################
###########################################
###########################################
matched.dtaB <- matched.dta[as.numeric(matched.dta$distance) < 0.99,]
matched.dtaB <- matched.dtaB[as.numeric(matched.dta$distance) > 0.01,]

tDta$LTDRoutcome <- tDta$postLTDR - tDta$preLTDR
matched.dta.ltdr <- merge(matched.dtaB, tDta[c("LTDRoutcome")], by="row.names")
matched.dta.ltdr$LTDRoutcome <- matched.dta.ltdr$LTDRoutcome / 10000
transOutcome <- list(rep(0,nrow(matched.dta.ltdr)))

for(i in 1:nrow(matched.dta.ltdr))
{
  if(matched.dta.ltdr$Treatment[i] == 1)
  {
    #Treated
    transOutcome[i] = matched.dta.ltdr$LTDRoutcome[i] * 
      (1 / matched.dta.ltdr$distance[i])
  }
  else
  {
    #Untreated
    transOutcome[i] = -1 * (matched.dta.ltdr$LTDRoutcome[i] * 
                              ((1-0) / (1 - matched.dta.ltdr$distance[i])))
  }
}
matched.dta.ltdr$transOutcome <- unlist(transOutcome)

###
tDtab <- tDta
tDtab$distance <- predict(pscore.Calc$model, tDtab, type="response")
tDtab <- tDtab[as.numeric(tDtab$distance) < 0.8,] #.9, -381 #.8 -379
tDtab <- tDtab[as.numeric(tDtab$distance) > 0.2,]
tDtab[is.na(tDtab$Treatment),] <- 1
transOutcome <- list(rep(0,nrow(tDta)))

for(i in 1:nrow(tDtab))
{
  if(tDtab$Treatment[i] == 1)
  {
    #Treated
    transOutcome[i] = tDtab$LTDRoutcome[i] * 
      (1 / tDtab$distance[i])
  }
  else
  {
    #Untreated
    transOutcome[i] = -1 * (tDtab$LTDRoutcome[i] * 
                              ((1-0) / (1 - tDtab$distance[i])))
  }
}
tDtab$transOutcome <- unlist(transOutcome)
###

alist <- list(eval=ctev, split=ctsplit, init=ctinit)
dbb = matched.dta.ltdr
k = 10
n = dim(dbb)[1]
crxvdata = dbb
crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
list = 1:k
m.split = 275

errset = list()



for (i in 1:k){
  errset[[i]] = list()
  trainingset <- subset(crxvdata, id %in% list[-i])
  sub.fit = rpart(cbind(trainingset$LTDRoutcome,trainingset$Treatment,
                        trainingset$distance,trainingset$transOutcome) ~
                    prehansen_absolute + prehansen +
                    preLTDR + preNTL +
                    distance_to_coast.na.mean+
                    dist_to_groads.na.mean+
                    dist_to_water.na.mean+
                    dist_to_gadm28_borders.na.mean+
                    wdpa_iucn.na.categorical_mix+
                    udel_precip_v4_01_yearly_max.2000.mean+
                    udel_air_temp_v4_01_yearly_max.2000.mean+
                    udel_air_temp_v4_01_yearly_mean.2000.mean+
                    udel_precip_v4_01_yearly_min.2000.mean+
                    modis_lst_day_yearly_mean.2001.mean+
                    srtm_slope_500m.na.mean+
                    srtm_elevation_500m.na.mean+
                    gpw_v3_density.2000.mean+
                    gpw_v3_count.2000.sum+
                    accessibility_map.na.mean,
                  trainingset,
                  control = rpart.control(cp = 0,minsplit = m.split),
                  method=alist)
  sub.fit.dm = data.matrix(sub.fit$frame)
  index = as.numeric(rownames(sub.fit$frame))
  removed_nodes = 0
  #fit1$frame$var = as.numeric(fit1$frame$var)
  removed_nodes = cross_validate(sub.fit.dm, index,removed_nodes)
  removed_nodes = removed_nodes[-1]
  #errset[i] = rep(0,length(removed_nodes))
  for(l in 1:length(removed_nodes)){
    error = 0
    sub.fit.pred = snip.rpart(sub.fit, removed_nodes[1:l])
    
    #Subset Fit
    testset <- subset(crxvdata, id %in% c(i))
    pt = predict(sub.fit.pred,testset,type = "matrix")
    y = data.frame(pt)
    val = data.matrix(y)
    idx = as.numeric(rownames(testset))
    dbidx = as.numeric(rownames(dbb))
    
    for(pid in 1:(dim(y)[1])){
      id = match(idx[pid],dbidx)
      error = error + (dbb$transOutcome[id] - val[pid])^2
    }
    
    if(error == 0){
      errset[[i]][l] = 1000000
    }
    else{
      errset[[i]][l] = error/k
    }
  }
}

#Identify the average error to depth ratio across all cross-validations
avg.index <- vector()
for(e in 1:length(errset))
{
  avg.index[e] <- which.min(errset[[e]])
}

#---------------
#Build Final Tree
#---------------
fit1 = rpart(cbind(crxvdata$LTDRoutcome,crxvdata$Treatment,
                   crxvdata$distance,crxvdata$transOutcome) ~
               prehansen_absolute + prehansen +
               preLTDR + preNTL +
               distance_to_coast.na.mean+
               dist_to_groads.na.mean+
               dist_to_water.na.mean+
               dist_to_gadm28_borders.na.mean+
               wdpa_iucn.na.categorical_mix+
               udel_precip_v4_01_yearly_max.2000.mean+
               udel_air_temp_v4_01_yearly_max.2000.mean+
               udel_air_temp_v4_01_yearly_mean.2000.mean+
               udel_precip_v4_01_yearly_min.2000.mean+
               modis_lst_day_yearly_mean.2001.mean+
               srtm_slope_500m.na.mean+
               srtm_elevation_500m.na.mean+
               gpw_v3_density.2000.mean+
               gpw_v3_count.2000.sum+
               accessibility_map.na.mean,
             crxvdata,
             control = rpart.control(cp = 0,minsplit = m.split),
             method=alist)
fit = data.matrix(fit1$frame)
index = as.numeric(rownames(fit1$frame))


removed_nodes = 0
removed_nodes = cross_validate(fit, index,removed_nodes)
removed_nodes = removed_nodes[-1]
pruned_nodes = removed_nodes[1:round(mean(avg.index))]
final.tree <- snip.rpart(fit1, pruned_nodes)

#Prep for output
print.tree <- final.tree
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen"] <- "Base Defor."
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_count.2000.sum"] <- "Absolute Population"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_precip_v4_01_yearly_max.2000.mean"] <- "Max Precip"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="air_pollution_o3.2000.mean"] <- "Ozone Concentration"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="preLTDR"] <- "Vegetation Density"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_density.2000.mean"] <- "Population Density"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_elevation_500m.na.mean"] <- "Elevation"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_slope_500m.na.mean"] <- "Slope"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_air_temp_v4_01_yearly_mean.2000.mean"] <- "Average Precip"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="modis_lst_day_yearly_mean.2001.mean"] <- "Land Surf. Temp."

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen_absolute"] <- "Init. Abs. Defor."

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="dist_to_gadm28_borders.na.mean"] <- "Dist. to Borders"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="accessibility_map.na.mean"] <- "Travel Time to Urban"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="distance_to_coast.na.mean"] <- "Dist. to Coast"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen"] <- "Base Defor."

png("./GEF_LTDR.png", width = 1280, height = 720, res=300)
rpart.plot(print.tree , cex=0.3, extra=1, branch=1, type=4, tweak=0.9, clip.right.labs=FALSE,
           box.col=c("pink", "palegreen")[findInterval(print.tree $frame$yval, v = c(-10000,0))],
           faclen=0,
           varlen=0
)
dev.off()


tDta$pred_LTDR <- predict(final.tree, tDta)

###########################################
###########################################
###########################################
###########################################
###########################################
#Tree - NTL
###########################################
###########################################
###########################################
###########################################
###########################################
matched.dtaB <- matched.dta[as.numeric(matched.dta$distance) < 0.99,]
matched.dtaB <- matched.dtaB[as.numeric(matched.dta$distance) > 0.01,]

tDta$NTLoutcome <- tDta$postNTL - tDta$preNTL
matched.dta.ltdr <- merge(matched.dtaB, tDta[c("NTLoutcome")], by="row.names")

transOutcome <- list(rep(0,nrow(matched.dta.ltdr)))

for(i in 1:nrow(matched.dta.ltdr))
{
  if(matched.dta.ltdr$Treatment[i] == 1)
  {
    #Treated
    transOutcome[i] = matched.dta.ltdr$NTLoutcome[i] * 
      (1 / matched.dta.ltdr$distance[i])
  }
  else
  {
    #Untreated
    transOutcome[i] = -1 * (matched.dta.ltdr$NTLoutcome[i] * 
                              ((1-0) / (1 - matched.dta.ltdr$distance[i])))
  }
}
matched.dta.ltdr$transOutcome <- unlist(transOutcome)

###
tDtab <- tDta
tDtab$distance <- predict(pscore.Calc$model, tDtab, type="response")
tDtab <- tDtab[as.numeric(tDtab$distance) < 0.99,] #.9, -381 #.8 -379
tDtab <- tDtab[as.numeric(tDtab$distance) > 0.2,]
tDtab[is.na(tDtab$Treatment),] <- 1
transOutcome <- list(rep(0,nrow(tDta)))

for(i in 1:nrow(tDtab))
{
  if(tDtab$Treatment[i] == 1)
  {
    #Treated
    transOutcome[i] = tDtab$NTLoutcome[i] * 
      (1 / tDtab$distance[i])
  }
  else
  {
    #Untreated
    transOutcome[i] = -1 * (tDtab$NTLoutcome[i] * 
                              ((1-0) / (1 - tDtab$distance[i])))
  }
}
tDtab$transOutcome <- unlist(transOutcome)
###

alist <- list(eval=ctev, split=ctsplit, init=ctinit)
dbb = matched.dta.ltdr
k = 10
n = dim(dbb)[1]
crxvdata = dbb
crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
list = 1:k
m.split = 250

errset = list()



for (i in 1:k){
  errset[[i]] = list()
  trainingset <- subset(crxvdata, id %in% list[-i])
  sub.fit = rpart(cbind(trainingset$NTLoutcome,trainingset$Treatment,
                        trainingset$distance,trainingset$transOutcome) ~
                    prehansen_absolute + prehansen +
                    preLTDR + preNTL +
                    dist_to_groads.na.mean+
                    dist_to_water.na.mean+
                    dist_to_gadm28_borders.na.mean+
                    wdpa_iucn.na.categorical_mix+
                    udel_precip_v4_01_yearly_max.2000.mean+
                    udel_air_temp_v4_01_yearly_max.2000.mean+
                    udel_air_temp_v4_01_yearly_mean.2000.mean+
                    udel_precip_v4_01_yearly_min.2000.mean+
                    modis_lst_day_yearly_mean.2001.mean+
                    srtm_slope_500m.na.mean+
                    srtm_elevation_500m.na.mean+
                    gpw_v3_density.2000.mean+
                    gpw_v3_count.2000.sum+
                    accessibility_map.na.mean,
                  trainingset,
                  control = rpart.control(cp = 0,minsplit = m.split),
                  method=alist)
  sub.fit.dm = data.matrix(sub.fit$frame)
  index = as.numeric(rownames(sub.fit$frame))
  removed_nodes = 0
  #fit1$frame$var = as.numeric(fit1$frame$var)
  removed_nodes = cross_validate(sub.fit.dm, index,removed_nodes)
  removed_nodes = removed_nodes[-1]
  #errset[i] = rep(0,length(removed_nodes))
  for(l in 1:length(removed_nodes)){
    error = 0
    sub.fit.pred = snip.rpart(sub.fit, removed_nodes[1:l])
    
    #Subset Fit
    testset <- subset(crxvdata, id %in% c(i))
    pt = predict(sub.fit.pred,testset,type = "matrix")
    y = data.frame(pt)
    val = data.matrix(y)
    idx = as.numeric(rownames(testset))
    dbidx = as.numeric(rownames(dbb))
    
    for(pid in 1:(dim(y)[1])){
      id = match(idx[pid],dbidx)
      error = error + (dbb$transOutcome[id] - val[pid])^2
    }
    
    if(error == 0){
      errset[[i]][l] = 1000000
    }
    else{
      errset[[i]][l] = error/k
    }
  }
}

#Identify the average error to depth ratio across all cross-validations
avg.index <- vector()
for(e in 1:length(errset))
{
  avg.index[e] <- which.min(errset[[e]])
}

#---------------
#Build Final Tree
#---------------
fit1 = rpart(cbind(crxvdata$NTLoutcome,crxvdata$Treatment,
                   crxvdata$distance,crxvdata$transOutcome) ~
               prehansen_absolute + prehansen +
               preLTDR + preNTL +
               dist_to_groads.na.mean+
               dist_to_water.na.mean+
               dist_to_gadm28_borders.na.mean+
               wdpa_iucn.na.categorical_mix+
               udel_precip_v4_01_yearly_max.2000.mean+
               udel_air_temp_v4_01_yearly_max.2000.mean+
               udel_air_temp_v4_01_yearly_mean.2000.mean+
               udel_precip_v4_01_yearly_min.2000.mean+
               modis_lst_day_yearly_mean.2001.mean+
               srtm_slope_500m.na.mean+
               srtm_elevation_500m.na.mean+
               gpw_v3_density.2000.mean+
               gpw_v3_count.2000.sum+
               accessibility_map.na.mean,
             crxvdata,
             control = rpart.control(cp = 0,minsplit = m.split),
             method=alist)
fit = data.matrix(fit1$frame)
index = as.numeric(rownames(fit1$frame))


removed_nodes = 0
removed_nodes = cross_validate(fit, index,removed_nodes)
removed_nodes = removed_nodes[-1]
pruned_nodes = removed_nodes[1:round(mean(avg.index))]
final.tree <- snip.rpart(fit1, pruned_nodes)

#Prep for output
print.tree <- final.tree
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen"] <- "Base Defor."
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_count.2000.sum"] <- "Absolute Population"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_precip_v4_01_yearly_max.2000.mean"] <- "Max Precip"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="air_pollution_o3.2000.mean"] <- "Ozone Concentration"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="preLTDR"] <- "Vegetation Density"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="gpw_v3_density.2000.mean"] <- "Population Density"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_elevation_500m.na.mean"] <- "Elevation"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="srtm_slope_500m.na.mean"] <- "Slope"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="udel_air_temp_v4_01_yearly_mean.2000.mean"] <- "Average Precip"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="modis_lst_day_yearly_mean.2001.mean"] <- "Land Surf. Temp."

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="prehansen_absolute"] <- "Init. Abs. Defor."

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="dist_to_gadm28_borders.na.mean"] <- "Dist. to Borders"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="accessibility_map.na.mean"] <- "Travel Time to Urban"

levels(print.tree $frame$var)[levels(print.tree $frame$var)=="distance_to_coast.na.mean"] <- "Dist. to Coast"


levels(print.tree $frame$var)[levels(print.tree $frame$var)=="dist_to_water.na.mean"] <- "Dist. to Water"
levels(print.tree $frame$var)[levels(print.tree $frame$var)=="dist_to_groads.na.mean"] <- "Dist. to Roads"

png("./GEF_NTL.png", width = 1280, height = 720, res=300)
rpart.plot(print.tree , cex=0.3, extra=1, branch=1, type=4, tweak=0.9, clip.right.labs=FALSE,
           box.col=c("pink", "palegreen")[findInterval(print.tree $frame$yval, v = c(-10000,0))],
           faclen=0,
           varlen=0
)
dev.off()


tDta$pred_NTL <- predict(final.tree, tDta)

write.csv(tDta, "map_results.csv")
