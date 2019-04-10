library(rpart)
library(doBy)
library(Rcpp)
library(MatchIt)
library(rpart.plot)
library(rgdal)

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
control_dta["Date_start_year"] <- runif(nrow(control_dta),1999,2013)


######################### Construct pre-trend and post-trend for each outcome variable - NTL
gefDta$ID <- seq.int(nrow(gefDta))

d <- data.frame('pre_year' = gefDta$Date_start_year - 1, ID = gefDta$ID, 'post_year' = gefDta$Date_start_year)
d$NTL_prefix = "gefDta$v4composites."
d$NTL_post = ".mean"
NTL_pre_cols <- c('NTL_prefix', 'pre_year', 'NTL_post')
NTL_post_cols <- c('NTL_prefix', 'post_year', 'NTL_post')

d$LTDR_prefix = "gefDta$ltdr_avhrr_yearly_ndvi."
d$LTDR_post = ".mean"
LTDR_pre_cols <- c('LTDR_prefix', 'pre_year', 'LTDR_post')
LTDR_post_cols <- c('LTDR_prefix', 'post_year', 'LTDR_post')

d$hansen_prefix = "gefDta$lossyear.na.categorical_"
d$hansen_post = "/gefDta$lossyear.na.categorical_count"
hansen_pre_cols <- c('hansen_prefix', 'pre_year', 'hansen_post')
hansen_post_cols <- c('hansen_prefix', 'post_year', 'hansen_post')
hansen_pre_absolute_cols <- c('hansen_prefix', 'pre_year')

d$pre_NTL_vars = do.call(paste, c(d[ , NTL_pre_cols], list(sep = "")))
d$post_NTL_vars = do.call(paste, c(d[ , NTL_post_cols], list(sep = "")))

d$pre_LTDR_vars = do.call(paste, c(d[ , LTDR_pre_cols], list(sep = "")))
d$post_LTDR_vars = do.call(paste, c(d[ , LTDR_post_cols], list(sep = "")))

d$pre_hansen_vars = do.call(paste, c(d[ , hansen_pre_cols], list(sep = "")))
d$post_hansen_vars = do.call(paste, c(d[ , hansen_post_cols], list(sep = ""))) 
d$pre_hansen_absolute_vars = do.call(paste, c(d[ , hansen_pre_absolute_cols], list(sep = "")))

d$preNTL <- eval(parse(text=d$pre_NTL_vars))
d$postNTL <- eval(parse(text=d$post_NTL_vars))

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



