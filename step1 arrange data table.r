break()
# for(nn in ls(1)) eval(parse(text=paste0("rm(",nn,")")))

##################
##Load required packages
##################
# library(dclone)
# library(R2jags)
# library(RColorBrewer)
# library(rstan)
# library(MCMCpack)
# library(RODBC)

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())


logit <- function(xx) log(xx/(1-xx))
ilogit <- function(xx) exp(xx)/(1+exp(xx))

paste1 <- function(xx){
  output <- xx[1]
  for(ii in 2:length(xx)) output <- paste0(output, xx[ii])
  return(as.character(output))
} 

extractSims <- function(jj, nn) c(jj[[1]][,nn], jj[[2]][,nn], jj[[3]][,nn])
sumNA <- function(xx) sum(xx, na.rm=T)
replaceNA <- function(xx,yy=-1) ifelse(is.na(xx),yy,xx)
replaceBlank <- function(xx,yy=-1) ifelse(length(xx)==0,yy,xx)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


apply0 <- function(tt,ii,ff, cond=F, kk=1){
  if(!cond) tmp_out <- apply(tt,ii,ff)
  if(cond){
    if(!is.null(kk)) tmp_out <- apply(tt,kk,ff)
    # if(is.null(kk)) tmp_out <- do.call(ff, list(tt))
    if(is.null(kk)) tmp_out <- tt
  } 
  return(tmp_out)
} 


all_stats <- function(xx, rr=8) round(c(mean(xx), sd(xx), quantile(xx, c(0.5,0.025, 0.975), na.rm=T), 
                                        Mode(xx), HPDinterval(as.mcmc(xx), prob=0.95)[1:2]),rr)


simplex <- function(xx) xx/sum(xx)


############################


##########
### set directories and folder paths
setwd("C:/Users/Mustelid/OneDrive for Business/Predation Rates/")
dir.folder <- "C:\\Users\\Mustelid\\OneDrive for Business\\"
dir_name <- "C:/Users/Mustelid/OneDrive for Business/Predation Rates/2021/cumulative/"
tag_folder = "G:/tag.data/"


surv.folder <- paste0(dir.folder,"One Off Survival\\")		
pred.folder <- paste0(dir.folder,"Predation Rates\\")

gpud.folder <- "G:\\2017 GPUD synthesis\\"



############################################
##!!! specify years and recapture sites
############################################

year_list = 2008:2021
n_years = length(year_list)


## list of interrogation sites for this analysis
all.int.sites <- c("MCN","JDA","BON","EST","ADULT")
n_segments = length(all.int.sites)

## create lookup table: individual passage routes --> interrogation site
obs_sites <- data.frame(site=as.character(c("B2J","BCC","BO1","BO2","BO3","BO4","JDJ","MC1","MC2","MCJ","TD1","TD2","TWX","PRA","RIA")),
                        abr=as.character(c("BON","BON","ADULT","ADULT","ADULT","ADULT","JDA","ADULT","ADULT","MCN","ADULT","ADULT","EST","ADULT","ADULT")))


################################





################################
### Import and assemble release, resight, and recover data
################################


##########
### released tags
##########
### read in tags for all years
### subset to steelhead, set release week
RIS_tags    <- read.table(paste0(tag_folder,"RIS_ConditionxYear_2008_2018.csv"), sep=',', fill=T, header=T, stringsAsFactors = F)


### add 2019
tmp <- read.table(paste0(tag_folder,"RIS_ConditionxYear_2019.csv"), sep=',', fill=T, header=T, stringsAsFactors = F)
names(tmp) = names(RIS_tags)
RIS_tags = rbind(RIS_tags, tmp)


### simplify
RIS_tags = RIS_tags[,c("Tag_Code","Species","Run","Rearing","Sample_Year","Sample_Date","Length")]
RIS_tags$Rearing = substr(RIS_tags$Rearing, 1, 1)

### add 2020
tmp <- read.table(paste0(tag_folder,"RIS SW 2020v2.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Species_Age=="ST")
tmp$Tag_Code = tmp$TagID
tmp$Rearing = tmp$Rear
tmp$Sample_Year = tmp$RIS_Year
tmp$Sample_Date = tmp$RIS_Date


RIS_tags = rbind(RIS_tags, tmp[,names(RIS_tags)])



### add 2021
tmp <- read.table(paste0(tag_folder,"RIS SW 2021v0.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Species_Age=="ST")
tmp$Tag_Code = tmp$TagID
tmp$Rearing = tmp$Rear
tmp$Sample_Year = tmp$RIS_Year
tmp$Sample_Date = tmp$RIS_Date


RIS_tags = rbind(RIS_tags, tmp[,names(RIS_tags)])
RIS_tags$release_day <- strptime(as.character(RIS_tags$Sample_Date), "%m/%d/%Y")$yday
RIS_tags$release_week <- strptime(as.character(RIS_tags$Sample_Date), "%m/%d/%Y")$yday%/%7+1
# head(RIS_tags)



##########
### interroagtions
##########
### pull in interogations from dams
### translate interrogation codes, set year, reduce to unique records (one int per dam) 

int.file.colnames = c("Tag.Code", "Site.Name", "First.Time.Value")

## 2008 - 2016
observation_history <- read.csv(paste0(gpud.folder,"data\\RIS UCR STHD Interrogation Summary 0.csv"))[,int.file.colnames]
for(ii in 1:5){ # ii <- 0
  tmp <- read.csv(paste0(gpud.folder,"data\\RIS UCR STHD Interrogation Summary ",ii,".csv"))
  observation_history <- rbind(observation_history, tmp[,int.file.colnames])
  rm(tmp)
}

## 2017
tmp <- read.csv(paste0(gpud.folder,"data\\2017 RIS STHD Interrogation Summary.csv"))
observation_history <- rbind(observation_history, tmp[,int.file.colnames])
rm(tmp)

## 2018
tmp <- read.table(paste0(tag_folder,"2018 RIS STHD Interrogation Summary.csv"), sep=',', fill=T, header=T)
observation_history <- rbind(observation_history, tmp[,int.file.colnames])
rm(tmp)


## 2019
tmp <- read.table(paste0(tag_folder,"2019 RIS STHD Interrogation Summary.csv"), sep=',', fill=T, header=T)
observation_history <- rbind(observation_history, tmp[,int.file.colnames])
rm(tmp)


## 2014 and 2015 returns
# tmp <- read.table(paste0(gpud.folder,"data\\2014_2015_returns.csv"), sep=',', fill=T, header=T)
# observation_history <- rbind(observation_history, tmp[,int.file.colnames])
# rm(tmp)

## all returns
tmp <- read.table(paste0(tag_folder,"\\UCR ST ADULT returns 2008_2021.csv"), sep=',', fill=T, header=T)
#### !!! this has juvenile observations as well as returns, some of these were not used in 2018's analysis 
tmp <- subset(tmp, substr(tmp$Site.Name,1,3)%in%subset(obs_sites, abr=="ADULT")$site)
observation_history <- rbind(observation_history, tmp[,int.file.colnames])
rm(tmp)

### clean and ready
observation_history     <- unique(observation_history)
observation_history$abr <- substr(observation_history$Site.Name,1,3)
# head(observation_history)


## lookup main site location, make sure its a juvenile, filter for unique records
observation_history$site <- obs_sites$abr[match(observation_history$abr, obs_sites$site)]
# head(observation_history)




##########
### xtab int histories by tag code to make final histories table
observation_history$day = strptime(as.character(observation_history$First.Time.Value), "%m/%d/%Y")$yday
tmp = aggregate(data.frame(day=observation_history$day), list(Tag.Code=observation_history$Tag.Code, site=observation_history$site), min)

t_hist <- as.data.frame(as.matrix(xtabs(tmp$day ~ tmp$Tag.Code + tmp$site, sparse=T)) )
t_hist <- t_hist[,c("MCN","JDA","BON","EST","ADULT")]
t_hist$ADULT = (t_hist$ADULT>0)*1
# table(t_hist$hist); head(t_hist)



## 2020
tmp <- read.table(paste0(tag_folder,"RIS SW 2020v2.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Species_Age=="ST")
tmp$MCN = replaceNA(strptime(as.character(tmp$MCJ_Date), "%m/%d/%Y")$yday, 0)
tmp$JDA = replaceNA(strptime(as.character(tmp$JDJ_Date), "%m/%d/%Y")$yday, 0)
tmp$BON = replaceNA(strptime(as.character(tmp$BON_Date), "%m/%d/%Y")$yday, 0)
tmp$EST = replaceNA(strptime(as.character(tmp$TWX_Date), "%m/%d/%Y")$yday, 0)
tmp$ADULT = 0
row.names(tmp) = tmp$TagID
tmp = tmp[,all.int.sites]

t_hist = rbind(t_hist, tmp)
rm(tmp)
# head(t_hist)



## 2021
tmp <- read.table(paste0(tag_folder,"RIS SW 2021v0.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Species_Age=="ST")
tmp$MCN = replaceNA(strptime(as.character(tmp$MCJ_Date), "%m/%d/%Y")$yday, 0)
tmp$JDA = replaceNA(strptime(as.character(tmp$JDJ_Date), "%m/%d/%Y")$yday, 0)
tmp$BON = replaceNA(strptime(as.character(tmp$BON_Date), "%m/%d/%Y")$yday, 0)
tmp$EST = replaceNA(strptime(as.character(tmp$TWX_Date), "%m/%d/%Y")$yday, 0)
tmp$ADULT = 0
row.names(tmp) = tmp$TagID
tmp = tmp[,all.int.sites]

t_hist = rbind(t_hist, tmp)
rm(tmp)




##########
### recoveries
##########
### pull in recoveries from colonies
### only use recoveries from same migration year
### translate to consistent colony and island names 

##########
## 2008 - 2017
tags_recovered <- read.csv(paste0(gpud.folder,"\\data\\RIS_Recoveries_2008_2017.csv"))
tags_recovered <- subset(tags_recovered, Recovery.Year==Migration.Year) 
tags_recovered <- tags_recovered[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]


##########
## ESI DCCO 2016
tmp <- read.csv(paste0(gpud.folder,"\\data\\2016ESIDCCO.csv"), stringsAsFactors = F)
tmp <- data.frame(Tag_ID=unique(tmp$Tag.Code), Bird.Colony="DCCO", Recovery.Site="ESI", Migration.Year=2016, Recovery.Year=2016)

tags_recovered <- unique(rbind(tags_recovered, tmp[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]))
rm(tmp)


##########
## 2018
tmp <- read.table(paste0(tag_folder,"2018_Recoveries.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Sample_Year==Recovery.Year)
tmp$Migration.Year = tmp$Sample_Year
tmp$Tag_ID = tmp$Tag_Code

tags_recovered <- unique(rbind(tags_recovered, tmp[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]))
rm(tmp)


## 2019
tmp <- read.table(paste0(tag_folder,"RIS_Recoveries_2019.csv"), sep=',', fill=T, header=T)
tmp <- subset(tmp, Sample_Year==Recovery.Year)
tmp$Migration.Year = tmp$Sample.Year
tmp$Tag_ID = tmp$Tag_Code

tags_recovered <- unique(rbind(tags_recovered, tmp[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]))
rm(tmp)



## 2020
tmp <- read.table(paste0(tag_folder,"RIS SW 2020v2.csv"), sep=',', fill=T, header=T, stringsAsFactors = F)
tmp <- subset(tmp, Species_Age=="ST")
tmp <- subset(tmp, !Recovery_Site=="")

names(tmp) = gsub("_",".",names(tmp))
tmp$Migration.Year = 2020
tmp$Tag_ID = tmp$TagID
tags_recovered <- unique(rbind(tags_recovered, tmp[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]))
rm(tmp)



## 2021
tmp <- read.table(paste0(tag_folder,"RIS SW 2021v0.csv"), sep=',', fill=T, header=T, stringsAsFactors = F)
tmp <- subset(tmp, Species_Age=="ST")
tmp <- subset(tmp, !Recovery_Site=="")
tmp <- subset(tmp, toupper(Recovery.Type)=="NESTING")
tmp$Recovery_Site <- ifelse(tmp$Recovery_Site=="AMB", ifelse(tmp$Recovery.Location=="Plot1_5", "AMB", "AMB6"), tmp$Recovery_Site)

names(tmp) = gsub("_",".",names(tmp))
tmp$Migration.Year = 2021
tmp$Tag_ID = tmp$TagID
tags_recovered <- unique(rbind(tags_recovered, tmp[,c("Tag_ID","Bird.Colony","Recovery.Site","Migration.Year")]))
rm(tmp)



tags_recovered$Recovery.Site = as.character(tags_recovered$Recovery.Site)
tags_recovered$Bird.Colony   = as.character(tags_recovered$Bird.Colony)
tags_recovered$Bird.Colony <- ifelse(tags_recovered$Bird.Colony%in%c("CAGU","RBGU"),"LAXX",tags_recovered$Bird.Colony)
tags_recovered$Bird.Colony <- ifelse(grepl("MIX",toupper(tags_recovered$Bird.Colony)), "MIX", tags_recovered$Bird.Colony)


### adjustments
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Recovery.Site%in%c("CBI") & tags_recovered$Bird.Colony=="CATE", "BLA", tags_recovered$Recovery.Site)

# 2016
tags_recovered$Bird.Colony <- ifelse(tags_recovered$Migration.Year==2016 & tags_recovered$Recovery.Site%in%c("ESI") & tags_recovered$Bird.Colony=="Unkown", "DCCO", tags_recovered$Bird.Colony)

# 2019
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Migration.Year==2019 & tags_recovered$Recovery.Site%in%c("LEN-NORTH"), "LEN", tags_recovered$Recovery.Site)

# 2020
tags_recovered$Bird.Colony <- ifelse(tags_recovered$Migration.Year==2020 & tags_recovered$Recovery.Site%in%c("LEN-NORTH"), "CATE", tags_recovered$Bird.Colony)
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Migration.Year==2020 & tags_recovered$Recovery.Site%in%c("LEN-NORTH"), "LEN", tags_recovered$Recovery.Site)
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Migration.Year==2020 & tags_recovered$Recovery.Site%in%c("ESI-SAT"), "SAT", tags_recovered$Recovery.Site)

# 2021
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Migration.Year==2021 & tags_recovered$Recovery.Site%in%c("LEN-SHOAL"), "LEN", tags_recovered$Recovery.Site)
tags_recovered$Recovery.Site <- ifelse(tags_recovered$Migration.Year==2021 & tags_recovered$Recovery.Site%in%c("ESI-SAT"), "SAT", tags_recovered$Recovery.Site)


tags_recovered$colony_island <-  paste0(tags_recovered$Bird.Colony, tags_recovered$Recovery.Site)
tags_recovered$colony_island = ifelse(tags_recovered$colony_island=="", "Not Recovered", tags_recovered$colony_island)

tags_recovered$colony_island = paste0(tags_recovered$Bird.Colony, tags_recovered$Recovery.Site)
tags_recovered$colony_island = ifelse(tags_recovered$Sample_Year==tags_recovered$Recovery.Year, tags_recovered$colony_island, "")
tags_recovered$colony_island = replaceNA(tags_recovered$colony_island, "")


# CBI adjustements
tCBI  <- read.table(paste0(tag_folder,"CBI islands by year.csv"), sep=',', fill=T, header=T, stringsAsFactors = F)

mm = match(paste0(tags_recovered$Migration.Year, tags_recovered$colony_island), paste0(tCBI$Year, tCBI$PredRatesIslandName))
tags_recovered$colony_island = ifelse(is.na(mm), tags_recovered$colony_island, tCBI$DetPosteriorName[mm])






##########
### assemble all data and distill down to the essentials; calculate segment of disappearance 
tData <- data.frame(pit_code=RIS_tags$Tag_Code, rear_type=RIS_tags$Rearing, length=RIS_tags$Length,
                    release_year=RIS_tags$Sample_Year, release_week=RIS_tags$release_week, release_day=RIS_tags$release_day,
                    colony_island=replaceNA(tags_recovered$colony_island[match(RIS_tags$Tag_Code, tags_recovered$Tag_ID)],"Not Recovered"),
                    RIS=1, t_hist[match(RIS_tags$Tag_Code,row.names(t_hist)), all.int.sites]
)
row.names(tData) <- NULL
for(cc in all.int.sites) tData[,cc] <- replaceNA(tData[,cc],0)
for(cc in all.int.sites) tData[,cc] = (tData[,cc]>0)*1
tData$int.disappear <- apply(tData[,c("release_week",all.int.sites)], 1, function(xx) max(which(rev(cumsum(rev(xx)))>0))-1)
# head(subset(tData, release_year==2016))
# table(tData$release_year, tData$colony_island)






################################
### construct colony_range
################################
##########
### Define forage ranges for each colony to be considered; also define colony_list 
### (this used to be inferred directly from the observed data but was prone to error and is now hard-coded)
colony_range <- data.frame(colony.island=c("CATEBKL", "CATEPTI", "CATENEP", "CATELEN",  "LAXXIS20",  
                                           "CATEBGI", "LAXXBGI", "DCCOFDI", 
                                           "CATECSI", "LAXXCSI", 
                                           "CATEBLA", "CATERKI", "CATELONG", "CATEMID",
                                           "MIXBGI", "AWPEBGI",
                                           "LAXXMRI", "LAXXANV", "LAXXSIX",
                                           "CATEESI",  "CATERCI", "CATESAT", "DCCOESI", "DCCOAMB", "DCCOAMB6"), 
                           min=c(rep(1,8), 
                                 rep(1,6), 
                                 rep(1,2),
                                 rep(2,3),
                                 rep(5,6)), 
                           max=c(rep(1,8),
                                 rep(2,6),
                                 rep(3,2),
                                 rep(3,3),
                                 rep(5,6)))
colony_range$colony.island <- as.character(colony_range$colony.island)
colony_list <- unique(colony_range$colony.island)
colony_list = c(colony_list[!grepl("MIX",colony_list)], colony_list[grepl("MIX",colony_list)])
colony_range = colony_range[match(colony_list,colony_range$colony.island), ]


n_colonies = length(colony_list)


### NA out colonies in tData not under explicit consideration in this analysis
tData$colony_island <- replaceNA(colony_list[match(tData$colony_island, colony_list)], "Not Recovered")
# table(tData$release_year, tData$colony_island)


################################
### construct det_table
################################

### read in detection files, assign week, distill essential data
det_table <- NULL
for(yy in year_list){   # yy <- 2012
  for(ii in colony_list){  # yy <- 2011; ii="CATECBI"
    tmp_colony <- substr(ii, 1,4)
    tmp_island <- substr(ii, 5,99)
    
    ######
    if(tmp_island=="AMB6") tmp_island = "AMB"
    ######
    
    save_list <- NULL
    try(load(paste0(pred.folder,"SE.data\\Posteriors\\",tmp_colony,tmp_island,yy,".Rdata")), silent=T)
    if(!is.null(save_list)){
      tmp <- data.frame(year=yy, colony_island=ii)
      det_table <- rbind(det_table, tmp)
    } 
  }
}
# det_table 
warnings()



#### !!!!!! ####
### replace colonies without any detection
tmp <- unique(paste0(det_table$year,det_table$colony))

# subset(tData, !paste0(release_year,colony_island)%in%tmp & !colony_island%in%"Not Recovered")
tData <- rbind(subset(tData, paste0(release_year,colony_island)%in%tmp),
               transform(subset(tData, !paste0(release_year,colony_island)%in%tmp), colony_island="Not Recovered"))


################################
### col pred
################################

#### identify active colonies by year
col_pred <- array(0, dim=c(n_years,n_segments,n_colonies+2))
for(yy in 1:n_years){
  # yy=1
  tmp_col_pred <- array(0, dim=c(n_segments,n_colonies))
  
  tmp_col_list <- colony_list[which(colony_list%in%subset(det_table, year==as.character(year_list)[yy])$colony_island)]
  # if(yy==3) tmp_col_list <- c(tmp_col_list, "CATEBLA")   #### no detection data available from 2010 but plenty otherwise
  
  for(cc in tmp_col_list){  # cc <- "LAXXMRI"
    mm <- match(cc, colony_list)
    tmp_col_pred[subset(colony_range, colony.island==cc)$min:min(subset(colony_range, colony.island==cc)$max,n_segments),mm] <- 1
  }
  tmp_col_pred <- cbind(1,tmp_col_pred,1)
  
  col_pred[yy,,] <- tmp_col_pred
  # print(col_pred[yy,,])
}
colonies_by_year <- (apply(col_pred[,,1:n_colonies+1], c(1,3), sum)>0)*1


n_only_colonies = apply(colonies_by_year %*% diag(!grepl("MIX",colony_list)*1), 1, sum)
n_mix_colonies  = apply(colonies_by_year %*% diag(grepl("MIX",colony_list)*1), 1, sum)


data_list = list(year_list = year_list,
                 n_years = n_years, 
                 
                 all.int.sites = all.int.sites,
                 n_segments = n_segments,
                 
                 tData=tData,
                 
                 det_table=det_table,
                 
                 colony_list=colony_list,
                 n_colonies=n_colonies,
                 n_only_colonies=n_only_colonies,
                 n_mix_colonies=n_mix_colonies,
                 colony_range=colony_range, 
                 
                 col_pred=col_pred, 
                 colonies_by_year=colonies_by_year)

save(data_list, file=paste0(dir_name,"data/tData.Rdata"))


















