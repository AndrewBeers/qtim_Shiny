# rm(list=ls())
#setwd("~/Documents/QIN/featureComp/submissions")
#setwd("~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final")
#source("~/SyncplicityFolders/QIN/featureComp/submissions/scripts/calculateCCC.R")

setwd("C:/Users/azb22/Documents/Scripting/Interval_Lung_QIN_Challenge/Public_ShinyApp/Interval_Lung_Challenge_ShinyApp")
# source("~/Syncplicity/QIN/featureComp/submissions/scripts/calculateCCC.R")

library(RSQLite)
library(sqldf)
library(ggplot2)
library(XLConnect)
library(reshape2)
library(plyr)
library(tidyr)
library(doBy)
library(corrplot)
# library(agRee)
library(data.table)


### features in individual tables 
temp = list.files(pattern="*.csv")
sqlitePath<-"QINLungSegmentationStatsDb_CUMCUpdate_3.sqlite"
db <- dbConnect(SQLite(), sqlitePath)

for(i in temp){
  print(i)
  a<- unlist(strsplit(i, ".", fixed=TRUE))
  fName<-a[1]
 # sqlitePath <-paste(fName,".sqlite",sep="")
  feat <-read.csv(i)
  dbWriteTable(conn = db, name = fName, value = feat, row.names = FALSE, overwrite = TRUE)
 # dbReadTable(db, "features")

  
}

# sql <- readLines("test.sql")
sql = 'CREATE VIEW "volumes" AS select PID, Indx, Outcome, UMC_Voxel_Integral as UMC, UCL_Voxel_Integral as UCL, UCLman_Voxel_Integral as UCL_manual, CMC_Voxel_Integral as CMC, CMCman_Voxel_Integral as CMC_manual, DFCI_Voxel_Integral as DFCI, MCC_Voxel_Integral as MCC, CMCman_diameter as CMCman_diameter from QINLung_VolumeData'
dbSendQuery(conn=db, sql)

sql = 'CREATE VIEW "volume_changes" AS select PID, Outcome, UMC_Voxel_Integral as UMC, UCL_Voxel_Integral as UCL, UCLman_Voxel_Integral as UCL_manual, CMC_Voxel_Integral as CMC, CMCman_Voxel_Integral as CMC_manual, DFCI_Voxel_Integral as DFCI, MCC_Voxel_Integral as MCC, CMCman_diameter as CMCman_diameter from QINLung_VolumeChangeData'
dbSendQuery(conn=db, sql)

sql = 'CREATE VIEW "volume_changes_percent" AS select PID, Outcome, UMC_Voxel_Integral as UMC, UCL_Voxel_Integral as UCL, UCLman_Voxel_Integral as UCL_manual, CMC_Voxel_Integral as CMC, CMCman_Voxel_Integral as CMC_manual, DFCI_Voxel_Integral as DFCI, MCC_Voxel_Integral as MCC, CMCman_diameter as CMCman_diameter from QINLung_VolumeChangeData_Percent'
dbSendQuery(conn=db, sql)

dbDisconnect(db)


# ### dictionary
# 
# 
# 
# ## dictionary
# #sqlitePath <-"~/Documents/QIN/featureComp/submissions/featComp.sqlite"
# setwd("~/SyncplicityFolders/QIN/featureComp/submissions/all/dictionary")
# dictFiles = list.files(pattern="*.csv")
# 
# sqlitePath <-"~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# 
# db <- dbConnect(SQLite(), sqlitePath)
# dbListTables(db)
# 
# for(i in dictFiles){
#   print(i)
#   a<- unlist(strsplit(i, ".", fixed=TRUE))
#   fName<-a[1]
#   # sqlitePath <-paste(fName,".sqlite",sep="")
#   dict <-read.csv(i)
#   dbWriteTable(conn = db, name = "featureDic", value = dict, row.names = FALSE, append=TRUE)
#   # dbReadTable(db, "features")
# }
# dbDisconnect(db)
# 
# 
# #### correletion coefficients of features########
# sqlitePath <-"~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# 
# db <- dbConnect(SQLite(), sqlitePath)
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD ccc_intra real")
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD ccc_inter real")
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD ccc_all real")
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD occ real")
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD ccc2 real")
# # results <- dbSendQuery(db,"ALTER TABLE featureDic ADD cccwscv real")
# 
# 
# 
# ##################calculate ccc by site ##############################
# sqlitePath <-"~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# db <- dbConnect(SQLite(), sqlitePath)
# tbls<- dbListTables(db)
# 
# for (i in tbls) {
#   sqlcmd1 <- paste("SELECT  count(uNoduleId) as count from ", i, sep="")
#   #print(sqlcmd1)
#   #data.annual <- dbGetQuery(db, sqlcmd)
# 
#   results <- tryCatch(
#   {
#     results <- dbSendQuery(db,sqlcmd1)
#     res <-fetch(results)
#     #print(res)
#     print(i)
#     sqlcmd2 <- paste("SELECT * from ", i, sep="")
#     # print(sqlcmd2)
#   #  sqlcmd2<- "SELECT * from pm_features"
#    # sqlcmd2<- "SELECT * from nStanford_features"
#     results2 <- dbSendQuery(db,sqlcmd2)
#     features <-fetch(results2)
#     #print(res2)
#     cn<-colnames(features)
#     cnLength<-length(cn)
#     numFeat <-cnLength-11
#     featNames<-cn[c(12:length(cn))]
#     d <- matrix(nrow=numFeat, ncol=6)
#     #for stanford
#     features[features == 0] <- NA
#     stds<-apply(features[,c(12:cnLength)], 2, sd)
#     for (j in 1:numFeat) {
#       thisFeat <- featNames[j]
#      # featRow <-fetch(results3)
#       #  print(featRow)
#       res1<-agreementCCC(features,featNames,j)
#       thisfeatCCC <- unlist(res1)
#       #ccc_intra,ccc_inter, ccc_all,occ, ccc2, cccwscv
#       ccc_intra<-thisfeatCCC[1]
#       if (!is.na(ccc_intra)) {
#         sqlcmd3 <- paste("Update featureDic set ccc_intra = ",ccc_intra," where featureName = \"",thisFeat,"\"", sep="")
#         results3 <- dbSendQuery(db,sqlcmd3)
#       }
#       ccc_inter<-thisfeatCCC[2]
#       if (!is.na(ccc_inter)) {
#         sqlcmd4 <- paste("Update featureDic set ccc_inter = ",ccc_inter," where featureName = \"",thisFeat,"\"", sep="")
#         results4 <- dbSendQuery(db,sqlcmd4)
#       }
# 
#       ccc_all<-thisfeatCCC[3]
#       if (!is.na(ccc_all)) {
#         sqlcmd5 <- paste("Update featureDic set ccc_all = ",ccc_all," where featureName = \"",thisFeat,"\"", sep="")
#         results5 <- dbSendQuery(db,sqlcmd5)
#       }
# 
#       occ<-thisfeatCCC[4]
#       if (!is.na(occ)) {
#         sqlcmd6 <- paste("Update featureDic set occ = ",occ," where featureName = \"",thisFeat,"\"", sep="")
#         results6 <- dbSendQuery(db,sqlcmd6)
#       }
# 
#       ccc2<-thisfeatCCC[5]
#       if (!is.na(ccc2)) {
#         sqlcmd7 <- paste("Update featureDic set ccc2 = ",ccc2," where featureName = \"",thisFeat,"\"", sep="")
#         results7 <- dbSendQuery(db,sqlcmd7)
#       }
# 
#       cccwscv<-thisfeatCCC[6]
#       if (!is.na(cccwscv)) {
#         sqlcmd8 <- paste("Update featureDic set cccwscv = ",cccwscv," where featureName = \"",thisFeat,"\"", sep="")
#         results8 <- dbSendQuery(db,sqlcmd8)
#       }
# 
# 
#     }
# 
#   },
#   error=function(cond) {
#     # message("Here's the original error message:")
#     message(cond)
#     cat(paste("no uid\n"))
#     # Choose a return value in case of error
#     #  return(NA)
#   }
#   )
# }
# 
# 
# ################ pairwise correlation coefficients#############
# #"cStanford_features"  "cumc_features"       "featureDic"          "iowa_features"       "nStanford_features"  "pairwiseCorrelation" "pm_features"
# # "stanford_features"   "ucla_features"       "umich_features"      "usf_features"
# 
# sqlitePath <-"~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# sqlitePath <-"~/Syncplicity/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# 
# db <- dbConnect(SQLite(), sqlitePath)
# 
# site1<-"pm_features"
# site1<-"cumc_features"
# sqlcmd1<- paste("SELECT * from ", site1, sep="")
# print(sqlcmd1)
# results1 <- dbSendQuery(db,sqlcmd1)
# features1 <-fetch(results1)
# cn1<-colnames(features1)
# numFeat1 <-length(cn1)-11
# featNames1<-cn1[c(12:length(cn1))]
# 
# featVals1<-features1[,c(12:length(cn1))]
# nfv1<-featVals1[,complete.cases(t(featVals1))]
# 
# bad <- sapply(nfv1, function(x) any(is.nan(x)))
# nfv1<-nfv1[,!bad]
# 
# featVals1.scale <- scale(nfv1[1:ncol(nfv1)],center=TRUE,scale=TRUE)
# corMatFeatScale1 <- cor(featVals1.scale)
# 
# corrplot(corMatFeatScale1, order = "hclust", method = "color")
# corrplot(corMatFeatScale1)
# 
# 
# 
# site2<-"ucla_features"
# site2<-"usf_features"
# sqlcmd2 <- paste("SELECT * from ", site2, sep="")
# # print(sqlcmd2)
# results2 <- dbSendQuery(db,sqlcmd2)
# features2 <-fetch(results2)
# cn2<-colnames(features2)
# numFeat2 <-length(cn2)-11
# featNames2<-cn1[c(12:length(cn2))]
# 
# featVals2<-features2[,c(12:length(cn2))]
# nfv2<-featVals2[,complete.cases(t(featVals2))]
# 
# bad <- sapply(nfv2, function(x) any(is.nan(x)))
# nfv2<-nfv2[,!bad]
# 
# featVals2.scale <- scale(nfv2[1:ncol(nfv2)],center=TRUE,scale=TRUE)
# corMatFeatScale2 <- cor(featVals2.scale)
# 
# corrplot(corMatFeatScale2, order = "hclust", method = "color")
# corrplot(corMatFeatScale2)
# 
# corMatFeatScale12 <- cor(featVals1.scale,featVals2.scale)
# 
# corrplot(corMatFeatScale12, method = "color",order = "hclust",)
# corrplot(corMatFeatScale12)
# 
# 
# feature1
# sqlcmd1<- paste("SELECT uNoduleID, cumc_feature_37  from ", site1, sep="")
# print(sqlcmd1)
# sqlcmd1<-"SELECT u.uNoduleID, c.cumc_feature_37 as cFeature, u.ucla_feature_2 as uFeature from cumc_features as c, ucla_features as u where u.uNoduleID=c.uNoduleID"
# results1 <- dbSendQuery(db,sqlcmd1)
# features1 <-fetch(results1)
# 
# g<-ggplot(features1, aes(x=))
# 
# p1 <- ggplot(features1, aes(x = cFeature, y = uFeature))
# 
# p1 + geom_point()