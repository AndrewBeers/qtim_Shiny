# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)

library(RSQLite)
library(ggplot2)
library(ROCR)
library(Agreement)
library(reshape2)
library(RColorBrewer)
#sqlitePath <<-"~/Syncplicity/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
#sqlitePath <<-"~/SyncplicityFolders/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
#sqlitePath <<-"~/Documents/QIN/featureComp/submissions/all/features/final/featuresDb.sqlite"
# sqlitePath <<-"/Users/kalpathy/Documents/QIN/intervalChallenge/shiny/dbBlandAltmanVolume/QINLungSegmentationDb.sqlite"
sqlitePath <<-"QINLungSegmentationStatsDb_CUMCUpdate_3.sqlite"

db <- dbConnect(SQLite(), sqlitePath)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
na.var <- function(x){
  var(x, na.rm=TRUE)
}

script <- "$('tbody tr td:nth-child(5)').each(function() {
var cellValue = $(this).text();
if (cellValue == 'Small') {
$(this).css('background-color', '#0c0');
}
else if (cellValue== 'Large') {
$(this).css('background-color', '#f00');
}
})"

shinyServer(function(input, output, session) {
  
  output$plot1 <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)

    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)    
    
    s1<-thisSite1
    s2<-thisSite2
    f<-paste(s2, " ~", s1)
    fit <- lm(f, data = vols)
    if (s1==s2){
      slope = 1
      pvalue = 0
    }
    else {
      slope = signif(fit$coef[[2]], 5)
      pvalue = signif(summary(fit)$coef[2,4], 5)
    }
    
    stat1 = as.numeric(unlist(vols[s1]))
    stat2 = as.numeric(unlist(vols[s2]))
    
    mean1 = sapply(vols[s1], mean, na.rm=TRUE)
    mean2 = sapply(vols[s2], mean, na.rm=TRUE)
    sd1 = sd(stat1, na.rm=TRUE)
    sd2 = sd(stat2, na.rm=TRUE)
    cor12 = cor(vols[s1], vols[s2],use='complete.obs')
    
    ccc = (2 * cor12 * sd1 * sd2) / (sd1^2 + sd2^2 + (mean1 - mean2)^2)
    
    g1<- ggplot(vols, aes_string(s1,s2, colour = "colorOutcome")) + geom_point() + labs(title = "\n Segmentation Volume, Site vs. Site") + guides(colour=FALSE) + geom_text(x = .3*colMax(vols[s1]), y= .7*colMax(vols[s2]), colour="black", aes(label=paste("CCC = ",signif(ccc, 5),'\n',"Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),'\n',"Intercept =",signif(fit$coef[[1]],5 ),'\n'," Slope =",slope,'\n'," P =",pvalue)))
    g1+ stat_smooth(method = "lm", col = "red") +theme_bw()
    
  })
  
  output$plot2 <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    s1<-thisSite1
    s2<-thisSite2
    vols$Difference<-vols[[s2]]-vols[[s1]]
    vols$Average_Volume<-(vols[[s2]]+vols[[s1]])/2
    g1<- ggplot(vols, aes(Average_Volume,Difference, colour = colorOutcome)) + geom_point()
    g1 + geom_hline(yintercept = 0)+theme_bw()+scale_y_continuous() + labs(title = "\n Segmentation Volume Differences vs. Average Segmentation Volumes") + guides(colour=guide_legend(title="Patient Outcome")) + theme(legend.position="bottom")
    
  })
  
  output$vol_histogram_sites <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    
    limits = input$SliderHistogramLimits
    bins = input$SliderHistogramBins
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols = vols[,c(thisSite1, thisSite2)]
    # vols = data.frame(cbind(vols[,thisSite1], vols[,thisSite2]))
    vols = melt(vols, measure.vars = c(thisSite1, thisSite2))
    
    if (input$histogram_mode == "density"){
    g1 = ggplot(vols, aes(x=value))+ geom_density(aes(colour=variable, fill=variable, alpha=.2)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "cumulative"){
    g1 = ggplot(vols, aes(x=value))+ stat_ecdf(aes(colour=variable)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "histogram"){
     g1 = ggplot(vols, aes(x=value))+ geom_histogram(data=subset(vols,variable == thisSite1),fill = "red", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins)) + geom_histogram(data=subset(vols,variable == thisSite2),fill = "blue", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins))
    }

    g1 + labs(title = "\n Segmentation Volume Histogram, Site 1 (Red) vs Site 2 (Blue)")
    
  })
  
  output$vol_histogram_status <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    
    limits = input$SliderHistogramLimits
    bins = input$SliderHistogramBins
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols = vols[,c(thisSite1, "Outcome")]

    if (input$histogram_mode == "density"){
      g1 = ggplot(vols, aes_string(x=thisSite1))+ geom_density(aes(colour=Outcome, fill=Outcome, alpha=.2)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "cumulative"){
      g1 = ggplot(vols, aes_string(x=thisSite1))+ stat_ecdf(aes(colour=Outcome)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "histogram"){
      g1 <- ggplot(vols, aes_string(x=thisSite1)) + geom_histogram(data=subset(vols,Outcome == 'NC'),fill = "blue", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins)) + geom_histogram(data=subset(vols,Outcome == 'CC'),fill = "red", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins))
    }
    
    g1 + labs(title = "\n Segmentation Volume Histogram, Site 1, NC (Blue) vs. CC (Red)")
    
  })
  
  output$vol_histogram_date <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    
    limits = input$SliderHistogramLimits
    bins = input$SliderHistogramBins
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols$Indx = factor(vols$Indx)
    vols = vols[,c(thisSite1, "Indx")]

    print(vols)
    
    if (input$histogram_mode == "density"){
      g1 = ggplot(vols, aes_string(x=thisSite1))+ geom_density(aes(colour=Indx, fill=Indx, alpha=.2)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "cumulative"){
      g1 = ggplot(vols, aes_string(x=thisSite1))+ stat_ecdf(aes(colour=Indx)) + scale_x_continuous(limits=c(limits[1], limits[2]))     
    }
    if (input$histogram_mode == "histogram"){
      g1 = ggplot(vols, aes_string(x=thisSite1))+ geom_histogram(data=subset(vols,Indx == 1),fill = "red", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins)) + geom_histogram(data=subset(vols,Indx == 2),fill = "blue", alpha=0.2, colour='black', breaks=seq(limits[1], limits[2], by = limits[2] / bins))
    }
    
    g1 + labs(title = "\n Segmentation Volume Histogram, Site 1, 1999 (Red) vs 2000 (Blue)")
    

  })
  
  output$changeplot1 <- renderPlot({

    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    print(baselinevols)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
      slider = input$SliderChangePercent
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
      slider = input$SliderChangeAbsolute
    }
    v <- dbSendQuery(db, sqlcmd)
    
    vols <-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols <- replace(vols, vols > slider[2], NA)
    vols <- replace(vols, vols < slider[1], NA)
    s1<-thisSite1
    s2<-thisSite2
    f<-paste(s2, " ~", s1)
    fit <- lm(f, data = vols)
    if (s1==s2){
      slope = 1
      pvalue = 0
    }
    else {
      slope = signif(fit$coef[[2]], 5)
      pvalue = signif(summary(fit)$coef[2,4], 5)
    }
    
    stat1 = as.numeric(unlist(vols[s1]))
    stat2 = as.numeric(unlist(vols[s2]))
    
    mean1 = sapply(vols[s1], mean, na.rm=TRUE)
    mean2 = sapply(vols[s2], mean, na.rm=TRUE)
    sd1 = sd(stat1, na.rm=TRUE)
    sd2 = sd(stat2, na.rm=TRUE)
    cor12 = cor(vols[s1], vols[s2],use='complete.obs')
    
    ccc = (2 * cor12 * sd1 * sd2) / (sd1^2 + sd2^2 + (mean1 - mean2)^2)
    
    g1<- ggplot(vols, aes_string(s1,s2, colour = "colorOutcome")) + geom_point() + labs(title = "\n Segmentation Volume Change, Site vs. Site") + guides(colour=FALSE) + geom_text(x = .3*colMax(vols[s1]), y= .7*colMax(vols[s2]), colour="black", aes(label=paste("CCC = ",signif(ccc, 5),'\n',"Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),'\n',"Intercept =",signif(fit$coef[[1]],5 ),'\n'," Slope =",slope,'\n'," P =",pvalue)))
    g1+ stat_smooth(method = "lm", col = "red") +theme_bw() + scale_x_continuous(limits=c(slider[1]*1.05, slider[2]*1.05))
    
  })
  
  output$changeplot2 <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
      slider = input$SliderChangePercent
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
      slider = input$SliderChangeAbsolute
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    vols$colorOutcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols <- replace(vols, vols > slider[2], NA)
    vols <- replace(vols, vols < slider[1], NA)
    s1<-thisSite1
    s2<-thisSite2
    vols$Difference<-vols[[s2]]-vols[[s1]]
    vols$Average_Volume_Change<-(vols[[s2]]+vols[[s1]])/2
    g1<- ggplot(vols, aes(Average_Volume_Change,Difference,colour=colorOutcome)) + geom_point()
    g1 + geom_hline(yintercept = 0)+theme_bw()+scale_x_continuous(limits=c(slider[1]*1.05, slider[2]*1.05)) + scale_y_continuous() + labs(title = "\n Volume Change Difference vs. Average Volume Change") + guides(colour=guide_legend(title="Patient Outcome")) + theme(legend.position="bottom")
    
  })
  
  output$ROCplot <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
      slider = input$SliderChangePercent2
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
      slider = input$SliderChangeAbsolute2
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    print(vols)
    s1<-thisSite1
    vols = vols[c("Outcome",s1)]
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    model = glm(Outcome~.,family=binomial(link='logit'),data=vols)
    p = predict(model, vols, type="response")
    pr = prediction(p, vols$Outcome)
    prf <- performance(pr, measure="tpr", x.measure="fpr")
    
    Nulldf = as.data.frame(c(0,1))
    colnames(Nulldf)[1] = "fpr"
    Nulldf$tpr = c(0,1)
    
    ROCdf = data.frame(prf@x.values, prf@y.values)
    colnames(ROCdf) = c('fpr','tpr')
    lo = loess(tpr~fpr, ROCdf)
    # polyfit = lm(tpr~poly(fpr,2,raw=TRUE), ROCdf)
    SmoothROCdf = data.frame(seq(0,1,.01))
    colnames(SmoothROCdf) = 'fpr'
    SmoothROCdf$tpr = predict(lo, newdata=SmoothROCdf)
    SmoothROCdf = replace(SmoothROCdf, SmoothROCdf < 0, NA)
    
    g1 = ggplot(ROCdf, aes(fpr, tpr)) + geom_point()
    g1 + theme_bw()+scale_y_continuous()+scale_x_continuous()+labs(title = "\n ROC Curve, Site 1", x="False Positive Rate", y='True Positive Rate') + geom_text(x = .9, y= .1, colour="black", aes(label=paste("AUC = ",signif(as.numeric(performance(pr, measure = "auc")@y.values),2),", n = ",sum(!is.na(vols[s1]))))) +
      geom_line(data = Nulldf, linetype=2, aes_string("fpr","tpr"))  +
      geom_line(data = SmoothROCdf, linetype=1, aes_string("fpr","tpr"))})
  
  output$LogisticPlot <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
      slider = input$SliderChangePercent2
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
      slider = input$SliderChangeAbsolute2
    }
    
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    vols = vols[c("Outcome",s1)]
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    print(vols)
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols$averages = NULL
    vols$CMCman_diameter = NULL
    
    vols <- replace(vols, vols > slider[2], NA)
    vols <- replace(vols, vols < slider[1], NA)
    
    model = glm(Outcome~.,family=binomial(link='logit'),data=vols)
    vols$colorOutcome = factor(vols$Outcome)
    logdf = as.data.frame(seq(from=min(vols[s1], na.rm=TRUE),to=max(vols[s1], na.rm=TRUE),by=.05))
    colnames(logdf)[1] = s1
    logdf$predictions = predict(model, logdf,type="resp") + 1
    
    g1 = ggplot(vols, aes_string(s1,'Outcome')) + geom_point(aes_string(colour = "colorOutcome"))
    g1 + theme_bw() + scale_x_continuous(limits=c(slider[1]*1.05, slider[2]*1.05)) + labs(title = "\n Logistic Regression, Site 1", x=paste(s1, 'Change')) + guides(colour=guide_legend(title="Patient Outcome")) + theme(legend.position="bottom") + 
      geom_line(data = logdf, linetype=2, aes_string(s1,'predictions'))})
  
  output$ROCCompareplot <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    selectedSites <- input$ROCsites
    # if (is.null(selectedSites))
    # return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    vols = vols[c("Outcome",selectedSites)]
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter_thresh = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter_thresh > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter_thresh < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter_thresh <= 8 | CMCman_diameter_thresh >= 8)
    
    vols$averages = NULL
    vols$CMCman_diameter_thresh = NULL
    
    print(vols)
    
    Nulldf = as.data.frame(c(0,.5,1))
    colnames(Nulldf)[1] = "fpr"
    Nulldf$tpr = c(0,.5,1)
    g1 = ggplot(Nulldf, aes(fpr, tpr)) + geom_line(linetype=2)
    
    colorsROC = c("UCL"="red","UCL_manual"="green","CMC"="blue","CMC_manual"="brown","DFCI"="orange","MCC"="purple","UMC"="darkgoldenrod","CMCman_diameter"="black")
    
    AUCText=""
    
    for (i in 2:dim(vols)[2]){
      model = glm(Outcome~vols[,i],family=binomial(link='logit'),data=vols)
      p = predict(model, vols, type="response")
      pr = prediction(p, vols$Outcome)
      prf <- performance(pr, measure="tpr", x.measure="fpr")
      ROCdf = data.frame(prf@x.values, prf@y.values)
      colnames(ROCdf) = c('fpr','tpr')
      AUC = signif(as.numeric(performance(pr, measure = "auc")@y.values),2)
      AUCText = paste(AUCText,colnames(vols)[i], " AUC = ",AUC,", n = ",sum(!is.na(vols[,i])),"\n")
      g1 = g1 + geom_line(data = ROCdf, linetype=1, aes(fpr,tpr), color=colorsROC[colnames(vols)[i]])
    }    
    
    g1 = g1 + geom_text(x = .8, y= .3, colour="black", aes(label=AUCText))
    
    g1 = g1 + scale_colour_manual(name = 'Site Key', values =c('red'='red','green'='green',"blue"="blue","orange"="orange","purple"="purple","darkgoldenrod"="darkgoldenrod","black"="black"), labels = c("UCLA","UCLA_Manual","CMC","DFCI","MCC","UMC","STAPLE"), guide_legend(override.aes=aes(fill=NA)))
    g1 + theme_bw()+scale_y_continuous()+scale_x_continuous()+labs(title = "\n ROC Curve, Site 1", x="False Positive Rate", y='True Positive Rate')
    
    
    # ROCdf = data.frame(prf@x.values, prf@y.values)
    # colnames(ROCdf) = c('fpr','tpr')
    # lo = loess(fpr~tpr, ROCdf)
    # SmoothROCdf = data.frame(predict(lo), ROCdf$tpr)
    # colnames(SmoothROCdf) = c('fpr','tpr')
    # 
    # g1 = ggplot(ROCdf, aes(fpr, tpr)) + geom_point()
    # g1 + theme_bw()+scale_y_continuous()+scale_x_continuous()+labs(title = "\n ROC Curve, Site 1", x="False Positive Rate", y='True Positive Rate') + geom_text(x = .9, y= .1, colour="black", aes(label=paste("AUC = ",signif(as.numeric(performance(pr, measure = "auc")@y.values),2)))) +
    #   geom_line(data = Nulldf, linetype=2, aes_string("fpr","tpr"))  +
    #   geom_line(data = SmoothROCdf, linetype=1, aes_string("fpr","tpr"))
  })
  
  output$AllDICE <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    sqlcmd <- paste("SELECT * from QINLung_DiceData",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    
    vols = vols[1:200,]
    
    vols$averages = baselinevols$averages
    vols$Outcome = baselinevols$Outcome
    vols$CMCman_diameter = baselinevols$CMCman_diameter
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols = vols[,c(-1,-17,-18, -19)]
    vols = vols[,!grepl("ucl", colnames(vols)) | grepl("man", colnames(vols))]
    vols$Average_DICE = rowMeans(vols, na.rm = TRUE)
    qplot(vols$Average_DICE, geom="histogram",
          binwidth = 1/input$HistogramBins,  
          main = "Average DICE Coeffecient, All Sites", 
          xlab = "DICE Coeffecient",  
          fill=I("blue"), 
          col=I("red"), 
          alpha=I(.2))
    # histogram = ggplot(vols, aes(x=Average_DICE) + geom_freqpoly(bins=input$binwidth, size=2, aes(y=..count..)) + ColorScale + ylab("count"))
    
  })
  
  output$TwoDICE <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
        
    sqlcmd <- paste("SELECT * from QINLung_DiceData",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    
    vols = vols[1:200,]
    
    vols$averages = baselinevols$averages
    vols$Outcome = baselinevols$Outcome
    vols$CMCman_diameter = baselinevols$CMCman_diameter
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    s1<-thisSite1
    s2 = thisSite2
    vols = vols[,-1]
    if (s1 == 'UCL_manual')
      s1 = 'UCLman'
    if (s2 == 'UCL_manual')
      s2 = 'UCLman'
    if (length(grep(paste(s1,'_',s2,sep=''), colnames(vols))) == 0)
      vols = vols[,grep(paste(s2,'_',s1,sep=''), colnames(vols))]
    else
      vols = vols[,grep(paste(s1,'_',s2,sep=''), colnames(vols))]
    vols = cbind(vols, vols)
    vols$Average_DICE = rowMeans(vols, na.rm = TRUE)
    
    qplot(vols$Average_DICE, geom="histogram",
          binwidth = 1/input$HistogramBins,  
          main = "Average DICE Coeffecient, Two Chosen Sites", 
          xlab = "DICE Coeffecient",  
          fill=I("red"), 
          col=I("blue"), 
          alpha=I(.2),
          xlim=c(-.01,1.01))
    # histogram = ggplot(vols, aes(x=Average_DICE) + geom_freqpoly(bins=input$binwidth, size=2, aes(y=..count..)) + ColorScale + ylab("count"))
    
  })
  
  output$vol_CCC <- renderPlot({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$averages = rowMeans(vols[4:9], na.rm = TRUE)
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols$averages = NULL
    
    relevant_vols = vols[, 5:ncol(vols)-1]
    
    cor_xy = cor(relevant_vols, use="pairwise.complete.obs")
    var_y = rep.row(apply(relevant_vols, 2, na.var), 7)
    var_x = rep.col(apply(relevant_vols, 2, na.var), 7)
    mean_y = rep.row(colMeans(relevant_vols, na.rm=TRUE),7)
    mean_x = rep.col(colMeans(relevant_vols, na.rm=TRUE),7)

    ccc_matrix = 2 * cor_xy * sqrt(var_x) * sqrt(var_y)
    ccc_matrix = ccc_matrix / (var_x + var_y + (mean_x - mean_y)^2)
    
    melted_ccc = melt(ccc_matrix)
    
    hm.palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')  
    
    ggplot(data = melted_ccc, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile()+labs(title = "Concordance Correlation Coeffecient, Absolute Volumes", x="Site", y='Site')+scale_fill_gradientn(colours = hm.palette(100), limits = c(0,1.01))+geom_text(aes(fill = melted_ccc$value, label = signif(melted_ccc$value,2)))
    
  })
  
  output$change_CCC <- renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    selectedSites <- input$ROCsites
    # if (is.null(selectedSites))
    # return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)

    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    vols$averages = NULL
    vols$CMCman_diameter = NULL
    
    vols = replace(vols, vols > 20, NA)
    vols = replace(vols, vols < -1, NA)
  
    relevant_vols = vols[, 3:ncol(vols)]
    print(relevant_vols)
    
    cor_xy = cor(relevant_vols, use="pairwise.complete.obs")
    var_y = rep.row(apply(relevant_vols, 2, na.var), 7)
    var_x = rep.col(apply(relevant_vols, 2, na.var), 7)
    mean_y = rep.row(colMeans(relevant_vols, na.rm=TRUE),7)
    mean_x = rep.col(colMeans(relevant_vols, na.rm=TRUE),7)
    
    ccc_matrix = 2 * cor_xy * sqrt(var_x) * sqrt(var_y)
    ccc_matrix = ccc_matrix / (var_x + var_y + (mean_x - mean_y)^2)
    
    melted_ccc = melt(ccc_matrix)
    
    hm.palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')  
    
    ggplot(data = melted_ccc, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile()+labs(title = "Concordance Correlation Coeffecient, % Volume Change", x="Site", y='Site')+scale_fill_gradientn(colours = hm.palette(100), limits = c(0,1.01))+geom_text(aes(fill = melted_ccc$value, label = signif(melted_ccc$value,2)))
    
  })
  
  output$SizeTable = renderDataTable({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    
    if (input$sizetable_measure == 'first')
      vols = vols[c(TRUE, FALSE),]
    if (input$sizetable_measure == 'second')
      vols = vols[c(FALSE, TRUE),]
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols$manual_averages = rowMeans(vols[,c(6,8)], na.rm = TRUE)
    vols$automatic_averages = rowMeans(vols[,c(4,7,9,10)], na.rm = TRUE)
    
    # baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    # patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    size_vols = vols[,(1:3)]
    size_vols$CUMC_Diameter = vols$CMCman_diameter
    size_vols$CUMC_Volume = vols$CMC_manual
    size_vols$UCLA_Volume = vols$UCL_manual
    size_vols$Ave_Volume_Manual = vols$manual_averages
    size_vols$Ave_Volume_Automatic = vols$automatic_averages
    
    output_vols = size_vols
    
    for (i in 4){
      output_vols[,i][size_vols[,i]<8] = 'Small'
      output_vols[,i][size_vols[,i]>=8] = 'Large'
    }
    
    
    for (i in 5:8){
      output_vols[,i][size_vols[,i]<268] = 'Small'
      output_vols[,i][size_vols[,i]>=268] = 'Large'
    }
  
    output_vols
      
  })
  
  output$SizeHeatmap = renderPlot({
    
  thisSite1 <- input$siteSelected_1
  thisSite2 <- input$siteSelected_2
  if (is.null(thisSite1))
    return(NULL)
  if (is.null(thisSite2))
    return(NULL)
  sqlcmd <- paste("SELECT * from volumes",  sep="")
  v <- dbSendQuery(db, sqlcmd)
  vols <-fetch(v)
  
  if (input$sizetable_measure == 'first')
    vols = vols[c(TRUE, FALSE),]
  if (input$sizetable_measure == 'second')
    vols = vols[c(FALSE, TRUE),]
  
  if (!('benign' %in% input$status))
    vols <- subset(vols, Outcome != 'NC')
  if (!('malignant' %in% input$status))
    vols <- subset(vols, Outcome != 'CC')
  
  vols$manual_averages = rowMeans(vols[,c(6,8)], na.rm = TRUE)
  vols$automatic_averages = rowMeans(vols[,c(4,7,9,10)], na.rm = TRUE)
  
  # baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
  # patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
  
  size_vols = vols[,(1:3)]
  size_vols$CUMC_Diameter = vols$CMCman_diameter
  size_vols$CUMC_Volume = vols$CMC_manual
  size_vols$UCLA_Volume = vols$UCL_manual
  size_vols$Ave_Volume_Manual = vols$manual_averages
  size_vols$Ave_Volume_Automatic = vols$automatic_averages
  
  output_vols = size_vols
  
  for (i in 4){
    output_vols[,i][size_vols[,i]<8] = "Small"
    output_vols[,i][size_vols[,i]>=8] = "Large"
  }
  
  
  for (i in 5:8){
    output_vols[,i][size_vols[,i]<268] = "Small"
    output_vols[,i][size_vols[,i]>=268] = "Large"
  }
  
  output_vols$PID = factor(with(output_vols, paste0(PID, Outcome, Indx)))
  output_vols$Indx = NULL
  output_vols$Outcome = NULL
  
  ggplot(output_vols, aes(x=variable, y=color, label=value, fill=as.factor(value))) + 
    geom_text(colour="black") +
    geom_tile(alpha=0.2)
  
  output_vols.m = melt(output_vols, id.vars="PID", variable_name="size")
  
  gplot = ggplot(output_vols.m, aes(x=variable, y=PID)) 
  gplot + geom_tile(aes(fill = factor(value)), colour="black") + theme(legend.position="bottom")
  
})
  
  output$OutcomeTable = renderDataTable({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    # vols = vols[c("Outcome",s1)]
    
    output_vols = vols[,c(1,2)]
    
    for (i in 3:9) {
    temp_vols = vols[c('Outcome', colnames(vols)[i])]
    # print(temp_vols)
    model = glm(Outcome~.,family=binomial(link='logit'),data=temp_vols)
    p = predict(model, temp_vols)
    pr = prediction(p, temp_vols$Outcome)
    prf = performance(pr, measure="tpr", x.measure="fpr")

    df <- data.frame(cut = prf@alpha.values[[1]], sens = prf@x.values[[1]], spec = prf@y.values[[1]])
    cutoff = df[which.min(((df$sens)^2 + (1-df$spec)^2)^.5), "cut"]
    
    print(df[which.min(((-df$sens)^2 + (1-df$spec)^2)^.5), "sens"])
    print(df[which.min(((-df$sens)^2 + (1-df$spec)^2)^.5), "spec"])
    # print(prf@alpha.values[[1]])
    # print(((-df$sens)^2 + (1-df$spec)^2)^.5)
    # print(df)
    print(colnames(temp_vols)[2])
    print(cutoff)
    
    # print(ROCdf)
    # print(dim(temp_vols))
    # print(pr@labels)
    
    output_vols[colnames(temp_vols)[2]] = temp_vols[,2]
    output_vols[,i][temp_vols[,2]<cutoff] = 'NC'
    output_vols[,i][temp_vols[,2]>=cutoff] = 'CC'
    # print(intersect(as.character(output_vols[,2]) %in% as.character(output_vols[,i])))
    
    matches = 0
    for (j in 1:dim(output_vols)[1]){
      if (!(is.na((as.character(output_vols[j,i]) == as.character(output_vols[j,2]))))){
      if (as.character(output_vols[j,i]) == as.character(output_vols[j,2])){
        matches = matches + 1
      }}
    }
    print(matches/dim(output_vols)[1])
    
    }
    
    output_vols
    
  })
  
  output$OutcomeHeatmap = renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    # vols = vols[c("Outcome",s1)]
    
    output_vols = vols[,c(1,2)]
    
    for (i in 3:9) {
      temp_vols = vols[c('Outcome', colnames(vols)[i])]
      # print(temp_vols)
      model = glm(Outcome~.,family=binomial(link='logit'),data=temp_vols)
      p = predict(model, temp_vols)
      pr = prediction(p, temp_vols$Outcome)
      prf = performance(pr, measure="tpr", x.measure="fpr")
      
      df <- data.frame(cut = prf@alpha.values[[1]], sens = prf@x.values[[1]], spec = prf@y.values[[1]])
      cutoff = df[which.min(((df$sens)^2 + (1-df$spec)^2)^.5), "cut"]
      
      output_vols[colnames(temp_vols)[2]] = temp_vols[,2]
      output_vols[,i][temp_vols[,2]<cutoff] = 'NC'
      output_vols[,i][temp_vols[,2]>=cutoff] = 'CC'
      # print(intersect(as.character(output_vols[,2]) %in% as.character(output_vols[,i])))
      
    }
    
    output_vols$PID = factor(output_vols$PID)
    
    output_vols.m = melt(output_vols, id.vars="PID", variable_name="outcome")
  
    
    gplot = ggplot(output_vols.m, aes(x=variable, y=PID)) 
    gplot + geom_tile(aes(fill = factor(value)), colour="black") + theme(legend.position="bottom")
    
  })
  
  output$ConfusionGroundTruth = renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    output_vols = vols[,c(1,2)]
    
    for (i in 3:9) {
      temp_vols = vols[c('Outcome', colnames(vols)[i])]
      # print(temp_vols)
      model = glm(Outcome~.,family=binomial(link='logit'),data=temp_vols)
      p = predict(model, temp_vols)
      pr = prediction(p, temp_vols$Outcome)
      prf = performance(pr, measure="tpr", x.measure="fpr")
      
      df <- data.frame(cut = prf@alpha.values[[1]], sens = prf@x.values[[1]], spec = prf@y.values[[1]])
      cutoff = df[which.min(((df$sens)^2 + (1-df$spec)^2)^.5), "cut"]
      
      output_vols[colnames(temp_vols)[2]] = temp_vols[,2]
      output_vols[,i][temp_vols[,2]<cutoff] = 'NC'
      output_vols[,i][temp_vols[,2]>=cutoff] = 'CC'
      
    }
    
    conf_matrix = data.frame(matrix(ncol = 3, nrow = 2))
    conf_matrix[] = 0
    conf_matrix[,1] = c('T','F')
    colnames(conf_matrix) = c('GroundTruth','T','F')
    for (j in 1:dim(output_vols)[1]){
      comp1 = as.character(output_vols[j,thisSite1])
      comp2 = as.character(output_vols[j,2])
      if (!(is.na((comp1) == comp2))){
        if (comp1 == comp2){
          if (comp1 == 'CC')
            conf_matrix[1,2] = conf_matrix[1,2] + 1
          else
            conf_matrix[2,3] = conf_matrix[2,3] + 1
        }
        else{
          if (comp1 == 'CC')
            conf_matrix[2,2] = conf_matrix[2,2] + 1
          else
            conf_matrix[1,3] = conf_matrix[1,3] + 1
        }
      }}
    
    conf_matrix$GroundTruth = factor(conf_matrix$GroundTruth)
    
    conf_matrix.m = melt(conf_matrix, id.vars="GroundTruth", variable_name="outcome")
    
    print(conf_matrix.m)
    print(str(conf_matrix.m))
    
    gplot = ggplot(conf_matrix.m, aes(x=variable, y=GroundTruth)) 
    gplot + geom_tile(colour="black", fill="white") + theme(legend.position="bottom") + geom_text(aes(label = round(value, 1))) + labs(title="Confusion Matrix for Site 1 and Ground Truth",x=thisSite1)
    
  })
  
  output$ConfusionTwoSites = renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)    
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    output_vols = vols[,c(1,2)]
    
    for (i in 3:9) {
      temp_vols = vols[c('Outcome', colnames(vols)[i])]
      # print(temp_vols)
      model = glm(Outcome~.,family=binomial(link='logit'),data=temp_vols)
      p = predict(model, temp_vols)
      pr = prediction(p, temp_vols$Outcome)
      prf = performance(pr, measure="tpr", x.measure="fpr")
      
      df <- data.frame(cut = prf@alpha.values[[1]], sens = prf@x.values[[1]], spec = prf@y.values[[1]])
      cutoff = df[which.min(((df$sens)^2 + (1-df$spec)^2)^.5), "cut"]
      
      output_vols[colnames(temp_vols)[2]] = temp_vols[,2]
      output_vols[,i][temp_vols[,2]<cutoff] = 'NC'
      output_vols[,i][temp_vols[,2]>=cutoff] = 'CC'
      
    }
    
    conf_matrix = data.frame(matrix(ncol = 3, nrow = 2))
    conf_matrix[] = 0
    conf_matrix[,1] = c('T','F')
    colnames(conf_matrix) = c('GroundTruth','T','F')
    for (j in 1:dim(output_vols)[1]){
      comp1 = as.character(output_vols[j,thisSite1])
      comp2 = as.character(output_vols[j,thisSite2])
      if (!(is.na((comp1) == comp2))){
        if (comp1 == comp2){
          if (comp1 == 'CC')
            conf_matrix[1,2] = conf_matrix[1,2] + 1
          else
            conf_matrix[2,3] = conf_matrix[2,3] + 1
        }
        else{
          if (comp1 == 'CC')
            conf_matrix[2,2] = conf_matrix[2,2] + 1
          else
            conf_matrix[1,3] = conf_matrix[1,3] + 1
        }
      }}
    
    conf_matrix$GroundTruth = factor(conf_matrix$GroundTruth)
    
    conf_matrix.m = melt(conf_matrix, id.vars="GroundTruth", variable_name="outcome")
    
    print(conf_matrix.m)
    print(str(conf_matrix.m))
    
    gplot = ggplot(conf_matrix.m, aes(x=variable, y=GroundTruth)) 
    gplot + geom_tile(colour="black", fill="white") + theme(legend.position="bottom") + geom_text(aes(label = round(value, 1))) + labs(title="Confusion Matrix for Site 1 and Site 2", x=thisSite1, y=thisSite2)
    
  })
  
  output$KappaMatrix = renderPlot({
    
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    baselinevols <-fetch(v)
    
    baselinevols$averages = rowMeans(baselinevols[4:10], na.rm = TRUE)
    patientaverages = aggregate(baselinevols[11], list(baselinevols$PID), mean)
    
    thisSite1 <- input$siteSelected_1
    if (is.null(thisSite1))
      return(NULL)
    if (input$mode == 'percent'){
      sqlcmd <- paste("SELECT * from volume_changes_percent",  sep="")
    }
    else{
      sqlcmd <- paste("SELECT * from volume_changes",  sep="")
    }
    v <- dbSendQuery(db, sqlcmd)
    vols <<-fetch(v)
    s1<-thisSite1
    
    vols$Outcome = factor(vols$Outcome)
    vols$averages = patientaverages$averages
    vols$CMCman_diameter = baselinevols$CMCman_diameter[c(TRUE, FALSE)]
    
    if (!('small' %in% input$size))
      vols <- subset(vols, CMCman_diameter > 8)
    if (!('large' %in% input$size))
      vols <- subset(vols, CMCman_diameter < 8)
    if (!('medium' %in% input$size))
      vols <- subset(vols, CMCman_diameter <= 8 | CMCman_diameter >= 8)
    
    # vols = vols[c("Outcome",s1)]
    
    output_vols = vols[,c(1,2)]
    
    for (i in 3:9) {
      temp_vols = vols[c('Outcome', colnames(vols)[i])]
      # print(temp_vols)
      model = glm(Outcome~.,family=binomial(link='logit'),data=temp_vols)
      p = predict(model, temp_vols)
      pr = prediction(p, temp_vols$Outcome)
      prf = performance(pr, measure="tpr", x.measure="fpr")
      
      df <- data.frame(cut = prf@alpha.values[[1]], sens = prf@x.values[[1]], spec = prf@y.values[[1]])
      cutoff = df[which.min(((df$sens)^2 + (1-df$spec)^2)^.5), "cut"]
      
      output_vols[colnames(temp_vols)[2]] = temp_vols[,2]
      output_vols[,i][temp_vols[,2]<cutoff] = 'NC'
      output_vols[,i][temp_vols[,2]>=cutoff] = 'CC'
      # print(intersect(as.character(output_vols[,2]) %in% as.character(output_vols[,i])))
      
    }
    
    output_vols$PID = factor(output_vols$PID)
    
    kappa_matrix = data.frame(matrix(ncol = dim(output_vols)[2], nrow = dim(output_vols)[2]-1))
    kappa_matrix[] = 0
    print(kappa_matrix)
    kappa_matrix[,1] = colnames(output_vols)[-1]
    colnames(kappa_matrix)[-1] = colnames(output_vols)[-1]
    print('kappa created!')
    for (i in 2:dim(output_vols)[2]){
      for (j in 2:dim(output_vols)[2]){
        conf_matrix = data.frame(matrix(ncol = 3, nrow = 2))
        conf_matrix[] = 0
        conf_matrix[,1] = c('T','F')
        colnames(conf_matrix) = c('GroundTruth','T','F')
        for (r in 1:dim(output_vols)[1]){
        comp1 = as.character(output_vols[r,i])
        comp2 = as.character(output_vols[r,j])
        if (!(is.na((comp1) == comp2))){
          if (comp1 == comp2){
            if (comp1 == 'CC')
              conf_matrix[1,2] = conf_matrix[1,2] + 1
            else
              conf_matrix[2,3] = conf_matrix[2,3] + 1
          }
          else{
            if (comp1 == 'CC')
              conf_matrix[2,2] = conf_matrix[2,2] + 1
            else
              conf_matrix[1,3] = conf_matrix[1,3] + 1
          }
        }}
        a = conf_matrix[1,2]
        b = conf_matrix[1,3]
        c = conf_matrix[2,2]
        d = conf_matrix[2,3]
        pe = (((a + b)*(a+c)) / (a+b+c+d) +  ((c + d)*(b+d)) / (a+b+c+d)) / (a+b+c+d)
        po = (a + d) / (a+b+c+d)
        kappa_matrix[i-1,j] = (po - pe) / (1 - pe)
      }}
        
    output_vols.m = melt(kappa_matrix, id.vars="X1", variable_name="outcome")
    
    hm.palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')  
    
    gplot = ggplot(output_vols.m, aes(x=variable, y=X1, fill=value))+scale_fill_gradientn(colours = hm.palette(100), limits = c(0,1.01)) 
    gplot + geom_tile(colour="black") + theme(legend.position="bottom") + labs(title="Cohen's Kappa Matrix", x="Site 1", y="Site 2") + geom_text(aes(label = round(value, 3)))
    
    
  })
  
  output$AllTable = renderDataTable({
    
    thisSite1 <- input$siteSelected_1
    thisSite2 <- input$siteSelected_2
    if (is.null(thisSite1))
      return(NULL)
    if (is.null(thisSite2))
      return(NULL)
    sqlcmd <- paste("SELECT * from volumes",  sep="")
    v <- dbSendQuery(db, sqlcmd)
    vols <-fetch(v)
    
    if (!('benign' %in% input$status))
      vols <- subset(vols, Outcome != 'NC')
    if (!('malignant' %in% input$status))
      vols <- subset(vols, Outcome != 'CC')
    
    vols

  })    
  
})