sites = c("UCLA"="UCL", "UCLA_Manual"="UCL_manual", "CMC"="CMC", "CMC_Manual" = "CMC_manual", "DFCI"="DFCI", "MCC"="MCC", "UMC"="UMC")
sites_colorcoded = c("UCLA, Red"="UCL", "UCLA_Manual, Green"="UCL_manual", "CMC, Blue"="CMC", "CMC_Manual, Brown" = "CMC_manual", "DFCI, Orange"="DFCI", "MCC, Purple"="MCC", "UMC, Yellow"="UMC", "CMC Manual Diameter, Black" = "CMCman_diameter")

shinyUI(
  fluidPage(
    titlePanel("Site Selection"),
    sidebarLayout(
      sidebarPanel(
        selectInput("siteSelected_1", label = h3("Select site 1"), 
                    choices = sites, selected="CMC"),
        selectInput("siteSelected_2", label = h3("Select site 2"), 
                    choices = sites, selected="MCC"),
        radioButtons("mode", "Change Measure:",c("Percent Difference" = "percent", "Absolute Difference" = "absolute")),
        checkboxGroupInput("size", "Nodule Size:",c("Above 8mm Diameter (268 mm^3)" = "large", "Below 8mm Diameter (268 mm^3)" = "small"), selected=c('large','small')),
        checkboxGroupInput("status", "Nodule Status:",c("Benign" = "benign", "Malignant" = "malignant"), selected=c('benign', 'malignant')),
        downloadButton('downloadAll', 'Download all Volumes')
      ), 
      mainPanel(
        tabsetPanel(type = "tabs", 
                    tabPanel("Volume Estimates",     
                             
                             plotOutput(outputId = "plot1", height = "400px"),
                             plotOutput(outputId = "plot2", height = "400px")
                    ),
                    tabPanel("Volume Histogram",     
                             
                             plotOutput(outputId = "vol_histogram_sites", height = "400px"),
                             plotOutput(outputId = "vol_histogram_status", height = "400px"),
                             plotOutput(outputId = "vol_histogram_date", height = "400px"),
                             fluidRow(column(10,radioButtons("histogram_mode", "Chart Type:",c("Density" = "density", "Cumulative Density" = "cumulative", "Histogram" = "histogram")))),
                             fluidRow(column(10,
                                             sliderInput("SliderHistogramLimits", "Histogram Limits", 
                                                         min = 0, max = 10000, value = c(0,2000), step= 50))),
                             fluidRow(column(10,
                                             sliderInput("SliderHistogramBins", "Histogram Bins", 
                                                         min = 0, max = 100, value = 20, step= 1)))
                    ),
                    tabPanel("Volume Change Estimates",
                             
                             plotOutput(outputId = "changeplot1", height = "350px"),
                             plotOutput(outputId = "changeplot2", height = "350px"),
                             fluidRow(column(8,
                                             sliderInput("SliderChangePercent", "Cutoff Values, Percent", 
                                                         min = -1, max = 20, value = c(-1,20), step= 0.1),
                                             sliderInput("SliderChangeAbsolute", "Cutoff Values, Absolute", 
                                                         min = -3000, max = 3000, value = c(-3000,3000), step= 100))
                                      
                             )),
                    tabPanel("ROC Curves",
                             
                             plotOutput(outputId = "ROCplot", height = "350px"),
                             plotOutput(outputId = "LogisticPlot", height = "350px"),
                             fluidRow(column(8,
                                             sliderInput("SliderChangePercent2", "Cutoff Values, Percent", 
                                                         min = -1, max = 20, value = c(-1,20), step= 0.1),
                                             sliderInput("SliderChangeAbsolute2", "Cutoff Values, Absolute", 
                                                         min = -3000, max = 3000, value = c(-3000,3000), step= 100))
                                      
                             )
                             
                    ),
                    tabPanel("ROC Curve Comparison",
                             
                             plotOutput(outputId = "ROCCompareplot", height = "350px"),
                             fluidRow(column(12,
                                             checkboxGroupInput("ROCsites", "Sites",
                                                                sites_colorcoded, selected=c("UCL_manual", "DFCI"))
                             ))),
                    tabPanel("DICE Histogram",
                             
                             plotOutput(outputId = "AllDICE", height = "350px"),
                             plotOutput(outputId = "TwoDICE", height = "350px"),
                             fluidRow(column(12,
                                             sliderInput("HistogramBins", "Number of Bins", 
                                                         min = 2, max = 150, value = 50, step= 1))
                                      
                             )
                             
                    ),
                    tabPanel("CCC Matrix",
                             
                             plotOutput(outputId = "vol_CCC", height = "400px"),
                             plotOutput(outputId = "change_CCC", height = "400px")
    
                    ),
                    tabPanel("Small/Large Category Comparison",
                             
                             fluidRow(column(10,radioButtons("sizetable_measure", "Size Measure",c("First Visit" = "first", "Second Visit" = "second", "Both" = "both")))),
                             plotOutput(outputId = "SizeHeatmap", height = "1200px"),
                             dataTableOutput('SizeTable')
                             
                    ),
                    tabPanel("Outcome Prediction",
                             
                             plotOutput(outputId = "OutcomeHeatmap", height = "1200px"),
                             dataTableOutput('OutcomeTable')
                             
                    ),
                    tabPanel("Confusion Matrices",
                             
                             plotOutput(outputId = "ConfusionGroundTruth", height = "200px"),
                             plotOutput(outputId = "ConfusionTwoSites", height = "200px"),
                             plotOutput(outputId = "KappaMatrix", height = "400px")
                             
                    ),                    
                    tabPanel("All Data",
                             
                             dataTableOutput('AllTable')
                             
                    )
        )
      ))))