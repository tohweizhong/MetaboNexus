library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("MetaboNexus pre-processing module"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    tags$h4("Raw data pre-processed here will be linked to MetaboNexus"),
    numericInput("cores",tags$h5("No. of cores detected (please change if necessary)"),as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"])),
    checkboxInput('choose_dir', tags$h5('Choose your working directory (where raw files are stored)'), FALSE),
    helpText("*Suitable file formats are .CDF, .mzXML and .mzData"),
    textOutput("Raw_data"),
    tableOutput("File_list"),
   # tableOutput("phenoData"),
#     checkboxInput("sample_class",tags$h5("Click here to edit sample class/phenotype (optional)"),FALSE),
    #uiOutput("sample_pheno"),
    tags$hr(),
    tags$h4(tags$u("Settings for peak detection")),
    helpText("Pre-processing is performed using the XCMS R package"),
    tags$br(),
    numericInput("minfraction",tags$h5("Peaks need to appear in at least x % of each sample class"),50),

#     textInput("snames","Sample names",),
#     textInput("sclass","Sample class",),
#     textInput("pheno","Phenotype Data",),
    

    selectInput("presets",tags$h5("Option A: Select preset or manual settings for your instrument"),
                list( "Manual settings"="manual",
                      "HPLC/Q-TOF"="HP_QTOF",
                      "HPLC/Q-TOF (high res)"="HP_QTOF.high",
                      "HPLC/Orbitrap"="HP_Orbi",
                      "UPLC/Q-TOF"="UP_QTOF",
                      "UPLC/Q-TOF (high res)"="UP_QTOF.high",
                      "UPLC/Orbitrap"="UP_Orbi"
                    )),
                                                                              
    htmlOutput("presets_details"),
    tags$hr(),
    tags$h5("Option B: Manual settings"),
#     selectInput("profmethod","Profiling method",list("bin"="bin",
#                                                      "binlin"="binlin",
#                                                      "binlinbase"="binlinbase",
#                                                      "intlin"="intlin")),
    selectInput("findPeaks",tags$h5("Choose findPeaks method"),list("matchedFilter"="matchedFilter",
                                                           "centWave"="centWave")),
    uiOutput("fwhm"),
    uiOutput("max"),
    uiOutput("snthresh"),
    uiOutput("step"),
    uiOutput("steps"),
    uiOutput("ppm"),
    uiOutput("peakwidth"),
    uiOutput("prefilter1"),
    uiOutput("prefilter2"),
    uiOutput("snthresh1"),
    uiOutput("integrate"),
    uiOutput("mzdiff"),
    uiOutput("fitgauss"),
    tags$hr(),
    uiOutput("retention"),
    uiOutput("retention_span")
    
#     checkboxInput("start",tags$h4("Start XCMS processing"),FALSE)
   
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Raw data processing (xset)",
               tableOutput("phenoData"),
#                checkboxInput("start",tags$h4("Start XCMS processing"),FALSE),
               
               htmlOutput("xset_processing"),
               htmlOutput("start_xcms"),
               tags$hr(),
               htmlOutput("xset_progress")
#                verbatimTextOutput("xset_progress2")
               ),
#       tabPanel("Match peaks across samples (group)",
#                htmlOutput("group")),
      tabPanel("Retention Time Correction",
               htmlOutput("retcor"),
               plotOutput("retcor_plot",height="100%")),
#       tabPanel("Re-match peaks across samples",
#                htmlOutput("Regroup")),
#       tabPanel("Fill peaks",
#                htmlOutput("Fillpeaks")),
      tabPanel("Peak Table",
               helpText("Please wait while peak table loads"),
               downloadButton("downloadPeaklist","Download Peak List"),
               verbatimTextOutput("temp"),
               selectInput("norm_method","Select normalization method",c("None"="None",
                                                                         "Internal Standard"="ISTD",
                                                                         "Quantile Normalization"="quantile")),
               uiOutput("norm_out"),
               checkboxInput("apply_method","Apply method selected",FALSE),
               tableOutput("Peaklist"))

      
      )
    
  )
))
