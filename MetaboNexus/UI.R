require(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel("MetaboNexus - an interactive platform for integrated metabolomics analysis"),
  sidebarPanel(

    imageOutput("MetaboNexusPNG",width="150px",height="120px"),
    tags$body(
      tags$h4(
        tags$strong(
          tags$em("Accelerating metabolite discovery")))),
    

    helpText(tags$strong("To start analysis, simply upload a file using the button below")),
          
    fileInput('file1', tags$h5(tags$strong('Step 1: Choose .csv or .txt file from local drive OR click button below to access preprocessed file')),
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    checkboxInput('tempfile',"Click here to access the preprocessed file (UPLC data may take a while)",FALSE),
    verbatimTextOutput('templocation'),
    
    checkboxInput('header', 'Header contains column names', TRUE),
    tags$hr(),
    radioButtons('sep', tags$h5(tags$strong('Step 2: Choose file type')),
                 c("Comma(.csv)"=',',
                   
                    "Tab(.txt/.tsv)"='\t'),
                 'Tab(.txt/.tsv)'),
    radioButtons('quote', 'Quote',
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 'Double Quote'),
    tags$head(tags$style(type="text/css",
                         "label.radio { display: inline-block; margin:0 10 0 0;  }",
                         ".radio input[type=\"radio\"] { float: none; }")),
    tags$hr(),
    tags$body(
      tags$h4("Transform & Transpose Data")),
    
    selectInput('log_10', 'Select log transform', list("log base 10"="log10",
                                               "log base e"="ln",
                                               "log base 2"="log2",
                                               "None"="None")),
    checkboxInput('transpose','Transpose data?',FALSE),
    
    tags$hr(),
    tags$body(
      tags$h4("Y Variable Annotation")),
 
    textInput("Y","Enter the column name containing Y variable",),
    checkboxInput('Supervised', 'Check to use Y variable', FALSE),
    tags$hr(),
    
    selectInput('scaling',tags$h4(tags$strong('Heatmap Data Scaling')),list("Row"="row",
                                                                            "Column"="column",
                                                                            "None"="none"
                                                 
                                                 ))
    
   
    
    
  ),
  
  # Where the panels start
  mainPanel(
    tabsetPanel(
      tabPanel("1a. Uploaded Data",
               checkboxInput('impute','Fill missing or NULL values?',TRUE),
               checkboxInput('zerovar','Remove variables with zero/near zero variance',TRUE),
               tags$hr(),
               uiOutput("table_status"),
               checkboxInput('full_table',"Show full table? (may take a while to load)",FALSE),
                
               #helpText(tags$strong("Remember to enter column name of Y variable")),
               tableOutput("RAW"),imageOutput("welcome",height="100%")),
      tabPanel("1b. Annotate sample class",
               
               sidebarPanel(
                 radioButtons('xcms_anno', tags$h4('Annotation of groupings/sample classes'),
                              c("Use XCMS annotation"='xcms_yes',"Use Y variable"="y_variable",
                                'Manual annotation'='manual_anno')),
                 uiOutput("class_confirmation"),
                 
                 uiOutput("manual_classes")
                 ),
               
               
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               
              
               uiOutput("annotate")),
      tabPanel("2. PCA",
               tags$body(
                 tags$h5("Principal Components Analysis")),
               checkboxInput('PCA_scale', 'Scale Data?', TRUE),
               checkboxInput("PCA_3D_scores","Display score plot in 3D?",FALSE),
               checkboxInput("PCA_3D_loads","Display loadings plot in 3D?",FALSE),
               tags$br(),
               tags$hr(),
               tags$body(
                 tags$h4("Plots")),
               selectInput('pca_scoreload', 'Select plot type to load', list("Score plot"="score",
                                                                            "Loadings plot"="loadings"
                                                                             )),
               checkboxInput("PCA_labels","Show data labels?",TRUE),
               checkboxInput("edit_points","Edit color of points?",FALSE),
               uiOutput("editpoints"),
               plotOutput("NIPALS",height="100%"),
               plotOutput("NIPALS_var",height="100%"),
               plotOutput("NIPALS_3D",height="100%"),
               plotOutput("NIPALS_3D1",height="100%")
               ),
      tabPanel("3. PCA Scree Plot",
               tags$body(
                 tags$h5("PCA Screeplot")),
               sidebarPanel(
                 sliderInput("PCAcomp", "Select number of components:", 
                           min = 2, max = 18, value = 3, step= 1),
                checkboxInput("cum","Display cumulative eigenvalues?",TRUE),
                checkboxInput('line', 'Display 80% cutoff line?', TRUE)),
               plotOutput("RAWScree", height="100%")
                
               ),

      tabPanel("4. PLS(-DA)",

               tags$head(tags$style(type="text/css",
                                    "label.radio { display: inline-block; margin:0 10 0 0;  }",
                                    ".radio input[type=\"radio\"] { float: none; }")),
               sidebarPanel(
                 helpText(tags$strong("Partial Least Squares (-Discriminant Analysis)")),
                 helpText("*Suitable for two-group analysis"),
                 sliderInput("ncomp", "Select number of components:", 
                             min = 2, max = 18, value = 3, step= 1),
                 checkboxInput('PLS_3D_scores','Display score plot in 3D?', FALSE),
                 checkboxInput('PLS_3D_loads','Display loadings plot in 3D?', FALSE)),

               
               sidebarPanel(
                 

                 helpText(tags$strong("For Q2 ")),
                 tableOutput("Q2")),

                sidebarPanel(
                  helpText(tags$strong("For all R2")),
                  tableOutput("expvar")),
                
                 
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$hr(),
               tags$body(
                 tags$h4("Score Plot")
               ),
               selectInput('pls_scoreload', 'Select plot type to load', list("Score plot"="score",
                                                                             "Loadings plot"="loadings"
               )),
               checkboxInput("PLS_labels","Show data labels?",TRUE),
               plotOutput("pls",height="100%"),
              
               
               tags$br(),
               plotOutput("pls3D"),
               plotOutput("pls_loadings3D")

               ),
      
               
      tabPanel("5. PLS(-DA) VIP Values",
               
               tags$h5("This page gives a quick overview of influential features"),
               tags$h5("Use the displayed VIP values to estimate suitable cut-off values for Univariate Analysis & Heatmap"),
               plotOutput("VIP_chart",height="100%"),
               helpText(tags$em("Only the top 30 variables are plotted")),
               tags$hr()

               ),
      
      tabPanel("6. PLS Model Diagnostics",
               sidebarPanel(
                 
                 selectInput('diag_type',tags$h5(tags$strong('Plot model diagnostics')),list("Predicted Y"="y.pred",
                                                                                             "PRESS"="PRESS",
                                                                                             "RSS"="RSS",
                                                                                             "R2"="R2",
                                                                                             "R2X"="R2X",
                                                                                             "R2Y"="R2Y",
                                                                                             "R2Xcum"="R2Xcum",
                                                                                             "R2Ycum"="R2Ycum",
                                                                                             "Q2"="Q2",
                                                                                             "Q2cum"="Q2cum",
                                                                                             #"Predicted Y"="y.pred",
                                                                                             "Residuals"="resid")),
                 checkboxInput('show_all','Show all plots?',FALSE)
                 ),
               plotOutput("mod_diag", height="100%")),
               
      
      tabPanel("7. Random Forest",
               sidebarPanel(
                 tags$body(tags$h5("Random Forest")),
                 helpText("*Suitable for analysis of two or more groups"),
                 textInput("trees","Enter no. of trees",500),
                 checkboxInput("tree_var_edit","Click to edit default no. of variables tried",FALSE),
#                  textInput("tree_var","Enter no. of variables",21)  
                 uiOutput("tree_var_in"),
                 helpText(tags$em("*For large datasets with >1000 variables, analysis might take a while...")),
                 checkboxInput("RF_labels","Show data labels?",TRUE)
                 ),
               tags$div(class="span1"),
               tags$div(class="span1"),
               
               sidebarPanel(
                 uiOutput("RF_stats")
               ),
               
               
               
               
               plotOutput("MDS",height="100%")),
      
      
      tabPanel("8. Random Forest VIP",
               downloadButton('downloadData2', 'Download'),
               tags$h5("This page gives a quick overview of influential features"),
               tags$h5("Use the displayed variable importance to estimate suitable cut-off values for Univariate Analysis & Heatmap"),
               plotOutput("RF_VIP",height="100%"),
              helpText(tags$em("Only the top 30 variables are plotted"))),
     
      tabPanel("9. Univariate analysis",
                              
               sidebarPanel(
                                  
                 tags$body(tags$h4(tags$strong("Select method for univariate analysis"))),
                 
                 tags$br(),
                 radioButtons('test_type', 'Select Test',
                              c("Kruskal-Wallis"='KW',
                                "ANOVA (1-way)"='anov',
                                "Mann-Whitney U test" = "man",
                                "Student's t-test" = "tt",
                                "None"="none"
                              ),
                              "None"),
                 checkboxInput('p.adjust','Apply False Discovery Rate?',FALSE),tags$br(),
                 helpText(tags$i("p-values will be displayed next to VIP values")),
                 helpText(tags$i("Large datasets may take a while to complete..."))),
               sidebarPanel(
                 tags$body(tags$h4(tags$strong("Select criteria for shortlisting features/peaks"))),
                 radioButtons('varimp_cutoff',tags$h5('Type of variable importance cut-off'),
                              c("None"="None","PLS"="PLS","Random Forest"="RF")),
                 uiOutput("varimp_input"),
                 numericInput("uni_pval",tags$h5("Enter cut-off for p-value"),0.05),
                 checkboxInput("uni_subset","Shortlist variables based on p-value",TRUE)
                 
                 ),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
               helpText("Choose a shortlisted feature and search for putative matches using the next tab"),
               tags$hr(),
               uiOutput("boxplot_control"),
               plotOutput("boxplots"),
               plotOutput("ROC"),
               tags$br(),
               downloadButton('downloadData4', 'Download shortlisted peak'),
               tableOutput("uni_var")
               ),
      
      tabPanel("10. Metabolite Search",
               tags$body(
                 tags$h4("Human Metabolome Database & Metlin Search"),
                 tags$h5("Search for metabolite spectra & pathway using m/z values or names")),
               tags$hr(),
               radioButtons("qType",tags$h5("Step 1: Search by"),
                            c("Variable Name"="var_mz", "m/z value"="mz",
                              "Metabolite name"="name")),
               
               uiOutput("var_mz1"),
               tags$hr(),
               radioButtons("ion_mode",tags$h5("Step 2: Select Adducts & Tolerance"),c("Positive"="pos", "Negative"="neg")),
               
               uiOutput("adduct_ctrl"),
               
               
               
               numericInput("tol","Tolerance (+/- m/z value)",0.3),
               tags$hr(),
               tags$h5("Step 3: Browse Results"),
               tags$h5(tags$em("Metlin results")),
               htmlOutput("metlin"),
               tags$h5(tags$em("HMDB results")),
               helpText("*Copy and paste the HMDB ID into the box to generate MS/MS spectra & pathways"),
               
               tableOutput("hits"),
               
               tags$hr(),
               
               textInput("hmdbid",tags$h4("Enter HMDB ID to display list of associated pathways and MS/MS spectra:"),),
               
               radioButtons("MassBank","Select reference database",c("HMDB"="HMDB","MassBank (referenced within HMDB)"="MassBank")),
               
               uiOutput("pwTable"),
               
               plotOutput("specPlot",height="100%")
               
      ),
      
      tabPanel("11. Heatmap",
               sidebarPanel(
                 tags$h5("Heatmap Settings"),
                 selectInput("heatmap.col","Choose heatmap colors",list("Green to Red"="greenred",
                                                                        "Red to Green"="redgreen",
                                                                        "Blue | White | Pink"="bluewhitepink",
                                                                        "Blue | White | Red"="bluewhitered",
                                                                        "Heat colors"="heatcol")),
                 selectInput("dendro_ctrls","Dendrogram control",list("Both"="both",
                                                                      "Row dendrogram"="row",
                                                                      "Column dendrogram"="col",
                                                                      "None"="none"
                 )),
                 selectInput("ordering","Re-ordering of elements",list("By row & column means"="both",
                                                                       "By row means"="rowmeans",
                                                                       "By column means"="colmeans",
                                                                       "No ordering"="no_order")),
                 radioButtons('varimp_cutoff1',tags$h5('Apply variable importance cut-off'),
                              c("None"="None","PLS"="PLS","Random Forest"="RF")),
                 uiOutput("varimp_input1"),
                 numericInput("uni_pval1","Enter cut-off for p-value",0.05),
                 checkboxInput("uni_subset1","Shortlist variables based on p-value",TRUE),
                 checkboxInput("heat_trans","Transpose heatmap",TRUE),
                 helpText(tags$em("* If FDR correction was applied, p-value will be modified accordingly"))
               ),
               
               
               plotOutput("QHeatmap",height="2500px")),
      
      tabPanel("12. Overall Results",
               downloadButton('downloadData3', 'Download'),
               tableOutput("overall")),
               
     
      
      
      tabPanel("Log",
               tags$body(
                 tags$h5("For reproducibility of analysis, we have compile all the settings used for this session"),
                 tags$h5("Click on download to obtain the list & share with your collaborators")),
               downloadButton('download_log', 'Download'),
               tableOutput("all.settings")),
      
      tabPanel("Help",
               tags$body(
                 tags$h4("How do I use MetaboNexus?"),
                 tags$em("More information can be found in the tutorial"),
                 tags$br(),
                 tags$br(),
                 tags$h5("1a. Pre-process your data using our pre-processing module (runs on XCMS):"),
                 tags$h5("OR"),
                 tags$h5("1b. Prepare a .csv/txt file with your data arranged as follows:"),
                 tags$p("  - Variables as columns"),
                 tags$p("  - Observations as rows"),
                 tags$br(),
                 tags$h5("2. Upload the file and toggle the checkboxes to your desired settings"),
                 tags$p("  - Header: Choose this to let first row of uploaded text to be column names"),
                 tags$p("  - Tranpose?: If your variables are vertically aligned, use this to flip it to horizontal alignment. Data analysis only works in this manner"),
                 tags$p("  - Apply log transformation?: This allows you to transform your data if it is not normally distributed"),
                 
                 tags$br(),
                 tags$h5("3. Annotating samples"),
                 tags$p(" - Peak data coming from our pre-processed module can be easily annotated in the 'Annotate sample class' tab"),
                 tags$p(" - Peak data loaded in as .csv/.tsv/.txt will require a Y variable (e.g. stated as 0 or 1) to represent sample class"),
                 tags$br(),
                 tags$h5("4. Performing Unsupervised Analysis"),
                 
                 tags$p("  - Principal components analysis is an unsupervised method for multivariate analysis and the NIPALS algorithm is used here"),
                 
                 tags$br(),
                 tags$h5("5. Performing Supervised Analysis with PLS(-DA)/Random Forest"),
                 
                 tags$h5("5.1 Partial Least Squares (-Discriminant Analysis)"),
                 tags$p(" - PLS is a regression method that finds relations between X (predictors/metabolites) and Y (responses/phenotype)"),
                 tags$p(" - When Y is binary (i.e. 0 and 1), PLS can be used to perform discriminant analysis (hence DA) on the two sample classes"),
                 tags$h5("5.2 Random Forest (RF)"),
                 tags$p(" - RF is an ensemble learning method for classification (and also regression)"),
                 tags$p(" - Here we use RF for classification purposes instead of regression since most biological applications are categorical"),
                 tags$p(" - The out-of-bag (OOB) error provides unbiased estimate of the classification error"),
                 tags$p(" - The no. of variables used for RF is the sqrt(no. of variables), which is automatically calculated"),
                 tags$br(),
                 tags$h5("6. Univariate analysis"),
                 tags$p(" - Here you can perform parametric and non-parametric tests to test for statistical significance"),
                 tags$p(" - Available tests are for 2 sample classes (t-test & Mann-Whitney U) and >2 sample classes (ANOVA & Kruskal-Wallis)"),
                 tags$p(" - Additionally, Receiver Operator Characteristic (ROC) will be automatically performed for 2 class samples"),
                 tags$h5("6.1 Shortlisting peaks"),
                 tags$p(" - Here you can shortlist peaks based on variable importance* (e.g. VIP values) and p-values in a simple fast manner"),
                 tags$p(" - The shortlisted peaks are downloadable for further analysis"),
                 tags$br(),
                 tags$h5("7. Metabolite Search"),
                 tags$p(" - You can search for putative identities by selecting the feature of interest from a drop-down menu"),
                 tags$p(" - And selecting the appropriate molecular adduct (e.g. [M+H]+, [M+Na]+)"),
                 tags$h5("7.1 Metabolite Search using 3 databases"),
                 tags$p(" - After selecting the feature & adduct, the m/z value will be matched to HMDB, MassBank and Metlin"),
                 tags$p(" - For Metlin, users can click on a customized hyperlink that brings them to Metlin to view metabolite matches"),
                 tags$p(" - For MassBank, HMDB has incorporated chemical information from MassBank hence by searching HMDB, one can match data from both databases "),
                 tags$h5("7.2 Tandem MS/MS spectra"),
                 tags$p(" - By entering the HMDB ID, users can display tandem MS/MS spectra from HMDB and MassBank"),
                 tags$p(" - For Metlin, the hyperlink would also contain tandem MS/MS spectra if available"),
                 tags$br(),
                 tags$h5("8. Heatmap Data Scaling"),
                 tags$p(" - The heatmap can be customized using the buttons made available in the heatmap tab"),
                 tags$p(" - To scale the data in heatmap, use the settings in the leftmost panel"),
                 tags$br(),
                 tags$h5("9. Log"),
                 tags$p(" - To ensure reproducibility of work, a log is generated and can be downloaded to keep track of settings used"),
                 tags$br(),
                 
                 
                 tags$br(),
                 tags$h5("* Feature selection through variable importance"),
                 tags$p("  - The variable importance is an indication of how influential a particular variable is for the model "),
                 tags$p("  - For PLS(-DA), we select important variables based on a cut-off of >1 in VIP"),
                 tags$p("  - For Random Forest, the top 30 important variables are displayed to guide in feature selection"),
                 tags$p("  - These variables are influential as they result in a decrease in accuracy or Gini coefficient when removed from the model"))
               ),
  
    
    tabPanel("About",
             tags$body(
               tags$h4("Our motivation"),
               tags$p("Metabolomics analysis typically require users to switch between programs to generate data and visuals"),
               tags$p("By bringing together these necessary tools and pipelining it, analysis is now made simpler and seamless for any user"),
               tags$br(),
               tags$h4("Sharing your data along with this program"),
               tags$p("The open source nature of the program makes it easier for collaborators to share and understand each other's methodology"),
               tags$p("Users can right click on images to save or click on 'Download' buttons to retrieve .csv tables. A log file is provided for reproducibility of analysis. "),
               tags$p("This may promote better understanding and planning among groups in future"),
               
               tags$h4("Referencing MetaboNexus"),
               tags$p("If you have used MetaboNexus in your research, do us a a favour and cite our paper:"),
               tags$a(href="http://link.springer.com/article/10.1007/s11306-014-0648-8",
                      "SM. Huang et al. MetaboNexus - an interactive platform for metabolomics analysis. Metabolomics 6 (2014): 1098 - 1093.",
                      target="_blank")
               
               
             )
           
    )
    
    ))
     
      
))