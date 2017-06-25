library(shiny)
library(xcms)
library(snow) #prompt user to install snow
library(preprocessCore)
# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
  
  location<-function(){
    if(input$choose_dir==TRUE){
      cdfpath<-choose.dir()
    }
  }
  Raw_data <- reactive({
    if(input$choose_dir==TRUE){
    location()
    }
  })
  output$Raw_data<-renderText({
    Raw_data()
  })
  
  #optional?
  output$sample_pheno<-renderUI({
    if(input$sample_class==TRUE){
      list(textInput("snames","Sample names",),
           textInput("sclass","Sample class",),
           textInput("pheno","Phenotype Data",))
    }

  })
  #==========================================================
  output$fwhm<-renderUI({
    if(input$findPeaks=="matchedFilter"){
      numericInput("fwhm","FWHM",30)
    } else NULL
  })
  
  output$max<-renderUI({
    if(input$findPeaks=="matchedFilter"){
      numericInput("max","Max no. of peaks per EIC",5)
    } else NULL
  })
  
  output$snthresh<-renderUI({
    if(input$findPeaks=="matchedFilter"){
      numericInput("snthresh","Signal-to-Noise threshold",10)
    } else NULL
  })
  
  output$step<-renderUI({
    if(input$findPeaks=="matchedFilter"){
      numericInput("step","Step-size",0.1)
    } else NULL
  })
  
  output$steps<-renderUI({
    if(input$findPeaks=="matchedFilter"){
      numericInput("steps","No. of steps to merge prior to filtration",2)
    } else NULL
  })
  #========================================================
  #centWave
  output$ppm<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("ppm","Set ppm",15)
    } else NULL
  })
  
  output$ppm<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("ppm","Set ppm",15)
    } else NULL
  })
  
  output$peakwidth<-renderUI({
    if(input$findPeaks=="centWave"){
      sliderInput("peakwidth", "Peakwidth:",
                   min = 1, max = 100, value = c(10,60))
    } else NULL
  })
  
  output$snthresh1<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("snthresh","Signal-to-Noise threshold",10)
    } else NULL
  })
  
  output$prefilter1<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("prefilter1","Contains at least N peaks",3)
    } else NULL
  })
  
  output$prefilter2<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("prefilter2","...of given Intensity",100)
    } else NULL
  })
  
  output$integrate<-renderUI({
    if(input$findPeaks=="centWave"){
      selectInput("integrate","Integration method based on",list("Mexican hat (robust to noise)"=1,
                                                                 "Real data (prone to noise)"=2))
    } else NULL
  })
  
  output$mzdiff<-renderUI({
    if(input$findPeaks=="centWave"){
      numericInput("mzdiff","Minimum difference in m/z for peaks with overlapping RT",-0.001)
    } else NULL
  })
  
  output$fitgauss<-renderUI({
    if(input$findPeaks=="centWave"){
      checkboxInput("fitgauss","Fit to Gaussian Model?",FALSE)
    } else NULL
  })
  


  output$retention<-renderUI({
    OUT<-list()
    OUT[[1]]<-tags$h5("Grouping parameters")
    OUT[[2]]<-numericInput("bandwidth","Select bandwidth (e.g. 2, 5, 10 or 30)",5)
    OUT[[3]]<-numericInput("mzwid","Width of overlapping m/z slices (e.g. 0.015, 0.025, 0.25)",0.25)
    OUT[[4]]<-selectInput("retcor_method",tags$h5("RT correction method"),list("Obiwarp"="obiwarp",
                                                                     "Loess"="loess",
                                                                     "Linear"="linear"))
    OUT
                                 
  })
  
  output$retention_span<-renderUI({
    if(input$retcor_method=="loess"|input$retcor_method=="linear")
      numericInput("rt_span","Enter retention time span",0.2)
  })
  
  output$File_list<-renderTable({
    if(length(Raw_data()!=0)){
      LIST<-data.frame(list.files(Raw_data()))
      colnames(LIST)[1]<-"Folders/File Names"
      LIST
    }
    else{return(NULL)}
  })
  
  phenoData<-reactive({
    if(input$choose_dir==TRUE){
      files<-Raw_data()
      #from XMCS
      filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                       "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
      filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), 
                           collapse = "|")
      if (is.null(files)) 
        files <- getwd()
      info <- file.info(files)
      listed <- list.files(files[info$isdir], pattern = filepattern, 
                           recursive = TRUE, full.names = TRUE)
      files <- c(files[!info$isdir], listed)
      phenoData<-xcms:::phenoDataFromPaths(files)
      phenoData<-t(phenoData)
    }
    
  })
  
  output$phenoData<-renderTable({
      
    pd<-data.frame(phenoData())
    #sink("XCMStemp/phenoData.txt")
    write.table(pd,"XCMStemp/phenoData.txt")
    print(pd)
    #sink()
      if(is.null(pd)){
        NULL
      }
      else
        pd
  })
  
  start_xcms<-function(){
    if(input$start==TRUE){
      
      cdfpath<-Raw_data()
      cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
      
      
      if(input$presets=="HP_QTOF"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=30,peakwidth=c(10,60),prefilter=c(3,1000),snthresh=10,integrate=1,noise=200)
      }
      
      else if(input$presets=="HP_QTOF.high"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=15,peakwidth=c(10,60),prefilter=c(3,1000),snthresh=10,integrate=1,noise=200)
      }
      
      else if(input$presets=="HP_Orbi"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=2.5,peakwidth=c(10,60),prefilter=c(3,5000),snthresh=10,integrate=1,noise=200)
      }
      
      else if(input$presets=="UP_QTOF"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=30,peakwidth=c(5,20),prefilter=c(3,1000),snthresh=10,integrate=1,noise=200) #HELLO
      }
      
      else if(input$presets=="UP_QTOF.high"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=15,peakwidth=c(5,20),prefilter=c(3,1000),snthresh=10,integrate=1,noise=200)
      }
      
      else if(input$presets=="UP_Orbi"){
        xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=2.5,peakwidth=c(5,20),prefilter=c(3,5000),snthresh=10,integrate=1,noise=200)
      }
      
      else if(input$presets=="manual"){
        if(input$findPeaks=="matchedFilter"){
          xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="matchedFilter",fwhm=as.numeric(input$fwhm),
                          max=as.numeric(input$max),snthresh=as.numeric(input$snthresh),
                          step=as.numeric(input$step),steps=as.numeric(input$steps))
          
        }
        else if(input$findPeaks=="centWave"){
          xset <- xcmsSet(cdffiles,nSlaves=input$cores,method="centWave",ppm=as.numeric(input$ppm),
                          peakwidth=as.numeric(input$peakwidth),snthresh=as.numeric(input$snthresh),
                          prefilter=c(as.numeric(input$prefilter1),as.numeric(input$prefilter2)),
                          integrate=as.numeric(input$integrate),mzdiff=as.numeric(input$mzdiff),fitgauss=input$fitgauss)
        }
      }
    }
    else NULL
  }
    
  output$presets_details<-renderText({
    if(input$presets=="manual"){
      Details<-paste("<h5><u><i>Please adjust peak detection settings below (Option B)</i></u></h5>")
      Details<-paste0("<div style='border:1px solid black'><h5><u><i>Please adjust peak detection settings below using Option B</i></u></h5></div>")
    }
    else{
      if(input$presets=="HP_QTOF"){
        ppm<-30;peakwidth<-"10 to 60";bw<-5;mzwid<-0.025;prefilter<-"c(3,1000) [at least 3 peaks of height=1000]"
      }
      else if(input$presets=="HP_QTOF.high"){
        ppm<-15;peakwidth<-"10 to 60";bw<-5;mzwid<-0.015;prefilter<-"c(3,1000) [at least 3 peaks of height=1000]"
      }
      else if(input$presets=="HP_Orbi"){
        ppm<-2.5;peakwidth<-"10 to 60";bw<-5;mzwid<-0.015;prefilter<-"c(3,5000) [at least 3 peaks of height=5000]"
      }
      else if(input$presets=="UP_QTOF"){
        ppm<-30;peakwidth<-"5 to 20";bw<-2;mzwid<-0.025;prefilter<-"c(3,1000) [at least 3 peaks of height=1000]"
      }
      else if(input$presets=="UP_QTOF.high"){
        ppm<-15;peakwidth<-"5 to 20";bw<-2;mzwid<-0.015;prefilter<-"c(3,1000) [at least 3 peaks of height=1000]"
      }
      else if(input$presets=="UP_Orbi"){
        ppm<-2.5;peakwidth<-"5 to 20";bw<-2;mzwid<-0.015;prefilter<-"c(3,5000) [at least 3 peaks of height=5000]"
      }
      
      Details1<-paste("<h5><u><i>Peak detection settings</u></i></h5>","\n<p>Method=centWave</p> \n<p>ppm=",ppm,"</p>\nPeak Width=",peakwidth,
                     "</p>\nPrefilter=",prefilter,"</p>\n<h5><u><i>Retention Time Correction</u></i></h5> \n<p>Method=Obiwarp</p> \n<p>bw=",bw,"</p>\n<p>mzwid=",mzwid,"</p>",
                     sep="")
      Details<-paste("<div style='border:1px solid black'>",Details1,"</div>")
    }
    Details
    
  })
  
  Processing<-reactive({
    if(input$start==TRUE)
    xset<-start_xcms()
    else NULL
      
  })
  #===========================================================
  Group<-reactive({
    minFRAC<-as.numeric(input$minfraction)/100
    if(input$presets=="HP_QTOF"){
      xset<-group(Processing(),bw=5,mzwid=0.025,minfrac=minFRAC)
    }
    else if(input$presets=="HP_QTOF.high"){
      xset<-group(Processing(),bw=5,mzwid=0.015,minfrac=minFRAC)
    }
    else if(input$presets=="HP_Orbi"){
      xset<-group(Processing(),bw=5,mzwid=0.015,minfrac=minFRAC) 
    }
    else if(input$presets=="UP_QTOF"){
      xset<-group(Processing(),bw=2,mzwid=0.025,minfrac=minFRAC) 
    }
    else if(input$presets=="UP_QTOF.high"){
      xset<-group(Processing(),bw=2,mzwid=0.015,minfrac=minFRAC)   
    }
    else if(input$presets=="UP_Orbi"){
      xset<-group(Processing(),bw=2,mzwid=0.015,minfrac=minFRAC)  
    }
    else if(input$presets=="manual"){
      xset<-group(Processing(),bw=as.numeric(input$bandwidth),mzwid=as.numeric(input$mzwid),minfrac=minFRAC)   #Please update and improve
    }
    xset
  })
  
  output$xset_processing<-renderUI({
    
    #     if(input$choose_dir==TRUE){
    #       PROGRESS<-paste("<h4>Adjust the settings in the left panel</h4><h4>Then click on 'Start XCMS processing' at the bottom to initialize</h4>",
    #                       "<br><h4>Thereafter you can monitor the progress in the command prompt (cmd.exe) window</h4><hr>")
    #     }
    #     else NULL
    list(tags$body(
      tags$h4("1. Choose the directory containing your files"),
      tags$h4("2. Adjust the settings in the left panel"),
      tags$h4("3. Then click on 'Start pre-processing'"),tags$br()),
         checkboxInput("start",tags$h4("Start pre-processing"),FALSE),
         tags$br(),
         tags$body(tags$h5("Important: Please do not change any settings after starting pre-processing"),
                   tags$h5("*Tip:You can monitor the progress in the command prompt (cmd.exe) window")))
    
  })
  
  output$xset_progress<-renderText({
    xset<-Processing()
    
    if(is.null(xset)) NULL
    else
    if(!is.null(xset)){
      sink("XCMStemp/progress.txt")
      print(xset)
      sink()
      READ<-readLines("XCMStemp/progress.txt")
    }
    else NULL
#     sink("XCMStemp/progress.txt")
#     print(xset)
#     sink()
#     READ<-readLines("XCMStemp/progress.txt")
    if(!is.null(xset)){
      for (i in 1:length(READ)){
        READ[i]<-paste("<p>",READ[i],"</p>",sep="")
      }
      READ[length(READ)+1]<-paste("<br><h4>Please click on the next tab to continue...(Retention Time Correction)</h4>")
      READ1<-c(paste("<h4>Processing complete!</h4><br>"),READ)
      READ1 
    }
    else NULL

    })
#   output$start_xcms<-renderUI({
#     if(input$choose_dir==TRUE)
#       checkboxInput("start",tags$h4("Start XCMS processing"),FALSE)
#     else NULL
#   })
#   

  
  output$group<-renderText({
    grouping<-Group()
    sink("XCMStemp/group.txt")
    print(grouping)
    sink()
    READ<-readLines("XCMStemp/group.txt")
    if(!is.null(grouping)){
      for (i in 1:length(READ)){
        READ[i]<-paste("<p>",READ[i],"</p>",sep="")
      }
      READ[length(READ)+1]<-paste("<br><h4>Please click on the next tab to continue...(Retention Time Correction)</h4>")
      READ1<-c(paste("<h4>Processing complete!</h4><br>"),READ)
      READ1
    }
  })
  #===========================================================
  Retcor<-reactive({
    if(input$presets=="HP_QTOF"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation'))
    }
    else if(input$presets=="HP_QTOF.high"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation'))
    }
    else if(input$presets=="HP_Orbi"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation'))
    }
    else if(input$presets=="UP_QTOF"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation'))
    }
    else if(input$presets=="UP_QTOF.high"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation'))
    }
    else if(input$presets=="UP_Orbi"){
      xset<-retcor(Group(),method='obiwarp',plottype=c('deviation')) 
    }
    else if(input$presets=="manual"){
      xset <- retcor(Group(), method=input$retcor_method, plottype = "deviation")
    }
   
  })
  
  output$retcor<-renderText({
    retcor<-Retcor()
    sink("XCMStemp/retcor.txt")
    print(retcor)
    sink()
    READ<-readLines("XCMStemp/retcor.txt")
    if(!is.null(retcor)){
      for (i in 1:length(READ)){
        READ[i]<-paste("<p>",READ[i],"</p>",sep="")
      }
      READ[length(READ)+1]<-paste("<br><h4>Please click on the next tab to retrieve peak table</h4>")
      READ1<-c(paste("<h4>Processing complete!</h4><br>"),READ)
      READ1
    }
    
  })
  
  output$retcor_plot<-renderPlot({
    retcor<-Retcor()
    plotrt(retcor)
  },width=1000,height=750)
  #===========================================================
  Regroup<-reactive({
    minFRAC<-as.numeric(input$minfraction)/100
    if(input$presets=="HP_QTOF"){
      set4<-group(Retcor(),bw=5,mzwid=0.025,minfrac=minFRAC)
    }
    else if(input$presets=="HP_QTOF.high"){
      set4<-group(Retcor(),bw=5,mzwid=0.015,minfrac=minFRAC)
    }
    else if(input$presets=="HP_Orbi"){
      set4<-group(Retcor(),bw=5,mzwid=0.015,minfrac=minFRAC) 
    }
    else if(input$presets=="UP_QTOF"){
      set4<-group(Retcor(),bw=2,mzwid=0.025,minfrac=minFRAC) 
    }
    else if(input$presets=="UP_QTOF.high"){
      set4<-group(Retcor(),bw=2,mzwid=0.015,minfrac=minFRAC)   
    }
    else if(input$presets=="UP_Orbi"){
      set4<-group(Retcor(),bw=2,mzwid=0.015,minfrac=minFRAC)  
    }
    else if(input$presets=="manual"){
      set4<-group(Retcor(),bw=as.numeric(input$bandwidth),mzwid=as.numeric(input$mzwid),minfrac=minFRAC) #changes made here 24th Dec
    }
    set4
  })
  
  output$Regroup<-renderText({
    regroup<-Regroup()
    sink("XCMStemp/regroup.txt")
    print(regroup)
    sink()
    READ<-readLines("XCMStemp/regroup.txt")
    if(!is.null(regroup)){
      for (i in 1:length(READ)){
        READ[i]<-paste("<p>",READ[i],"</p>",sep="")
      }
      READ[length(READ)+1]<-paste("<br><h4>Please click on the next tab to continue...(Fill peaks)</h4>")
      READ1<-c(paste("<h4>Processing complete!</h4><br>"),READ)
      READ1
    }
  })
  #===========================================================
  
  Fillpeaks<-reactive({
    set5<-fillPeaks(Regroup())
  })
  
  output$Fillpeaks<-renderText({
    fillpeaks<-Fillpeaks()
    sink("XCMStemp/fillpeaks.txt")
    print(fillpeaks)
    sink()
    READ<-readLines("XCMStemp/fillpeaks.txt")
    if(!is.null(fillpeaks)){
      for (i in 1:length(READ)){
        READ[i]<-paste("<p>",READ[i],"</p>",sep="")
      }
      READ[length(READ)+1]<-paste("<br><h4>Click on next tab to access/download peak table</h4>")
      READ1<-c(paste("<h4>Processing complete!</h4><br>"),READ)
      READ1
    }
  })
  #===========================================================
  Peaklist<-reactive({
    peaklist<-peakTable(Fillpeaks(),filebase="Result")
  })
  
  
  Peaklist1<-reactive({
    #peaklist<-peakTable(Fillpeaks(),filebase="Result")
    peaklist<-Peaklist()
    if(input$apply_method==FALSE){
      return(peaklist)
    }
    else if(input$apply_method==TRUE){
      if(input$norm_method=="None"){
        return(peaklist)
      }
      else if(input$norm_method=="ISTD"){
        rel_path<-"XCMStemp/phenoData.txt"
        classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
        classinfo<-t(classinfo)
        unique_classes<-unique(classinfo)
        num_classes<-length(unique_classes)
        #Find out range of peak table to use
        columns<-as.numeric(7+num_classes)
        columns1<-columns+1
        meta<-data.frame(peaklist[,c(1:columns)])  #note that this might not hold for multiple groups
        Peaktable1<-peaklist[,c(columns1:ncol(peaklist))]
        
#         cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
#         Normed3<-snow::parApply(cl,Peaktable1,MARGIN=1,"/",Peaktable1[as.numeric(input$ISTD_list),])
#         stopCluster(cl)
#         Normed3<-data.frame(matrix(unlist(Normed3),nrow=length(Normed3),byrow=T))
#         names(Normed3)<-names(Peaktable1)
        PeakMAT<-as.matrix(Peaktable1)
        cl<-makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
        PeakMAT<-snow::parApply(cl,PeakMAT,MARGIN=1,"/",PeakMAT[as.numeric(input$ISTD_list),])
        stopCluster(cl)
        PeakMAT<-data.frame(t(PeakMAT))
        colnames(PeakMAT)<-colnames(Peaktable1)
        rownames(PeakMAT)<-rownames(Peaktable1)
        PeakMAT<-cbind(meta,PeakMAT)
        rownames(PeakMAT)<-rep(1:nrow(PeakMAT))
        return(PeakMAT)
      }
      else if(input$norm_method=="quantile"){
        rel_path<-"XCMStemp/phenoData.txt"
        classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
        classinfo<-t(classinfo)
        unique_classes<-unique(classinfo)
        num_classes<-length(unique_classes)
        #Find out range of peak table to use
        columns<-as.numeric(7+num_classes)
        columns1<-columns+1
        meta<-data.frame(peaklist[,c(1:columns)])  #note that this might not hold for multiple groups
        Peaktable1<-peaklist[,c(columns1:ncol(peaklist))]
        Quant.Norm<-normalize.quantiles(as.matrix(Peaktable1))
        colnames(Quant.Norm)<-colnames(Peaktable1)
        Quant.Norm<-data.frame(Quant.Norm)
        Out<-cbind(meta,Quant.Norm)
        return(Out)
      }
      
    }
    
    
  })
  
  Peaklist2<-reactive({ #not in use
    peaklist<-Peaklist()
    if(input$apply_ISTD){
      
      rel_path<-"XCMStemp/phenoData.txt"
      classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)
      #Find out range of peak table to use
      columns<-as.numeric(7+num_classes)
      columns1<-columns+1
      meta<-data.frame(peaklist[,c(1:columns)])  #note that this might not hold for multiple groups
      Peaktable1<-peaklist[,c(columns1:ncol(peaklist))]
      #Need inspiration for this
#       Normed<-apply(Peaktable1,MARGIN=1,"/",Peaktable1[as.numeric(input$ISTD_list),])
#       Normed1<-do.call(rbind.data.frame, Normed)
#       Normed1<-cbind(meta,Normed1)
#       rownames(Normed1)<-rep(1:nrow(Normed1))
#       return(Normed1)
      
      cl<-makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
      Normed3<-snow::parApply(cl,Peaktable1,MARGIN=1,"/",Peaktable1[as.numeric(input$ISTD_list),])
      stopCluster(cl)
      Normed3<-data.frame(matrix(unlist(Normed3),nrow=length(Normed3),byrow=T))
      
      names(Normed3)<-names(Peaktable1)
      Normed3<-cbind(meta,Normed3)
      rownames(Normed3)<-rep(1:nrow(Normed3))
      return(Normed3)
    }
      else{
        peaklist
    }
  })
  
  Peaklist3<-reactive({
    peaklist<-Peaklist()
    if(input$apply_quantile){
      rel_path<-"XCMStemp/phenoData.txt"
      classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)
      #Find out range of peak table to use
      columns<-as.numeric(7+num_classes)
      columns1<-columns+1
      meta<-data.frame(peaklist[,c(1:columns)])  #note that this might not hold for multiple groups
      Peaktable1<-peaklist[,c(columns1:ncol(peaklist))]
      Quant.Norm<-normalize.quantiles(as.matrix(Peaktable1))
      colnames(Quant.Norm)<-colnames(Peaktable1)
      Quant.Norm<-data.frame(Quant.Norm)
      Out<-cbind(meta,Quant.Norm)
      return(Out)
    }
    else peaklist
  })
  output$Peaklist<-renderTable({
#     if(input$norm_method=="None" |input$norm_method=="ISTD"|input$norm_method=="quantile"){
#       Peaklist()
#     }
#     else if(input$norm_method=="ISTD"){
#       Peaklist2()
#     }
#     else if(input$norm_method=="quantile"){
#       Peaklist3()
#     }
    Peaklist1()

  })
  
  output$ISTD_details<-renderText({
    X<-Peaklist()
    mz<-X[input$ISTD_list,1]
    mz<-round(mz,digits=3)
    rt<-X[input$ISTD_list,4]
    out<-paste("m/z=",mz,"& RT =",rt)
  })
  
  output$quantile_out<-renderText({
    if(input$apply_method==TRUE)
    X<-"Quantile normalization performed"
    else
      X<-"Please check 'Apply method selected' to proceed"
  })
  
  output$norm_out<-renderUI({
    if(input$norm_method=="None"){
      NULL
    }
    else if(input$norm_method=="ISTD"){
      list(
        selectInput("ISTD_list","Choose peak number containing internal standard",rownames(Peaklist())),
        verbatimTextOutput("ISTD_details")
        #checkboxInput("apply_ISTD","Apply selected internal standard?",FALSE)
        )
    }
    else if(input$norm_method=="quantile"){
      list(
        verbatimTextOutput("quantile_out")
        #checkboxInput("apply_quantile","Apply normalization?",FALSE)
      )
    }
  })
  #===========================================================
  output$downloadPeaklist<-downloadHandler(
    filename = function() { 
      type<-switch(input$norm_method, 
                   None="No.Norm",
                   ISTD="Internal.Std",
                   quantile="Quantile.Norm")
      paste("Peaklist",type,Sys.time(), '.tsv', sep='') },
    content = function(file) {
      if(input$norm_method=="None" & input$apply_method==TRUE){
        write.table(Peaklist1(), file, sep="\t",col.names=NA)
      }
      else if(input$norm_method=="ISTD" & input$apply_method==TRUE){
        write.table(Peaklist1(), file, sep="\t",col.names=NA)
      }else if(input$norm_method=="quantile" & input$apply_method==TRUE){
        write.table(Peaklist1(), file, sep="\t",col.names=NA)
      }
      
    }
    )
  output$temp<-renderPrint({
    path<-"XCMStemp"
    DIR<-dir(getwd())
    TEMPDIR<-grep("XCMStemp",DIR)
    if(!DIR[TEMPDIR]%in%"XCMStemp")
    dir.create(path)
#     tmppath<-paste(tempfile("XCMS_Peaklist"),'.tsv',sep='')
    tmppath<-paste(getwd(),"/XCMStemp/temp.tsv",sep="")
    if(input$norm_method=="None" & input$apply_method==TRUE){
      write.table(Peaklist1(),"XCMStemp/temp.tsv",sep="\t")
    }
    else if(input$norm_method=="ISTD" & input$apply_method==TRUE){ #used to be apply_ISTD==TRUE
      write.table(Peaklist1(),"XCMStemp/temp.tsv",sep="\t")
    }
    else if(input$norm_method=="quantile" & input$apply_method==TRUE){
      write.table(Peaklist1(),"XCMStemp/temp.tsv",sep="\t")
    }
#     else{
#       NULL
#     }
    
#     on.exit(unlink("XCMStemp/temp.tsv"))
#     if(input$norm_method=="None" & input$apply_method==TRUE){
#       cat("No normalization method applied \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
#     }
#     else if(input$norm_method!="None" & input$apply_method==TRUE){
#       if(input$norm_method=="ISTD" & input$apply_method==FALSE){
#         cat("Select peak belonging to internal standard and check the apply box")
#       }
#       else if(input$norm_method=="ISTD" & input$apply_method==TRUE){
#         cat("Internal Standard normalization applied. \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
#       }
#       else if(input$norm_method=="quantile" & input$apply_method==FALSE){
#         cat("Check the apply box to complete quantile normalization")
#       }
#       else if(input$norm_method=="quantile" & input$apply_method==TRUE){
#         cat("Quantile normalization applied. \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
#       }
#       else{
#         cat("Please apply your normalization method by checking the box")
#       }
#     }
# #     cat("Temporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
    if(input$apply_method==FALSE){
        if(input$norm_method=="None"){
          cat("Check the apply box to complete the process")
        }
        else if(input$norm_method=="ISTD"){
          cat("Select peak belonging to internal standard and check the apply box")
        }
        else if(input$norm_method=="quantile"){
          cat("Check the apply box to complete quantile normalization")
        }
      }
    else if(input$apply_method==TRUE){
      if(input$norm_method=="None"){
        cat("No normalization method applied \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
      }
      else if(input$norm_method=="ISTD"){
        cat("Internal Standard normalization applied. \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
      }
      else if(input$norm_method=="quantile"){
        cat("Quantile normalization applied. \nTemporary file written to",tmppath,"\nYou may now launch MetaboNexus to retrieve and work on this peaklist")
      }
    }
    
  })

})
