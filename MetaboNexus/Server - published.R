require(shiny)
require(gplots)
require(rgl)
require(randomForest)
require(plsdepot)
require(XML)
require(ROCR)
require(snow)

source("House_Blend.R")

Welcome<-"Images/Welcome.gif"
MetaboNexusPNG<-"Images/MetaboNexus.png"
Reminder<-"Images/reminder.png"

#The following lines are for HMDB purposes
HMDB<-read.csv("Database/Full HMDB with status origin1.csv",stringsAsFactors=FALSE)
HMDB<-data.frame(HMDB,stringsAsFactors=FALSE)
HMDB_xml_path<-paste("hmdb_spectra_xml/",sep="")
LC_spectra<-read.delim("Database/Full_LC_spectra.txt",stringsAsFactors=F)
HMDB_pwList<-read.delim("Database/allPWs.txt",stringsAsFactors=F)
HMDB_pwList<-data.frame(HMDB_pwList)
colnames(HMDB_pwList)<-c("ID","Names","SMPDB","KEGG")
#for MSMS spectra plotting
HMDB_specList<-list.files("~/hmdb_spectra_xml")
HMDB_specList<-data.frame(HMDB_specList)
colnames(HMDB_specList)<-"specList"
#=========================================================================

shinyServer(function(input, output) {
  
  output$welcome<-renderImage({
    if(is.null(Data())){
      
      list(src = Welcome, alt = NULL)
    }
    else
    {
       list(src= Reminder, alt="Please annotate your sample classes using the next tab")
      
    }
  }, deleteFile = FALSE)
  
  output$MetaboNexusPNG<-renderImage({
      
      list(src = MetaboNexusPNG, alt = NULL)
   
  },deleteFile = FALSE)
  

  output$templocation<-renderText({
    if(input$tempfile==TRUE){
      rel_path<-paste(DIR,"XCMS GUI/XCMStemp/",sep="")
      FILELIST<-list.files(rel_path)
      TEMP<-grep("temp",FILELIST)
      
      if(FILELIST[TEMP]%in%"temp.tsv"){
        print("XCMS data found! \nNow only displaying peaks --->")
      }
          
      else stop("Temp file not found")
    }
    else NULL
  })
  

  output$Peaklist<-renderTable({
     Data()[[1]]
  })
  
  Data <- reactive({
    
    if(input$tempfile==TRUE){
      rel_path<-paste(DIR,"XCMS GUI/XCMStemp/phenoData.txt",sep="")
      classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)

    }

    
    if(input$tempfile==TRUE){
      rel_path2<-paste(DIR,"XCMS GUI/XCMStemp/temp.tsv",sep="")
      df.raw<-read.table(rel_path2)
      df.raw<-data.frame(df.raw,stringsAsFactors=FALSE)
      columns<-as.numeric(7+num_classes)
      columns1<-columns+1
      meta<-data.frame(df.raw[,c(1:columns)])  
      df.raw<-df.raw[,c(columns1:ncol(df.raw))]
      df.raw<-t(df.raw)
      df.raw<-data.frame(df.raw)

      
    } 
    else {
      inFile <- input$file1
      if (is.null(inFile)) return(NULL)
      
      df.raw <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      # Play with df.raw 
      row.names(df.raw)<-df.raw[,1]
      df.raw<-df.raw[,-1]
      MZ_position<-grep("mz",colnames(df.raw))
      RT_position<-grep("RT",colnames(df.raw))
      if(length(MZ_position)%in%0){
        meta<-NULL
      }
      else if(length(MZ_position)>0){
        mz<-as.numeric(df.raw[,MZ_position])
        mz<-mz[!is.na(mz)]
        
        rt<-as.numeric(df.raw[,RT_position])
        rt<-rt[!is.na(rt)]
        meta<-cbind(mz,rt)
        df.raw<-df.raw[,-c(MZ_position,RT_position)] 
      }

      invisible(df.raw)
    }

    
    
    
    if(input$transpose){
      df.raw<-t(df.raw)
      df.raw<-data.frame(df.raw)
      if(input$Supervised){
        Y<-df.raw[,input$Y]
        Y<-Y[!is.na(Y)]

        df.raw<-df.raw[,!names(df.raw) %in% input$Y]
        myList<-list(x=df.raw,meta,y=Y)
        invisible(myList)
        ###Imputation############################################
        if(input$impute){
          myList$x<-impute(myList$x)
        }
        ###NEARZEROVAR###########################################
        if(input$zerovar){
          nzv<-nearZeroVar(myList$x)
          if(length(nzv)==0){
            myList$x<-myList$x
          }
          else{
            myList[[2]]<-myList[[2]][-nzv,]
            myList$x<-myList$x[,-nzv]
            
          }
          
        }
        #option for log-transformation

          myList$x<-switch(input$log_10,
                           log10=parLOG10(myList$x),
                           ln=parLN(myList$x),
                           log2=parLOG2(myList$x),
                           None=myList$x)
                           

        invisible(myList)
       
      }
      else{#not supervised
        #imputation##########################################
        if(input$impute){
          df.raw<-impute(df.raw)
        }
        #NEARZEROVAR########################################
        if(input$zerovar){
          nzv<-nearZeroVar(df.raw)
          if(length(nzv)==0){
            df.raw<-df.raw
          }
          else{
            meta<-meta[-nzv,] #new addition
            df.raw<-df.raw[,-nzv]
          }
          
        }
        
        
        #option for log-transformation

        df.raw<-switch(input$log_10,
                         log10=parLOG10(df.raw),
                         ln=parLN(df.raw),
                         log2=parLOG2(df.raw),
                         None=df.raw)
        invisible(df.raw)
        
        metaList<-list(df.raw,meta) #NOTE
      }
    }
    #not transposed
    else{
      df.raw
      if(input$Supervised){
        Y<-df.raw[,input$Y]

        df.raw<-df.raw[,!names(df.raw) %in% input$Y]
        myList<-list(x=df.raw,meta,y=Y) #NOTE
        invisible(myList)
        #IMPUTATION#######################################
        if(input$impute){
          myList$x<-impute(myList$x)
        }
        if(input$zerovar){
          nzv<-nearZeroVar(myList$x)
          if(length(nzv)==0){
            myList$x<-myList$x
          }
          else{
            myList[[2]]<-myList[[2]][-nzv,]
            myList$x<-myList$x[,-nzv]
          }
        }        
        #option for log-transformation

        myList$x<-switch(input$log_10,
                         log10=parLOG10(myList$x),
                         ln=parLN(myList$x),
                         log2=parLOG2(myList$x),
                         None=myList$x)
        invisible(myList)
        
      }
      else{#not supervised
        #Imputation
        if(input$impute){
          df.raw<-impute(df.raw)
        }
        if(input$zerovar){
          nzv<-nearZeroVar(df.raw)
          if(length(nzv)==0){
            df.raw<-df.raw
          }
          else{
            meta<-meta[-nzv,] 
            df.raw<-df.raw[,-nzv]
          }
        }
        #option for log-transformation

        df.raw<-switch(input$log_10,
                         log10=parLOG10(df.raw),
                         ln=parLN(df.raw),
                         log2=parLOG2(df.raw),
                         None=df.raw)

        metaList<-list(df.raw,meta)
      }
    }
	print(str(df.raw))
  }) #end of Data()
  #Important note:
  #Data()[[1]]= df.raw; Data()[[2]]=meta; Data()[[3]]=annotations(from XCMS/Y)
  output$annotate<-renderUI({
    if(input$xcms_anno=="xcms_yes"|input$xcms_anno=="y_variable"){
      NULL
    }
    else{
      ROWNAMES<-rownames(Data()[[1]])
      OUTPUT<-list()
      INPUT<-NULL
      classes<-NULL
      for(i in 1:input$classes){
        tmp<-paste("Group",LETTERS[i])
        classes[[i]]<-tmp
      }
      #New stuff======================
      classes[[i+1]]<-"QC"
      classes[[i+2]]<-"Exclude"
      #===============================
      for(i in 1:length(ROWNAMES)){
        INPUT<-selectInput(ROWNAMES[i],ROWNAMES[i],classes)
        OUTPUT[[i]]<-list(INPUT)
        
      }
      OUTPUT
    }
    
  })
  
  output$manual_classes<-renderUI({
    if(input$xcms_anno=="manual_anno"){
      list(tags$body(tags$h4("Please select the correct groupings/classes below")),
      numericInput("classes",tags$h5("Select no. of classes"),2))
    }
    else if(input$xcms_anno=="y_variable"){
      tags$body(tags$h5("Please ensure you have checked the box for Y variable annotation"),
                tags$h5("See 'Y Variable Annotation' box on the left panel"))
           
    }
    else
      NULL
  })
  #===============================================================================
  #For Changing XMCS labels
  #===============================================================================
  
  labels<-reactive({
    if(input$xcms_anno=="xcms_yes"){
      classinfo<-read.table(paste(DIR,"XCMS GUI/XCMStemp/phenoData.txt",sep=""), header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)
      name_list<-list()
      for(i in 1:num_classes){
        name_list[[i]]<-rownames(classinfo)[which(classinfo==unique_classes[i])]
      }
    
      classinfo1<-classinfo
      for(i in 1:length(unique_classes)){
        classinfo1[which(classinfo==unique_classes[i])]<-i-1
      }
      classinfo1<-as.integer(classinfo1)
      if(input$label_override==TRUE){

        anno_list<-paste("input$class",1:num_classes,sep="")
        Y<-anno_list
        Y1<-NULL
        for(i in 1:length(Y)){
          tmptext<-Y[i]
          tmp<-eval(parse(text=tmptext))
          Y1<-c(Y1,tmp)
        }
        Y2<-NULL
        convert<-NULL
        #new stuff==================================
        for(i in 1:length(Y1)){
          convert<-which(LETTERS==substring(Y1[i],7))
          if(length(convert)==0){
            convert<-switch(Y1[i],"QC"=999,"Exclude"=1000)
          }
          Y2<-c(Y2,convert)
        }
        #===========================================
        Y2<-as.integer(Y2-1)
        Y2
        
      }
    }
    else NULL
  })
  
  output$class_confirmation<-renderUI({
    if(input$xcms_anno=="xcms_yes"){
      classinfo<-read.table(paste(DIR,"XCMS GUI/XCMStemp/phenoData.txt",sep=""), header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)
      name_list<-list()
      for(i in 1:num_classes){
        name_list[[i]]<-rownames(classinfo)[which(classinfo==unique_classes[i])]
      }
      
      classinfo1<-classinfo
      for(i in 1:length(unique_classes)){
        classinfo1[which(classinfo==unique_classes[i])]<-i-1
      }
      classinfo1<-as.integer(classinfo1)
      OUT<-list()
      
      class_num<-length(unique(classinfo))
      class<-NULL
      classes<-NULL
      #to generate alphabets
      for(i in 1:class_num){
        tmp<-paste("Group",LETTERS[i])
        classes[[i]]<-tmp
      }
      #new stuff=======================
      classes[[i+1]]<-"QC"
      classes[[i+2]]<-"Exclude"
      #================================
      for (i in 1:class_num){
        class<-paste("class",i,sep="")
        OUT[[i]]<-selectInput(class,paste("Class",unique(classinfo)[i], "will be labeled as"),classes)
      }
      OUT[[1+i]]<-checkboxInput("label_override",tags$h5("Apply class labelings"),FALSE)
      OUT
    }
    else NULL
  })
  #================================================================
  #For sample class annotation
  Annotation<-reactive({
    if(input$xcms_anno=="manual_anno"){
      anno_list<-NULL
      ROWNAMES<-rownames(Data()[[1]])
      for(i in 1:length(ROWNAMES)){
        tmp<-paste("input$",ROWNAMES[i],sep="")
        anno_list<-c(anno_list,tmp)
      }
      #     Y<-data.frame(anno_list)
      Y<-anno_list
      Y1<-NULL
      for(i in 1:length(Y)){
        tmptext<-Y[i]
        tmp<-eval(parse(text=tmptext))
        Y1<-c(Y1,tmp)
      }
      Y2<-NULL
      convert<-NULL
      
      for(i in 1:length(Y1)){
        convert<-which(LETTERS==substring(Y1[i],7))
        if(length(convert)==0){
          convert<-switch(Y1[i],"QC"=999,"Exclude"=1000)
        }
        
        Y2<-c(Y2,convert)
      }
      Y2<-as.integer(Y2-1)
      Y2
    }
    
    
    if(input$xcms_anno=="xcms_yes"){
      classinfo<-read.table(paste(DIR,"XCMS GUI/XCMStemp/phenoData.txt",sep=""), header=T, quote="\"",stringsAsFactors=F)
      classinfo<-t(classinfo)
      unique_classes<-unique(classinfo)
      num_classes<-length(unique_classes)
      name_list<-list()
      for(i in 1:num_classes){
        name_list[[i]]<-rownames(classinfo)[which(classinfo==unique_classes[i])]
      }

      classinfo1<-classinfo
      for(i in 1:length(unique_classes)){
        classinfo1[which(classinfo==unique_classes[i])]<-i-1
      }
      classinfo1<-as.integer(classinfo1)
      classinfo2<-NULL
      if(input$label_override==TRUE){
        positions<-list()
        for (i in 1:num_classes){
          j<-i-1
          positions[[i]]<-which(classinfo1==j)
        }
        for (i in 1:num_classes){
          classinfo1[positions[[i]]]<-labels()[i]
        }
      }
      else{
        classinfo1
      }
      classinfo1
    }
    else if(input$xcms_anno=="manual_anno"){
      Y2
    }
  })
  

  output$annotate2<-renderTable({
    Y<-data.frame(Annotation())
    Y
  })
  #==================================================
  output$meta<-renderTable({
   data.frame(Data()[[2]])
  })

  
  output$RAW <- renderTable({
    if (is.null(Data())){return()}

    else{
      if(input$full_table==FALSE){
        if(ncol(Data()[[1]])<500)
          columns<-ncol(Data()[[1]])
        else
          columns<-500
        data.frame(Data()[[1]][c(1:nrow(Data()[[1]])),c(1:columns)])
      }
      else{
        data.frame(Data()[[1]])
      }
    } 
  },digits=3)#end of output$RAW
  
  
  output$table_status<-renderUI({
    if(is.null(Data()[[1]])){return(NULL)}
    else{
      if(input$full_table==FALSE){
        if(ncol(Data()[[1]])<500)
          columns<-ncol(Data()[[1]])
        else
          columns<-500
        OUT<-paste("Now displaying the observations and the first",columns, "features/peaks only")
        helpText(OUT)
      }
      else{
        helpText("Now displaying the full table")
      }
    }
   
        
  })
  
  output$RAW_Y<-renderTable({
    if (is.null(input$file1)) { return() }
    
  })
  
  #====================================================================
  #For PCA
  #Uses {plsdepot}
  #====================================================================
  
  #Mother function of NIPALS PCA
  
  NIPALS<-reactive({
    if(length(Data())==3){
      y<-Data()[[3]]
      NIPALS<-nipals(Data()[[1]],comps=input$PCAcomp,scaled=input$PCA_scale)
    }
    else{
      y<-Annotation()
      QC<-which(y==998)
      Exclude<-which(y==999)
      if(length(Exclude)!=0){
        y<-y[-Exclude]
        NIPALS<-nipals(Data()[[1]][-Exclude,],comps=input$PCAcomp,scaled=input$PCA_scale)
      }
      else{
        NIPALS<-nipals(Data()[[1]],comps=input$PCAcomp,scaled=input$PCA_scale)
      }
      
      
    }
     
    NIPALS
            
  })
  #=====================================================================
  #Function to edit point of a plot
  #This will affect all score plots & proximity matrix plots
  output$editpoints<-renderUI({
    if(input$edit_points==TRUE){
      if(length(Data())==3){
        y<-Data()[[3]]
      }
      else{
        y<-Annotation()
      }
      
      unique_y<-sort(unique(y))
      label_y<-paste("label_y",unique_y,sep="")
      colorlist<-c("Black","Red","Green","Blue","Light Blue","Magenta","Yellow")
      outUI<-list()
      for (i in 1:length(unique_y)){
        if(unique_y[i]==998){
          outUI[[i]]<-selectInput(label_y[i],"Group labeled as QC",colorlist)
        }
        else{
          outUI[[i]]<-selectInput(label_y[i],paste("Group labeled as",unique_y[i],"or",LETTERS[unique_y[i]+1]),colorlist)
        }
        
      }
      
      outUI[[i+1]]<-helpText("*Changes made here will be reflected in other plots as well")
      outUI
      
      }
    else{
      NULL
    }
      
    
    
    
  })
  
  EditPoints<-reactive({
    if(input$edit_points==TRUE){
      if(length(Data())==3){
        y<-Data()[[3]]
      }
      else{
        y<-Annotation()
      }
      
      unique_y<-unique(y)
      newcolor<-NULL
      newinput<-NULL
      for (i in 1:length(unique_y)){
        newinput<-paste("input$label_y",unique_y[i],sep="")
        newcolor[i]<-eval(parse(text=newinput))
      }
      
      y_color<-as.numeric(y)
      selectedcolor<-NULL
      index<-list()
      for (i in 1:length(unique_y)){
        index[[i]]<-which(y_color==unique_y[i])
      }
      
      for (i in 1:length(unique_y)){
        selectedcolor<-newcolor[i]
        y_color[index[[i]]]<-switch(selectedcolor,   "Black"=1,
                                    "Red"=2,
                                    "Green"=3,
                                    "Blue"=4,
                                    "Light Blue"=5,
                                    "Magenta"=6,
                                    "Yellow"=7)
        
      }
      
      as.numeric(y_color)
        
      

      
    }
  })
  
  #Child functions from NIPALS PCA
  ### PCA plot
  
  output$NIPALS<-renderPlot({
    PCval<-round(NIPALS()$values[,2],digits=1)
    if(input$pca_scoreload=="score"){
      if(length(Data())==3){
        y<-Data()[[3]]
      }
      else{
        y<-Annotation()
        QC<-which(y==998)
        Exclude<-which(y==999)
        if(length(Exclude)!=0){
          y<-y[-Exclude]
          
        }
        else{
          y<-Annotation()
        }

        Unique<-sort(unique(y))
        max_y<-max(Unique[-length(Unique)])
        y[which(y==998)]<-max_y+1
      }
      
      if(input$PCA_labels==TRUE){
        if(input$edit_points==TRUE){
          plot(NIPALS(),what="observations",show.names=TRUE,cex=1.5,col.labels=as.numeric(EditPoints()),
               xlab=paste("PC1 (",PCval[1],"%)",sep=""),ylab=paste("PC2 (",PCval[2],"%)",sep=""),
               main="Principal Components Analysis")
        }
        else{ #no edit of points
          plot(NIPALS(),what="observations",show.names=TRUE,cex=1.5,col.labels=as.numeric(y)+1,
               xlab=paste("PC1 (",PCval[1],"%)",sep=""),ylab=paste("PC2 (",PCval[2],"%)",sep=""),
               main="Principal Components Analysis")
        }
      }
      else{#no labels, just points
        if(input$edit_points==TRUE){
          plot(NIPALS(),what="observations",cex=2,col.points=as.numeric(EditPoints()),pch=as.numeric(y)+15,
               xlab=paste("PC1 (",PCval[1],"%)",sep=""),ylab=paste("PC2 (",PCval[2],"%)",sep=""),
               main="Principal Components Analysis")
        }
        else{#no edit of points
          plot(NIPALS(),what="observations",cex=2,col.points=as.numeric(y)+1,pch=as.numeric(y)+15,
               xlab=paste("PC1 (",PCval[1],"%)",sep=""),ylab=paste("PC2 (",PCval[2],"%)",sep=""),
               main="Principal Components Analysis")
        }
      }
     
    }
    else if(input$pca_scoreload=="loadings"){
      plot(NIPALS(),what="variables",show.names=TRUE)
    }
    
  },width=1024,height=768)
  
  ##PCA plot in 3D
  output$NIPALS_3D<-renderPlot({
     if(input$PCA_3D_scores==TRUE){
       if(length(Data())==3){
         y<-Data()[[3]]
       }
       else{
         y<-Annotation()
         QC<-which(y==998)
         Exclude<-which(y==999)
         if(length(Exclude)!=0){
           y<-y[-Exclude]
         }
         else{
           y<-Annotation()
         }
         
         if(length(y)==0){
           y<-as.numeric(Annotation())
         }
       }
       EditPoints<-EditPoints()[-c(QC,Exclude)]
       if(input$edit_points==TRUE){
         plot3d_obs(NIPALS(),col.score=as.numeric(EditPoints)-1)
       }
       else{
         plot3d_obs(NIPALS(),col.score=y)
       }
       
      }
  },width=1024,height=768)
  
  ##PCA loadings  
  output$NIPALS_var<-renderPlot({
    if(input$pca_scoreload=="loadings"){
      plot(NIPALS(),what="variables",show.names=TRUE)
    }
     
  },width=1024,height=768)
  
  ##PCA loadings in 3D
  
  output$NIPALS_3D1<-renderPlot({
    if(input$PCA_3D_loads==TRUE){
        plot3d_var(NIPALS())
      }      
  },width=1024,height=768)
  
  ##Scree Plot for PCA
  output$RAWScree <- renderPlot({
    scree(NIPALS(),line=input$line,cum=input$cum)
  },width=1024,height=768)
  
  
  #===================================================================
  #For PLS
  #Uses {plsdepot}
  #===================================================================
  #Mother function for PLS analysis
  PLS<-reactive({
    x<-Data()[[1]]
    x<-data.frame(x)
    if(length(Data())==2){
      y<-as.numeric(Annotation())
      QC<-which(y==998)
      Exclude<-which(y==999)
      y<-y[-c(QC,Exclude)]
      if(length(y)==0){
        y<-as.numeric(Annotation())
        x<-data.frame(x)
      }
      else{
        x<-x[-c(QC,Exclude),]
      }
      
    }
    else{
      y<-Data()[[3]]
      x<-Data()[[1]]
    }
    PLS<-mod_plsreg1(x,y,comps=input$ncomp,crosval=TRUE)
  })
  #==================================================================
  ##Score plot for PLS
  output$pls<-renderPlot({
    if(input$pls_scoreload=="score"){
      if(length(Data())==3){
        y<-Data()[[3]]
      }
      else{
        y<-Annotation()
        QC<-which(y==998)
        Exclude<-which(y==999)
        y<-y[-c(QC,Exclude)]
        if(length(y)==0){
          y<-as.numeric(Annotation())
        }
        else{
          y
        }
      }
      if(length(Data())==2){
        EditPoints<-EditPoints()[-c(QC,Exclude)]
      }
      else{
        EditPoints<-EditPoints()
      }
      
      if(input$PLS_labels==TRUE){
        if(input$edit_points==TRUE){
          plot(PLS(),what="observations",col.points=as.numeric(EditPoints),col.xlabels=as.numeric(EditPoints),
               show.names=TRUE,cex=1.5,xlab="t1",ylab="t2",main="Partial Least Squares")
        }
        else{
          plot(PLS(),what="observations",col.points=as.numeric(y)+1,col.xlabels=as.numeric(y)+1,
               show.names=TRUE,cex=1.5,xlab="t1",ylab="t2",main="Partial Least Squares")
        }
      
      }
      else{#no labels, just points
        if(input$edit_points==TRUE){
          plot(PLS(),what="observations",col.points=as.numeric(EditPoints),col.xlabels=as.numeric(EditPoints),
               pch=as.numeric(y)+15,show.names=FALSE,cex=2,xlab="t1",ylab="t2",main="Partial Least Squares")
        }
        else{
          plot(PLS(),what="observations",col.points=as.numeric(y)+1,col.xlabels=as.numeric(y)+1,
               pch=as.numeric(y)+15,show.names=FALSE,cex=2,xlab="t1",ylab="t2",main="Partial Least Squares")
        }
        
      }
  
    }
    else if(input$pls_scoreload=="loadings"){
      plot(PLS(),what="variables",
           show.names=TRUE)
    }
    
  },width=1024,height=768)
  
  ##Score plot for PLS in 3D
  output$pls3D<-renderPlot({

    if(input$PLS_3D_scores==TRUE){
      if(length(Data())==3){
        y<-Data()[[3]]
      }
      else{
        y<-Annotation()
        QC<-which(y==998)
        Exclude<-which(y==999)
        y<-y[-c(QC,Exclude)]
        if(length(y)==0){
          y<-as.numeric(Annotation())
        }
        else{
          y
        }
      }
      if(length(Data())==2){
        EditPoints<-EditPoints()[-c(QC,Exclude)]
      }
      else{
        EditPoints<-EditPoints()
      }
      
      if(input$edit_points==TRUE){
        pls3d_obs(PLS(),col.score=as.numeric(EditPoints)-1) 
      }
      else{
        pls3d_obs(PLS(),col.score=y) 
      }
       
    }
             
  })
  
  ##Loadings for PLS
  output$pls_loadings<-renderPlot({
      plot(PLS(),what="variables",
         show.names=TRUE)
  },width=1024,height=768)
  
  ##Loadings for PLS in 3D
  output$pls_loadings3D<-renderPlot({
    if(input$PLS_3D_loads==TRUE){
       pls3d_var(PLS())
    }
  })
  

  
  #====================================================================
  #Q2 & R2 for PLS main tab
  #====================================================================
  output$Q2<-renderTable({
      
    PLS<-PLS()
    Q2<-PLS$Q2
    Q2<-as.data.frame(Q2)
  })
  
  output$R2<-renderTable({
   
    PLS<-PLS()
    R2<-PLS$R2
    R2<-as.data.frame(R2)
  })
  
  output$expvar<-renderTable({
    
    PLS<-PLS()
    expvar<-PLS$expvar
    expvar<-as.data.frame(expvar)
  })
 
  #====================================================================
  #Model diagnostics tab
  #====================================================================
  #From custom script
  output$mod_diag<-renderPlot({

    PLS<-PLS()
    if(input$show_all==FALSE){
      
      if(input$diag_type=="PRESS"){
        plot_R2Q2(PLS$Q2,type="PRESS")
        
      }
      else if(input$diag_type=="RSS"){
        plot_R2Q2(PLS$Q2,type="RSS")
      }
      else if(input$diag_type=="R2"){
        plot_R2Q2(PLS$R2,type="R2")
      }
      else if(input$diag_type=="R2X"){
        plot_R2Q2(PLS$expvar,type="R2X")##
      }
      else if(input$diag_type=="R2Y"){
        plot_R2Q2(PLS$expvar,type="R2Y")##
      }
      else if(input$diag_type=="R2Xcum"){
        plot_R2Q2(PLS$expvar,type="R2Xcum")##
      }
      else if(input$diag_type=="R2Ycum"){
        plot_R2Q2(PLS$expvar,type="R2Ycum")
      }
      else if(input$diag_type=="Q2"){
        plot_R2Q2(PLS$Q2,type="Q2")
      }
      else if(input$diag_type=="Q2cum"){
        plot_R2Q2(PLS$Q2,type="Q2cum")
      }
      else if(input$diag_type=="y.pred"){
        plot_R2Q2(PLS$y.pred,type="y.pred")
      }
      else if(input$diag_type=="resid"){
        plot_R2Q2(PLS$resid,type="resid")
      }
    }
    else{
      par(mfrow=c(4,3))
      plot_R2Q2(PLS$Q2,type="PRESS")
      plot_R2Q2(PLS$Q2,type="RSS")
      plot_R2Q2(PLS$R2,type="R2")
      plot_R2Q2(PLS$expvar,type="R2X")
      plot_R2Q2(PLS$expvar,type="R2Y")
      plot_R2Q2(PLS$expvar,type="R2Xcum")
      plot_R2Q2(PLS$expvar,type="R2Ycum")
      plot_R2Q2(PLS$Q2,type="Q2")
      plot_R2Q2(PLS$Q2,type="Q2cum")
      plot_R2Q2(PLS$y.pred,type="y.pred")
      plot_R2Q2(PLS$resid,type="resid")
    }
  },width=1024,height=768)
  #====================================================================
  #VIP value (PLS)
  #====================================================================
  
  #Mother function
  #Generates VIP values per component and supplies a p-value for each VIP value
  #Four types of statistical tests
  #To provide output for VIP download
  VIPlist<-reactive({
    PLS<-PLS()
    VIP<-PLS$VIP
    VIP
  })
  #============================================================================
  #Download VIP value list
  output$downloadData1 <- downloadHandler(
    
    filename = function() { paste("PLSDA_VIP", '.csv', sep='') },
    content = function(file) {
      write.csv(VIPlist(), file)
    }
  )
  
  #Plots VIP values from last component (cumulative)
  output$VIP_chart<-renderPlot({
    VIP<-VIPlist()
    VIP<-sort(VIP[,ncol(VIP)],decreasing=TRUE)
    VIP<-VIP[c(1:30)]
    dotchart(rev(VIP[c(1:30)]),main="Variable Importance in the Projection for PLS(-DA)",xlab="VIP value")
  },width=750,height=750)
  
  #Renders VIP values from all components with p-values
  output$pls_vip<-renderTable({
    VIPlist()
  },digits=3)
  
  #====================================================================
  #uiRender for setting VarImp cut-off
  #====================================================================
  output$varimp_input<-renderUI({
    if(input$varimp_cutoff=="None"){
      NULL
    }
    else if(input$varimp_cutoff=="PLS"){
      numericInput("PLS_cutoff","Enter VIP cut-off value",1)
    }
    else if(input$varimp_cutoff=="RF"){
      
      OUTPUT<-list()
      OUTPUT<-list(radioButtons("RF_varType",tags$h5("Select type of variable importance"),
                                c("Mean decrease in accuracy"="acc","Mean decrease in Gini"="gini")),
                   numericInput("RF_cutoff1","Enter RF cut-off value",0.15)
           
        )
    }
  })
  #============================================================================
  #FOR HEATMAP CONTROL
  #============================================================================
  output$varimp_input1<-renderUI({
    if(input$varimp_cutoff1=="None"){
      NULL
    }
    else if(input$varimp_cutoff1=="PLS"){
      numericInput("PLS_cutoff2","Enter VIP cut-off value",1)
    }
    else if(input$varimp_cutoff1=="RF"){
      
      OUTPUT<-list()
      OUTPUT<-list(radioButtons("RF_varType1",tags$h5("Select type of variable importance"),
                                c("Mean decrease in accuracy"="acc","Mean decrease in Gini"="gini")),
                   numericInput("RF_cutoff2","Enter RF cut-off value",0.15)
                   
      )
    }
  })
  #====================================================================
  #Univariate analysis
  #====================================================================
  Univariate<-reactive({
    x<-Data()[[1]]
    if(length(Data())==2){
      y<-as.numeric(Annotation())
      QC<-which(y==998)
      Exclude<-which(y==999)
      y<-y[-c(QC,Exclude)]

      if(length(y)==0){
        y<-as.numeric(Annotation())
      }
      else{
        x<-x[-c(QC,Exclude),]
      }

      
    }
    else{
      y<-Data()[[3]]
      x<-Data()[[1]]
    }
    df<-cbind(x,y)
    if(input$test_type=="none"){
      stop("Please select a hypothesis testing method. For >2 groups: Kruskal-Wallis/ANOVA | For 2 groups: Mann-Whitney/t-test")
    }
    
    else if(input$test_type=="KW"){
      #perform kruskal-wallis rank sum test
      KWp<-NULL #KW p-values

      cat(">>> Initializing Univariate test (Kruskal-Wallis)...\n")
      temp_time<-system.time({
      cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
      KWp<-parCapply(cl,x,function(x,y){
        KWt<-kruskal.test(x~y)
        KWp<-KWt[["p.value"]]
      },y=y)
      names(KWp)<-NULL
      stopCluster(cl)
      })
      
      
      Pval<-data.frame(colnames(x),KWp)
      colnames(Pval)[1]<-"Variable"
      colnames(Pval)[2]<-"KW p-value"
      cat(paste(">>> Univariate test (Kruskal-Wallis) completed in", round(temp_time[3],1),"seconds...\n"))
      Pval


    }
    else if(input$test_type=="anov"){
      #perform anova
      cat(">>> Initializing Univariate test (ANOVA)...\n")
      temp_time<-system.time({
      cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
      AOVp<-parCapply(cl,x,function(x,y){
        AOVt<-aov(x~y)
        AOVp<-summary(AOVt)[[1]][[5]][[1]]
      },y=y)
      stopCluster(cl)
      })
      names(AOVp)<-NULL
      
      
      Pval<-data.frame(colnames(x),AOVp)
      colnames(Pval)[1]<-"Variable"
      colnames(Pval)[2]<-"ANOVA p-value"
      cat(paste(">>> Univariate test (ANOVA) completed in", round(temp_time[3],1),"seconds...\n"))
      Pval

    }
    else if(input$test_type=="man"){
      #Mann-Whitney U tests/2-samples Wilcoxon tests
      if(length(unique(y))>2)
        #more than 2 treatment groups
        stop("Please use Kruskal-Wallis test or ANOVA for multiple treatments")
      else{
        cat(">>> Initializing Univariate test (Mann-Whitney U)...\n")
        temp_time<-system.time({
        cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
        MANp<-parCapply(cl,x,function(x,y){
          MANt<-wilcox.test(x~y)
          MANp<-MANt$p.value
        },y)
        stopCluster(cl)
        })
        names(MANp)<-NULL
        
        Pval<-data.frame(colnames(x),MANp)
        colnames(Pval)[1]<-"Variable"
        colnames(Pval)[2]<-"Mann-Whitney p-value"
        cat(paste(">>> Univariate test (Mann-Whitney U) completed in", round(temp_time[3],1),"seconds...\n"))
        Pval

      }
    }
    else if(input$test_type=="tt"){
      #t-tests
      if(length(unique(y))>2)
        #more than 2 treatment groups
        stop("Please use Kruskal-Wallis test or ANOVA for multiple treatments")
      else{
        cat(">>> Initializing Univariate test (Student's t-test)...\n")
        temp_time<-system.time({
        cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
        Tp<-parCapply(cl,x,function(x,y){
          Tt<-t.test(x~y)
          Tp<-Tt$p.value
        },y)
        stopCluster(cl)
        })
        names(Tp)<-NULL
        Pval<-data.frame(colnames(x),Tp)
        colnames(Pval)[1]<-"Variable"
        colnames(Pval)[2]<-"t-test p-value"
        cat(paste(">>> Univariate test (Student's t-test) completed in", round(temp_time[3],1),"seconds!<<<\n"))
        Pval

      }
    }
    
    if(input$p.adjust){
      tmp<-p.adjust(Pval[,2],method="fdr",n=length(Pval[,2]))
      Pval<-cbind(Pval,tmp)
      colnames(Pval)[3]<-"FDR-adjusted p-values"
      cat(">>> False discovery rate correction completed!<<<\n")
    }
    rm(x)
    cat(">>> Univariate test 100% completed! <<<\n")
    Pval
    

  })
  
  
  ###================================
  #ROC AUC plot
  output$ROC<-renderPlot({
    X<-Data()[[1]]
    if(length(Data())==2){
      Y<-as.numeric(Annotation())
      QC<-which(Y==998)
      Exclude<-which(Y==999)
      Y<-Y[-c(QC,Exclude)]
      if(length(Y)==0){
        Y<-as.numeric(Annotation())
        
      }
      else{
        X<-X[-c(QC,Exclude),]
      }
      
    }
    else{
      Y<-Data()[[3]]
      X<-Data()[[1]]
    }
    if(length(unique(Y))!=2){
      stop("Receiver Operating Characteristic (ROC) is unavailable for more than two sample classes, please try again with only two classes")
    }

    #ROC-AUC
    if(is.null(AUCval())) stop
    
    if(input$test_type=="none"){
      stop("Please select a statistical test before proceeding")
    }
    
        
    pred<-prediction(X[,input$predictor],Y)
    perf<-performance(pred, "tpr", "fpr")
    auc<-performance(pred, "auc")
    auc<-unlist(slot(auc, "y.values"))
    auc<-round(auc,digits=3)
    auc<-paste("AUC = ",auc,sep="")
    plot(perf, colorize = TRUE, main = paste("ROC-AUC for",input$predictor))
    legend(0,1,auc,border="white",cex=0.9,bty='n')
  },width=600)
  ###===============================================
  ###ROC1 function
  ROC1<-function(X,Y){
    pred<-ROCR::prediction(X,Y)
    auc<-ROCR::performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    auc<-round(auc,digits=4)
    return(auc)
  }
  
  
  #ROC AUC values
  AUCval<-reactive({
    X<-Data()[[1]]
    if(length(Data())==2){
      Y<-as.numeric(Annotation())
      QC<-which(Y==998)
      Exclude<-which(Y==999)
      Y<-Y[-c(QC,Exclude)]
      if(length(Y)==0){
        Y<-as.numeric(Annotation())
        
      }
      else{
        X<-X[-c(QC,Exclude),]
      }
      
    }
    else{
      Y<-Data()[[3]]
      X<-Data()[[1]]
    }
    if(length(unique(Y))!=2){
      stop("Receiver Operating Characteristic (ROC) is unavailable for more than two sample classes, please try again with only two classes")
    }
    cat(">>> Initializing ROC calculations...\n")
    temp_time<-system.time({
    cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
    allAUC<-parApply(cl,X,2,ROC1,Y=Y)
    stopCluster(cl)
    })
    cat(paste(">>> ROC calculations completed in",round(temp_time[3],1),"seconds! <<<\n"))
    return(allAUC)
  })
  
  ###=============================================== 
  
  output$uni_var<-renderTable({
    Univariate_table()
    
  },digits=4)
  
  Univariate_table<-reactive({
    if(is.null(Univariate())){
      stop("Please select a statistical test before proceeding")
    }
    else{
      Uni<-Univariate()
    }
    
    meta<-Data()[[2]]
    col_pval<-colnames(Uni)[-1]

    if(length(Data())==2){
      Y<-as.numeric(Annotation())
    }
    else{
      Y<-Data()[[3]]
    }
    
    if(length(unique(Y))==2){
      if(is.null(meta)){
        Uni<-data.frame(Uni[,1],AUCval(),Uni[,-1])
        colnames(Uni)[2:ncol(Uni)]<-"ROC AUC"
        colnames(Uni)[3:ncol(Uni)]<-col_pval
      }
      else{
        Uni<-data.frame(Uni[,1],meta,AUCval(),Uni[,-1])
        colnames(Uni)[(2+ncol(meta)):ncol(Uni)]<-"ROC AUC"
        colnames(Uni)[(3+ncol(meta)):ncol(Uni)]<-col_pval
      }
      
    }
    else{
      if(is.null(meta)){
        Uni<-data.frame(Uni[,1],Uni[,-1])
        colnames(Uni)[2:ncol(Uni)]<-col_pval
      }
      else{
        Uni<-data.frame(Uni[,1],meta,Uni[,-1])
        colnames(Uni)[(2+ncol(meta)):ncol(Uni)]<-col_pval
      }
    
    } 
    
    colnames(Uni)[1]<-"Peaks"
    
    
    if(input$uni_subset==FALSE){
      Uni2<-Uni
      rownames(Uni2)<-c(1:nrow(Uni2))
      Uni2
    }
      
      
    else if(input$uni_subset==TRUE & input$varimp_cutoff=="None"){
      Uni2<-Uni[which(Uni[,ncol(Uni)]<=as.numeric(input$uni_pval)),]
      
      if(length(Uni2[,1])==0){ 
        stop("No variables were found to match selection criteria. Consider reducing stringency of p-value")
      }
      else{
        if(length(unique(Y))==2){
          ROC<-Uni2[,"ROC AUC"]
          Uni2<-Uni2[,-which(colnames(Uni2)=="ROC AUC")]
          Uni2<-cbind(Uni2,ROC)
          colnames(Uni2)[ncol(Uni2)]<-"ROC AUC"
          rownames(Uni2)<-c(1:nrow(Uni2))
          Uni2
        }
        else{
          rownames(Uni2)<-c(1:nrow(Uni2))
          Uni2
        }
       
      }
      
    }
    
    else if(input$varimp_cutoff!="None" & input$uni_subset==TRUE){
      
      if(input$varimp_cutoff=="PLS"){
        VIP<-VIPlist()
        VIP<-VIP[,ncol(VIP)]
        Uni2<-cbind(Uni,VIP)
        Uni2<-Uni2[which(Uni2[,"VIP"]>=as.numeric(input$PLS_cutoff)),]
        Uni2<-Uni2[which(Uni2[,ncol(Uni2)-1]<=as.numeric(input$uni_pval)),]
        
        if(length(rownames(Uni2))%in%0){
          stop("No variables found to match PLS & p-value criteria, please consider reducing stringency")
        }
        else{
          if(length(unique(Y))==2){
            ROC<-Uni2[,"ROC AUC"]
            Uni2<-Uni2[,-which(colnames(Uni2)=="ROC AUC")]
            Uni2<-cbind(Uni2,ROC)
            colnames(Uni2)[ncol(Uni2)]<-"ROC AUC"
            rownames(Uni2)<-c(1:nrow(Uni2))
            Uni2
          }
          else{
            rownames(Uni2)<-c(1:nrow(Uni2))
            Uni2
          }
        
        }
        
      }
      else if(input$varimp_cutoff=="RF"){
        RF<-RF()
        Importance<-importance(RF)
        MeanGini<-ncol(Importance)
        MeanAcc<-MeanGini-1
        
        if(input$RF_varType=="acc"){
          Uni2<-cbind(Uni,Importance[,MeanAcc])
          colnames(Uni2)[ncol(Uni2)]<-"Decrease in \nMean Accuracy"
        }
        else if(input$RF_varType=="gini"){
          Uni2<-cbind(Uni,Importance[,MeanGini])
          colnames(Uni2)[ncol(Uni2)]<-"Decrease in \nMean Gini"
        }
        
        Uni3<-Uni2[which(Uni2[,ncol(Uni2)]>=as.numeric(input$RF_cutoff1)),]
        Uni3<-Uni3[which(Uni3[,ncol(Uni3)-1]<=as.numeric(input$uni_pval)),]
        
        if(length(rownames(Uni3))%in%0){
          stop("No variables found to match RF & p-value criteria, please consider reducing stringency")
        }
        else{
          if(length(unique(Y))==2){
            ROC<-Uni2[,"ROC AUC"]
            Uni2<-Uni2[,-which(colnames(Uni2)=="ROC AUC")]
            Uni2<-cbind(Uni2,ROC)
            colnames(Uni2)[ncol(Uni2)]<-"ROC AUC"
            rownames(Uni2)<-c(1:nrow(Uni2))
            Uni2
          }
          else{
            rownames(Uni2)<-c(1:nrow(Uni2))
            Uni2
          }
   
        }
        
      }
      
    }
  })
  
  output$downloadData4 <- downloadHandler(
    
    filename = function() { paste("Shortlisted peaks", '.csv', sep='') },
    content = function(file) {
      write.csv(Univariate_table(), file) #need edit!
    }
  )
  #====================================================================
  #Boxplots
  #====================================================================
  output$boxplot_control<-renderUI({
    if(is.null(Univariate())){
      stop("Please select a statistical test before proceeding")
    }
    else{
      Uni<-Univariate()
    }
    cat(">>> Filtering based on assigned p-value and/or variable importance...\n")        
    
      if(input$uni_subset==FALSE)
        selectInput("predictor","Please select feature/variable to plot",colnames(Data()[[1]]))
      
      else if(input$uni_subset==TRUE & input$varimp_cutoff=="None"){
        Uni2<-Uni[which(Uni[,ncol(Uni)]<=as.numeric(input$uni_pval)),]
              
        if(length(Uni2[,1])==0){ 
          stop("No variables were found to match selection criteria. Consider reducing stringency of p-value")
        }
        else{
                    
          ##=========ADDED BY WZ========================
          peakNames<-as.vector(Uni2[,1])
          #=============================================
          cat(">>> Filtering complete! <<<\n")
          selectInput("predictor","Please select shortlisted feature/variable to plot",peakNames)
        }
        
      }
        
    else if(input$varimp_cutoff!="None" & input$uni_subset==TRUE){
      
      if(input$varimp_cutoff=="PLS"){
        VIP<-VIPlist()
        VIP<-VIP[,ncol(VIP)]
        Uni2<-cbind(Uni,VIP)
        Uni2<-Uni2[which(Uni2[,"VIP"]>=as.numeric(input$PLS_cutoff)),]
        Uni2<-Uni2[which(Uni2[,ncol(Uni2)-1]<=as.numeric(input$uni_pval)),]
        
        if(length(rownames(Uni2))%in%0){
          stop("No variables found to match PLS & p-value criteria, please consider reducing stringency")
        }
        else{
          peakNames<-as.vector(Uni2[,1])
          cat(">>> Filtering complete! <<<\n")
          selectInput("predictor","Please select shortlisted PLS feature/variable to plot",peakNames)
        }
        
      }
      else if(input$varimp_cutoff=="RF"){
        RF<-RF()
        Importance<-importance(RF)
        MeanGini<-ncol(Importance)
        MeanAcc<-MeanGini-1
        
        if(input$RF_varType=="acc"){
          Uni2<-cbind(Uni,Importance[,MeanAcc])
        }
        else if(input$RF_varType=="gini"){
          Uni2<-cbind(Uni,Importance[,MeanGini])
        }
        
        Uni3<-Uni2[which(Uni2[,ncol(Uni2)]>=as.numeric(input$RF_cutoff1)),]
        Uni3<-Uni3[which(Uni3[,ncol(Uni3)-1]<=as.numeric(input$uni_pval)),]
        
        if(length(rownames(Uni3))%in%0){
          stop("No variables found to match RF & p-value criteria, please consider reducing stringency")
        }
        else{
          peakNames<-as.vector(Uni3[,1])
          cat(">>> Filtering complete! <<<\n")
          selectInput("predictor","Please select shortlisted RF feature/variable to plot",peakNames)
        }
  
      }
      
    }
    
    
  })
  
  
  
  #=========================================================================================
  output$boxplots<-renderPlot({
    if(input$test_type=="none"){
      stop("Please select a statistical test before proceeding")
    }
    
      X<-Data()[[1]]
      if(length(Data())==2){
        Y<-as.numeric(Annotation())
        QC<-which(Y==998)
        Exclude<-which(Y==999)
        Y<-Y[-c(QC,Exclude)]
        if(length(Y)==0){
          Y<-as.numeric(Annotation())
          
        }
        else{
          X<-X[-c(QC,Exclude),]
        }
        
      }
      else{
        if(length(Data()[[3]])%in%0)
          stop("Please ensure that Y variable/Annotation is entered/declared")
        else
        Y<-Data()[[3]]
        X<-Data()[[1]]
      }
      meta<-data.frame(Data()[[2]])
      if(length(meta)==0){
        meta<-NULL
        MZ_chosen<-NULL
        RT_chosen<-NULL
      }
      else{
        MZ<-grep("mz",colnames(meta))
        MZ<-meta[,MZ[1]]
        MZ_chosen<-which(colnames(X)==input$predictor)
        MZ_chosen<-round(MZ[MZ_chosen],digits=3)
        
        RT<-grep("rt",colnames(meta))
        RT<-meta[,RT[1]]
        RT_chosen<-which(colnames(X)==input$predictor)
        RT_chosen<-round(RT[RT_chosen],digits=3)
      }
      
 
  
      predictor<-input$predictor
      if(substring(predictor,0,1)=="X"){
        predictor<-substring(predictor,2)
      }
      else{
        predictor
      }
    
      ###Work on including mz & RT time over here
      testtype<-switch(input$test_type,KW="Kruskal-Wallis",anov="ANOVA",man="Mann-Whitney",tt="t-test")
      pvalue<-round(Univariate()[which(Univariate()[,1]==input$predictor),2],digits=6)
      if(input$xcms_anno=="xcms_yes"){
        rel_path<-paste(DIR,"XCMS GUI/XCMStemp/phenoData.txt",sep="")
        classinfo<-read.table(rel_path, header=T, quote="\"",stringsAsFactors=F)
        classinfo<-t(classinfo)
        unique_classes<-unique(classinfo)
        num_classes<-length(unique_classes)
        labels<-labels()
        NAMES<-sort(labels)
        for(i in 1:length(labels)){
          NAMES[which(NAMES==labels[i])]<-unique_classes[i]
        }
        NAMES
      }
      else if(input$xcms_anno=="manual_anno"){
        num_classes<-1:input$classes
        NAMES<-paste("Group",LETTERS[num_classes])
      }
      else if(input$xcms_anno!="manual_anno" & input$xcms_anno!="xcms_yes"){
        num_classes<-1:length(unique(Y))
        NAMES<-sort(unique(Y))
        NAMES<-paste("Group",NAMES)
      }
   
      
      if(input$p.adjust){
      FDR<-round(Univariate()[which(Univariate()[,1]==input$predictor),3],digits=6)
      }
      
      boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                              MZ_chosen,
                                                              "& RT=",
                                                              RT_chosen,
                                                              "\nTest:",testtype,"\np-value=",pvalue))
      if(input$p.adjust & input$varimp_cutoff=="None")
      boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,
                                                              "with m/z=",MZ_chosen,"& RT=",RT_chosen,
                                                              "\nTest:",testtype,"\np-value=",pvalue,"   ",
                                                              "FDR adj. p-value=",FDR))
   
      else if(input$varimp_cutoff=="PLS"){
        VIP<-VIPlist()
        VIP<-VIP[,ncol(VIP)]
        VIPval<-VIP[which(names(VIP)==input$predictor)]
        boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                                MZ_chosen,
                                                                "& RT=",RT_chosen,
                                                                "\nTest:",testtype,"   ","VIP value=",round(VIPval,digits=3),"\np-value=",pvalue))
        if(input$p.adjust)
          boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying feature",input$predictor,"with m/z=",
                                                                  MZ_chosen,"& RT=",RT_chosen,
                                                                  "\nTest:",testtype,"VIP value=",round(VIPval,digits=3),"\np-value=",pvalue,"   ",
                                                                  "FDR adj. p-value=",FDR))
      }
      
      else if(input$varimp_cutoff=="RF"){
        RF<-RF()
        Importance<-importance(RF)
 
        MeanAcc<-ncol(Importance)-1
        MeanGini<-ncol(Importance)
        Acc_val<-Importance[which(rownames(Importance)==input$predictor),MeanAcc]
        Gini_val<-Importance[which(rownames(Importance)==input$predictor),MeanGini]
        
        if(input$RF_varType=="acc"){
          boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                                  MZ_chosen,
                                                                  "& RT=",RT_chosen,
                                                                  "\nTest:",testtype,"   ","Mean Decrease in Accuracy=",round(Acc_val,digits=3),"\np-value=",pvalue))
          if(input$p.adjust)
            boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                                MZ_chosen,
                                                                "& RT=",RT_chosen,
                                                                "\nTest:",testtype,"   ","Mean Decrease in Accuracy=",round(Acc_val,digits=3),"\np-value=",pvalue,
                                                                "   ","FDR adj. p-value=",FDR))
        }
        else if(input$RF_varType=="gini"){
          boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                                MZ_chosen,
                                                                "& RT=",RT_chosen,
                                                                "\nTest:",testtype,"   ","Mean Decrease in Gini=",round(Gini_val,digits=3),"\np-value=",pvalue))
          if(input$p.adjust)
            boxplot(X[,input$predictor]~Y,names=NAMES,sub=paste("Displaying Feature",input$predictor,"with m/z=",
                                                                MZ_chosen,
                                                                "& RT=",RT_chosen,
                                                                "\nTest:",testtype,"   ","Mean Decrease in Gini=",round(Gini_val,digits=3),"\np-value=",pvalue,
                                                                "   ","FDR adj. p-value=",FDR))
        }
      }
    

      
    
  },width=600)
  
  #====================================================================
  #Overall Results
  #====================================================================
  
  Overall<-reactive({
    PLS<-VIPlist()
    PLS<-PLS[,ncol(PLS)]
    RF<-RF()
    Peaks<-colnames(Data()[[1]])
    if(length(Data())==2){
      Y<-as.numeric(Annotation())
      QC<-which(Y==998)
      Exclude<-which(Y==999)
      Y<-Y[-c(QC,Exclude)]
      if(length(Y)==0){
        Y<-as.numeric(Annotation())
      }
      else{
        Y
      }
      
    }
    else{
      Y<-Data()[[3]]
    }
    
      
    meta<-Data()[[2]]
    if(length(meta)%in%0){
      meta<-NULL
    }
    
    RFImportance<-importance(RF)
    MeanAcc<-ncol(RFImportance)-1
    MeanGini<-ncol(RFImportance)
    RFImportance<-RFImportance[,c(MeanAcc,MeanGini)]
    Univariate<-Univariate()
        
    
    if(!is.null(meta)){
      if(length(unique(Y))==2){
        Overall<-data.frame(Peaks,meta,PLS,RFImportance,Univariate[,-1])
        rownames(Overall)<-rep(1:nrow(Overall))
        if(ncol(Univariate)==2){
          colnames(Overall)[ncol(Overall)]<-c(colnames(Univariate())[ncol(Univariate())])
        }
        
        Overall<-cbind(Overall,AUCval())
        colnames(Overall)[ncol(Overall)]<-"ROC AUC"
        rownames(Overall)<-rep(1:nrow(Overall))
        

      }
      else{
        Overall<-data.frame(Peaks,meta,PLS,RFImportance,Univariate[,-1])
        rownames(Overall)<-NULL
        if(ncol(Univariate)==2){
          colnames(Overall)[ncol(Overall)]<-c(colnames(Univariate())[ncol(Univariate())])
        }

      }
      
      VIP_pos<-grep("PLS",colnames(Overall))
      colnames(Overall)[VIP_pos]<-"PLS VIP"
      rownames(Overall)<-rep(1:nrow(Overall))
      Overall
    }
          
    else{
      if(length(unique(Y))==2){
        Overall<-data.frame(Peaks,PLS,RFImportance,Univariate[,-1])
        rownames(Overall)<-NULL
        colnames(Overall)[c(1,2)]<-c("Peaks","PLS VIP")
        if(ncol(Univariate)==2){
          colnames(Overall)[ncol(Overall)]<-c(colnames(Univariate())[ncol(Univariate())])
        }
        Overall<-cbind(Overall,AUCval())
        colnames(Overall)[ncol(Overall)]<-"ROC AUC"
        
      }
      else{
        Overall<-data.frame(Peaks,PLS,RFImportance,Univariate[,-1])
        rownames(Overall)<-NULL
        colnames(Overall)[c(1,2)]<-c("Peaks","PLS VIP")
        if(ncol(Univariate)==2){
          colnames(Overall)[ncol(Overall)]<-c(colnames(Univariate())[ncol(Univariate())])
        }
        
      }
      rownames(Overall)<-rep(1:nrow(Overall))

      Overall
    }
   
  })
  
  output$overall<-renderTable({
    Overall()
  })
  
  output$downloadData3 <- downloadHandler(
    
    filename = function() { paste("Overall Results", '.csv', sep='') },
    content = function(file) {
      write.csv(Overall(), file)
    }
  )
  
  #====================================================================
  #Random Forest
  #====================================================================
  #Mother function
  ###Random Forest
  output$tree_var_in <- renderUI({
    if(input$tree_var_edit==TRUE){
      var<-sqrt(ncol(Data()[[1]]))
      var<-round(var,digits=0)
      numericInput("tree_var","No. of variables tried", var)
    }
    else{
      NULL
    }
  
    
  })
  
  output$tree_var1<-renderText({
    X<-ncol(Data()[[1]])
    X<-sqrt(X)
    print(X)
  })
  

  RF<-reactive({
    set.seed(40)
    X<-Data()[[1]]
    X<-as.data.frame(X)
    if(length(Data())==2){
      Y<-as.numeric(Annotation())
      QC<-which(Y==998)
      Exclude<-which(Y==999)
      Y<-Y[-c(QC,Exclude)]
      if(length(Y)==0){
        Y<-as.numeric(Annotation())
      }
      else{
        X<-X[-c(QC,Exclude),]
      }
    }
    else{
      Y<-Data()[[3]]
      X<-as.data.frame(Data()[[1]])
    }
    
    Y<-as.factor(Y)
    if(input$tree_var_edit==FALSE){
      RF<-randomForest(X,Y,ntree=as.numeric(input$trees),importance=TRUE,proximity=TRUE,na.action=na.roughfix)
    }
    else{
      RF<-randomForest(X,Y,ntree=as.numeric(input$trees),mtry=as.numeric(input$tree_var),importance=TRUE,proximity=TRUE,na.action=na.roughfix)
    }
    
          
  })
  
  output$RF<-renderPrint({
    
    print(RF())
    
  })
  
  output$RF_confusion<-renderTable({
    RF()$confusion
  })
  
  #==========================
  output$RF_stats<-reactive({
    OOB<-function(x){round(x$err.rate[x$ntree,"OOB"] * 100, digits = 2)}
    conf<-RF()$confusion
    ntree<-RF()$ntree
    mtry<-RF()$mtry
    nclass<-RF()$forest$nclass
    oobErr<-OOB(RF())
    
    htmlCode<-paste(
      "<table border=0 width=300>
        <tr>
          <th align=left colspan=",nclass+5,">Type: Classification</th>
        </tr>
        <tr>
          <th align=left colspan=",nclass+5,">No. of trees: ",ntree,"</th>
        </tr>
        <tr>
          <th align=left colspan=",nclass+5,">No. of variables tried: ",mtry,"</th>
        </tr>
        <tr>
          <th align=left colspan=",nclass+5,">OOB error: ",oobErr, "% </th>
        </tr>
        <tr>
          <th align=left colspan=",nclass+5,"><hr></th>
        </tr>
        <tr>
          <th align=left colspan=",nclass+5,">Confusion matrix</th>
        </tr>
        <tr>
          <td> </td>",sep="")
    
    conf<-data.frame(conf)
    #html is row-major
    for(i in 0:(nclass-1)){htmlCode<-paste(htmlCode,"<td align=center>",i,"</td>",sep="")}
    htmlCode<-paste(htmlCode,"<td align=center>Class error</td></tr>",sep="")
                    #<tr>",sep="")
    
    for(i in 0:(nclass-1)){#row number
      htmlCode<-paste(htmlCode,
            "<tr><td align=center>",i,"</td>",sep="")
      for(j in 0:(nclass)){#column number
        htmlCode<-paste(htmlCode,
            "<td align=center>",conf[i+1,j+1],"</td>",sep="")
      }
      htmlCode<-paste(htmlCode,"</tr>",sep="")
    }
    htmlCode<-paste(htmlCode,"</table>",sep="")
    
    
  })
  
  
  
  output$RF_VIP<-renderPlot({

    RF<-RF()
    rn<-round(importance(RF),2)
    rn[order(rn[,3],decreasing=TRUE)]
    VIP<-varImpPlot(RF,main="Variable Importance for Random Forest")
    
  }, width=750,height=750)
  
  output$MDS<-renderPlot({
    RF<-RF()
    if(length(Data())==2){
      kay<-Annotation()
      QC<-which(kay==998)
      Exclude<-which(kay==999)
      kay<-kay[-c(QC,Exclude)]
      if(length(kay)==0){
        kay<-as.numeric(Annotation())
      }
      else{
        kay
      }
      
    }
    else{
      kay<-Data()[[3]]
    }
    prox <- cmdscale(as.dist(1-RF$proximity), eig=TRUE, k=length(unique(kay)))
    x <- prox$points[,1]
    y <- prox$points[,2]
    if(length(Data())==2){
      coloring<-kay
    }
    else{
      coloring<-Data()[[3]]
    }
    if(length(Data())==2){
      EditPoints<-EditPoints()[-c(QC,Exclude)]
    }
    else{
      EditPoints<-EditPoints()
    }
    
    plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Proximity Matrix",sub="Multi-dimensional scaling", type="n")
    if(input$RF_labels==TRUE){
      if(input$edit_points==TRUE){
        text(x, y, labels = rownames(prox$points), col=as.numeric(EditPoints), cex=1.5)
      }
      else{
        text(x, y, labels = rownames(prox$points), col=as.numeric(coloring)+1, cex=1.5)
      }
      
    }
    else{
      if(input$edit_points==TRUE){
        points(x, y,col=as.numeric(EditPoints), pch=as.numeric(EditPoints)+15, cex=2)
      }
      else{
        points(x, y,col=as.numeric(coloring)+1, pch=as.numeric(coloring)+15, cex=2)
      }
     
    }
    
    
  },width=1024,height=768)
  
  output$RF_VIPval<-renderTable({
    RF<-RF()
    Importance<-importance(RF)
    
  })
  
  output$downloadData2 <- downloadHandler(
    
    filename = function() { paste("Random_Forest_VIP", '.csv', sep='') },
    content = function(file) {
      write.csv(importance(RF()), file)
    }
  )
  
  
  
  #====================================================================
  #HMDB search
  #====================================================================
  output$adduct_ctrl<-renderUI({
    if(input$ion_mode=="pos"){
      selectInput("adduct","Choose Adduct",
                  list("None (Non-ionized)"=0,
                       "M+H"=1.007825,
                       "M+Na"=22.989769,
                       "M+K"=38.963706,
                       "M+NH4"=14.0067+4*1.007825,
                       "M+H-2H2O"=1.007825-(2*18.0153),
                       "M+H-H2O"=1.007825-18.0153,
                       "M-H2O+NH4"=(14.003074+4*1.007825)-18.0153,
                       "M+CH3OH+H"=32.02621+1.007825,
                       "M+ACN+H"=41.02655+1.007825,
                       "M+2Na-H"=2*22.989769-1.007825,
                       "M+ACN+Na"=41.02655+22.989769,
                       "M+Li"=7.016004,
                       "M+2H"=2*1.007825,
                       "M+H+Na"=1.007825+22.989769,
                       "M+2Na"=2*22.989769,
                       "M+3H"=3*1.007825,
                       "M+2H+Na"=2*1.007825+22.989769,
                       "M+2Na+H"=2*22.989769+1.007825))
    }
    else{
      selectInput("adduct","Choose Adduct",
                  list("None (Non-ionized)"=0,
                       "M-H"=-1.007825,
                       "M-H2O-H"=-18.0153-1.007825,
                       "M+F"=18.998403,
                       "M+Na-2H"=22.989769-2*1.007825,
                       "M+Cl"=34.968852,
                       "M+K-2H"=38.963706-2*1.007825,
                       "M+FA-H"=46.005478-1.007825,
                       "M+CH3COO"=60.02113-1.007825,
                       "M-2H"=-2*1.007825,
                       "M-3H"=-3*1.007825))
    }
  })
  
  output$meta_info<-renderText({
    if(input$qType=="var_mz"){
      if(substring(input$var_mz2,0,1)=="X"){
        INPUT<-substring(input$var_mz2,2)
      }
      else{
        INPUT<-input$var_mz2
      }

      X<-Data()[[1]]            
      meta<-data.frame(Data()[[2]])
      if(length(meta)==0){
        meta<-NULL
        MZ_chosen<-NULL
        RT_chosen<-NULL
      }
      else{
        MZ<-grep("mz",colnames(meta))
        MZ<-meta[,MZ[1]]
        MZ_chosen<-which(colnames(X)==input$var_mz2)
        MZ_chosen<-round(MZ[MZ_chosen],digits=3)

        
        RT<-grep("rt",colnames(meta))
        RT<-meta[,RT[1]]
        RT_chosen<-which(colnames(X)==input$var_mz2)
        RT_chosen<-round(RT[RT_chosen],digits=3)

      }

      print(paste("m/z value =",MZ_chosen,"\nRetention Time=",RT_chosen))
    }
    else NULL
  })
  output$var_mz1<-renderUI({
    if(input$qType=="var_mz"){

      list(selectInput("var_mz2",tags$h5("Select variable to search for matches"),colnames(Data()[[1]])),
           verbatimTextOutput("meta_info"))
    }
    else{
      textInput("query","Enter query:",)
    }
  })
  
  output$metlin<-renderText({
    q<-input$query
    if(input$qType=="var_mz"){
      X<-Data()[[1]]            
      meta<-data.frame(Data()[[2]])
      if(length(meta)==0){
        meta<-NULL
        MZ_chosen<-NULL
        RT_chosen<-NULL
      }
      else{
        MZ<-grep("mz",colnames(meta))
        MZ<-meta[,MZ[1]]
        MZ_chosen<-which(colnames(X)==input$var_mz2)
        MZ_chosen<-round(MZ[MZ_chosen],digits=3)
 
        
        RT<-grep("rt",colnames(meta))
        RT<-meta[,RT[1]]
        RT_chosen<-which(colnames(X)==input$var_mz2)
        RT_chosen<-round(RT[RT_chosen],digits=3)
     
      }

      q<-MZ_chosen
      
    }
    if(q=="") return("No data was selected")
    #check if query is cpd name or mz value
    else if(input$qType=="name"){
      return("Metlin works only for a given m/z value. Please search by Variable Name or m/z value")
    }
    else if(input$qType=="mz"||input$qType=="var_mz"){
      
      #consider adduct
      ad<-as.numeric(input$adduct)
      q<-as.numeric(q)
      q<-q-ad
      
      
      #use tolerance to improve mz search
      tol<-input$tol
      mass_min<-q-tol
      mass_max<-q+tol
      metlink<-paste("http://metlin.scripps.edu/metabo_list.php?mass_min=",mass_min,"&mass_max=",mass_max,sep="")
      metlink1<-paste("<a href=",metlink," target=_blank> Click here to access Metlin results for chosen adduct & tolerance</a>",sep="")

    }
  })
  
  #purpose of function: output$hits renders a table of search results, based on query
  output$hits<-renderTable({
    q<-input$query
    if(input$qType=="var_mz"){
      X<-Data()[[1]]            
      meta<-data.frame(Data()[[2]])
      if(length(meta)==0){
        meta<-NULL
        MZ_chosen<-NULL
        RT_chosen<-NULL
      }
      else{
        MZ<-grep("mz",colnames(meta))
        MZ<-meta[,MZ[1]]
        MZ_chosen<-which(colnames(X)==input$var_mz2)
        MZ_chosen<-round(MZ[MZ_chosen],digits=3)
 
        
        RT<-grep("rt",colnames(meta))
        RT<-meta[,RT[1]]
        RT_chosen<-which(colnames(X)==input$var_mz2)
        RT_chosen<-round(RT[RT_chosen],digits=3)

      }
      q<-MZ_chosen

    }
    if(q=="") return()
    #check if query is cpd name or mz value
    else if(input$qType=="name"){
      #use agrep to search from local database using cpd names
      allHits<-agrep(q,HMDB$Name,max.distance=0.01)
      if(length(allHits)==0){
        tmp<-"No results found"
        tmp<-data.frame(tmp)
        colnames(tmp)<-"Error"
        tmp
      }
      else{
        df<-HMDB[allHits,1:ncol(HMDB)]
        data.frame(df)
        df
      }
    }
    else if(input$qType=="mz"||input$qType=="var_mz"){

      #consider adduct
      ad<-as.numeric(input$adduct)
      q<-as.numeric(q)
      q<-q-ad
      
      
      #use tolerance to improve mz search
      tol<-input$tol
      allHits<-which(abs(HMDB[,4]-q)<tol)
      df<-HMDB[allHits,1:ncol(HMDB)]

      
      if(length(allHits)==0||is.na(q)){
        tmp<-"No results found"
        tmp<-data.frame(tmp)
        colnames(tmp)<-"Error"
        tmp
      }
      else{
                df<-HMDB[allHits,1:ncol(HMDB)]
                df<-data.frame(df)
                df

      }
    }
    

  })#end of output$hits
  
  
  
  
  
  
  
  
  Pathways<-reactive({
    id<-input$hmdbid
    if(id==""){return()}
    else{
      #search for pathway names and ids
      pwNames<-HMDB_pwList$Names[which(id==HMDB_pwList$ID)]
      KID<-HMDB_pwList$KEGG[which(id==HMDB_pwList$ID)]
      SID<-HMDB_pwList$SMPDB[which(id==HMDB_pwList$ID)]
      
      #helper functions
      KPaste<-function(x){
        if(!is.na(x)) paste("http://www.genome.jp/kegg/pathway/map/",x,".html",sep="")
        else x<-NA}
      SPaste<-function(x){
        if(!is.na(x)) paste("http://pathman.smpdb.ca/pathways/",x,"/pathway?reset=true",sep="")
        else x<-NA}
      
      #formulating the URL
      KURL<-unlist(lapply(KID,KPaste))
      SURL<-unlist(lapply(SID,SPaste))
      
      KID[which(is.na(KID))]<-"Pathway not available on KEGG"
      SID[which(is.na(SID))]<-"Pathway not available on SMPDB"
      
      #store in a data frame
      tb<-cbind(pwNames,SID,SURL,KID,KURL)
      tb<-data.frame(tb,stringsAsFactors=FALSE)
      colnames(tb)<-c("Pathway","SMPDB","SMPDB_URL","KEGG","KEGG_URL")
      tb
    }
  })#end of Pathways()
  
  output$pwTable<-reactive({
    id<-input$hmdbid
    tb<-Pathways()
    tb<-data.frame(tb,stringsAsFactors=F)
    
    
    if(id==""||is.null(id)||is.na(id)||is.na(tb[1,1])){return()}
    else{
      htmlCode<-"<table border=1> 
      <tr>
      <th>Pathway</th>
      <th>SMPDB</th>
      <th>KEGG</th>
      </tr>" #open table tag in html
      
      #write html script in rows (html is row-major)
      for(i in 1:nrow(tb)){
        htmlCode<-paste(htmlCode,
                        "<tr>
                        <td>",tb[i,1],"</td>",sep="")
        if(!is.na(tb[i,3])){
          htmlCode<-paste(htmlCode,
                          "<td align=center><a href=",tb[i,3]," target=_blank>",tb[i,2],"</a></td>",sep="")
        }
        else{
          htmlCode<-paste(htmlCode,
                          "<td align=center>",tb[i,2],"</td>",sep="")
        }
        if(!is.na(tb[i,5])){
          htmlCode<-paste(htmlCode,
                          "<td align=center><a href=",tb[i,5]," target=_blank>",tb[i,4],"</a></td>",sep="")
        }
        else{
          htmlCode<-paste(htmlCode,
                          "<td align=center>",tb[i,4],"</td>",sep="")
        }
        htmlCode<-paste(htmlCode,"</tr>",sep="")
      }
      htmlCode<-paste(htmlCode,"</table>",sep="")
      htmlCode
  }
  })#end of pwTable()
  
  

  
  
  output$specPlot<-renderPlot({
    
    id<-input$hmdbid
    if(id==""){return()}
    else{

      cpdName<-NULL
      cpdName<-which(HMDB[,1]==id) #made to optimize speed
      cpdName<-HMDB[cpdName,2]
      
      fname<-paste(id,"_ms_ms_spectrum_",sep="")
      #use grep to search using file name
      specs<-grep(fname,LC_spectra[,1])
      if(length(specs)==0){
        stop("LC Spectra not found in both HMDB and MassBank")
      }
      specsdf<-data.frame(LC_spectra[specs,],stringsAsFactors=FALSE)
      nrow_specs<-nrow(specsdf)
      #use par to plot 3 spectra: low, med, high energy
    
      if(input$MassBank=="MassBank"){
        List<-grep("MassBank",specsdf[,4])
        specsdf<-specsdf[List,1]
        specsdf<-data.frame(specsdf)
        if(nrow(specsdf)>8){
          specsdf<-specsdf[c(1:8),1]
          specsdf<-data.frame(specsdf)
        }
        else specsdf
      }
      else{
        List<-grep("HMDB",specsdf[,4])
        specsdf<-specsdf[List,1]
        specsdf<-data.frame(specsdf)
       
      }

      if(input$MassBank=="MassBank"){
        if(ceiling(nrow(specsdf)/2)>1){
          par(mfrow=c(ceiling(nrow(specsdf)/2),2))
        }
        else if(ceiling(nrow(specsdf)/2)==1)
          par(mfrow=c(ceiling(nrow(specsdf)/2),1))
        else if(ceiling(nrow(specsdf)/2)==0)
          stop("No matching spectra found in MassBank")

       
      }
      else{
        par(mfrow=c(3,1))
      }
      
      for(i in 1:nrow(specsdf)){
        fpath<-paste(HMDB_xml_path,specsdf[i,1],sep="")
        if(fpath=="hmdb_spectra_xml/NA"){
          alternateDB<-switch(input$MassBank,HMDB="MassBank", MassBank="HMDB")
          stop(paste("Invalid ID/Spectra not found in ", input$MassBank, ", try the other reference database (",alternateDB,") instead",sep=""))
        }
          
        tree<-xmlTreeParse(fpath)
        root<-xmlRoot(tree)
        
        #only ms-ms spectra are plotted
        #check using energy levels
        energy<-root[["collision-energy-level"]]
        enLevel<-xmlValue(energy)
        if(length(enLevel)==0) {
          
          MassBank<-xmlValue(root[["references"]][[1]][[1]])
          Instrument<-xmlValue(root[["instrument-type"]])
          Instrument<-substring(Instrument,0,2)
          
          if(input$MassBank=="MassBank"){
            
            ionization<-xmlValue(root[["ionization-mode"]])
            voltage<-xmlValue(root[["collision-energy-voltage"]])
            peaks<-root[["ms-ms-peaks"]]
            mz<-NULL
            inten<-NULL
            for(i in 1:length(peaks)){
              #extract m/z values and respective intentsity
              inten<-c(inten,as.numeric(unlist(peaks[[i]][[2]][[1]])[2]))
              mz<-c(mz,as.numeric(unlist(peaks[[i]][[3]][[1]])[2]))
            }
            df<-data.frame(cbind(mz,inten))
            colnames(df)<-c("m/z","intensity")
            #plot
            
            plot(df, type="h",col="red",main=paste(cpdName, "MS/MS -", enLevel,"energy (",voltage,"V)", "\nin ", ionization,"mode"),
                 cex.lab=2.0, cex.axis=1.5, cex.main=2.2) #removed sub=fpath
            for(j in 1:nrow(df)){
              text(df[j,1],df[j,2],labels=df[j,1],cex=1.5)
            }
          }
          else{
            #warning is thrown in console, not browser
            simpleWarning("Spectra not found in local database.") #edited
          }
          
        }
        else{
          #extract peaks
          ionization<-xmlValue(root[["ionization-mode"]])
          voltage<-xmlValue(root[["collision-energy-voltage"]])
          peaks<-root[["ms-ms-peaks"]]
          mz<-NULL
          inten<-NULL
          for(i in 1:length(peaks)){
            #extract m/z values and respective intentsity
            inten<-c(inten,as.numeric(unlist(peaks[[i]][[2]][[1]])[2]))
            mz<-c(mz,as.numeric(unlist(peaks[[i]][[3]][[1]])[2]))
          }
          df<-data.frame(cbind(mz,inten))
          colnames(df)<-c("m/z","intensity")
          #plot
          
          plot(df, type="h",col="red",main=paste(cpdName, "MS/MS -", enLevel,"energy (",voltage,"V)", "in ", ionization,"mode"),
               cex.lab=2.0, cex.axis=1.5, cex.main=2.5)
          for(j in 1:nrow(df)){
            text(df[j,1],df[j,2],labels=df[j,1],cex=1.5)
          }
        }
      }
    } 
  },width=1000,height=1000)#end of output$specPlot
  
  #====================================================================
  #Heatmap tab
  #====================================================================
  ###Heatmap Panel
  Heatmap<-reactive({

    Data<-Data()[[1]]
    if(length(Data())==2){
      Y<-as.numeric(Annotation())
      QC<-which(Y==998)
      Exclude<-which(Y==999)
      Y<-Y[-c(QC,Exclude)]
      if(length(Y)==0){
        Y<-as.numeric(Annotation())
        
      }
      else{
        Data<-Data[-c(QC,Exclude),]
      }
      
    }
    else{
      Y<-Data()[[3]]
    }
    
    Uni<-Univariate()
    if(input$uni_subset1==FALSE){
      Data_out<-Data
    }
    else if(input$uni_subset1==TRUE & input$varimp_cutoff1=="None"){
      Shortlist<-which(Uni[,ncol(Uni)]<=as.numeric(input$uni_pval1))
      Data_out<-Data[,Shortlist]
    }
    else if(input$uni_subset1==TRUE & input$varimp_cutoff1=="PLS"){
      VIP<-VIPlist()
      VIP<-VIP[,ncol(VIP)]
      Uni2<-cbind(Uni,VIP)
      Shortlist<-which(Uni2[,ncol(Uni2)]>=as.numeric(input$PLS_cutoff2) & Uni2[,ncol(Uni2)-1]<=as.numeric(input$uni_pval1))
      Data_out<-Data[,Shortlist]
    }
    else if(input$uni_subset1==TRUE & input$varimp_cutoff1=="RF"){

      RF<-RF()
      RF<-importance(RF)
      RF<-RF[,c(ncol(RF)-1,ncol(RF))]
      Uni2<-cbind(Uni,RF)
      if(input$RF_varType1=="acc"){
        Shortlist<-which(Uni2[,ncol(Uni2)-1]>=as.numeric(input$RF_cutoff2) & Uni2[,ncol(Uni2)-2]<=as.numeric(input$uni_pval1))
      }
      else{
        Shortlist<-which(Uni2[,ncol(Uni2)]>=as.numeric(input$RF_cutoff2) & Uni2[,ncol(Uni2)-2]<=as.numeric(input$uni_pval1))
      }
      
      Data_out<-Data[,Shortlist]
      
    }
 
  })
  
  
  output$QHeatmap<-renderPlot({
    cat(">>> Creating Heatmap...\n")
    dendro<-input$dendro_ctrls
    select.col<-switch(input$heatmap.col,
                       greenred=greenred(75),
                       redgreen=redgreen(75),
                       bluewhitepink=cm.colors(75),
                       bluewhitered=bluered(75),
                       heatcol=heat.colors(75))
    ordering<-switch(input$ordering,rowmeans=eval(parse(text=c("ROW<-TRUE","COL<-FALSE"))),
                     colmeans=eval(parse(text=c("COL<-TRUE","ROW<-FALSE"))),
                     no_order=eval(parse(text=c("ROW<-FALSE","COL<-FALSE"))),
                     both=eval(parse(text=c("ROW<-TRUE","COL<-TRUE"))))
    if(input$heat_trans==TRUE){
      heatmap.2(as.matrix(t(Heatmap())),Rowv=ROW, Colv=COL, dendrogram=dendro, 
                col=select.col,key=TRUE,keysize=0.5,scale=input$scaling,trace="none")
    }
    else{
      heatmap.2(as.matrix(Heatmap()),Rowv=ROW, Colv=COL, dendrogram=dendro,
                col=select.col,key=TRUE,keysize=0.5,scale=input$scaling,trace="none")
    }
    
  },res=81)


  all.settings<-reactive({

    Control_list<-NULL
    if(is.null(input$file1[[1]])){
      filename<-paste("Name of file analyzed: Pre-processed file")
    }
    else{
      filename<-paste("Name of file analyzed:",input$file1[[1]])
    }
    
    session<-paste("Current time of session:",date())
    header<-paste("Select header as column name:",input$header)
    if(input$sep=='\t'){
      sep<-"Tab"
    }
    if(input$sep==","){
      sep<-"Comma (.csv)"
    }
    if(input$sep==";"){
      sep<-"Semicolon"
    }
    separator<-paste("Type of separator used:",sep)
    if(input$quote==''){
      quote<-"None"
    }
    if(input$quote=="'"){
      quote<-"Single Quote"
    }
    if(input$quote=='"'){
      quote<-"Double Quote"
    }
    quote<-paste("Quote used:",quote)
    transpose<-paste("Data was transposed:",input$transpose)
    supervised<-paste("Supervised analysis was performed:",input$Supervised)
#     Y_vec<-paste("Name of Y variable used:",input$Y)
    impute<-paste("Imputation for zero/NA values was performed:",input$impute)
    zerovar<-paste("Removal of variables with zero/near zero variance was performed:",input$zerovar)
    log_10<-paste("log10 transformation was performed:",input$log_10)
    PCAcomp<-paste("No. of PCA components used:",input$PCAcomp)
    PCA_scale<-paste("Scaling of data for PCA:",input$PCA_scale)
    PCA_3D_scores<-paste("3D view of PCA scores:",input$PCA_3D_scores)  #(might be changed)
    PCA_3D_loads<-paste("3D view of PCA loadings:",input$PCA_3D_loads) 
    line<-paste("80% line displayed for PCA Scree plot:",input$line) #(80% line)
    cum<-paste("Scree plot plotted as cumulative:",input$cum) #(cumulative format)
    ncomp<-paste("No. of components used for PLS(-DA):",input$ncomp) #(ncomp for PLS)
    PLS_3D_scores<-paste("3D view of PLS(-DA) scores:",input$PLS_3D_scores) #(might be changed)
    PLS_3D_loads<-paste("3D view of PLS(-DA) loads:",input$PLS_3D_loads)
    show_all<-paste("Show all diagnostic plots:",input$show_all)
    diag_type<-paste("Type of diagnostic plot last viewed:",input$diag_type) #(diagnostic plot type, multiple options)
    if(input$p.adjust==TRUE){
      padjust<-"with FDR correction"
    }
    else{
      padjust<-"with no FDR correction"
    }
    
    
    if(input$test_type=='KW'){
      test<-"Kruskal-Wallis"
    }
    if(input$test_type=="anov"){
      test<-"ANOVA (1-way)"
    }
    if(input$test_type=="man"){
      test<-"Mann-Whitney U"
    }
    if(input$test_type=="tt"){
      test<-"Student's t-test"
    }
    test_type<-paste("Type of univariate statistical test performed:",test, padjust) #(type of stats analysis)
    #p.adjust<-paste("Apply False Discovery Rate:",input$p.adjust)
    trees<-paste("No. of trees used for Random Forest:",input$trees)
    if(input$tree_var_edit==FALSE){
      tree_var<-paste("No. of variables used for each tree in Random Forest:",round(sqrt(ncol(Data()[[1]])),0))
    }
    else{
      tree_var<-paste("No. of variables used for each tree in Random Forest:",input$tree_var)
    }
    
    query<-paste("Value entered for MS search:",input$query)  #(MS search)
    qType<-paste("Type of query:",input$qType) #(query type : name or mz)
    adduct<-paste("Type of adduct:",input$adduct)
    tol<-paste("Tolerance used:",input$tol,"+/- mz")
    hmdbid<-paste("HMDB ID last entered:",input$hmdbid)
    MassBank<-paste("MS database reference viewed:",input$MassBank)
    Heatmap_cutofftype<-switch(input$varimp_cutoff1,"None"="None","PLS"="VIP value","RF"="Random Forest")
    
    if(Heatmap_cutofftype=="Random Forest"){
      RF_type<-switch(input$RF_varType1,"acc"="Mean Decrease in Accuracy","gini"="Mean Decrease in Gini")
      heatmap_subset<-paste("Heatmap was subsetted using:", RF_type, ">=", input$RF_cutoff2, "and p-value =<", input$uni_pval1 , padjust )
    }
    else if(Heatmap_cutofftype=="VIP value"){
      heatmap_subset<-paste("Heatmap was subsetted using: VIP value",  ">=", input$PLS_cutoff2, "and p-value =<", input$uni_pval1, padjust )
    }
    else if(Heatmap_cutofftype=="None"){
      heatmap_subset<-paste("Heatmap was subsetted using:", "p-value =<", input$uni_pval1 ,padjust)
    }
    
    heat_trans<-paste("Heatmap was tranposed:",input$heat_trans) #(tranpose)
    scaling<-paste("Type of scaling used for heatmap:",input$scaling) #(heatmap scaling)
    Control_list<-rbind(filename,session,header,separator,quote,transpose,supervised,impute,zerovar,log_10,PCAcomp,PCA_scale,
                             line,cum,ncomp,trees,tree_var,test_type,
                             heatmap_subset,heat_trans,scaling)
    Control_list<-data.frame(Control_list,stringsAsFactors=F)
    
    
    tmp<-c(rep("Session",2),rep("Data",3),rep("Data manipulation",5),rep("PCA",4),
           rep("PLS(-DA)",1),rep("Random Forest",2),rep("Univariate analysis",1),rep("Heatmap",3))
    Control_list<-cbind(tmp,Control_list)
    
    
    colnames(Control_list)<-c("Operation Type","Settings")
    
    rownames(Control_list)<-NULL
    Control_list
  })
  
  output$download_log <- downloadHandler(
    
    filename = function() { paste("Log Session ",date(), '.txt', sep='') },
    content = function(file) {
      write.table(all.settings(), file)
    }
  )
  
  output$all.settings<-renderTable({
    x<-all.settings()
  })
  
  output$sessiontime<-renderText({
    x<-date()
  })
  gc(TRUE)
  gcinfo(FALSE)

})

