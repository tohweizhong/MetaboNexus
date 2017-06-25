@echo off

:: A batch script for running R scripts in M$ XP
:: Just attach your R script to the end of this batch file (below the
:: last ":::" comment line) and run this script.

:: A portion of this file was taken from Gabor Grothendieck's
:: batchfiles http://cran.r-project.org/contrib/extra/batchfiles/
:: Copyright by Gabor Grothendieck 2005, GPL v2

setlocal
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo For first time installation, the following packages will be installed
echo ===================================================
echo 1."shiny" (For GUI)
echo 2."xcms"  (For raw data processing)
echo 3."snow"  (For parallel processing)
echo This program was built and tested 
echo Using R version 3.0.1 (2013-05-16) -- "Good Sport"
echo And shiny 0.7.0
echo From Windows XP to Windows 8
echo ===================================================
echo Best results when Firefox/Chrome is used as browser
echo ===================================================
echo Important Note: This batch file will only work 
echo if you have installed R with its registry key
echo Hence it is recommended that you do so prior to 
echo using this software
echo ===================================================
PAUSE

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Find R.
:: use environment variable R_HOME if defined
:: else current folder if bin\rcmd.exe exists
:: else most current R as determined by registry entry
:: else error, not found.
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set relativeup=%~dp0
set R_HOME=%relativeup%\R-3.0.1\
set cmdpath=%R_HOME%\bin\R.exe
set thisfile=%~f0
set relative=%relativeup%\XCMS GUI\

:: Run R.  Make it parse the commands at the end of this file
echo x=readLines(Sys.getenv('thisfile'));setwd(Sys.getenv('relative'));eval(parse(text=x[-(1:grep('Put R code below',x)[2])]))  | "%cmdpath%" --vanilla
PAUSE
endlocal
exit /b

:::::::::::::::::::  Put R code below this line :::::::::::::::::::::

if(!require("shiny")){
  install.packages("shiny",repos="http://cran.us.r-project.org")
}
if(!require("xcms")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("xcms")
}
if(!require("snow")){
  install.packages("snow",repos="http://cran.us.r-project.org")
}
if(!require("preprocessCore")){
  source("http://bioconductor.org/biocLite.R")
    biocLite("preprocessCore")
}


library(shiny)
runApp(port=8112L,launch.browser=TRUE)