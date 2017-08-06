WD <- '/Users/zhuo/Desktop/'
system(paste0("mkdir -p ", "'",WD, "'"))
setwd(WD)

if(F){
	install.packages("formatR")	
}

packageName <- "MetaSparseKmeans"
devtools::create(packageName) 

## licenses
setwd(WD)
setwd(packageName)

## put some code in the R package

## make the code neat
formatR::tidy_dir("R")

## prepare .Rd refer to the slides

## create .Rd R documentation
devtools::document() 

## check the package
devtools::check()

## build the package
devtools::build()


## install the package
remove.packages("MetaSparseKmeans")
devtools::install()

install.packages("../MetaSparseKmeans_0.0.3.tar.gz",repos=NULL,type="source")


setwd(WD)
devtools::use_vignette(name = "MetaSparseKmeans", pkg = "MetaSparseKmeans/")

browseVignettes("MetaSparseKmeans")

## use the package
library(MetaSparseKmeans)
S = list(t(S1),t(S2), t(S2))
?MetaSparseKmeans

if(F){
  browseVignettes()
}

help(package="MISKmeans")

if(F){
  nstart = 20
  ntrial = 1
  maxiter = 20
  lambda = 1/2
  sampleSizeAdjust = FALSE
  wsPre = NULL
  silence = FALSE
  K<- 3
  wbounds <- 10
  awbound <- 10

  S = list(t(S1),t(S2), t(S2))
  res = MetaSparseKmeans(x=S,K=3,wbounds=10,lambda=0.5)
}

