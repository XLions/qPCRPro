
Get_File <- function(ZipAddress){
  library(stringr)
  names <- unzip(ZipAddress, list = TRUE)$Name
  names <- names[-which(str_detect(names, 'MACOSX'))]
  Output <- list()
  Output[[1]] <- read.csv(unzip(ZipAddress, 
                                names[which(str_detect(names,'Amplification'))]))[,-1]
  Output[[2]] <- read.csv(unzip(ZipAddress, 
                                names[which(str_detect(names,'Cq Results'))]))[,-1]
  Output[[3]] <- read.csv(unzip(ZipAddress, 
                                names[which(str_detect(names,'Derivative'))]))[,-1]
  Output[[4]] <- read.csv(unzip(ZipAddress, 
                                names[which(str_detect(names,'RFU'))]))[,-1]
  names(Output) <- c('AmplificationResults',
                     'CqResults',
                     'MeltCurve_DerivativeResults',
                     'MeltCurve_RFUResults')
  return(Output)
}
