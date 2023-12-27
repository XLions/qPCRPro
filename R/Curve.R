#Covert Data to Paint Format
ConvertPaintFormatCurveData <- function(RawCurveData){
  Output <- data.frame(matrix(nrow = (nrow(RawCurveData) * ncol(RawCurveData)),
                              ncol = 3))
  colnames(Output) <- c('Rownames', 'Value', 'Colnames')
  Output$Rownames <- rep(rownames(RawCurveData), ncol(RawCurveData))
  Output$Colnames <- rep(colnames(RawCurveData),
                         each = nrow(RawCurveData))
  for (i in 1 : ncol(Format_MeltCurve_RFU)) {
    Output$Value[((i-1) * nrow(RawCurveData)+1) : ((i) * nrow(RawCurveData))] <-
      as.numeric(RawCurveData[,i])
  }
  return(Output)
}

#Some Process Funtcions
##Transfer Well Names (with or without 0, like 'C01' or 'C1')
Col0X_to_ColX<-function(Col0X){
  library(stringr)
  first_letter<-str_sub(Col0X,1,1)
  second_letter<-str_sub(Col0X,2,2)
  third_letter<-str_sub(Col0X,3,3)
  second_letter<-str_replace_all(second_letter,'0','')
  return(paste0(first_letter,second_letter,third_letter))
}
##Get Well Locations and Co-responsive Target and Sample from Cq File
Get_Locs_Targets_Samples<-function(CqResults){
  library(stringr)
  Output<-CqResults[,which(colnames(CqResults) %in%
                             c('Well','Target','Sample'))]
  Output$Well<-Col0X_to_ColX(Output$Well)
  return(Output)
}

#Select Sample and Target
SelectTargetSample <- function(targets,
                               treatment,
                               Format_MeltCurve_RFU,
                               Locs_Targets_Samples){
  selected_wells <- Locs_Targets_Samples$Well[which(
    Locs_Targets_Samples$Target %in% as.character(targets) &
      Locs_Targets_Samples$Sample %in% as.character(treatment)
  )]
  Output <- Format_MeltCurve_RFU[,which(
    colnames(Format_MeltCurve_RFU) %in% selected_wells)]
  return(Output)
}

#MeltCurve_RFU
Plot_MeltCurve_RFU <- function(Format_MeltCurve_RFU,
                               Locs_Targets_Samples,
                               Group){
  library(ggplot2)
  #Define Curve Plot Style Theme
  CurveTheme <- {ggplot2::theme_bw()+
      theme(axis.title = element_text(face = 'bold'),
            plot.title = element_text(face = 'bold',
                                      hjust = 0.5),
            legend.title = element_text(face = 'bold',
                                        hjust = 0.5))}
  PaintFormatCurveData <- ConvertPaintFormatCurveData(Format_MeltCurve_RFU)
  PaintFormatCurveData <- merge(x = PaintFormatCurveData,
                                y = Locs_Targets_Samples,
                                by.x = 'Colnames',
                                by.y = 'Well')
  PaintFormatCurveData$Rownames <- as.numeric(PaintFormatCurveData$Rownames)
  Plot_Basic <- ggplot(data=PaintFormatCurveData,
                       mapping=aes(x = Rownames,
                                   y = Value,
                                   group = Colnames))
  if (Group == 'none'){
    Plot_Advanced <- Plot_Basic+
      geom_line(color = '#008B00')+
      labs(x = 'Temperature',
           y = 'RFU',
           title = 'Melt Curve')+
      CurveTheme
  } else if (Group == 'target'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Target))+
      labs(x = 'Temperature',
           y = 'RFU',
           title = 'Melt Curve')+
      CurveTheme
  } else if (Group == 'treatment'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Sample))+
      labs(x = 'Temperature',
           y = 'RFU',
           title = 'Melt Curve')+
      CurveTheme
  }
  return(Plot_Advanced)
}
#MeltCurve_Derivative
Plot_MeltCurve_Derivative <- function(Format_MeltCurve_Derivative,
                               Locs_Targets_Samples,
                               Group){
  library(ggplot2)
  #Define Curve Plot Style Theme
  CurveTheme <- {ggplot2::theme_bw()+
      theme(axis.title = element_text(face = 'bold'),
            plot.title = element_text(face = 'bold',
                                      hjust = 0.5),
            legend.title = element_text(face = 'bold',
                                        hjust = 0.5))}
  PaintFormatCurveData <- ConvertPaintFormatCurveData(Format_MeltCurve_Derivative)
  PaintFormatCurveData <- merge(x = PaintFormatCurveData,
                                y = Locs_Targets_Samples,
                                by.x = 'Colnames',
                                by.y = 'Well')
  PaintFormatCurveData$Rownames <- as.numeric(PaintFormatCurveData$Rownames)
  Plot_Basic <- ggplot(data=PaintFormatCurveData,
                       mapping=aes(x = Rownames,
                                   y = Value,
                                   group = Colnames))
  if (Group == 'none'){
    Plot_Advanced <- Plot_Basic+
      geom_line(color = '#008B00')+
      labs(x = 'Temperature',
           y = '-d(RFU)/dT',
           title = 'Melt Peak')+
      CurveTheme
  } else if (Group == 'target'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Target))+
      labs(x = 'Temperature',
           y = '-d(RFU)/dT',
           title = 'Melt Peak')+
      CurveTheme
  } else if (Group == 'treatment'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Sample))+
      labs(x = 'Temperature',
           y = '-d(RFU)/dT',
           title = 'Melt Peak')+
      CurveTheme
  }
  return(Plot_Advanced)
}
#Amplification_Curve
Plot_Amplification_Curve <- function(Format_Amplification,
                                      Locs_Targets_Samples,
                                      Group){
  library(ggplot2)
  #Define Curve Plot Style Theme
  CurveTheme <- {ggplot2::theme_bw()+
      theme(axis.title = element_text(face = 'bold'),
            plot.title = element_text(face = 'bold',
                                      hjust = 0.5),
            legend.title = element_text(face = 'bold',
                                        hjust = 0.5))}
  PaintFormatCurveData <- ConvertPaintFormatCurveData(Format_Amplification)
  PaintFormatCurveData <- merge(x = PaintFormatCurveData,
                                y = Locs_Targets_Samples,
                                by.x = 'Colnames',
                                by.y = 'Well')
  PaintFormatCurveData$Rownames <- as.numeric(PaintFormatCurveData$Rownames)
  Plot_Basic <- ggplot(data=PaintFormatCurveData,
                       mapping=aes(x = Rownames,
                                   y = Value,
                                   group = Colnames))
  if (Group == 'none'){
    Plot_Advanced <- Plot_Basic+
      geom_line(color = '#008B00')+
      labs(x = 'Cycle',
           y = 'RFU',
           title = 'Amplification')+
      CurveTheme
  } else if (Group == 'target'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Target))+
      labs(x = 'Cycle',
           y = 'RFU',
           title = 'Amplification')+
      CurveTheme
  } else if (Group == 'treatment'){
    Plot_Advanced <- Plot_Basic+
      geom_line(mapping = aes(color = Sample))+
      labs(x = 'Cycle',
           y = 'RFU',
           title = 'Amplification')+
      CurveTheme
  }
  return(Plot_Advanced)
}
