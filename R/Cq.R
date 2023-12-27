#Calculte Mean and Sd
Cq_Mean_Sd <- function(CqResults){
  Output<-data.frame(Target = rep(unique(CqResults$Target),
                                  each = length(unique(CqResults$Sample))),
                     Sample = rep(unique(CqResults$Sample),
                                  length(unique(CqResults$Target))),
                     Mean = NA,
                     Sd = NA)
  for (i in 1:nrow(Output)) {
    values <- CqResults$Cq[which(CqResults$Target == Output$Target[i] &
                                   CqResults$Sample == Output$Sample[i])]
    Output$Mean[i] <- mean(values)
    Output$Sd[i] <- sd(values, na.rm = T)
  }
  return(Output)
}

#Cq Theme
CqTheme <- function(){
  library(ggplot2)
  output<-
    theme_classic() +
    theme(axis.text = element_text(face = 'bold', size = 10),
                   axis.title = element_text(face = 'bold', size = 14),
                   legend.title = element_text(face = 'bold', hjust = 0.5),
                   axis.line = element_line(linewidth = 1.5),
                   plot.title = element_text(face = 'bold', hjust = 0.5))
  return(output)
  }

#Cq Value Col Plot
Plot_Col_UnRef <- function(Cq_Mean_Sd_Frame){
  library(ggplot2)
  Plot <- ggplot(data = Cq_Mean_Sd_Frame,
               mapping = aes(x = Target,
                             y = Mean,
                             fill = Sample)) +
    geom_col(position = position_dodge(0.85),
             width = 0.75) +
    scale_fill_manual(values = c('#696969',
                                 '#CFCFCF',
                                 '#828282',
                                 '#E8E8E8',
                                 '#4F4F4F'))+
    geom_errorbar(mapping = aes(ymin = Mean - Sd,
                                ymax = Mean + Sd),
                  position = position_dodge(0.85),
                  width = .3) +
    labs(x = NULL,
         y = 'Cq') +
    CqTheme() +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, max(Cq_Mean_Sd_Frame$Mean)*1.05))
  return(Plot)
}

#Cq DeltaCq
Plot_Delta_Cq <- function(Cq_Mean_Sd_Frame,
                         Ref_target){
  Other_target <-
    unique(Cq_Mean_Sd_Frame$Target)[-which(unique(Cq_Mean_Sd_Frame$Target) ==
                                                   Ref_target)]
  Output <- data.frame(Target = rep(Other_target,
                                    each = length(unique(Cq_Mean_Sd_Frame$Sample))),
                       Sample = rep(unique(Cq_Mean_Sd_Frame$Sample),
                                           length(Other_target)),
                       DeltaCq = NA,
                       Sd = NA)
  for (i in 1:nrow(Output)) {
    Output$DeltaCq[i] <-
      (Cq_Mean_Sd_Frame$Mean[which(Cq_Mean_Sd_Frame$Target == Output$Target[i] &
                                     Cq_Mean_Sd_Frame$Sample == Output$Sample[i])]) -
      (Cq_Mean_Sd_Frame$Mean[which(Cq_Mean_Sd_Frame$Target == Ref_target &
                                     Cq_Mean_Sd_Frame$Sample == Output$Sample[i])])
    Output$Sd[i] <-
      (Cq_Mean_Sd_Frame$Sd[which(Cq_Mean_Sd_Frame$Target == Output$Target[i] &
                                 Cq_Mean_Sd_Frame$Sample == Output$Sample[i])])
  }
  library(ggplot2)
  Plot <- ggplot(data = Output,
                 mapping = aes(x = Target,
                               y = DeltaCq,
                               fill = Sample)) +
    geom_col(position = position_dodge(0.85),
             width = 0.75) +
    scale_fill_manual(values = c('#696969',
                                 '#CFCFCF',
                                 '#828282',
                                 '#E8E8E8',
                                 '#4F4F4F'))+
    geom_errorbar(mapping = aes(ymin = DeltaCq - Sd,
                                ymax = DeltaCq + Sd),
                  position = position_dodge(0.85),
                  width = .3) +
    labs(x = NULL,
         y = 'Delta Cq') +
    CqTheme() +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, max(Output$DeltaCq)*1.05))
  return(Plot)
}

#Fold Change - 2^(-DeltaDeltaCq)
Plot_Fold_Change <- function(Cq_Mean_Sd_Frame,
                         Ref_target,
                         Control_sample){
  Other_target <-
    unique(Cq_Mean_Sd_Frame$Target)[-which(unique(Cq_Mean_Sd_Frame$Target) ==
                                             Ref_target)]
  Output <- data.frame(Target = rep(Other_target,
                                    each = length(unique(Cq_Mean_Sd_Frame$Sample))),
                       Sample = rep(unique(Cq_Mean_Sd_Frame$Sample),
                                    length(Other_target)),
                       DeltaCq = NA,
                       Sd = NA,
                       FC=NA,
                       Sd_FC=NA)
  for (i in 1:nrow(Output)) {
    Output$DeltaCq[i] <-
      (Cq_Mean_Sd_Frame$Mean[which(Cq_Mean_Sd_Frame$Target == Output$Target[i] &
                                     Cq_Mean_Sd_Frame$Sample == Output$Sample[i])]) -
      (Cq_Mean_Sd_Frame$Mean[which(Cq_Mean_Sd_Frame$Target == Ref_target &
                                     Cq_Mean_Sd_Frame$Sample == Output$Sample[i])])
    Output$Sd[i] <-
      (Cq_Mean_Sd_Frame$Sd[which(Cq_Mean_Sd_Frame$Target == Output$Target[i] &
                                   Cq_Mean_Sd_Frame$Sample == Output$Sample[i])])
  }
  for (i in 1:nrow(Output)) {
    Control_value <- Output$DeltaCq[which(Output$Target == Output$Target[i] &
                                            Output$Sample == Control_sample)]
    Output$FC[i] <- 2^(-(Output$DeltaCq[i]-Control_value))
    Output$Sd_FC[i] <- (Output$Sd[i]/Control_value)
  }

  library(ggplot2)
  Plot <- ggplot(data = Output,
                 mapping = aes(x = Target,
                               y = FC,
                               fill = Sample)) +
    geom_col(position = position_dodge(0.85),
             width = 0.75) +
    scale_fill_manual(values = c('#696969',
                                 '#CFCFCF',
                                 '#828282',
                                 '#E8E8E8',
                                 '#4F4F4F'))+
    geom_errorbar(mapping = aes(ymin = FC - Sd_FC,
                                ymax = FC + Sd_FC),
                  position = position_dodge(0.85),
                  width = .3) +
    labs(x = NULL,
         y = 'Fold Change') +
    CqTheme() +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, max(Output$FC)*1.05))
  return(Plot)
}

#Plate Visualization
Plate96_Cq <- function(CqResults){
  library(ggplot2)
  library(ggplate)
  return(plate_plot(
    data = CqResults,
    position = Well,
    value = Cq,
    label = round(Cq,2),
    plate_size = 96,
    label_size = 2.5,
    limits = c(0, NA)
  ))
}

#Cq PCA
CqPCA <- function(CqResults){
  PCA_data <- data.frame(matrix())
  df_pca <- prcomp(PCA_data[which(colnames(PCA_data) %in%
                                    c('Cq'))])
  df_pcs <-data.frame(df_pca$x)
}
