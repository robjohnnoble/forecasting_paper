#!/usr/bin/Rscript
# script for creating plots and summary data from a batch of simulations (first argument)
library(demonanalysis)
library(ggplot2)
require(data.table)
require(dplyr)
library(reshape)
library(tidyr)
library(readr)
library(RColorBrewer)
library(gridExtra)
library("lazyeval")
library(RVAideMemoire)
library(hexbin)
library(ggmuller)# function of demonanalysis
library(ggrepel)
library(moments) #skewness
library(arsenal)
library(gridExtra)
########################## Functions ##########################
#  function to add aes and aes_string 
`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

#function to compute diversity
inv_Simpson_index <- function(p) 1/sum(p * p)
inv_Simpson_index2 <- function(n, N) N * (N - 1)/sum(n * (n - 1))

#function to rescale the generation between -1 and 1 with 0 being the time were the treatment is given
ComputeScaledGeneration<-function(x, GenerationBeforeTTT, MaxGeneration){
  
  y<-vector(mode="numeric", length=length(x))
  
  for(i in seq_along(x)){
    
    if(x[i]<=GenerationBeforeTTT){
      y[i]<-((x[i]-GenerationBeforeTTT)/GenerationBeforeTTT)
    }else{
      y[i]<-((x[i]-GenerationBeforeTTT)/(MaxGeneration-GenerationBeforeTTT))
    }
  }
  
  return(y)
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
########################## Setting the folders ########################## 
subfolder_name<-"Batch_TTT_100Seeds"

n_cores <- 4

input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files
output_dir_plots_paper<- paste0(output_dir_plots, "/PlotsPaper")

#all_statuses(input_dir, summary = TRUE) # all should be "Exit code 0"
# Exit code 0 Exit code 4
#        1121          79

## First create the plot and data directories

ifelse(!dir.exists("plots"), dir.create("plots"), FALSE)
ifelse(!dir.exists("data"), dir.create("data"), FALSE)

## then create the subfolders (no recursive folder creation)

ifelse(!dir.exists(output_dir_plots), dir.create(output_dir_plots), FALSE)

ifelse(!dir.exists(output_dir_data), dir.create(output_dir_data), FALSE)


ifelse(!dir.exists(output_dir_plots_paper), dir.create(output_dir_plots_paper), FALSE)
########################## Variables ########################## 
KValue<-2^c(6, 9,12)
SValue<-c(0.05, 0.1, 0.15, 0.2)
MuValue<-c()

SampleSize<-c("1","4", "1Big", "4Big")
DiversityMeasure<-c("DriverDiversity",  "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
LabelsDiversityMeasure<-c("Clonal diversity", "Clonal diversity at boundary", "Mean cell division rate", "Mean number of drivers per cell")
LabelsDiversityMeasure_NoMeanDriverNb<-c("Clonal diversity", "Clonal diversity at boundary", "Mean cell division rate")
NewLabelsDiversityMeasure_NoMeanDriverNb<-c("Clonal diversity", "Clonal diversity at edge", "Mean cell division rate")

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))  

DiversityMeasureWithSweep<-c(DiversityMeasure, "mean_autocor")

DiversityMeasureWithSweepAndMuDriverBirth<-c(DiversityMeasureWithSweep, "mu_driver_birth")
LabelsDiversityMeasureWithSweep<-c(LabelsDiversityMeasure, "Mean clonal turnover")
LabelsDiversityMeasureWithSweepAndMuDriverBirth<-c(LabelsDiversityMeasure, "Mean clonal turnover", "Driver mutation rate")
SignificanceLevel=0.05

############## Data creation ############## 
data <- all_output_modified(input_dir, include_diversities = TRUE, max_generation = FALSE, n_cores = n_cores, ExitCode4=TRUE) # combined data for a batch of simulations, excluding diversity columns
fwrite(data, file = paste0(output_dir_data, "/dataWithExitCode4.csv"), row.names = FALSE)
print("Wrote dataWithExitCode4.csv to file")


############ Definition of the simulation ############

#creation of simulationDefinition
SimulationsDefinition<-as.data.frame(data)
SimulationsDefinition<-unique(SimulationsDefinition[, which(colnames(SimulationsDefinition) %in% c("K", "mu_driver_birth", "s_driver_birth", "seed"))])
SimulationsDefinition<-SimulationsDefinition[with(SimulationsDefinition, order(mu_driver_birth,s_driver_birth, seed)), ]
SimulationsDefinition<-SimulationsDefinition %>% mutate(SimulationNumber = seq_along(SimulationsDefinition[, 1]))

fwrite(SimulationsDefinition, file = paste0(output_dir_data, "/SimulationsDefinition.csv"), row.names = FALSE)

############## Add interesting columns ##############

data<-mutate(data,
                         JustAfterTTT_tmp=(c(0, diff(data$Treated))), 
                         JustBeforeTTT_tmp=c(diff(data$Treated), 0)
) %>% mutate(JustAfterTTT= ifelse(JustAfterTTT_tmp==1, 1, 0), 
             JustBeforeTTT=ifelse(JustBeforeTTT_tmp==1, 1, 0))
data$JustAfterTTT_tmp<-NULL
data$JustBeforeTTT_tmp<-NULL

data<-merge(data, SimulationsDefinition, by= c("K", "mu_driver_birth", "s_driver_birth", "seed"))


############## rescale the generation to center it at the TTT time ##############

for(i in unique(data$SimulationNumber)){
  
  tmpdata<-subset(data, data$SimulationNumber==i)
  GenerationTime<-tmpdata[which(tmpdata$JustBeforeTTT==1), "Generation"]
  
  
  # If length(GenerationTime)==0 we have a simulation with Exit code 4, 
  # i.e for which we have no information after the TTT because the tumor was killed
  # so we don't need to rescale the generations
  if( length(GenerationTime)==0){
    
    tmpdata<-tmpdata %>% mutate(ScaledGeneration = Generation, 
                                ExitCode4 =TRUE)
    
  }else{
    
    MaxGeneration<-max(tmpdata$Generation)
    
    tmpdata<-tmpdata %>% mutate(ScaledGeneration = ComputeScaledGeneration(Generation, 
                                                                           GenerationBeforeTTT=GenerationTime, 
                                                                           MaxGeneration=MaxGeneration), 
                                ExitCode4=FALSE)
    remove(MaxGeneration)
  }
  
  
  if(i==unique(data$SimulationNumber)[1]){
    #if(i==){
    data_rescaledGeneration<-tmpdata
  }else{
    data_rescaledGeneration<-rbind(data_rescaledGeneration,tmpdata )
  }
  
  remove(tmpdata)
  
  remove(GenerationTime)
  
}

fwrite(data_rescaledGeneration, file = paste0(output_dir_data, "/data_rescaledGeneration_WithExitCode4.csv"), row.names = FALSE)
############## Data loading ##############
SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)

dataWithExitCode4<-read_csv(paste0(output_dir_data, "/dataWithExitCode4.csv"), guess_max = 1E4)
data_rescaledGeneration_WithExitCode4<-read_csv( paste0(output_dir_data, "/data_rescaledGeneration_WithExitCode4.csv"),guess_max = 1E4)

# correlation analysis before TTT
wait_cor_summary_BeforeTTT<-read_csv( paste0(output_dir_data, "/waitingTime_cor_summary_dataBeforeTTT.csv"), guess_max = 1E4)
df_summary_BeforeTTT<-read_csv( paste0(output_dir_data, "/Summary_dataBeforeTTT.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_BeforeTTT<-read_csv( paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataBeforeTTT.csv"), guess_max = 1E4)
ExactAverageGrowthRate_cor_summary_BeforeTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataBeforeTTT.csv"), guess_max = 1E4)


Summary_allFinalSize_BeforeTTT<-read_csv( paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_BeforeTTT.csv"))
AverageGrowthRate_cor_summary_FinalSize_BeforeTTT<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_BeforeTTT.csv"))
ExactAverageGrowthRate_cor_summary__FinalSize_BeforeTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_BeforeTTT.csv"))

# correlation analysis after TTT
wait_cor_summary_AfterTTT<-read_csv( paste0(output_dir_data, "/waitingTime_cor_summary_dataAfterTTT.csv"), guess_max = 1E4)
# df_summary_AfterTTT<-read_csv( paste0(output_dir_data, "/Summary_dataAfterTTT.csv"), guess_max = 1E4)
# AverageGrowthRate_cor_summary_AfterTTT<-read_csv( paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataAfterTTT.csv"), guess_max = 1E4)
# ExactAverageGrowthRate_cor_summary_AfterTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataAfterTTT.csv"), guess_max = 1E4)
# 
# 
# Summary_allFinalSize_AfterTTT<-read_csv( paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_AfterTTT.csv"), guess_max = 1E4)
# AverageGrowthRate_cor_summary_FinalSize_AfterTTT<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT.csv"), guess_max = 1E4)
# ExactAverageGrowthRate_cor_summary_FinalSize_AfterTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT.csv"), guess_max = 1E4)

df_summary_AfterTTT<-read_csv( paste0(output_dir_data, "/Summary_dataAfterTTT_corrected.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_AfterTTT<-read_csv( paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataAfterTTT_corrected.csv"), guess_max = 1E4)
ExactAverageGrowthRate_cor_summary_AfterTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataAfterTTT_corrected.csv"), guess_max = 1E4)


Summary_allFinalSize_AfterTTT<-read_csv( paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_AfterTTT_corrected.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_FinalSize_AfterTTT<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT_corrected.csv"), guess_max = 1E4)
ExactAverageGrowthRate_cor_summary_FinalSize_AfterTTT<-read_csv(paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT_corrected.csv"), guess_max = 1E4)

############## Split data set ##############
#### we would like to compare the growth rate before and after the TTT
#### For this purpose, we could split the data set in 2 half, before and after the treatment and simply apply the function we had previously to compute the correlations.
#### Let's see if this works

dataBeforeTTT<-subset(data_rescaledGeneration_WithExitCode4, 
                      data_rescaledGeneration_WithExitCode4$Treated==0 & data_rescaledGeneration_WithExitCode4$K==512 & ! data_rescaledGeneration_WithExitCode4$s_driver_birth==0.15 )

dataAfterTTT<-subset(data_rescaledGeneration_WithExitCode4,
                     data_rescaledGeneration_WithExitCode4$Treated==1& data_rescaledGeneration_WithExitCode4$K==512 & ! data_rescaledGeneration_WithExitCode4$s_driver_birth==0.15)

remove(data_rescaledGeneration_WithExitCode4)

####### Before TTT : Correlation with different Start_site ####### 
start_size_range<-c(10000,30000,60000,  90000, 125000, 250000, 375000, 500000, 625000, 750000 )
gap_range<-0.1
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
min_count<-5 #only 20 seeds per parameter sets

df<-add_columns(dataBeforeTTT, 18) 
df<-add_relative_time(df, start_size=125000, 18)
df_summary<-get_summary(df, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

df_summary<-mutate(df_summary,
                   AverageGrowthRate = (final_size-start_size) / waiting_time, 
                   ExactAverageGrowthRate = (final_size-NumCells) / waiting_time)

cols_list<-colnames(df_summary)[grepl("DriverDiversity",colnames(df_summary))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10),
                      "DriverDiversityFrom2RandomSamples",
                      "DriverDiversityFrom2BigRandomSamples", 
                      paste0("CellsWith", 0:10 , "Drivers"))

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]
cols_list<-cols_list[!grepl("Depth0", cols_list)]
cols_list<-cols_list[!grepl("2Samples", cols_list)]
cols_list<-cols_list[!grepl("2BigSamples", cols_list)]
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers", "mean_autocor")
cols_list<-cols_list[! cols_list %in% "mu_driver_birth" ]


#the following lines generates summary dataframe of correlations with "waiting_time"

AverageGrowthRate_cor_summary<-get_cor_summary_RandomMu_TreatmentCase(df_summary,"AverageGrowthRate",  cols_list, num_parameters = num_parameters, min_count = min_count)


ExactAverageGrowthRate_cor_summary<-get_cor_summary_RandomMu_TreatmentCase(df_summary,"ExactAverageGrowthRate",  cols_list, num_parameters = num_parameters, min_count = min_count)

wait_cor_summary <-get_cor_summary_RandomMu_TreatmentCase(df_summary,"waiting_time",  cols_list, num_parameters = num_parameters, min_count = min_count)


# testwait_cor_summary <- get_wait_cor_summary_RandomMu_TreatmentCase(df_summary, cols_list, 
#                                                                     num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
# 
# testAverageGrowthRate_cor_summary<-get_AverageGrowthRate_cor_summary_RandomMu_TreatmentCase(df_summary, cols_list, 
#                                                                                             num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
# 
# testExactAverageGrowthRate_cor_summary<-get_ExactAverageGrowthRate_cor_summary_RandomMu_TreatmentCase(df_summary, cols_list, 
#                                                                                                       num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
# 
# 
# 
# comparedf(testwait_cor_summary,wait_cor_summary )
# comparedf(testAverageGrowthRate_cor_summary,AverageGrowthRate_cor_summary )
# comparedf(testExactAverageGrowthRate_cor_summary, ExactAverageGrowthRate_cor_summary)

fwrite(wait_cor_summary, file = paste0(output_dir_data, "/waitingTime_cor_summary_dataBeforeTTT.csv"), row.names = FALSE)
fwrite(df_summary, file = paste0(output_dir_data, "/Summary_dataBeforeTTT.csv"), row.names = FALSE)
fwrite(AverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataBeforeTTT.csv"), row.names = FALSE)
fwrite(ExactAverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataBeforeTTT.csv"), row.names = FALSE)

####### Before TTT : Correlation with different Final_size ####### 

start_size_range<-c(125000)
gap_range<-0.1
final_size <-c(250000, 375000, 500000, 625000, 750000,875000,  1E6) # waiting time is measured until tumour reaches this NumCells value
min_count<-20

Summary_allFinalSize<-vector(mode="list", length= length(final_size))

for(FinSizeIt in seq_along(final_size) ){
  
  Summary_allFinalSize[[FinSizeIt]] <- get_summary(df, start_size_range, gap_range, final_size[FinSizeIt], num_parameters = num_parameters) %>% mutate(FinalSize =final_size[FinSizeIt])
  
}
Summary_allFinalSize<-do.call(rbind, Summary_allFinalSize)

Summary_allFinalSize<-mutate(Summary_allFinalSize, 
                             AverageGrowthRate = (FinalSize-start_size) / waiting_time, 
                             ExactAverageGrowthRate = (FinalSize-NumCells) / waiting_time)

fwrite(Summary_allFinalSize, file = paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_BeforeTTT.csv"))

AverageGrowthRate_cor_summary_CombinedMutationRate<-get_cor_summary_RandomMu_VaryingFinalSize(Summary_allFinalSize,"AverageGrowthRate", 
                                                                                              col_names_list=DiversityMeasure,
                                                                                              num_parameters = num_parameters, min_count = min_count)

ExactAverageGrowthRate_cor_summary_CombinedMutationRate<-get_cor_summary_RandomMu_VaryingFinalSize(Summary_allFinalSize,"ExactAverageGrowthRate", 
                                                                                                   col_names_list=DiversityMeasure,
                                                                                                   num_parameters = num_parameters, min_count = min_count)


fwrite(AverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_BeforeTTT.csv"))

fwrite(ExactAverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_BeforeTTT.csv"))

########## Plots Before TTT Average growth time ########## 

####### Supp Fig4a ####### 

tmpSummary<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs Average growth rate \n Measurement before TTT"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

#### Fig4a comparison with ExactGrowthRate  ####### 
#actually useless because values of correlation are almost exactly the same.

tmpSummary1<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)
tmpSummary1<-mutate(tmpSummary1, 
                    type= "AverageGrowthRate")


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary1<-mutate(tmpSummary1,
                     tmpVar =(tmpSummary1[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary1[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary1)[which(colnames(tmpSummary1)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary1)[which(colnames(tmpSummary1)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

tmpSummary1$mean_AverageGrowthRate<-NULL

tmpSummary2<-subset(ExactAverageGrowthRate_cor_summary_BeforeTTT, !ExactAverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)
tmpSummary2<-mutate(tmpSummary2, 
                    type= "ExactAverageGrowthRate")

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary2<-mutate(tmpSummary2,
                      tmpVar =(tmpSummary2[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                      tmpVar2 =(tmpSummary2[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary2)[which(colnames(tmpSummary2)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary2)[which(colnames(tmpSummary2)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

tmpSummary2$mean_ExactAverageGrowthRate<-NULL

tmpSummary<-rbind(tmpSummary1, tmpSummary2)

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = type, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = type , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = type, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = type , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = type , color="c", shape ="solid"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = type, color="a", linetype=type))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = type , color="b", linetype=type))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = type , color="c", linetype=type))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="type",
                        values = c("AverageGrowthRate"="solid", "ExactAverageGrowthRate"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs Estimation of growth rate \n Measurement before TTT"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/ComparisonGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()



####### Supp Fig4c ####### 

tmpSummary<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(subset(tmpSummary, tmpSummary$K==512))
  
  for(SizeIt in seq_along(SampleSize)){
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    
    
    #if( ! SizeIt ==1){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    #}
    
                       
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
   
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="solid", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( . ~ s_driver_birth,
                labeller =label_bquote( cols= s== .(s_driver_birth)))+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("Correlation coefficient:\n Clonal diversity vs Average growth rate \n Before TTT"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_Correlations_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_TreatmentCase_BeforeTTT.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

##### Supp Fig3d  #####     

tmpSummary<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(subset(tmpSummary, 
                  tmpSummary$K==512)
)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  #then all points which are  significant
  #geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  #geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                   labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  # scale_size_manual(2)+
  facet_grid( . ~ s_driver_birth,labeller =label_bquote( cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu_NoS015_BeforeTTT.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

#### Supp Fig3b : Forecasting with fixeed start size and varying final size ########


tmpSummary<-subset(AverageGrowthRate_cor_summary_FinalSize_BeforeTTT, 
                   (AverageGrowthRate_cor_summary_FinalSize_BeforeTTT$K==512 & ! AverageGrowthRate_cor_summary_FinalSize_BeforeTTT$s_driver_birth==0.15))


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="endpoint size", breaks =c(250000, 375000, 500000, 625000,  750000, 875000, 1E6))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate \n Before TTT"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000_RandomMu_3DiversityMeasures_K512_BeforeTTT.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


#######After TTT :  Correlation with different Start_site  ####### 
start_size_range<-c(10000,30000,60000,  90000, 125000, 250000, 375000, 500000, 625000, 750000 )
gap_range<-0.1
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
min_count<-5 #only 20 seeds per parameter sets

# remove 2 simulations which reached 1e6 just after the TTT (no selection of cells/ all cells were selected together)
df<-subset(dataAfterTTT, 
           !dataAfterTTT$SimulationNumber  %in% c(1026, 1164))

# # We also remove simulation with ExitCode4 => NO ONE SHOULD EXIST IN df !!
# df<-subset(df, 
#            df$ExitCode4==FALSE)



df<-add_columns(df, 18) 
df<-add_relative_time(df, start_size=125000, 18)


df_summary<-get_summary(df, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

df_summary<-mutate(df_summary,
                   AverageGrowthRate = (final_size-start_size) / waiting_time, 
                   ExactAverageGrowthRate = (final_size-NumCells) / waiting_time)

cols_list<-colnames(df_summary)[grepl("DriverDiversity",colnames(df_summary))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10),
                      "DriverDiversityFrom2RandomSamples",
                      "DriverDiversityFrom2BigRandomSamples", 
                      paste0("CellsWith", 0:10 , "Drivers"))

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]
cols_list<-cols_list[!grepl("Depth0", cols_list)]
cols_list<-cols_list[!grepl("2Samples", cols_list)]
cols_list<-cols_list[!grepl("2BigSamples", cols_list)]
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers", "mean_autocor")
cols_list<-cols_list[! cols_list %in% "mu_driver_birth" ]


#the following lines generates summary dataframe of correlations with "waiting_time"

AverageGrowthRate_cor_summary<-get_cor_summary_RandomMu_TreatmentCase(df_summary,"AverageGrowthRate",  cols_list, num_parameters = num_parameters, min_count = min_count)


ExactAverageGrowthRate_cor_summary<-get_cor_summary_RandomMu_TreatmentCase(df_summary,"ExactAverageGrowthRate",  cols_list, num_parameters = num_parameters, min_count = min_count)

wait_cor_summary <-get_cor_summary_RandomMu_TreatmentCase(df_summary,"waiting_time",  cols_list, num_parameters = num_parameters, min_count = min_count)


# fwrite(wait_cor_summary, file = paste0(output_dir_data, "/waitingTime_cor_summary_dataAfterTTT.csv"), row.names = FALSE)
# fwrite(df_summary, file = paste0(output_dir_data, "/Summary_dataAfterTTT.csv"), row.names = FALSE)
# fwrite(AverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataAfterTTT.csv"), row.names = FALSE)
# fwrite(ExactAverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataAfterTTT.csv"), row.names = FALSE)

fwrite(wait_cor_summary, file = paste0(output_dir_data, "/waitingTime_cor_summary_dataAfterTTT_corrected.csv"), row.names = FALSE)
fwrite(df_summary, file = paste0(output_dir_data, "/Summary_dataAfterTTT_corrected.csv"), row.names = FALSE)
fwrite(AverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_dataAfterTTT_corrected.csv"), row.names = FALSE)
fwrite(ExactAverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_dataAfterTTT_corrected.csv"), row.names = FALSE)


#######AfterTTT :  Correlation with different Final_size ####### 

start_size_range<-c(125000)
gap_range<-0.1
final_size <-c(250000, 375000, 500000, 625000, 750000,875000,  1E6) # waiting time is measured until tumour reaches this NumCells value
min_count<-20

Summary_allFinalSize<-vector(mode="list", length= length(final_size))

# remove 2 simulations which reached 1e6 just after the TTT (no selection of cells/ all cells were selected together)
df<-subset(dataAfterTTT, 
           !dataAfterTTT$SimulationNumber  %in% c(1026, 1164))

df<-add_columns(df, 18) 
df<-add_relative_time(df, start_size=125000, 18)


for(FinSizeIt in seq_along(final_size) ){
  
  Summary_allFinalSize[[FinSizeIt]] <- get_summary(df, start_size_range, gap_range, final_size[FinSizeIt], num_parameters = num_parameters) %>% mutate(FinalSize =final_size[FinSizeIt])
  
}
Summary_allFinalSize<-do.call(rbind, Summary_allFinalSize)

Summary_allFinalSize<-mutate(Summary_allFinalSize, 
                             AverageGrowthRate = (FinalSize-start_size) / waiting_time, 
                             ExactAverageGrowthRate = (FinalSize-NumCells) / waiting_time)

# fwrite(Summary_allFinalSize, file = paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_AfterTTT.csv"))
fwrite(Summary_allFinalSize, file = paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000_AfterTTT_corrected.csv"))

Summary_allFinalSize_AfterTTT<-Summary_allFinalSize
remove(Summary_allFinalSize)
Summary_allFinalSize_AfterTTT<-subset(Summary_allFinalSize_AfterTTT, ! Summary_allFinalSize_AfterTTT$waiting_time==0)
cols_list<-DiversityMeasure

AverageGrowthRate_cor_summary_CombinedMutationRate<-get_cor_summary_RandomMu_VaryingFinalSize(Summary_allFinalSize_AfterTTT,"AverageGrowthRate", 
                                                                                              cols_list,
                                                                                              num_parameters = num_parameters, min_count = 5)

ExactAverageGrowthRate_cor_summary_CombinedMutationRate<-get_cor_summary_RandomMu_VaryingFinalSize(Summary_allFinalSize_AfterTTT,"ExactAverageGrowthRate", 
                                                                                                   cols_list,
                                                                                                   num_parameters = num_parameters, min_count = 5)



# fwrite(AverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT.csv"))
# 
# fwrite(ExactAverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT.csv"))

fwrite(AverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT_corrected.csv"))

fwrite(ExactAverageGrowthRate_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/ExactAverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_RandomMu_AfterTTT_corrected.csv"))


##### search for the bug##### 
setdiff( (which(is.na(df$GrowthRate))+1), which(!is.na(df$GrowthRate) & is.na(lag(df$GrowthRate, 1))))
unique(df[setdiff( (which(is.na(df$GrowthRate))+1), which(!is.na(df$GrowthRate) & is.na(lag(df$GrowthRate, 1)))), "SimulationNumber"])

test<-subset(df, df$SimulationNumber %in% c(1026, 1032, 1164, 1166))

test<-subset(df, df$SimulationNumber %in% c(1026)) # 2 Nan out if 2 rows =>

test<-subset(df, df$SimulationNumber %in% c(1032)) # 1 Na out of 209 rows
test<-subset(df, df$SimulationNumber %in% c(1164)) # 2 Nan out if 2 rows
test<-subset(df, df$SimulationNumber %in% c(1166)) # 1 Na out of 226 rows

data_rescaledGeneration_WithExitCode4<-read_csv( paste0(output_dir_data, "/data_rescaledGeneration_WithExitCode4.csv"),guess_max = 1E4)
test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1026 & data_rescaledGeneration_WithExitCode4$Treated==1)
test2[, c("NumCells", "Treated")]

test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1026)
as.data.frame(test2[, c("NumCells", "Passengers")])


test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1164 & data_rescaledGeneration_WithExitCode4$Treated==1)
test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1164)
as.data.frame(test2[, c("NumCells", "Passengers", "Treated")])


test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1032 & data_rescaledGeneration_WithExitCode4$Treated==1)
test2<-subset(data_rescaledGeneration_WithExitCode4, data_rescaledGeneration_WithExitCode4$SimulationNumber == 1032)
as.data.frame(test2[, c("NumCells", "Passengers", "Treated")])


#20 simulation with exit code 4
unique(dataBeforeTTT[which(dataBeforeTTT$ExitCode4==TRUE), which(colnames(dataBeforeTTT) %in% c("K", "seed", "mu_driver_birth", "s_driver_birth"))])
# after TTT, remaining 280 simulations without exit code 4
unique(dataAfterTTT[which(dataAfterTTT$ExitCode4==TRUE), which(colnames(dataAfterTTT) %in% c("K", "seed", "mu_driver_birth", "s_driver_birth"))])


test<-subset(dataAfterTTT, 
       !dataAfterTTT$SimulationNumber  %in% c(1026, 1164))

test<-add_columns(test, 18) 
test<-add_relative_time(test, start_size=125000, 18)
unique(test[, which(colnames(test) %in% c("K", "seed", "mu_driver_birth", "s_driver_birth"))])

test_summary<-get_summary(test, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
#in debugging get_summary
AllSimulationNumber<-summary[which(summary$start_size==375000), "SimulationNumber"]
setdiff(AllSimulationNumber, summary[which(summary$start_size==500000), "SimulationNumber"] ) #61
setdiff(AllSimulationNumber, summary[which(summary$start_size==625000), "SimulationNumber"] ) # 61
setdiff(AllSimulationNumber, summary[which(summary$start_size==750000), "SimulationNumber"] ) # 61

as.data.frame(test[which(test$SimulationNumber==61), "NumCells"]) # largest NumCells is 448171 < 500000
as.data.frame(unique(test[which(test$SimulationNumber==61),c(1:18)]))

test2<-subset(dataAfterTTT, 
             !dataAfterTTT$SimulationNumber  %in% c(1026, 1164, 61))

test2<-add_columns(test2, 18) 
test2<-add_relative_time(test2, start_size=125000, 18)
unique(test2[, which(colnames(test2) %in% c("K", "seed", "mu_driver_birth", "s_driver_birth"))])

test_summary<-get_summary(test2, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

# we get differnent summary
all.equal(df_summary[which(df_summary$SimulationNumber==41), ],
          df_summary_BeforeTTT[which(df_summary_BeforeTTT$SimulationNumber==41), ] )

####### Supp Fig4a ####### 

tmpSummary<-subset(AverageGrowthRate_cor_summary_AfterTTT, !AverageGrowthRate_cor_summary_AfterTTT$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
#   #first all points which are not significant
#   geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
#   #then all points which are  significant
#   geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
#   geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
#   # dashed line with all points
#   geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
#   #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
#   scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
#   #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
#   scale_linetype_manual(name="",
#                         values = c("dashed"="solid", "solid"="solid"))+
#   scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#                      values = c("dashed"=17, "solid"=16),
#                      labels=c("dashed"="FALSE", "solid"="TRUE"))+
#   # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
#   #                     values = c("solid"=16),
#   #                     labels=c("solid"="TRUE"))+
#   scale_size_manual(2)+
#   guides(linetype=FALSE)+
#   facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
#   theme_bw(base_size = 16)+
#   ylim(-1, 1)+
#   geom_hline(aes(yintercept=0), linetype="dashed")+
#   #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
#   scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
#   ylab(paste0("correlation coefficient:\n Predictor variables vs Average growth rate \n Measurement after TTT"))+
#   ggtitle("")+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)), 
#   #        shape= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)))+
#   guides(linetype=FALSE,
#          color= guide_legend(title.theme = element_text(size = 14),
#                              label.theme = element_text(size=14)), 
#          shape= guide_legend(title.theme = element_text(size = 14),
#                              label.theme = element_text(size=12)))+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
#         aspect.ratio=1)


# png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_AfterTTT.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()


gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs Average growth rate \n Measurement after TTT"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)
png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_AfterTTT_corrected.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


#### Supp Fig3b : Forecasting with fixeed start size and varying final size ########


tmpSummary<-subset(AverageGrowthRate_cor_summary_FinalSize_AfterTTT, 
                   (AverageGrowthRate_cor_summary_FinalSize_AfterTTT$K==512 & ! AverageGrowthRate_cor_summary_FinalSize_AfterTTT$s_driver_birth==0.15))


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="endpoint size", breaks =c(250000, 375000, 500000, 625000,  750000, 875000, 1E6))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate \n After TTT"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000_RandomMu_3DiversityMeasures_K512_AfterTTT.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Comparison before and after##### 

###### Fig 4a

tmpSummaryBefore<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


SampleSize<-c("1","4", "1Big", "4Big")
DepthIt<-1
Denomination<-c("a", "b", "c", "d")

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummaryBefore<-mutate(tmpSummaryBefore,
                     tmpVar =(tmpSummaryBefore[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummaryBefore[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummaryBefore)[which(colnames(tmpSummaryBefore)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummaryBefore)[which(colnames(tmpSummaryBefore)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}
tmpSummaryBefore<-mutate(tmpSummaryBefore, 
                         Time="Before")


tmpSummaryAfter<-subset(AverageGrowthRate_cor_summary_AfterTTT, !AverageGrowthRate_cor_summary_AfterTTT$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummaryAfter<-mutate(tmpSummaryAfter,
                           tmpVar =(tmpSummaryAfter[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                           tmpVar2 =(tmpSummaryAfter[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummaryAfter)[which(colnames(tmpSummaryAfter)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummaryAfter)[which(colnames(tmpSummaryAfter)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}
tmpSummaryAfter<-mutate(tmpSummaryAfter, 
                         Time="After")

tmpSummary<-rbind(tmpSummaryBefore, tmpSummaryAfter)


gp <- list()

gg<-ggplot(tmpSummary)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), 
                     labels=c("a"="Clonal diversity", "b"="Clonal diversity at edge", "c"="Mean cell division rate"))+
  #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("Before"="solid", "After"="dashed"), 
                        labels = c("Before"="Before TTT", "After"="After TTT"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  #guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables \n vs Average growth rate"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)

gp[[1]] <- gg+
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = Time, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = Time, color="a", shape ="solid"), size=2)+
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = Time, color="a", linetype=Time))
  
gp[[2]] <- gg+
  geom_point(aes(x= start_size, y =ifelse( DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = Time, color="b", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = Time, color="b", shape ="solid"), size=2)+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = Time, color="b", linetype=Time))

gp[[3]] <- gg+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = Time, color="c", shape ="solid"), size=2)+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = Time, color="c", linetype=Time))
  



gg1<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="Before"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="Before"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="After"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), 
                     labels=c("a"="Clonal diversity", "b"="Clonal diversity at edge", "c"="Mean cell division rate"))+
  scale_linetype_manual(name="",
                        values = c("Before"="solid", "After"="dashed"), 
                        labels = c("Before"="Before TTT", "After"="After TTT"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  scale_size_manual(2)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables \n vs Average growth rate"))+
  ggtitle("")+
  guides(linetype=guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


mylegend<-g_legend(gg1)


png(paste0(output_dir_plots_paper, "/ComparisonBeforeAfter_AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase.png"),  width = 1000, height = 1000, res = 100)
ggF<-grid.arrange(arrangeGrob(gp[[1]]+ theme(legend.position="none"), 
                              gp[[2]]+ theme(legend.position="none"), 
                              gp[[3]]+ theme(legend.position="none"),
                              mylegend,
                              layout_matrix=rbind(c(1,1,4),c(2,2,4), c(3,3,4))))
dev.off()





#### Supp Fig3b : Forecasting with fixeed start size and varying final size ########

tmpSummaryBefore<-subset(AverageGrowthRate_cor_summary_FinalSize_BeforeTTT, 
                        (AverageGrowthRate_cor_summary_FinalSize_BeforeTTT$K==512 & ! AverageGrowthRate_cor_summary_FinalSize_BeforeTTT$s_driver_birth==0.15))


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummaryBefore<-mutate(tmpSummaryBefore,
                          tmpVar =(tmpSummaryBefore[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                          tmpVar2 =(tmpSummaryBefore[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummaryBefore)[which(colnames(tmpSummaryBefore)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummaryBefore)[which(colnames(tmpSummaryBefore)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

tmpSummaryBefore<-mutate(tmpSummaryBefore, 
                        Time="Before")

tmpSummaryAfter<-subset(AverageGrowthRate_cor_summary_FinalSize_AfterTTT, 
                   (AverageGrowthRate_cor_summary_FinalSize_AfterTTT$K==512 & ! AverageGrowthRate_cor_summary_FinalSize_AfterTTT$s_driver_birth==0.15))


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummaryAfter<-mutate(tmpSummaryAfter,
                     tmpVar =(tmpSummaryAfter[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummaryAfter[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummaryAfter)[which(colnames(tmpSummaryAfter)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummaryAfter)[which(colnames(tmpSummaryAfter)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

tmpSummaryAfter<-mutate(tmpSummaryAfter, 
                        Time="After")


tmpSummary<-rbind(tmpSummaryBefore, tmpSummaryAfter)


gp <- list()


gg<-ggplot(tmpSummary)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), 
                     labels=c("a"="Clonal diversity", "b"="Clonal diversity at edge", "c"="Mean cell division rate"))+
  #scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("Before"="solid", "After"="dashed"), 
                        labels = c("Before"="Before TTT", "After"="After TTT"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  #guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="endpoint size", breaks =c(250000, 375000, 500000, 625000,  750000, 875000, 1E6))+
  ylab(paste0("correlation coefficient:\n Predictor variables \n vs average growth rate"))+
  ggtitle("")+
  guides(linetype=guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


gp[[1]] <- gg+
  geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = Time, color="a", shape ="solid"), size=2)+
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = Time, color="a", linetype=Time))

gp[[2]] <- gg+
  geom_point(aes(x= FinalSize, y =ifelse( DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = Time, color="b", shape ="solid"), size=2)+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = Time, color="b", linetype=Time))


gp[[3]] <- gg+
  geom_point(aes(x= FinalSize, y =ifelse( MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = Time, color="c", shape ="solid"), size=2)+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = Time, color="c", linetype=Time))


gg1<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="Before"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="Before"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="After"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), 
                     labels=c("a"="Clonal diversity", "b"="Clonal diversity at edge", "c"="Mean cell division rate"))+
  scale_linetype_manual(name="",
                        values = c("Before"="solid", "After"="dashed"), 
                        labels = c("Before"="Before TTT", "After"="After TTT"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  scale_size_manual(2)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables \n vs Average growth rate"))+
  ggtitle("")+
  guides(linetype=guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


mylegend<-g_legend(gg1)


png(paste0(output_dir_plots_paper, "/ComparisonBeforeAfter_AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000_RandomMu_3DiversityMeasures_K512.png"),  width = 1000, height = 1000, res = 100)
ggF<-grid.arrange(arrangeGrob(gp[[1]]+ theme(legend.position="none"), 
                              gp[[2]]+ theme(legend.position="none"), 
                              gp[[3]]+ theme(legend.position="none"),
                              mylegend,
                              layout_matrix=rbind(c(1,1,4),c(2,2,4), c(3,3,4))))
dev.off()

##### get data for Kapplan Meyer plot ####
DataKM<-select( df_summary_AfterTTT, c("K", 
                                       "mu_driver_birth",
                                       "s_driver_birth",
                                       "waiting_time", 
                                       "DriverDiversity", "ExitCode4"))


fwrite(DataKM, file = paste0(output_dir_data, "/DataKM.csv"), row.names = FALSE)

NewDataKM<-select( df_summary_AfterTTT, c("K", 
                                          "mu_driver_birth",
                                          "s_driver_birth",
                                          "waiting_time", 
                                          "DriverDiversity", "ExitCode4", "start_size", "seed", "MeanBirthRate"))

DiversityJustBeforeTTT<-dataBeforeTTT[which(dataBeforeTTT$JustBeforeTTT==1),  
                     which(colnames(dataBeforeTTT) %in% c("K","mu_driver_birth","s_driver_birth","ExitCode4", "DriverDiversity", "MeanBirthRate"))]
NewDataKM<-merge(NewDataKM, DiversityJustBeforeTTT, by=c("K","mu_driver_birth","s_driver_birth","ExitCode4"))

colnames(NewDataKM)[which(colnames(NewDataKM)=="DriverDiversity.x")]<-"DriverDiversity_df_summary_AfterTTT"
colnames(NewDataKM)[which(colnames(NewDataKM)=="DriverDiversity.y")]<-"DriverDiversityJustBeforeTTT"

colnames(NewDataKM)[which(colnames(NewDataKM)=="MeanBirthRate.x")]<-"MeanBirthRate_df_summary_AfterTTT"
colnames(NewDataKM)[which(colnames(NewDataKM)=="MeanBirthRate.y")]<-"MeanBirthRateJustBeforeTTT"


DiversityJustAfterTTT<-dataAfterTTT[which(dataAfterTTT$JustAfterTTT==1),  
                                      which(colnames(dataAfterTTT) %in% c("K","mu_driver_birth","s_driver_birth","ExitCode4", "DriverDiversity", "MeanBirthRate"))]
NewDataKM<-merge(NewDataKM, DiversityJustAfterTTT, by=c("K","mu_driver_birth","s_driver_birth","ExitCode4"))
colnames(NewDataKM)[which(colnames(NewDataKM)=="DriverDiversity")]<-"DriverDiversityJustAfterTTT"
colnames(NewDataKM)[which(colnames(NewDataKM)=="MeanBirthRate")]<-"MeanBirthRateJustAfterTTT"

fwrite(NewDataKM, file = paste0(output_dir_data, "/NewDataKM.csv"), row.names = FALSE)

##### Change the shape of the points #####



tmpSummary<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", color ="a"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", color ="b"),shape=21, size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", color ="solid"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", color ="solid"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c", color ="solid"),shape=21, size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
   #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  scale_fill_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "solid"="black"),
                     labels=c("a"="", "b"="", "c"="", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs Average growth rate \n Measurement before TTT"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         fill= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShape.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", color ="dashed"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", color ="dashed"),shape=21, size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", color ="a"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", color ="b"),shape=21, size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c", color ="c"),shape=21, size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  scale_fill_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "dashed"="black"),
                     labels=c("a"="", "b"="", "c"="","dashed"="FALSE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs Average growth rate \n Measurement before TTT"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         fill= guide_legend(title.theme = element_text(size = 14),
                            label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeOpposit.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Rob version : Change the shape of the points #####



tmpSummary<-subset(AverageGrowthRate_cor_summary_BeforeTTT, !AverageGrowthRate_cor_summary_BeforeTTT$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


thin_width = 0.1
thick_width = 2
point_size = 2



gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c"), linetype="solid")+
  
  # #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  #then all points which are  significant
  geom_point(aes(x= start_size, 
                 y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b",  stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c",  stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
  
  scale_size_manual(point_size)+
  #guides(linetype=FALSE)+
  guides(linetype=FALSE,
         color=FALSE,
         fill = guide_legend(override.aes = list(stroke = 0))) + 
  continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
                   palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                   breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs progression-free survival \n Measurement before Treatment"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 14),
  #                            label.theme = element_text(size=14)), 
  #        fill= guide_legend(title.theme = element_text(size = 14),
  #                           label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeRob.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


### other way
gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c"), linetype="solid")+
  
  # #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", stroke=DriverDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", stroke=DriverEdgeDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  #then all points which are  significant
  geom_point(aes(x= start_size, 
                 y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b",  stroke=DriverEdgeDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width), shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c",  stroke=MeanBirthRate_Unsignificant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
  scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
  
  scale_size_manual(point_size)+
  #guides(linetype=FALSE)+
  guides(linetype=FALSE,
         color=FALSE,
         fill = guide_legend(override.aes = list(stroke = 0))) + 
  continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
                   palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                   breaks = c(thin_width, thick_width), labels = c("yes", "no"))+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variables vs progression-free survival \n Measurement before Treatment"))+
  ggtitle("")+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)), 
  #        shape= guide_legend(title.theme = element_text(size = 12),
  #                            label.theme = element_text(size=10)))+
  # guides(linetype=FALSE,
  #        color= guide_legend(title.theme = element_text(size = 14),
  #                            label.theme = element_text(size=14)), 
  #        fill= guide_legend(title.theme = element_text(size = 14),
  #                           label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeRobOpposite.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()





# 
# gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
#   # dashed line with all points
#   geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
#   
#   #first all points which are not significant
#   geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", stroke=DriverDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", stroke=DriverEdgeDiversity_Unsignificant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
#   #then all points which are  significant
#   geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a",  stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b",  stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c",  stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
#   #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
#   scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
#   scale_linetype_manual(name="",
#                         values = c("dashed"="solid", "solid"="solid"))+
#   # scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#   #                    values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "dashed"="black"),
#   #                    labels=c("a"="", "b"="", "c"="","dashed"="FALSE"))+
#   scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#                      values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
#   # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
#   #                     values = c("solid"=16),
#   #                     labels=c("solid"="TRUE"))+
#   #scale_size_manual(2)+
#   #guides(linetype=FALSE)+
#   guides(linetype=FALSE,
#          color=FALSE,
#          fill = guide_legend(override.aes = list(stroke = 0))) + 
#   continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
#                    palette = function(x){scales::rescale(x, c(thin_width, thick_width))},
#                    breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
#   facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
#   theme_bw(base_size = 16)+
#   ylim(-1, 1)+
#   geom_hline(aes(yintercept=0), linetype="dashed")+
#   #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
#   scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
#   ylab(paste0("correlation coefficient:\n Predictor variables vs progression-free survival \n Measurement before Treatment"))+
#   ggtitle("")+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)), 
#   #        shape= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)))+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 14),
#   #                            label.theme = element_text(size=14)), 
#   #        fill= guide_legend(title.theme = element_text(size = 14),
#   #                           label.theme = element_text(size=12)))+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
#         aspect.ratio=1)
# 
# 
# png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeRob.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()
# 
# 
# gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
#   # dashed line with all points
#   geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
#   
#   #first all points which are not significant
#   geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", stroke=thin_width),shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", stroke=thin_width),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
#   #then all points which are  significant
#   geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a",  stroke=thick_width), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b",  stroke=thick_width), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c",  stroke=thick_width),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
#   #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
#   scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
#   scale_linetype_manual(name="",
#                         values = c("dashed"="solid", "solid"="solid"))+
#   # scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#   #                    values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "dashed"="black"),
#   #                    labels=c("a"="", "b"="", "c"="","dashed"="FALSE"))+
#   scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#                      values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
#   # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
#   #                     values = c("solid"=16),
#   #                     labels=c("solid"="TRUE"))+
#   #scale_size_manual(2)+
#   #guides(linetype=FALSE)+
#   guides(linetype=FALSE,
#          color=FALSE,
#          fill = guide_legend(override.aes = list(stroke = 0))) + 
#   # continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
#   #                  palette = function(x){scales::rescale(x, c(thin_width, thick_width))},
#   #                  breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
#   discrete_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
#                    palette = function(x){if(x){thin_width}else{thick_width}},
#                    breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
#   facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
#   theme_bw(base_size = 16)+
#   ylim(-1, 1)+
#   geom_hline(aes(yintercept=0), linetype="dashed")+
#   #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
#   scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
#   ylab(paste0("correlation coefficient:\n Predictor variables vs progression-free survival \n Measurement before Treatment"))+
#   ggtitle("")+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)), 
#   #        shape= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)))+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 14),
#   #                            label.theme = element_text(size=14)), 
#   #        fill= guide_legend(title.theme = element_text(size = 14),
#   #                           label.theme = element_text(size=12)))+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
#         aspect.ratio=1)
# 
# 
# png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeRob.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()
# 
# 
# gg<-ggplot(subset(tmpSummary,tmpSummary$K==512) )+
#   # dashed line with all points
#   geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
#   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
#   
#   #first all points which are not significant
#   geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a", stroke=FALSE),shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b", stroke=FALSE),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
#   #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
#   #then all points which are  significant
#   geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, fill="a",  stroke=TRUE), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , fill="b",  stroke=TRUE), shape=21, size=point_size)+
#   geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , fill="c",  stroke=TRUE),shape=21, size=point_size)+
#   #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
#   #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
#   scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
#   scale_linetype_manual(name="",
#                         values = c("dashed"="solid", "solid"="solid"))+
#   # scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#   #                    values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "dashed"="black"),
#   #                    labels=c("a"="", "b"="", "c"="","dashed"="FALSE"))+
#   scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
#                      values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
#   # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
#   #                     values = c("solid"=16),
#   #                     labels=c("solid"="TRUE"))+
#   #scale_size_manual(2)+
#   #guides(linetype=FALSE)+
#   guides(linetype=FALSE,
#          color=FALSE,
#          fill = guide_legend(override.aes = list(stroke = 0))) + 
#   # continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
#   #                  palette = function(x){scales::rescale(x, c(thin_width, thick_width))},
#   #                  breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
#   discrete_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
#                  palette = function(x){if(x){thin_width}else{thick_width}},
#                  breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
#   facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))+
#   theme_bw(base_size = 16)+
#   ylim(-1, 1)+
#   geom_hline(aes(yintercept=0), linetype="dashed")+
#   #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
#   scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
#   ylab(paste0("correlation coefficient:\n Predictor variables vs progression-free survival \n Measurement before Treatment"))+
#   ggtitle("")+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)), 
#   #        shape= guide_legend(title.theme = element_text(size = 12),
#   #                            label.theme = element_text(size=10)))+
#   # guides(linetype=FALSE,
#   #        color= guide_legend(title.theme = element_text(size = 14),
#   #                            label.theme = element_text(size=14)), 
#   #        fill= guide_legend(title.theme = element_text(size = 14),
#   #                           label.theme = element_text(size=12)))+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
#         aspect.ratio=1)
# 
# 
# png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_DiversityMeasuresNoDrivers_TreatmentCase_BeforeTTT_testShapeRob.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()

#######
