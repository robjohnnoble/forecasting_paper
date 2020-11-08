#!/usr/bin/Rscript
# script for creating plots and summary data from a batch of simulations (first argument)
library(demonanalysis)
require(data.table)
require(dplyr)
library(readr)


args <- commandArgs(trailingOnly = TRUE)
subfolder_name <- "Batch_RandomS_ForecastingPaper_Working"
n_cores <- 4

input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files

#all_statuses(input_dir, summary = TRUE) # all should be "Exit code 0"

if (!file.exists(output_dir_plots)){
  dir.create(output_dir_plots)
}
if (!file.exists(output_dir_data)){
  dir.create(output_dir_data)
}

# create_plots_batch(input_dir, output_dir = output_dir_plots, type = "chart", max_size = 10000, generation = NA)
# create_plots_batch(input_dir, output_dir = output_dir_plots, type = "plot")

data_RandomS <- all_output(input_dir, include_diversities = TRUE, max_generation = FALSE, n_cores = n_cores) # combined data for a batch of simulations, excluding diversity columns
unique(data_RandomS$s_driver_birth)
# test<-data_RandomMu[which(data_RandomMu$mu_driver_birth==1.472317e-05), which(colnames(data_RandomMu) %in% c("K", "s_driver_birth")) ]
# test<-data_RandomMu[which(data_RandomMu$mu_driver_birth==1.472317e-05), ]
fwrite(data_RandomS, file = paste0(output_dir_data, "/data_BatchRandomS.csv"), row.names = FALSE)
print("Wrote data_BatchRandomS.csv to file")


### Compute summary and correlation ### 

data <- read_csv(paste0(output_dir_data, "/data_BatchRandomS.csv"), guess_max = 1E4)

#creation of simulationDefinition
SimulationsDefinition_RandomS<-as.data.frame(unique(data[, which(colnames(data) %in% c("K", "mu_driver_birth", "s_driver_birth", "seed"))]))
SimulationsDefinition_RandomS<-SimulationsDefinition_RandomS[with(SimulationsDefinition_RandomS, order(mu_driver_birth,s_driver_birth, seed)), ]
SimulationsDefinition_RandomS<-SimulationsDefinition_RandomS %>% mutate(SimulationNumber = seq_along(SimulationsDefinition_RandomS[, 1]), 
                                                                          SimulationType =  ifelse(mu_driver_birth== 1e-06, 1, 
                                                                                                   ifelse(mu_driver_birth==1e-05, 2,3 )))

fwrite(SimulationsDefinition_RandomS, file = paste0(output_dir_data, "/SimulationsDefinition_RandomS.csv"), row.names = FALSE)

#SimulationsDefinition_RandomS<-read_csv(paste0(output_dir_data, "/SimulationsDefinition_RandomS.csv"), guess_max = 1E4)

data <- add_columns(data, num_parameters = num_parameters)
data <- add_relative_time(data, start_size = 250000, num_parameters = num_parameters)

data<-merge(data,SimulationsDefinition_RandomS, by= c("K", "mu_driver_birth", "s_driver_birth", "seed"))

#save the modification
fwrite(data, file = paste0(output_dir_data, "/data_BatchRandomSAugmentedK512.csv"), row.names = FALSE)
print(paste0("Wrote ", output_dir_data, "/data_BatchRandomSAugmentedK512.csv.csv to file"))

#data<-read_csv(paste0(output_dir_data, "/data_BatchRandomSAugmentedK512.csv"), guess_max = 1E4)

start_size_range<-c(10000,30000,60000,  90000,  125000, 250000, 375000, 500000, 625000, 750000 )
gap_range<-0.1
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
min_count<-20

summary <- get_summary(data, start_size_range, gap_range, final_size, num_parameters = num_parameters)
summary<- summary %>% mutate(AverageGrowthRate = (final_size-start_size) / waiting_time,
                             InverseWaitingTime = 1/ waiting_time)
fwrite(summary, file = paste0(output_dir_data, "/Summary_NewCorrelation_K512_RandomS.csv"), row.names = FALSE)

#summary<-read_csv(paste0(output_dir_data, "/Summary_NewCorrelation_K512_RandomS.csv"), guess_max = 1E4)


cols_list<-colnames(data)[grepl("DriverDiversity",colnames(data))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10),
                      "DriverDiversityFrom2RandomSamples",
                      "DriverDiversityFrom2BigRandomSamples")

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]
cols_list<-cols_list[!grepl("Depth0", cols_list)]
cols_list<-cols_list[!grepl("2Samples", cols_list)]
cols_list<-cols_list[!grepl("2BigSamples", cols_list)]
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers", "mean_autocor")



waitingTime_cor_summary_CombinedFitnessEffectK512<-get_wait_cor_summary_CombinedFitnessEffect(summary, cols_list,
                                                                                            num_parameters = num_parameters, min_count = min_count)
fwrite(waitingTime_cor_summary_CombinedFitnessEffectK512, file = paste0(output_dir_data, "/waitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"))


AverageGrowthRate_cor_summary_CombinedFitnessEffectK512<-get_AverageGrowthRate_cor_summary_CombinedFitnessEffect(summary, cols_list,
                                                                                                               num_parameters = num_parameters, min_count = min_count)
fwrite(AverageGrowthRate_cor_summary_CombinedFitnessEffectK512, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"))


InverseWaitingTime_cor_summary_CombinedFitnessEffectK512<-get_InverseWaitingTime_cor_summary_CombinedFitnessEffect(summary, cols_list,
                                                                                                                 num_parameters = num_parameters, min_count = min_count)
fwrite(InverseWaitingTime_cor_summary_CombinedFitnessEffectK512, file = paste0(output_dir_data, "/InverseWaitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"))


waitingTime_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0(output_dir_data, "/waitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)
InverseWaitingTime_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0(output_dir_data, "/InverseWaitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)


# dataset for kaplan meier plots
summary<-read_csv(paste0(output_dir_data, "/Summary_NewCorrelation_K512_RandomS.csv"), guess_max = 1E4)

SubsetRandomSdata<-summary[which(summary$mu_driver_birth==1e-06), which(colnames(summary) %in% c("DriverDiversity", "s_driver_birth",  "waiting_time", "MeanBirthRate"))]
fwrite(SubsetRandomSdata, file = paste0(output_dir_data, "/SubsetRandomSdata_ForKaplanMeierPlots.csv"))
