library(demonanalysis)
library(readr)
library(dplyr)
library(data.table)
library(moments)
library(gridExtra)
library(ggplot2)

input_dir <- "/all_results/Batch_TTT_100Seeds" # folder containing results of the batch
input_dir_data <- "/data/Batch_TTT_100Seeds" # folder containing data files
output_dir_data <- "/forecasting_analysis"
local_dir <- "/Data/ForecastingJune2020"

num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)

print("About to read")

df_all <- read_csv(paste0(input_dir_data, "/data_rescaledGeneration_WithExitCode4.csv"), guess_max = 1E4)

df_summary_list <- list()

for(iter in 1:20) {
  seedstart <- (iter - 1) * 5
  seedend <- iter * 5
  df <- filter(df_all, seed %% 1000 >= seedstart, seed %% 1000 < seedend)
  
  print("iter:")
  print(iter)
  print("data dimensions:")
  print(dim(df))
  # print(df[1:12, ])
  # 
  # print("new data dimensions:")
  # print(dim(df))
  
  df <- add_columns(df, 18) 
  # df<-add_relative_time(df, start_size=125000, 18)
  
  # df <- filter(df, JustBeforeTTT + JustAfterTTT > 0)
  
  # df[which(df$JustBeforeTTT == 1), "NumCells"] <- 1
  # df[which(df$JustAfterTTT == 1), "NumCells"] <- 2
  # df[which(df$JustBeforeTTT == 1), "Generation"] <- df[which(df$JustBeforeTTT == 1), "Generation"] * 0.9999
  # df[which(df$JustBeforeTTT == 1), "maxgen"] <- df[which(df$JustAfterTTT == 1), "maxgen"]
  
  # fwrite(df, file = paste0(output_dir_data, "/df_added_columns.csv"), row.names = FALSE)
  
  start_size_range<-c(-2, -1)
  gap_range<-0.1
  final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
  min_count<-5 #only 20 seeds per parameter sets
  
  df_summary_list[[iter]] <- get_summary(df, start_size_range, gap_range, final_size, num_parameters = num_parameters)
  
  df_summary_list[[iter]] <- mutate(df_summary_list[[iter]],
                     AverageGrowthRate = (final_size-start_size) / waiting_time,
                     ExactAverageGrowthRate = (final_size-NumCells) / waiting_time)
  
  remove(df)
}

df_summary <- bind_rows(df_summary_list)

fwrite(df_summary, file = paste0(output_dir_data, "/df_summary.csv"), row.names = FALSE)

print("All done")





