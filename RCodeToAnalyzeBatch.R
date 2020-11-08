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


########################## Functions ##########################
inv_Simpson_index <- function(p) 1/sum(p * p)
inv_Simpson_index2 <- function(n, N) N * (N - 1)/sum(n * (n - 1))


########################## Modified function ########################## 
plot_corr_waiting_time_versus_start_size_ModifColor <- function(df, col_name = "DriverDiversity", output_filename = "corr_waiting_time_versus_start_size", file_type = "png", output_dir = NA, title="") {
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir, output_filename, ".png"), width = 500, height = 500, res = 100)
    else pdf(paste0(output_dir, output_filename, ".pdf"), width = 5, height = 5)
  }
  
  if(substr(col_name, 1, 4) != "Cor_") col_name <- paste0("Cor_", col_name)
  col_name_suffix <- substr(col_name, 5, nchar(col_name))
  
  start_size_range <- unique(df$start_size)
  
  Param_range <- sort(unique(df$SimulationType))
  VecParam<-vector(mode="character", length=length(Param_range))
  # 
  # for(SimType in seq_along(Param_range)){
  #   #tmp<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("K", "mu_driver_birth","s_driver_birth" ))]
  #   KVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("K"))]
  #   MuVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("mu_driver_birth"))]
  #   SVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("s_driver_birth"))]
  # 
  #   VecParam[SimType]<-paste0("K",KVal, "M", MuVal, "S", SVal )
  #   
  # }
  # 
  
  for(SimType in seq_along(Param_range)){
    #tmp<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("K", "mu_driver_birth","s_driver_birth" ))]
    #KVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("K"))]
    MuVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("mu_driver_birth"))]
    SVal<-df[which(df$SimulationType==Param_range[SimType])[1], which(colnames(df) %in% c("s_driver_birth"))]
    
    VecParam[SimType]<-paste0("M", MuVal, "S", SVal )
    
  }
  
  
  
  #cols <- rainbow(length(Param_range))
  cols<-c(brewer.pal(8, "YlOrRd")[5:8],brewer.pal(8, "BuGn")[5:8],brewer.pal(8, "Blues")[5:8])
  
  
  
  par(mfrow=c(1, 1))
  par(mar=c(4, 4, 2, 2))
  for(j in 1:length(Param_range)) {
    Param_val <- Param_range[j]
    if(j == 1) {
      plot(1, type = "n", xlim = c(0, max(start_size_range)), ylim = c(-1, 1),
           main = title, xlab = "tumour size at measurement", ylab = paste0("correlation coefficient:\n", col_name_suffix, " vs waiting time"))
      legend("topleft", VecParam, title = "Parameters", ncol = 3, lwd = 2, 
             xpd = TRUE, horiz = FALSE, inset = c(0, 0), bty = "n", lty = 1, col = cols, cex=0.8)
    }
    df_filtered <- filter(df, SimulationType == Param_val)
    lines(df_filtered[[col_name]] ~ df_filtered$start_size, 
          col = cols[j], lwd = 2)
    abline(h = 0, untf = FALSE, lty = 3)
  }
  
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
}

#add a title to the default function
plot_corr_waiting_time_versus_start_size_Title <- function(df, col_name = "DriverDiversity", output_filename = "corr_waiting_time_versus_start_size", file_type = "png", output_dir = NA, title="") {
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir, output_filename, ".png"), width = 500, height = 500, res = 100)
    else pdf(paste0(output_dir, output_filename, ".pdf"), width = 5, height = 5)
  }
  
  if(substr(col_name, 1, 4) != "Cor_") col_name <- paste0("Cor_", col_name)
  col_name_suffix <- substr(col_name, 5, nchar(col_name))
  
  start_size_range <- unique(df$start_size)
  K_range <- unique(df$K)
  
  cols <- rainbow(length(K_range))
  
  par(mfrow=c(1, 1))
  par(mar=c(4, 4, 2, 2))
  for(j in 1:length(K_range)) {
    K_val <- K_range[j]
    if(j == 1) {
      plot(1, type = "n", xlim = c(0, max(start_size_range)), ylim = c(-1, 1),
           main = title, xlab = "tumour size at measurement", ylab = paste0("correlation coefficient:\n", col_name_suffix, " vs waiting time"))
      legend("topleft", as.character(K_range), title = "K", ncol = 3, lwd = 2, 
             xpd = TRUE, horiz = FALSE, inset = c(0, 0), bty = "n", lty = 1, col = cols)
    }
    df_filtered <- filter(df, K == K_val)
    lines(df_filtered[[col_name]] ~ df_filtered$start_size, 
          col = cols[j], lwd = 2)
    abline(h = 0, untf = FALSE, lty = 3)
  }
  
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
}

# plot 3 graphs
plot_corr_waiting_time_versus_start_size_3Graphs <- function(df, col_name = "DriverDiversity", output_filename = "corr_waiting_time_versus_start_size", file_type = "png", output_dir = NA, title="", ParamToVary=NA) {
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir, output_filename, ".png"), width = 500, height = 500, res = 100)
    else pdf(paste0(output_dir, output_filename, ".pdf"), width = 5, height = 5)
  }
  
  if(substr(col_name, 1, 4) != "Cor_") col_name <- paste0("Cor_", col_name)
  col_name_suffix <- substr(col_name, 5, nchar(col_name))
  
  start_size_range <- unique(df$start_size)
  K_range <- unique(df$K)
  
  cols <- rainbow(length(K_range))
  if(is.na(ParamToVary)){
    par(mfrow=c(1, 1))
    par(mar=c(4, 4, 2, 2))
    for(j in 1:length(K_range)) {
      K_val <- K_range[j]
      if(j == 1) {
        plot(1, type = "n", xlim = c(0, max(start_size_range)), ylim = c(-1, 1),
             main = title, xlab = "tumour size at measurement", ylab = paste0("correlation coefficient:\n", col_name_suffix, " vs waiting time"))
        legend("topleft", as.character(K_range), title = "K", ncol = 3, lwd = 2, 
               xpd = TRUE, horiz = FALSE, inset = c(0, 0), bty = "n", lty = 1, col = cols)
      }
      df_filtered <- filter(df, K == K_val)
      lines(df_filtered[[col_name]] ~ df_filtered$start_size, 
            col = cols[j], lwd = 2)
      abline(h = 0, untf = FALSE, lty = 3)
    }
  }else{
    
    #ParamToVary_range<-unique(df[, which(colnames(df)==ParamToVary)])
    
    ParamToVary_range<-sort(unique(df[[ParamToVary]]))
    
    par(mfrow=c(1,length(ParamToVary_range)))
    
    for(PIndex in seq_along(ParamToVary_range)){
      
      #par(mar=c(4, 4, 2, 2))
      for(j in 1:length(K_range)) {
        K_val <- K_range[j]
        if(j == 1) {
          
          plot(1, type = "n", xlim = c(0, max(start_size_range)), ylim = c(-1, 1),
               main = paste0(ParamToVary_range[PIndex]), xlab = "tumour size at measurement")
          #ylab = paste0("correlation coefficient:\n", col_name_suffix, " vs waiting time")
          legend("topleft", as.character(K_range), title = "K", ncol = 3, lwd = 2,
                 xpd = TRUE, horiz = FALSE, inset = c(0, 0), bty = "n", lty = 1, col = cols)
        }
        
        df_filtered <- filter(df, K == K_val )
        
        colnames(df_filtered)[which(colnames(df_filtered)==ParamToVary)]<-"ParamToVary"
        df_filtered<-filter(df_filtered,ParamToVary== ParamToVary_range[PIndex])
        
        lines(df_filtered[[col_name]] ~ df_filtered$start_size, 
              col = cols[j], lwd = 2)
        abline(h = 0, untf = FALSE, lty = 3)
      }
      
    }
  }
  
  
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
}



get_summary_modified <- function(data, start_size_range, gap_range, final_size, num_parameters) {
  summary <- data.frame()
  data <- data %>% group_by_at(1:num_parameters)
  start_size_index=0
  for(start_size in start_size_range) {
    start_size_index=start_size_index+1
    gap_index=0
    for(gap in gap_range) {
      gap_index=gap_index+1
      if(start_size < final_size) {
        new_summary1 <- data %>% 
          filter(NumCells >= start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(start_time = min(gen_adj, na.rm = TRUE))
        new_summary2 <- data %>% 
          filter(NumCells >= final_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(end_time = min(gen_adj, na.rm = TRUE))
        new_summary12 <- merge(new_summary1, new_summary2, all.x = TRUE)
        new_summary12 <- new_summary12 %>% 
          mutate(waiting_time = end_time - start_time)
      }
      else {
        new_summary12 <- data %>% 
          filter(NumCells >= start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(waiting_time = NA, start_time = NA, end_time = NA)
      }
      new_summary3 <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        filter(gen_adj < min(gen_adj, na.rm = TRUE) + gap) %>% 
        summarise(outcome = max(NumCells, na.rm = TRUE))
      new_summary3a <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        summarise(outcome = max(NumCells, na.rm = TRUE))
      new_summary3$outcome <- ifelse(new_summary3$outcome == new_summary3a$outcome, NA, new_summary3$outcome)
      new_summary4 <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        filter(gen_adj == min(gen_adj, na.rm = TRUE)) %>%
        mutate(gap = gap, start_size = start_size)
      summary <- rbind(summary, merge(merge(new_summary12, new_summary3, all.x = TRUE), new_summary4, all.x = TRUE))
    
      print(paste("start size :" , start_size, " gap :", gap, " dim(summary)[1]: ",dim(summary)[1] ))
      print(paste("expected dim(summary)[1]:", ((100*1200)*(start_size_index-1) + (1200* gap_index)) ))
      }
  }
  summary <- summary %>% 
    mutate(Ratio = DriverEdgeDiversity / DriverDiversity)
  
  # check number of rows:
  count1 <- dim(summary)[1]
  count2 <- sum(count_seeds(summary, num_parameters)) * length(start_size_range) * length(gap_range)
  if(count1 != count2) stop(paste0("Row count (", count1, ") is not as expected (", count2, ")."))
  
  # report number of replicates per parameter set:
  print("Number of seeds:", quote = FALSE)
  print(count_seeds(summary, num_parameters))
  
  return(summary)
}


########################## Setting the folders ########################## 
subfolder_name <- "Batch_ForecastingPaper_bis"

n_cores <- 4


input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files

all_statuses(input_dir, summary = TRUE) # all should be "Exit code 0"

## First create the plot and data directories

ifelse(!dir.exists("plots"), dir.create("plots"), FALSE)
ifelse(!dir.exists("data"), dir.create("data"), FALSE)

## then create the subfolders (no recursive folder creation)

ifelse(!dir.exists(output_dir_plots), dir.create(output_dir_plots), FALSE)

ifelse(!dir.exists(output_dir_data), dir.create(output_dir_data), FALSE)



############## ############## ############## ############## 
##############        forecasting plots      ############## 
############## ############## ############## ############## 
############## We need first to create the data  ############## 
#as they are quite big we need to create separate output files



############## Manual binding ##############
include_diversities = TRUE
numcells = NA
vaf_cut_off = NA
generation = NA
max_generation = FALSE
df_type = "output"
each_df <- function(x, res) {
  full_dir <- make_dir(input_dir, pars, x)
  msg <- final_error_message(full_dir)
  print(paste0(full_dir, " ", msg), quote = FALSE)
  if(!identical(msg, character(0))) if(msg == "Exit code 0") return(combine_dfs(full_dir, include_diversities, 
                                                                                df_type, max_generation, vaf_cut_off, generation, numcells, num_parameters))
  return(data.frame())
}

df_type_list <- c("output", "driver_genotype_properties", "genotype_properties", 
                  "allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts", "driver_phylo",
                  "diversities")
stopifnot(df_type %in% df_type_list)

pars_and_values <- parameter_names_and_values(input_dir)
if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
pars <- pars_and_values$name
final_values <- pars_and_values$final_value

num_parameters <- count_parameters(input_dir)

vecs <- mapply(seq, 0,final_values, SIMPLIFY = FALSE)
tmp <- do.call(expand.grid, vecs)
tmp<-tmp[with(tmp, order(Var1, Var2, Var3, Var4)),]

for(Kval in unique(tmp$Var1)){
  
  print(paste("Saving data for K=", Kval))
  tmp_Kval<-tmp[which(tmp$Var1==Kval ), ]
  results<-c()
  
  for (i in seq_along(tmp_Kval[, 1])){
    results<-rbind(results, each_df(tmp_Kval[i,], res))
  }
  
  fwrite(results, file = paste0(output_dir_data, "/datatestK", Kval, ".csv"), row.names = FALSE)
  

}


## Manual binding for driver_genotype_properties
df_type = "driver_genotype_properties"

each_df <- function(x, res) {
  full_dir <- make_dir(input_dir, pars, x)
  msg <- final_error_message(full_dir)
  print(paste0(full_dir, " ", msg), quote = FALSE)
  if(!identical(msg, character(0))) if(msg == "Exit code 0") return(combine_dfs(full_dir, include_diversities, 
                                                                                df_type, max_generation, vaf_cut_off, generation, numcells, num_parameters))
  return(data.frame())
}
df_type_list <- c("output", "driver_genotype_properties", "genotype_properties", 
                  "allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts", "driver_phylo",
                  "diversities")
stopifnot(df_type %in% df_type_list)

pars_and_values <- parameter_names_and_values(input_dir)
if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
pars <- pars_and_values$name
final_values <- pars_and_values$final_value

num_parameters <- count_parameters(input_dir)

vecs <- mapply(seq, 0,final_values, SIMPLIFY = FALSE)
tmp <- do.call(expand.grid, vecs)
tmp<-tmp[with(tmp, order(Var1, Var2, Var3, Var4)),]

for(Kval in unique(tmp$Var1)){
  
  print(paste("Saving data for K=", Kval))
  tmp_Kval<-tmp[which(tmp$Var1==Kval ), ]
  results<-c()
  
  for (i in seq_along(tmp_Kval[, 1])){
    results<-rbind(results, each_df(tmp_Kval[i,], res))
  }
  
  fwrite(results, file = paste0(output_dir_data, "/driver_genotype_propertiesK", Kval, ".csv"), row.names = FALSE)
  
  
}

############## If already created just need to read the data from an existing csv file ############## 
data<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  data[[Kiter]] <- read_csv(paste0(output_dir_data, "/datatestK", Kval, ".csv"), guess_max = 1E4)
}


driver_genotype_properties<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  driver_genotype_properties[[Kiter]] <- read_csv(paste0(output_dir_data, "/driver_genotype_propertiesK", Kval, ".csv"), guess_max = 1E4)
}



############## Define the simulation number ############## 

# SimulationsPath<-tmp
# colnames(SimulationsPath)<-c("K","mu_driver_birth","s_driver_birth", "seed")
# SimulationsPath<-as.data.frame(SimulationsPath[ with(SimulationsPath, order(K, mu_driver_birth,s_driver_birth,seed)),])
# SimulationsPath<-mutate(SimulationsPath, SimulationNb =c(1:length(SimulationsPath[,1])))

# SimulationsDefinition<-SimulationsPath
# SimulationsDefinition$K<-ifelse(SimulationsDefinition$K==0, 64, 
#                                 ifelse(SimulationsDefinition$K==1, 512, 4096))
# 
# SimulationsDefinition$mu_driver_birth<-ifelse(SimulationsDefinition$mu_driver_birth==0,1e-6, 
#                                               ifelse(SimulationsDefinitio$mu_driver_birth==1, 1e-5, 1e-4))
# 
# SimulationsDefinition$s_driver_birth<-ifelse(SimulationsDefinition$s_driver_birth==0,0.05, 
#                                               ifelse(SimulationsDefinition$s_driver_birth==1, 0.1,
#                                                      ifelse(SimulationsDefinition$s_driver_birth==2, 0.15, 0.2)))


SimulationsDefinition<-c()
for(Kiter in seq_along(unique(tmp$Var1))){
  SimulationsDefinition<-rbind(SimulationsDefinition,as.data.frame(unique(data[[Kiter]][, c("seed","K","s_driver_birth","mu_driver_birth")])) )
  #SimulationsDefinition<-rbind(SimulationsDefinition,as.data.frame(unique(driver_genotype_properties[[Kiter]][, c("seed","K","s_driver_birth","mu_driver_birth")])) )
  
  }
SimulationsDefinition<-as.data.frame(SimulationsDefinition[ with(SimulationsDefinition, order(K, mu_driver_birth,s_driver_birth,seed)),])
SimulationsDefinition<-mutate(SimulationsDefinition, SimulationNb =c(1:length(SimulationsDefinition[,1])))


#define the " type " of the simulation, st the 100 seeds run with the same parameters belong to the same type

SimulationsDefinition<-mutate(SimulationsDefinition,
                              SimulationType=0)
SimulationTypeInit=0
K_init=SimulationsDefinition$K[1]
s_driver_birth_init=SimulationsDefinition$s_driver_birth[1]
mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[1]
iter=1
type=0
while(iter < dim(SimulationsDefinition)[1]){
  while( (SimulationsDefinition$K[iter] == K_init) &
         (SimulationsDefinition$s_driver_birth[iter]==s_driver_birth_init) &
         (SimulationsDefinition$mu_driver_birth[iter]==mu_driver_birth_init) & 
         (iter < dim(SimulationsDefinition)[1] )
  ){
    SimulationsDefinition$SimulationType[iter]=type
    iter=iter+1
    
  }
  type=type+1
  K_init=SimulationsDefinition$K[iter+1]
  s_driver_birth_init=SimulationsDefinition$s_driver_birth[iter+1]
  mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[iter+1]
  
}

fwrite(SimulationsDefinition, file = paste0(output_dir_data, "/SimulationsDefinition.csv"), row.names = FALSE)
print("Wrote SimulationsDefinition.csv to file")

SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)



############## modification of data ############## 
# add_columns and add_relative_time : see datacollation.R
# for(Kiter in seq_along(unique(tmp$Var1))){
for(Kiter in c(2:3)){
  
  data[[Kiter]] <- add_columns(data[[Kiter]], num_parameters = num_parameters)
  data[[Kiter]] <- add_relative_time(data[[Kiter]], start_size = 50000, num_parameters = num_parameters)
  
  data[[Kiter]]<-mutate(data[[Kiter]],
               SimulationNb =0,
               SimulationType=0)
  
  
  
  for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
    
    rowIndex<-which( (data[[Kiter]]$K== SimulationsDefinition$K[i]) &
                       (data[[Kiter]]$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (data[[Kiter]]$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
                       (data[[Kiter]]$seed==SimulationsDefinition$seed[i]))
    data[[Kiter]][rowIndex, 
         which(colnames(data[[Kiter]])=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
    data[[Kiter]][rowIndex, 
         which(colnames(data[[Kiter]])=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    
  }
  
  #as it takes a while, save the modification
  Kval<-unique(tmp$Var1)[Kiter]
  
  fwrite(data[[Kiter]], file = paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), row.names = FALSE)
  print(paste0("Wrote ", output_dir_data, "/dataAugmentedK", Kval, ".csv to file"))
  
  
  
}

############## If already created just need to read the data from an existing csv file ############## 
data<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  data[[Kiter]] <- as.data.frame(read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), guess_max = 1E4))
}


############## Some required variables ############## 

min_count <- 20

#column of the seed 
seed_index <- which(colnames(data) == "seed")
# column of tall other paramaters
pars_without_seed <- (1:num_parameters)[which((1:num_parameters) != seed_index)]



############## Plot waiting time correlation for each K separately ############## 

start_size_range <- c(125000,250000,375000,500000, 750000) # NumCells at time of initial measurement for forecasting
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value

output_dir_plots_WaitCorrelation<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation), dir.create(output_dir_plots_WaitCorrelation), FALSE)


SimulationsDefinition<-as.data.frame(SimulationsDefinition)

for(Kiter in c(1:3)){
  
  
  #for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){

# # for(i in seq_along(unique(SimulationsDefinition$SimulationNb))){
#   
#   SimulationNbOfInterest<-unique(SimulationsDefinition$SimulationNb)[i]
#   print(paste("Simulation Nb : ", SimulationNbOfInterest))
#   
#   KOfInterest<-SimulationsDefinition[which(SimulationsDefinition$SimulationNb==SimulationNbOfInterest), 
#                                              which(colnames(SimulationsDefinition)=="K")]
#   
#   MuOfInterest<-SimulationsDefinition[which(SimulationsDefinition$SimulationNb==SimulationNbOfInterest), 
#                                      which(colnames(SimulationsDefinition)=="mu_driver_birth")]
#   SOfInterest<-SimulationsDefinition[which(SimulationsDefinition$SimulationNb==SimulationNbOfInterest), 
#                                       which(colnames(SimulationsDefinition)=="s_driver_birth")]
#   
#   #First filter the data such that only the 100 seeds of interest remain
#   FilteredData<- data[[Kiter]] %>% subset(SimulationNb  %in% SimulationNbOfInterest)
  
  FilteredData<- data[[Kiter]]
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  
  depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))
  
  
  # the following line generates summary dataframe of correlations with "outcome"
  cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
  
  #the following lines generates summary dataframe of correlations with "waiting_time"
  wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                           num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                                 num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
  
  cor_summary<-mutate(cor_summary, SimulationType=0)
  wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)
  
  for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
    
    rowIndex<-which( (cor_summary$K== SimulationsDefinition$K[i]) &
                       (cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    
    cor_summary[rowIndex, 
                  which(colnames(cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
                       (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    wait_cor_summary[rowIndex, 
                which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    
  }
  
  
  
  col=cols_list[[1]]
  Kval<-unique(FilteredData$K)
  plot_corr_waiting_time_versus_start_size_ModifColor(wait_cor_summary, 
                                           output_filename = paste0(col, "_wait_correlations_K",Kval ), output_dir = output_dir_plots_WaitCorrelation, title=paste0("K=", Kval))

  #}
   
}

#special case of Kiter=3
start_size_range <- c(125000,250000,375000,500000, 750000)
Kiter=3
FilteredData<- data[[Kiter]] 
col_names<- colnames(FilteredData)[1:18]
#FilteredData<-FilteredData %>% group_by(col_names) %>% filter(max(NumCells) >= max(start_size_range)) %>% ungroup()
#SimulationsToRemove<-FilteredData%>% filter(K==4096, mu_driver_birth==1e-6, s_driver_birth==0.05, seed %in% c(52, 1038, 7024, 38009))

FilteredData<-FilteredData %>%filter( ! (K==4096 & mu_driver_birth==1e-6 & s_driver_birth==0.05 & seed %in% c(52, 1038, 7024, 38009)))

summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

###search for bug
test<-data %>% filter(NumCells >=50000, !is.na(DriverDiversity),
                      !is.na(NumClones)) %>% summarise(end_time = min(gen_adj,
                                                                      na.rm = TRUE))
new_summary1<-data %>% filter(NumCells >=750000, !is.na(DriverDiversity),
                              !is.na(NumClones)) %>% summarise(end_time = min(gen_adj,
                                                                              na.rm = TRUE))

as.data.frame(setdiff(test[, !colnames(test) %in% c("start_time", "end_time")], new_summary1[, !colnames(new_summary1) %in% c("start_time", "end_time")]))

Pbdataset<-data %>% filter(K==4096, mu_driver_birth==1e-6, s_driver_birth==0.05, seed %in% c(52, 1038, 7024, 38009))
Pbdataset<-Pbdataset %>% filter(NumCells >= 750000)


############ Combine correlation df ############
Allcor_summary<-vector(mode="list", length=3)
All_wait_cor_summary<-vector(mode="list", length=3)
Allsummary<-vector(mode="list", length=3)
Alldepth_wait_cor_summary<-vector(mode="list", length=3)

for(Kiter in c(1:3)){
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- data[[Kiter]]
  
  if(Kiter==3){
    #those simulations don't achieve a final size larger thant 750000
    FilteredData<-FilteredData %>%filter( ! (K==4096 & mu_driver_birth==1e-6 & s_driver_birth==0.05 & seed %in% c(52, 1038, 7024, 38009)))
  
  }
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  
  depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))
  
  
  # the following line generates summary dataframe of correlations with "outcome"
  cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
  
  #the following lines generates summary dataframe of correlations with "waiting_time"
  wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                           num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                                 num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
  
  cor_summary<-mutate(cor_summary, SimulationType=0)
  wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)
  
  for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
    
    rowIndex<-which( (cor_summary$K== SimulationsDefinition$K[i]) &
                       (cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    
    cor_summary[rowIndex, 
                which(colnames(cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
                       (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    wait_cor_summary[rowIndex, 
                     which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    
  }
  
  Allcor_summary[[Kiter]]<-cor_summary
  All_wait_cor_summary[[Kiter]]<-wait_cor_summary
  Allsummary[[Kiter]]<-summary
  Alldepth_wait_cor_summary[[Kiter]]<-Alldepth_wait_cor_summary
}


Allcor_summary<-do.call(rbind, Allcor_summary)
All_wait_cor_summary<-do.call(rbind,All_wait_cor_summary)
Allsummary<-do.call(rbind,Allsummary)
Alldepth_wait_cor_summary<-do.call(rbind,Alldepth_wait_cor_summary)

fwrite(Allcor_summary, file = paste0(output_dir_data, "/Allcor_summary.csv"), row.names = FALSE)
fwrite(All_wait_cor_summary, file = paste0(output_dir_data, "/All_wait_cor_summary.csv"), row.names = FALSE)
fwrite(Allsummary, file = paste0(output_dir_data, "/Allsummary.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary.csv"), row.names = FALSE)


print(paste0("Wrote ", output_dir_data, "/driver_genotype_propertiesAugmentedK", Kval, ".csv to file"))


Allcor_summary<-read_csv(paste0(output_dir_data, "/Allcor_summary.csv"), guess_max = 1E4)
All_wait_cor_summary<-read_csv(paste0(output_dir_data, "/All_wait_cor_summary.csv"), guess_max = 1E4)
Allsummary<-read_csv(paste0(output_dir_data, "/Allsummary.csv"), guess_max = 1E4)

# Allcor_summary<-rbind(Allcor_summary, cor_summary)
# All_wait_cor_summary<-rbind(All_wait_cor_summary, wait_cor_summary)
# Allsummary<-rbind(Allsummary, summary )




########### waiting time correlation ########### 
#Let's redo the plots of waiting time correlation this time fixing mu and s but varying K
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)

col=cols_list[[1]]

for(MuIndex in seq_along(MuValues)){
  for(SIndex in seq_along(SValues)){
    tmpSummary<-All_wait_cor_summary[which(All_wait_cor_summary$mu_driver_birth==MuValues[MuIndex] & All_wait_cor_summary$s_driver_birth==SValues[SIndex]),]
    
    plot_corr_waiting_time_versus_start_size_Title(tmpSummary, 
                                                        output_filename = paste0(col, "_wait_correlations_Mu",MuValues[MuIndex], "_S",SValues[SIndex]  ), output_dir = output_dir_plots_WaitCorrelation,
                                                        title=paste0("Mu : ",MuValues[MuIndex], ", S :",SValues[SIndex]))
    
  }
}

for(SIndex in seq_along(SValues)){
  
  tmpSummary<-All_wait_cor_summary[which(All_wait_cor_summary$s_driver_birth==SValues[SIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(K), group = K ))+
    facet_grid(. ~ mu_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity vs waiting time"))+
    ggtitle(paste0("fitness effect = ",SValues[SIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))+
    scale_color_discrete(name="K")
  
  png(paste0(output_dir_plots_WaitCorrelation, paste0("/", col, "_wait_correlations_S",SValues[SIndex]  ), "ggplot.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
}

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary[which(All_wait_cor_summary$mu_driver_birth==MuValues[MuIndex]),]
  
    gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(K), group = K ))+
    facet_grid(. ~ s_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))+
    scale_color_discrete(name="K")
  
  png(paste0(output_dir_plots_WaitCorrelation, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "ggplot.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}


########### waiting time correlation with other start size ###########

start_size_range <- c(500,1000,2000,5000*2^(0:4))

Allcor_summary_smallerStartSize<-vector(mode="list", length=3)
All_wait_cor_summary_smallerStartSize<-vector(mode="list", length=3)
Allsummary_smallerStartSize<-vector(mode="list", length=3)
Alldepth_wait_cor_summary_smallerStartSize<-vector(mode="list", length=3)

for(Kiter in c(1:3)){
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- data[[Kiter]]
  
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  
  depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                       paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))
  
  
  # the following line generates summary dataframe of correlations with "outcome"
  cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
  
  #the following lines generates summary dataframe of correlations with "waiting_time"
  wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                           num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                                 num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
  
  cor_summary<-mutate(cor_summary, SimulationType=0)
  wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)
  
  for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
    
    rowIndex<-which( (cor_summary$K== SimulationsDefinition$K[i]) &
                       (cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    
    cor_summary[rowIndex, 
                which(colnames(cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
                       (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
    wait_cor_summary[rowIndex, 
                     which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    
  }
  
  Allcor_summary_smallerStartSize[[Kiter]]<-cor_summary
  All_wait_cor_summary_smallerStartSize[[Kiter]]<-wait_cor_summary
  Allsummary_smallerStartSize[[Kiter]]<-summary
  Alldepth_wait_cor_summary_smallerStartSize[[Kiter]]<-Alldepth_wait_cor_summary
}


Allcor_summary_smallerStartSize<-do.call(rbind, Allcor_summary_smallerStartSize)
All_wait_cor_summary_smallerStartSize<-do.call(rbind,All_wait_cor_summary_smallerStartSize)
Allsummary_smallerStartSize<-do.call(rbind,Allsummary_smallerStartSize)
Alldepth_wait_cor_summary_smallerStartSize<-do.call(rbind,Alldepth_wait_cor_summary_smallerStartSize)

fwrite(Allcor_summary_smallerStartSize, file = paste0(output_dir_data, "/Allcor_summary_smallerStartSize.csv"), row.names = FALSE)
fwrite(All_wait_cor_summary_smallerStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_smallerStartSize.csv"), row.names = FALSE)
fwrite(Allsummary_smallerStartSize, file = paste0(output_dir_data, "/Allsummary_smallerStartSize.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary_smallerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_smallerStartSize.csv"), row.names = FALSE)






############Combine driver_genotype_properties ############



df <- fread(paste0(mypath, name, "/data_driver_genotype_properties.csv"))
# par_names <- colnames(df)[1:23]
par_names <- colnames(df)[1:18]
sum_df <- group_by_at(df, par_names) %>% mutate(Diversity = inv_Simpson_index(Population/sum(Population)),
                                                Diversity2 = inv_Simpson_index2(Population, sum(Population)),
                                                Diversity_Descendants = inv_Simpson_index(Descendants/sum(Descendants))) %>%
  slice(1) %>% ungroup()









########################## combine outputs ########################## 

for(Kiter in c(2:3)){
  
  driver_genotype_properties[[Kiter]]<-mutate(driver_genotype_properties[[Kiter]],
                                              SimulationNb =0,
                                              SimulationType=0)
  
  
  
  for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
    
    rowIndex<-which( (driver_genotype_properties[[Kiter]]$K== SimulationsDefinition$K[i]) &
                       (driver_genotype_properties[[Kiter]]$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                       (driver_genotype_properties[[Kiter]]$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
                       (driver_genotype_properties[[Kiter]]$seed==SimulationsDefinition$seed[i]))
    driver_genotype_properties[[Kiter]][rowIndex, 
                                        which(colnames(driver_genotype_properties[[Kiter]])=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
    driver_genotype_properties[[Kiter]][rowIndex, 
                                        which(colnames(driver_genotype_properties[[Kiter]])=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
    
    
  }
  
  par_names <- colnames(driver_genotype_properties[[Kiter]])[1:18]
  driver_genotype_properties[[Kiter]] <- group_by_at(driver_genotype_properties[[Kiter]], par_names) %>% mutate(Diversity = inv_Simpson_index(Population/sum(Population)),
                                                  Diversity2 = inv_Simpson_index2(Population, sum(Population)),
                                                  Diversity_Descendants = inv_Simpson_index(Descendants/sum(Descendants))) %>%
    slice(1) %>% ungroup()
  
  
  #as it takes a while, save the modification
  Kval<-unique(tmp$Var1)[Kiter]
  
  fwrite(driver_genotype_properties[[Kiter]], file = paste0(output_dir_data, "/driver_genotype_propertiesAugmentedK", Kval, ".csv"), row.names = FALSE)
  print(paste0("Wrote ", output_dir_data, "/driver_genotype_propertiesAugmentedK", Kval, ".csv to file"))
  
  
  
}

############## If already created just need to read the data from an existing csv file ############## 
data<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  data[[Kiter]] <- as.data.frame(read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), guess_max = 1E4))
}



######################################## ALL THE FOLLOWING IS NOT USED  YET ######################################## 
AllSimulationNbOfInterest<-SimulationsDefinition[which( (SimulationsDefinition$K==K_OfInterest) &
                                                          (SimulationsDefinition$s_driver_birth==s_driver_birth_OfInterest) &
                                                          (SimulationsDefinition$mu_driver_birth==mu_driver_birth_OfInterest)), 
                                                 which(colnames(SimulationsDefinition)=="SimulationNb")]

FilteredData<- data %>% subset(SimulationNb  %in% AllSimulationNbOfInterest)

summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

cols_list <- c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))

#check that all required columns beginning by e.g "DriverDiversityFrom1SamplesAtDepth"are selected
#colnames(data)[grep("DriverDiversityFrom1SamplesAtDepth", colnames(data))]

# the following line generates summary dataframe of correlations with "outcome"
cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)

#the following lines generates summary dataframe of correlations with "waiting_time"
wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                         num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells

depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                               num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols

plot_corr_waiting_time_versus_start_size(wait_cor_summary, 
                                         output_filename = paste0(col, "_wait_correlations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)
                                         





## combine output_driver_genotype_propoerties
CombinedDriverGenotypeProp<-all_output(input_dir, include_diversities = FALSE, df_type = "driver_genotype_properties",
                                       max_generation = FALSE, vaf_cut_off = NA, generation = NA,
                                       numcells = NA, n_cores = n_cores)

edgesCombinedDriverGenotypeProp <- filter(CombinedDriverGenotypeProp, Descendants >= 1e6 / 100) %>% select(Parent, Identity, Population, seed,K, s_driver_birth,mu_driver_birth)


edgesCombinedDriverGenotypeProp_df<-data.frame(edgesCombinedDriverGenotypeProp)


##########Give a number to each simulation 
# SimulationsDefinition<-expand.grid(unique(edgesCombinedDriverGenotypeProp_df$K), 
#                                    unique(edgesCombinedDriverGenotypeProp_df$s_driver_birth),
#                                    unique(edgesCombinedDriverGenotypeProp_df$mu_driver_birth))
# 
# colnames(SimulationsDefinition)<-c("K", "s_driver_birth","mu_driver_birth")
# 
# SimulationsDefinition<-mutate(SimulationsDefinition, Simulation=rownames(SimulationsDefinition))
# edgesCombinedDriverGenotypeProp_df<-mutate(edgesCombinedDriverGenotypeProp_df,
#                                            SimulationNb =0)
# for(i in c(1:dim(SimulationsDefinition)[1]) ){
#   
#   edgesCombinedDriverGenotypeProp_df[which( (edgesCombinedDriverGenotypeProp_df$K== SimulationsDefinition$K[i]) &
#                                               (edgesCombinedDriverGenotypeProp_df$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
#                                               (edgesCombinedDriverGenotypeProp_df$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
#                                               (edgesCombinedDriverGenotypeProp_df$seed==SimulationsDefinition$seed[i])), 
#                                      which(colnames(edgesCombinedDriverGenotypeProp_df)=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
#   
#   
# }
# 
# SeedSimulationsDefinition<-expand.grid(SimulationsDefinition$Simulation, unique(edgesCombinedDriverGenotypeProp_df$seed))
# 
# colnames(SeedSimulationsDefinition)<-c("Simulation", "seed")
# 
# SeedSimulationsDefinition<-mutate(SeedSimulationsDefinition, SeedSimulationNb=rownames(SeedSimulationsDefinition))
# 
# 
# edgesCombinedDriverGenotypeProp_df<-mutate(edgesCombinedDriverGenotypeProp_df,
#                                            SeedSimulationNb =0)
# 
# for(i in c(1:dim(SeedSimulationsDefinition)[1]) ){
#   
#   edgesCombinedDriverGenotypeProp_df[which( (edgesCombinedDriverGenotypeProp_df$SimulationNb== SeedSimulationsDefinition$Simulation[i]) &
#                                               (edgesCombinedDriverGenotypeProp_df$seed== SeedSimulationsDefinition$seed[i]) ), 
#                                      which(colnames(edgesCombinedDriverGenotypeProp_df)=="SeedSimulationNb")]<-SeedSimulationsDefinition$SeedSimulationNb[i]
# }



# CompleteSimulationsDefinition<-merge(SeedSimulationsDefinition,SimulationsDefinition, by="Simulation", all.x=TRUE)


SimulationsDefinition<- as.data.frame(unique(CombinedDriverGenotypeProp[, which(colnames(CombinedDriverGenotypeProp) %in% c("seed","K","s_driver_birth","mu_driver_birth"))]))
SimulationsDefinition<-as.data.frame(SimulationsDefinition[ with(SimulationsDefinition, order(K, mu_driver_birth,s_driver_birth,seed)),])
SimulationsDefinition<-mutate(SimulationsDefinition, SimulationNb =c(1:length(SimulationsDefinition[,1])))


#define the " type " of the simulation, st the 100 seeds run with the same parameters belong to the same type

SimulationsDefinition<-mutate(SimulationsDefinition,
                              SimulationType=0)
SimulationTypeInit=0
K_init=SimulationsDefinition$K[1]
s_driver_birth_init=SimulationsDefinition$s_driver_birth[1]
mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[1]
iter=1
type=0
while(iter < dim(SimulationsDefinition)[1]){
  while( (SimulationsDefinition$K[iter] == K_init) &
         (SimulationsDefinition$s_driver_birth[iter]==s_driver_birth_init) &
         (SimulationsDefinition$mu_driver_birth[iter]==mu_driver_birth_init)
  ){
    SimulationsDefinition$SimulationType[iter]=type
    iter=iter+1
    
  }
  type=type+1
  K_init=SimulationsDefinition$K[iter+1]
  s_driver_birth_init=SimulationsDefinition$s_driver_birth[iter+1]
  mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[iter+1]
  
}








edgesCombinedDriverGenotypeProp_df<-mutate(edgesCombinedDriverGenotypeProp_df,
                                           SimulationNb =0)
for(i in c(1:dim(SimulationsDefinition)[1]) ){
  
  edgesCombinedDriverGenotypeProp_df[which( (edgesCombinedDriverGenotypeProp_df$K== SimulationsDefinition$K[i]) &
                                              (edgesCombinedDriverGenotypeProp_df$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                                              (edgesCombinedDriverGenotypeProp_df$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
                                              (edgesCombinedDriverGenotypeProp_df$seed==SimulationsDefinition$seed[i])), 
                                     which(colnames(edgesCombinedDriverGenotypeProp_df)=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
  
  
}

# create a list of the metrics for each seeds
# MetricCombinedDriverGenotypeProp<-by(edgesCombinedDriverGenotypeProp_df,
#                                      INDICES=edgesCombinedDriverGenotypeProp_df$SeedSimulationNb,
#                                      FUN=metrics)
MetricCombinedDriverGenotypeProp<-by(edgesCombinedDriverGenotypeProp_df,
                                     INDICES=edgesCombinedDriverGenotypeProp_df$SimulationNb,
                                     FUN=metrics)


MetricCombinedDriverGenotypeProp_df<-as.data.frame(do.call(rbind, MetricCombinedDriverGenotypeProp))

MetricCombinedDriverGenotypeProp_df<-mutate(MetricCombinedDriverGenotypeProp_df, 
                                            SeedSimulationNb= rownames(MetricCombinedDriverGenotypeProp_df))

#### add the caracteristic of the simulation

MetricCombinedDriverGenotypeProp_df<-mutate(MetricCombinedDriverGenotypeProp_df,
                                            K=0,
                                            s_driver_birth=0,
                                            mu_driver_birth=0
)

for(i in seq_along(unique(MetricCombinedDriverGenotypeProp_df$SeedSimulationNb))){
  
  SeedSimulation<-MetricCombinedDriverGenotypeProp_df$SeedSimulationNb[i]
  #get the associated simulation number
  SimulatioNb<-SeedSimulationsDefinition[which(SeedSimulationsDefinition$SeedSimulationNb==SeedSimulation), which(colnames(SeedSimulationsDefinition)=="Simulation")]
  rowInSimulationsDefinition<-which(SimulationsDefinition$Simulation==SimulatioNb)
  
  MetricCombinedDriverGenotypeProp_df$K[i]=SimulationsDefinition$K[rowInSimulationsDefinition]  
  MetricCombinedDriverGenotypeProp_df$s_driver_birth[i]=SimulationsDefinition$s_driver_birth[rowInSimulationsDefinition]
  MetricCombinedDriverGenotypeProp_df$mu_driver_birth[i]=SimulationsDefinition$mu_driver_birth[rowInSimulationsDefinition]
  
  
}


########################## Metrics plot ########################## 

line_df <- data.frame(x = seq(0.01, 20, length = 100))
line_df$y1 <- 2^(1/3*(line_df$x-1))
line_df$y2 <- 2^(line_df$x-1)
line_df$y3 <- 2^(3*(line_df$x-1))
line_df$y4 <- 1

# x1 and y1 => tree with two nodes
# x2 and y1 => linear tree with three nodes, with root node extinct
# x3 and y2 => branched tree with three nodes, where leaf nodes have equal size
x <- seq(0, 1, length = 500)
curve_df <- data.frame(x1 = 1 + x, 
                       x2 = 2 + x, 
                       x3 = 2 - x)
curve_df$y1 <- 1 / (x^2 + (1 - x)^2)
alpha <- 1/2
curve_df$y2 <- 1 / (x^2 + alpha^2 * (1 - x)^2 + (1 - alpha)^2 * (1 - x)^2)
r <- seq(1, 10, length = 500)
k <- 2
D_func <- function(r) {
  vec <- 0:(r-1)
  top <- sum(2^vec * (4 * (vec + 1) - 1))
  bot <- 2 * (2^r - 1)
  return(top / bot)
}
my_func <- function(n) 1 / ((n %% 1 - 1)^2 + (n %% 1)^2)
z_func <- function(n, r) {
  q <- (n - r)/(1 - r)
  return(1 / (q^2 + (1 - q)^2))
}
curve_df$v1 <- (k^r - 1) / (k - 1)
curve_df$u1 <- ((r * k^(r + 1) - r * k^r - k^r + 1) / (k - 1)^2) / curve_df$v1
curve_df$u2 <- 3 * curve_df$u1
curve_df$zmax <- 1 / (2 - curve_df$x1)^2
curve_df$zmin <- sapply(curve_df$x1, z_func, r = 2)
curve_df$linex <- seq(1, 50, length = 500)
curve_df$liney <- 1
curve_df$bumps <- sapply(curve_df$linex, my_func)
# alpha <- 0
# curve_df$x4 <- x + (3 - alpha) * (1 - x)
# curve_df$y4 <- 1 / (x^2 + alpha^2 * (1 - x)^2 + (1 - alpha)^2 * (1 - x)^2)
# alpha <- 1
# curve_df$x5 <- x + (3 - alpha) * (1 - x)
# curve_df$y5 <- 1 / (x^2 + alpha^2 * (1 - x)^2 + (1 - alpha)^2 * (1 - x)^2)

g1 <- ggplot() +
  geom_point(data=MetricCombinedDriverGenotypeProp_df,
             aes(x=n, y=D), alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


########### adding color in functions of the parameters

g1 <- ggplot(data=MetricCombinedDriverGenotypeProp_df) +
  geom_point(aes(x=n, y=D, color=as.factor(K)),  alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_colorByK.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


g1 <- ggplot(data=MetricCombinedDriverGenotypeProp_df) +
  geom_point(aes(x=n, y=D, color=as.factor(s_driver_birth)),  alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_colorByFitness.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot(data=MetricCombinedDriverGenotypeProp_df) +
  geom_point(aes(x=n, y=D, color=as.factor(mu_driver_birth)),  alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_colorByDriverMutRate.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


########## Now let's fix K and color by fitness effect
for( Kind in seq_along(unique(MetricCombinedDriverGenotypeProp_df$K)) ){
  
  KValue=unique(MetricCombinedDriverGenotypeProp_df$K)[Kind]
  
  g1 <- ggplot(data=subset(MetricCombinedDriverGenotypeProp_df,
                           MetricCombinedDriverGenotypeProp_df$K==KValue)
  ) +
    geom_point(aes(x=n, y=D, color=as.factor(s_driver_birth)),  alpha=0.5) +
    scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
    #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
    scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
    theme_classic() +
    # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
    geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
    geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
    geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey") +
    ggtitle(paste("Effect of fitness effect for deme size ", KValue))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots, "/MetricPlot_K",KValue,  "colorByFitnessEffect.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
  
}

########## Now let's fix mu_driver_birth (which seems to have the largest influence) and color by fitness effect and then by K

for( Mu_ind in seq_along(unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)) ){
  
  MuValue=unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)[Mu_ind]
  
  g1 <- ggplot(data=subset(MetricCombinedDriverGenotypeProp_df,
                           MetricCombinedDriverGenotypeProp_df$mu_driver_birth==MuValue)
  ) +
    geom_point(aes(x=n, y=D, color=as.factor(s_driver_birth)),  alpha=0.5) +
    scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
    #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
    scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
    theme_classic() +
    # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
    geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
    geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
    geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey") +
    ggtitle(paste("Effect of fitness effect for fixed mu_driver_birth = ", MuValue))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots, "/MetricPlot_mu",MuValue,  "colorByFitnessEffect.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
  
}



for( Mu_ind in seq_along(unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)) ){
  
  MuValue=unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)[Mu_ind]
  
  g1 <- ggplot(data=subset(MetricCombinedDriverGenotypeProp_df,
                           (MetricCombinedDriverGenotypeProp_df$mu_driver_birth==MuValue) & !(MetricCombinedDriverGenotypeProp_df$s_driver_birth==0.14))
  ) +
    geom_point(aes(x=n, y=D, color=as.factor(s_driver_birth)),  alpha=0.5) +
    scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
    #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
    scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
    theme_classic() +
    # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
    geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
    geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
    geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey") +
    ggtitle(paste("Effect of fitness effect for fixed mu_driver_birth = ", MuValue))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots, "/MetricPlot_mu",MuValue,  "colorByTWOFitnessEffect.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
  
}

for( Mu_ind in seq_along(unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)) ){
  
  MuValue=unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)[Mu_ind]
  
  g1 <- ggplot(data=subset(MetricCombinedDriverGenotypeProp_df,
                           MetricCombinedDriverGenotypeProp_df$mu_driver_birth==MuValue)
  ) +
    geom_point(aes(x=n, y=D, color=as.factor(K)),  alpha=0.5) +
    scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
    #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
    scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
    theme_classic() +
    # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
    geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
    geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
    geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey") +
    ggtitle(paste("Effect of deme size for fixed mu_driver_birth = ", MuValue))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots, "/MetricPlot_mu",MuValue,  "colorByDemeSize.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
  
}

for( Mu_ind in seq_along(unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)) ){
  
  MuValue=unique(MetricCombinedDriverGenotypeProp_df$mu_driver_birth)[Mu_ind]
  
  g1 <- ggplot(data=subset(MetricCombinedDriverGenotypeProp_df,
                           (MetricCombinedDriverGenotypeProp_df$mu_driver_birth==MuValue) & ! (MetricCombinedDriverGenotypeProp_df$K==128) )) +
    geom_point(aes(x=n, y=D, color=as.factor(K)),  alpha=0.5) +
    scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
    #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
    scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
    theme_classic() +
    # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
    # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
    #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
    geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
    geom_line(aes(x = linex, y = liney), data=curve_df, lty = 2, color = "black") +
    geom_line(aes(x = linex, y = bumps), data=curve_df, lty = 2, color = "grey") +
    ggtitle(paste("Effect of deme size for fixed mu_driver_birth = ", MuValue))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots, "/MetricPlot_mu",MuValue,  "colorByTWODemeSize.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
  
}

### Seach for parameters more responsible for points in the lower part of the graph

ParamOfInterest<-MetricCombinedDriverGenotypeProp_df %>% filter(D<=2)


#test the selection
g1 <- ggplot() +
  geom_point(data=ParamOfInterest,
             aes(x=n, y=D), alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_ForDSmallerThan2.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

###
g1 <- ggplot() +
  geom_point(data=ParamOfInterest,
             aes(x=n, y=D,color=as.factor(K)), alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_ForDSmallerThan2_ColorByK.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot() +
  geom_point(data=ParamOfInterest,
             aes(x=n, y=D,color=as.factor(mu_driver_birth)), alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_ForDSmallerThan2_ColorByMu.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot() +
  geom_point(data=ParamOfInterest,
             aes(x=n, y=D,color=as.factor(s_driver_birth)), alpha=0.5) +
  scale_y_log10(limits = c(1, 40), name = "Diversity of driver mutations") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 14), name = "Mean driver mutations per cell") +
  theme_classic() +
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  geom_line(data=curve_df, aes(x = x1, y = zmax), lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey")


png(paste0(output_dir_plots, "/MetricPlot_ForDSmallerThan2_ColorByFitness.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

###
# would need to k^look at the percentage compared to the number of seed with each of these paramters.
g1<-ggplot(ParamOfInterest)+
  geom_histogram(aes(x=log2(K)), binwidth = 1)

png(paste0(output_dir_plots, "/HistogramKValues_ForDSmallerThan2.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


############## ############## ############## ############## 
##############        forecasting plots      ############## 
############## ############## ############## ############## 
############## We need first to create the data  ############## 

## combined data for a batch of simulations, including diversity columns
data<- all_output(input_dir, include_diversities = TRUE, max_generation = FALSE, n_cores = n_cores) 
fwrite(data, file = paste0(output_dir_data, "/data.csv"), row.names = FALSE)
print("Wrote data.csv to file")


############## If already created just need to read the data from an existing csv file ############## 

data <- read_csv(paste0(output_dir_data, "/data.csv"), guess_max = 1E4)


############## Define the simulation number ############## 

SimulationsDefinition<- as.data.frame(unique(data[, c("seed","K","s_driver_birth","mu_driver_birth")]))
SimulationsDefinition<-as.data.frame(SimulationsDefinition[ with(SimulationsDefinition, order(K, mu_driver_birth,s_driver_birth,seed)),])
SimulationsDefinition<-mutate(SimulationsDefinition, SimulationNb =c(1:length(SimulationsDefinition[,1])))


#define the " type " of the simulation, st the 100 seeds run with the same parameters belong to the same type

SimulationsDefinition<-mutate(SimulationsDefinition,
                              SimulationType=0)
SimulationTypeInit=0
K_init=SimulationsDefinition$K[1]
s_driver_birth_init=SimulationsDefinition$s_driver_birth[1]
mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[1]
iter=1
type=0
while(iter < dim(SimulationsDefinition)[1]){
  while( (SimulationsDefinition$K[iter] == K_init) &
         (SimulationsDefinition$s_driver_birth[iter]==s_driver_birth_init) &
         (SimulationsDefinition$mu_driver_birth[iter]==mu_driver_birth_init) & 
         (iter < dim(SimulationsDefinition)[1] )
  ){
    SimulationsDefinition$SimulationType[iter]=type
    iter=iter+1
    
  }
  type=type+1
  K_init=SimulationsDefinition$K[iter+1]
  s_driver_birth_init=SimulationsDefinition$s_driver_birth[iter+1]
  mu_driver_birth_init=SimulationsDefinition$mu_driver_birth[iter+1]
  
}

fwrite(SimulationsDefinition, file = paste0(output_dir_data, "/SimulationsDefinition.csv"), row.names = FALSE)
print("Wrote SimulationsDefinition.csv to file")

SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)


############## Some required variables ############## 

min_count <- 5

#column of the seed 
seed_index <- which(colnames(data) == "seed")
# column of tall other paramaters
pars_without_seed <- (1:num_parameters)[which((1:num_parameters) != seed_index)]


############## modification of data ############## 
# add_columns and add_relative_time : see datacollation.R
data <- add_columns(data, num_parameters = num_parameters)
data <- add_relative_time(data, start_size = 50000, num_parameters = num_parameters)

data<-mutate(data,
             SimulationNb =0,
             SimulationType=0)

for(i in c(1:dim(SimulationsDefinition)[1]) ){
  
  rowIndex<-which( (data$K== SimulationsDefinition$K[i]) &
                     (data$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                     (data$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
                     (data$seed==SimulationsDefinition$seed[i]))
  data[rowIndex, 
       which(colnames(data)=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
  data[rowIndex, 
       which(colnames(data)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
  
  
}

#as it takes a while, save the modification
fwrite(data, file = paste0(output_dir_data, "/dataAugmented.csv"), row.names = FALSE)
print("Wrote dataAugmented.csv to file")


data <- read_csv(paste0(output_dir_data,"/dataAugmented.csv"), guess_max = 1E4)
data<-as.data.frame(data)


############## Some variables ############## 
start_size_range <- c(500,1000,2000,5000*2^(0:4)) # NumCells at time of initial measurement for forecasting
gap_range <- (1:10)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value

############## As a first test of the functions, we filter the data to have only some seeds, as the summary function is quite long ############## 

# to filter the data, we select 10 random seeds from each of the 27 different simulations.
NbRandomSeeds =10

set.seed(42)

SelectedSimulations_df<-by(SimulationsDefinition, INDICES=SimulationsDefinition$SimulationType, FUN=sample_n, NbRandomSeeds)

SelectedSimulations_df<-as.data.frame(do.call(rbind,SelectedSimulations_df))

SelectedSimulations<-unique(SelectedSimulations_df$SimulationNb)


FilteredData<- data %>% subset(SimulationNb  %in% SelectedSimulations)

############## Get summary of the data ############## 
# summary data for each simulation, for each combination of gap and final_size
summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
write.csv(summary, paste0(output_dir_data, "/summaryFilteredData.csv"), row.names = FALSE)

# summary data for each simulation, for each combination of gap and final_size
summary100K <- get_summary(FilteredData, start_size_range, gap_range, 1e5, num_parameters = num_parameters)
write.csv(summary100K, paste0(output_dir_data, "/summary100KFilteredData.csv"), row.names = FALSE)

# if need to read the summary
summary <- read_csv(paste0(output_dir_data, "/summaryFilteredData.csv"), guess_max = 1E4)


############## Correlation computation ############## 

# cols_list <- c("DriverDiversity", "DriverEdgeDiversity", "quad_div", "MeanBirthRate", "Drivers")
cols_list <- c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")


depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))

#check that all required columns beginning by e.g "DriverDiversityFrom1SamplesAtDepth"are selected
#colnames(data)[grep("DriverDiversityFrom1SamplesAtDepth", colnames(data))]

# the following line generates summary dataframe of correlations with "outcome"
cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)

#the following lines generates summary dataframe of correlations with "waiting_time"
wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                         num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells

depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                               num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols


###### Plots###### 
plot_curves_faceted(data, num_parameters, x_var = "Generation", y_var = "NumCells", 
                    output_filename = "curves", output_dir = output_dir_plots)

plot_curves_Simplified_faceted(data, num_parameters, x_var = "Generation", y_var = "NumCells", 
                               output_filename = "curves_Simplified_faceted", output_dir = output_dir_plots, paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))


plot_curves_Simplified_faceted(data, num_parameters, x_var = "new_time", y_var = "NumCells", 
                               output_filename = "trajectories_Simplified_faceted", 
                               output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
                               paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))

plot_curves_Simplified_faceted(data, num_parameters, x_var = "new_time", y_var = "NumCells", 
                               output_filename = "trajectories_Simplified_faceted", 
                               output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
                               paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))

plot_curves_Simplified_faceted(data, num_parameters, x_var = "Generation", y_var = "NumCells", 
                               output_filename = "trajectories_Simplified_faceted_Generation", 
                               output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
                               paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))


plot_curves_Simplified_faceted(data, num_parameters, x_var = "Generation", y_var = "NumCells", 
                               output_filename = "trajectories_Simplified_faceted_Generation_test", 
                               output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
                               paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))




i=0
j=0
plot_trajectories_by_diversity_bySimulationType(data, file_type = "png", 
                                                output_filename = paste0("trajectories_by_diversity_type", i, "_edge", j, "_SimulationType"), output_dir = output_dir_plots)

## search for the pb
testData<-FilteredData[which(FilteredData$SimulationType==0),]

plot_curves_Simplified_faceted(testData, num_parameters, x_var = "Generation", y_var = "NumCells", 
                               output_filename = "trajectories_Simplified_faceted_Generation_test", 
                               output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
                               paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))

plot_curves_faceted(data, num_parameters,x_var = "Generation", y_var = "NumCells", 
                    output_filename ="trajectories_Simplified_faceted_Generation_test", 
                    output_dir = output_dir_plots, line_col = "div0", alpha = 0.5)


######## plot simulation per simulation
SimulationTypeName<-sort(unique(FilteredData$SimulationType))

for(SimIter in seq_along(SimulationTypeName)) {
  
  SimType<-SimulationTypeName[SimIter]
  
  #tmpData<-FilteredData[which(FilteredData$SimulationType==SimType),]
  tmpData<-data[which(data$SimulationType==SimType),]
  
  
  # plot_curves_Simplified_faceted(tmpData, num_parameters, x_var = "Generation", y_var = "NumCells", 
  #                                output_filename = paste0("trajectories_Simplified_faceted_Generation_SimType",SimType), 
  #                                output_dir = output_dir_plots, line_col = "div0", alpha = 0.5,
  #                                paramOfinterest =c("seed","K","s_driver_birth","mu_driver_birth"))
  
  plot_curves_faceted(tmpData, num_parameters,x_var = "Generation", y_var = "DriverDiversity", 
                      output_filename =paste0("DriverDiversityVsGeneration_SimType",SimType), 
                      output_dir = output_dir_plots, line_col = "div0", alpha = 0.5)
  
}


#From the previous sets of plots, it seems that the following sets of parameters leads to very distinguishable curves depending on thhe diversity
# s=0.2, mu= 10^-6 K=128
# s=0.2, mu= 9*10^-6 K=256

#The following sets also lead to acceptable predictable situation
# s=0.2, mu= 9*10^-6 K=128
# s=0.1, mu= 9*10^-6 K=128

#Let's look at the scatter plots for those sets

###################
# First investigated set : s=0.2, mu= 9*10^-6 K=256
s_driver_birth_OfInterest<-0.2
mu_driver_birth_OfInterest<-0.000009
K_OfInterest<-256
###################

# ###################
# # Second investigated set : s=0.2, mu= 9*10^-6 K=64
# s_driver_birth_OfInterest<-0.2
# mu_driver_birth_OfInterest<-0.000009
# K_OfInterest<-64
# ###################

###################
# # Third investigated set : s=0.15, mu= 5*10^-6 K=64
# s_driver_birth_OfInterest<-0.15
# mu_driver_birth_OfInterest<-0.000005
# K_OfInterest<-64
###################

# ###################
# # 4th investigated set : s=0.15, mu= 10^-6 K=256
s_driver_birth_OfInterest<-0.15
mu_driver_birth_OfInterest<-0.000001
K_OfInterest<-64
# ###################


start_size_range <- c(125000,250000,375000,500000, 750000) # NumCells at time of initial measurement for forecasting
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value


SimulationsDefinition<-as.data.frame(SimulationsDefinition)
AllSimulationNbOfInterest<-SimulationsDefinition[which( (SimulationsDefinition$K==K_OfInterest) &
                                                          (SimulationsDefinition$s_driver_birth==s_driver_birth_OfInterest) &
                                                          (SimulationsDefinition$mu_driver_birth==mu_driver_birth_OfInterest)), 
                                                 which(colnames(SimulationsDefinition)=="SimulationNb")]

FilteredData<- data %>% subset(SimulationNb  %in% AllSimulationNbOfInterest)

summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 

cols_list <- c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")


depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                     paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))

#check that all required columns beginning by e.g "DriverDiversityFrom1SamplesAtDepth"are selected
#colnames(data)[grep("DriverDiversityFrom1SamplesAtDepth", colnames(data))]

# the following line generates summary dataframe of correlations with "outcome"
cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)

#the following lines generates summary dataframe of correlations with "waiting_time"
wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                         num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells

depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                               num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols

####plots
#define directory
output_dir_plots_SpecialCase<-paste0(output_dir_plots, "/s",s_driver_birth_OfInterest, "_mu", mu_driver_birth_OfInterest, "_K", K_OfInterest )
ifelse(!dir.exists(output_dir_plots_SpecialCase), dir.create(output_dir_plots_SpecialCase), FALSE)

#plot correlations
for(col in cols_list) {
  # plot_corr_outcome_versus_period(cor_summary, col_name = col, 
  #                                 output_filename = paste0(col, "_correlations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)
  # plot_corr_waiting_time_versus_start_size(wait_cor_summary, 
  #                                          output_filename = paste0(col, "_wait_correlations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)
  plot_corr_outcome_versus_period(cor_summary, col_name = col, 
                                  output_filename = paste0(col, "_LongRangecorrelations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)
  plot_corr_waiting_time_versus_start_size(wait_cor_summary, col_name = col, 
                                           output_filename = paste0(col, "_wait_LongRangecorrelations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)
  
}

plot_corr_outcome_versus_period(cor_summary, col_name = cols_list[1], 
                                output_filename = paste0(cols_list[1], "_correlations_type0_edge0"), output_dir = output_dir_plots_SpecialCase)

#in the present case, we only have migration_type==0 and migration_edge_only==0 and K=128

i=0
j=0
plot_corr_waiting_time_versus_depth(filter(depth_wait_cor_summary, migration_type == i, migration_edge_only == j, K == K_OfInterest), "DriverDiversityFrom1SamplesAtDepth", 
                                    output_filename = paste0("depth_wait_correlations_K", K_OfInterest, "_type", i, "_edge", j), output_dir = output_dir_plots_SpecialCase)
plot_corr_waiting_time_versus_depth(filter(depth_wait_cor_summary, migration_type == i, migration_edge_only == j, K == K_OfInterest), "DriverDiversityFrom4BigSamplesAtDepth", 
                                    output_filename = paste0("4BigSamples_depth_wait_correlations_K", K_OfInterest, "_type", i, "_edge", j), output_dir = output_dir_plots_SpecialCase)

plot_curves_faceted(FilteredData, num_parameters,x_var = "Generation", y_var = "NumCells", 
                    output_filename ="trajectories_faceted_Generation", 
                    output_dir = output_dir_plots_SpecialCase, line_col = "div0", alpha = 0.5)
plot_curves_faceted(FilteredData, num_parameters,x_var = "new_time", y_var = "NumCells", 
                    output_filename ="trajectories_faceted_new_time", 
                    output_dir = output_dir_plots_SpecialCase, line_col = "div0", alpha = 0.5)

# # try the filtered version : first filtered the data for each simulation to select only driver genotypes that accounted for at least 1% of the total population.
# 
# CombinedDriverGenotypeProp<-all_output(input_dir, include_diversities = FALSE, df_type = "driver_genotype_properties",
#                                        max_generation = FALSE, vaf_cut_off = NA, generation = NA,
#                                        numcells = NA, n_cores = n_cores)
# 
# CombinedDriverGenotypeProp<-filter(CombinedDriverGenotypeProp, Population > 1e6 / 100 )
# CombinedDriverGenotypeProp<-mutate(CombinedDriverGenotypeProp,
#                                    Diversity=nv_Simpson_index(Population))


# Similar plot as the paper : evolution of driver diversity with the number of generation
#First only select 20 seed per simulations type 
NbRandomSeeds =20
SimulationTypeName<-sort(unique(data$SimulationType))
for(SimIter in seq_along(SimulationTypeName)) {
  
  SimType<-SimulationTypeName[SimIter]
  
  set.seed(42)
  tmpSimulationsDefinition<-SimulationsDefinition[which(SimulationsDefinition$SimulationType==SimType), ]
  
  tmpSimulationsDefinition<-tmpSimulationsDefinition[sample(nrow(tmpSimulationsDefinition), NbRandomSeeds), ]
  
  SelectedSeed_df<-tmpSimulationsDefinition$seed
  
  tmpData<-data[which(data$SimulationType==SimType & data$seed %in% SelectedSeed_df ),]
  
  
  plot_curves_faceted(tmpData, num_parameters,x_var = "Generation", y_var = "DriverDiversity",
                      output_filename =paste0("DriverDiversityVsGeneration_20Seeds_SimType",SimType),
                      output_dir = output_dir_plots, alpha = 0.5)
  
  plot_curves_faceted(tmpData, num_parameters,x_var = "NumCells", y_var = "DriverDiversity", 
                      output_filename =paste0("DriverDiversityVsNumCells_20Seeds_SimType",SimType), 
                      output_dir = output_dir_plots, alpha = 0.5)
  
}


s_driver_birth_OfInterest<-0.15
mu_driver_birth_OfInterest<-0.000001
K_OfInterest<-64
SimulationsDefinition[which( (SimulationsDefinition$K==K_OfInterest) &
                               (SimulationsDefinition$s_driver_birth==s_driver_birth_OfInterest) &
                               (SimulationsDefinition$mu_driver_birth==mu_driver_birth_OfInterest))[1],]


cor_summary[which(cor_summary$start_size==125000), which(colnames(cor_summary) %in% c( "Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))]
cor_summary[which(cor_summary$start_size==250000), which(colnames(cor_summary) %in% c( "Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))]
cor_summary[which(cor_summary$start_size==375000), which(colnames(cor_summary) %in% c( "Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))]
cor_summary[which(cor_summary$start_size==500000), which(colnames(cor_summary) %in% c( "Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))]
cor_summary[which(cor_summary$start_size==750000), which(colnames(cor_summary) %in% c( "Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))]





NbRandomSeeds =20
SimulationTypeName<-sort(unique(data$SimulationType))

ListPlots_DriverDiversityVsGeneration_20Seeds<-vector("list", length(SimulationTypeName))
ListPlots_DriverDiversityVsNumCells_20Seeds<-vector("list", length(SimulationTypeName))


for(SimIter in seq_along(SimulationTypeName)) {
  
  SimType<-SimulationTypeName[SimIter]
  
  set.seed(42)
  tmpSimulationsDefinition<-SimulationsDefinition[which(SimulationsDefinition$SimulationType==SimType), ]
  
  tmpSimulationsDefinition<-tmpSimulationsDefinition[sample(nrow(tmpSimulationsDefinition), NbRandomSeeds), ]
  
  SelectedSeed_df<-tmpSimulationsDefinition$seed
  
  tmpData<-data[which(data$SimulationType==SimType & data$seed %in% SelectedSeed_df ),]
  
  
  ListPlots_DriverDiversityVsGeneration_20Seeds[[SimIter]]<-get_curves_faceted(tmpData, num_parameters,x_var = "Generation", y_var = "DriverDiversity",
                                                                               output_filename =paste0("DriverDiversityVsGeneration_20Seeds_SimType",SimType),
                                                                               output_dir = output_dir_plots, alpha = 0.5)
  
  ListPlots_DriverDiversityVsNumCells_20Seeds[[SimIter]]<-get_curves_faceted(tmpData, num_parameters,x_var = "NumCells", y_var = "DriverDiversity", 
                                                                             output_filename =paste0("DriverDiversityVsNumCells_20Seeds_SimType",SimType), 
                                                                             output_dir = output_dir_plots, alpha = 0.5)
  
}

png(paste0(output_dir_plots, "/AllDriverDiversityVsGeneration_20Seeds", ".png"), width = 1000, height = 1000, res = 100)

plot_grid(plotlist =ListPlots_DriverDiversityVsGeneration_20Seeds,
          nrow=3,
          ncol=9)

dev.off()



################# Plots correlations to waiting time for all simulations  ################# 


s_driver_birth_vector<-c(0.1,0.15, 0.2)
mu_driver_birth_vector<-c(0.000001, 0.000005, 0.000009)
K_vector<-2^c(6,7,8)

start_size_range <- c(125000,250000,375000,500000, 750000) # NumCells at time of initial measurement for forecasting
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
min_count<-20


png(paste0(output_dir_plots, "/AllDriverDiversityCorrelationWaitingTime", ".png"), width = 1000, height = 1000, res = 100)

par(mfrow=c(3,9))

for(K_OfInterest in K_vector ){
  for(mu_driver_birth_OfInterest in mu_driver_birth_vector){
    for(s_driver_birth_OfInterest in s_driver_birth_vector){
      
      SimulationsDefinition<-as.data.frame(SimulationsDefinition)
      AllSimulationNbOfInterest<-SimulationsDefinition[which( (SimulationsDefinition$K==K_OfInterest) &
                                                                (SimulationsDefinition$s_driver_birth==s_driver_birth_OfInterest) &
                                                                (SimulationsDefinition$mu_driver_birth==mu_driver_birth_OfInterest)), 
                                                       which(colnames(SimulationsDefinition)=="SimulationNb")]
      
      FilteredData<- data %>% subset(SimulationNb  %in% AllSimulationNbOfInterest)
      
      summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
      
      cols_list <- c("DriverDiversity")
      
      
      depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
                           paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
                           paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
                           paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))
      
      #check that all required columns beginning by e.g "DriverDiversityFrom1SamplesAtDepth"are selected
      #colnames(data)[grep("DriverDiversityFrom1SamplesAtDepth", colnames(data))]
      
      # the following line generates summary dataframe of correlations with "outcome"
      cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
      
      #the following lines generates summary dataframe of correlations with "waiting_time"
      wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
                                               num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
      
      depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
                                                     num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
      
      get_corr_waiting_time_versus_start_size(wait_cor_summary, col_name = cols_list[1], 
                                              output_filename = paste0(cols_list[1], "_wait_LongRangecorrelations_type0_edge0"), output_dir = output_dir_plots)
      
      
    }
  }
}

dev.off()




