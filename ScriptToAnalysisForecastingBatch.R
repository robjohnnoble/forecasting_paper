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
library(ggmuller)
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

corCI_low<-function(Var1, Var2){
  # cor.result<-cor.test(Var1,Var2,  method = "spearman", alternative = "two.sided", exact=FALSE)
  # if(is.null(cor.result$conf.int)){
  #   return(NA)
  # }else{
  #   return(cor.result$conf.int[1])
  # }
  
  cor.result<-spearman.ci(Var1, Var2, nrep = 500, conf.level = 0.95)
  if(is.null(cor.result$conf.int)){
    return(NA)
  }else{
    return(cor.result$conf.int[1])
  }

}

corCI_high<-function(Var1, Var2){
  #cor.result<-cor.test(Var1,Var2,  method = "spearman", alternative = "two.sided", exact=FALSE)
  cor.result<-spearman.ci(Var1, Var2, nrep = 500, conf.level = 0.95)
  if(is.null(cor.result$conf.int)){
    return(NA)
  }else{
    return(cor.result$conf.int[2])    
  }
  
}

cor_pval<-function(Var1, Var2){
  cor.result<-cor.test(Var1,Var2,  method = "spearman", alternative = "two.sided", exact=FALSE)
  if(is.null(cor.result$p.value)){
    return(NA)
  }else{
    return(cor.result$p.value)    
  }
  
}

find_correlations_modified <- function(summary, factor1, factor2, result_name, min_count) {
  summary %>% 
    mutate_(variance = interp(~var(var1, na.rm = TRUE), var1 = as.name(factor2))) %>% 
    filter(variance > 0) %>% # to avoid warnings when all values of factor2 are identical
    mutate_(count1 = interp(~length(var1), var1 = as.name(factor1)), 
            count2 = interp(~length(var2), var2 = as.name(factor2))) %>% 
    filter(count1 >= min_count, count2 >= min_count) %>% 
    summarise_(temp_name = interp(~cor(var1, var2, method = "spearman", use="na.or.complete"), var1 = as.name(factor1), var2 = as.name(factor2)),
               temp_name_pval = interp(~cor_pval(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2))) %>%
               #temp_name_cilow = interp(~corCI_low(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2)),
               #temp_name_cihigh = interp(~corCI_high(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2))) %>%
  #rename_(.dots = setNames(c("temp_name", "temp_name_pval","temp_name_cilow", "temp_name_cihigh") , paste0(result_name,  c("", "_pVal", "_CI_low", "_CI_high")) ))
  rename_(.dots = setNames(c("temp_name", "temp_name_pval") , paste0(result_name,  c("", "_pVal")) ))
  
}

get_wait_cor_summary_modified <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[col_nums != which(colnames(summary) == "seed")]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(waiting_time)) %>% 
    filter(var(waiting_time) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)){
    print(col_names_list[i])
    cor_summary_list[[i]] <- find_correlations_modified(summary, "waiting_time", col_names_list[i], result_names_list[i], min_count)
    
  }
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}


get_wait_cor_summary_CombinedMutationRate <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[which(! col_nums %in% (which(colnames(summary) %in% c("seed", "mu_driver_birth"))))]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(waiting_time)) %>% 
    filter(var(waiting_time) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations_modified(summary, "waiting_time", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}


get_wait_cor_summary_CombinedFitnessEffect <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[which(! col_nums %in% (which(colnames(summary) %in% c("seed", "s_driver_birth"))))]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(waiting_time)) %>% 
    filter(var(waiting_time) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations_modified(summary, "waiting_time", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}

#get the correlation with the AverageGrowthRate variable
get_AverageGrowthRate_cor_summary_modified <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[col_nums != which(colnames(summary) == "seed")]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(AverageGrowthRate)) %>% 
    filter(!is.infinite(AverageGrowthRate)) %>% 
    filter(var(AverageGrowthRate) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_AverageGrowthRate=mean(AverageGrowthRate),
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)){
    print(col_names_list[i])
    cor_summary_list[[i]] <- find_correlations_modified(summary, "AverageGrowthRate", col_names_list[i], result_names_list[i], min_count)
    
  }
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}




get_AverageGrowthRate_cor_summary_CombinedMutationRate <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[which(! col_nums %in% (which(colnames(summary) %in% c("seed", "mu_driver_birth"))))]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(AverageGrowthRate)) %>% 
    filter(!is.infinite(AverageGrowthRate)) %>% 
    filter(var(AverageGrowthRate) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_AverageGrowthRate=mean(AverageGrowthRate),
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations_modified(summary, "AverageGrowthRate", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}


get_AverageGrowthRate_cor_summary_CombinedFitnessEffect <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[which(! col_nums %in% (which(colnames(summary) %in% c("seed", "s_driver_birth"))))]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(AverageGrowthRate)) %>% 
    filter(!is.infinite(AverageGrowthRate)) %>% 
    filter(var(AverageGrowthRate) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_AverageGrowthRate=mean(AverageGrowthRate), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations_modified(summary, "AverageGrowthRate", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}


get_AverageGrowthRate_cor_summary_CombinedMutationRate_VaryingFinalSize <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "FinalSize"))
  col_nums <- col_nums[which(! col_nums %in% (which(colnames(summary) %in% c("seed", "mu_driver_birth"))))]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(AverageGrowthRate)) %>% 
    filter(!is.infinite(AverageGrowthRate)) %>% 
    filter(var(AverageGrowthRate) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_AverageGrowthRate=mean(AverageGrowthRate),
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations_modified(summary, "AverageGrowthRate", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, FinalSize)
  
  return(cor_summary)
}


get_MeanBirthRate_correlation <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[col_nums != which(colnames(summary) == "seed")]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(MeanBirthRate)) %>% 
    filter(var(MeanBirthRate) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_MeanBirthRate = mean(MeanBirthRate), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)){
    print(col_names_list[i])
    cor_summary_list[[i]] <- find_correlations_modified(summary, "MeanBirthRate", col_names_list[i], result_names_list[i], min_count)
    
  }
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}

########################## Setting the folders ########################## 
subfolder_name <- "Batch_ForecastingPaper_bis"

n_cores <- 4


input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files

##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)

########################## SimulationsDefinition##########################

SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)

SimulationsDefinition<-as.data.frame(SimulationsDefinition)


########### All summary ########
#summary contains the variable waiting_time
Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)

########################## Usefull for density plots ##########################
colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))



############## ############## ############## ############## 
##############        forecasting plots      ############## 
############## ############## ############## ############## 
############## We need first to create the data  ############## 
#as they are quite big we need to create separate output files



##############tmp def##############
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

############## SimulationDef
SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)
SimulationsDefinition<-as.data.frame(SimulationsDefinition)

############## Augmented data
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

############ Larger start size, summary computation, Combine correlation df ############
# Allcor_summary<-vector(mode="list", length=3)
# All_wait_cor_summary<-vector(mode="list", length=3)
# Allsummary<-vector(mode="list", length=3)
# Alldepth_wait_cor_summary<-vector(mode="list", length=3)
# 
# for(Kiter in c(1:3)){
#   
#   print(paste("working on Kiter=",Kiter ))
#   
#   FilteredData<- data[[Kiter]]
#   
#   if(Kiter==3){
#     #those simulations don't achieve a final size larger thant 750000
#     FilteredData<-FilteredData %>%filter( ! (K==4096 & mu_driver_birth==1e-6 & s_driver_birth==0.05 & seed %in% c(52, 1038, 7024, 38009)))
#     
#   }
#   summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
#   
#   cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
#   
#   depth_cols_list <- c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
#                        paste0("DriverDiversityFrom4SamplesAtDepth", 0:10), 
#                        paste0("DriverDiversityFrom1BigSamplesAtDepth", 0:10), 
#                        paste0("DriverDiversityFrom4BigSamplesAtDepth", 0:10))
#   
#   
#   # the following line generates summary dataframe of correlations with "outcome"
#   cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
#   
#   #the following lines generates summary dataframe of correlations with "waiting_time"
#   wait_cor_summary <- get_wait_cor_summary(summary, cols_list, 
#                                            num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
#   
#   depth_wait_cor_summary <- get_wait_cor_summary(summary, depth_cols_list, 
#                                                  num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
#   
#   cor_summary<-mutate(cor_summary, SimulationType=0)
#   wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)
#   
#   for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
#     
#     rowIndex<-which( (cor_summary$K== SimulationsDefinition$K[i]) &
#                        (cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
#                        (cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
#     
#     cor_summary[rowIndex, 
#                 which(colnames(cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
#     
#     rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
#                        (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
#                        (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
#     wait_cor_summary[rowIndex, 
#                      which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
#     
#     
#   }
#   
#   Allcor_summary[[Kiter]]<-cor_summary
#   All_wait_cor_summary[[Kiter]]<-wait_cor_summary
#   Allsummary[[Kiter]]<-summary
#   Alldepth_wait_cor_summary[[Kiter]]<-Alldepth_wait_cor_summary
# }
# 
# 
# Allcor_summary<-do.call(rbind, Allcor_summary)
# All_wait_cor_summary<-do.call(rbind,All_wait_cor_summary)
# Allsummary<-do.call(rbind,Allsummary)
# Alldepth_wait_cor_summary<-do.call(rbind,Alldepth_wait_cor_summary)
# 
# fwrite(Allcor_summary, file = paste0(output_dir_data, "/Allcor_summary.csv"), row.names = FALSE)
# fwrite(All_wait_cor_summary, file = paste0(output_dir_data, "/All_wait_cor_summary.csv"), row.names = FALSE)
# fwrite(Allsummary, file = paste0(output_dir_data, "/Allsummary.csv"), row.names = FALSE)
# fwrite(Alldepth_wait_cor_summary, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary.csv"), row.names = FALSE)

########### Smaller start size summary computation ###########

#start_size_range <- c(500,1000,2000,5000*2^(0:4))
start_size_range <- seq(10000,85000, by=15000)

# Allcor_summary_smallerStartSize<-vector(mode="list", length=3)
# All_wait_cor_summary_smallerStartSize<-vector(mode="list", length=3)
# Allsummary_smallerStartSize<-vector(mode="list", length=3)
# Alldepth_wait_cor_summary_smallerStartSize<-vector(mode="list", length=3)
# 
#for(Kiter in c(1:3)){
for(Kiter in c(2:3)){
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- data[[Kiter]]
  
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  print("summary computed")
  
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
  
  fwrite(cor_summary, file = paste0(output_dir_data, "/cor_summary_smallerStartSizeK", Kiter,  ".csv"), row.names = FALSE)
  fwrite(wait_cor_summary, file = paste0(output_dir_data, "/wait_cor_summary_smallerStartSizeK", Kiter,".csv"), row.names = FALSE)
  fwrite(summary, file = paste0(output_dir_data, "/summary_smallerStartSizeK", Kiter,".csv"), row.names = FALSE)
  fwrite(depth_wait_cor_summary, file = paste0(output_dir_data, "/depth_wait_cor_summary_smallerStartSizeK", Kiter,".csv"), row.names = FALSE)
  
  
  
  # Allcor_summary_smallerStartSize[[Kiter]]<-cor_summary
  # All_wait_cor_summary_smallerStartSize[[Kiter]]<-wait_cor_summary
  # Allsummary_smallerStartSize[[Kiter]]<-summary
  # Alldepth_wait_cor_summary_smallerStartSize[[Kiter]]<-Alldepth_wait_cor_summary
}


# Allcor_summary_smallerStartSize<-do.call(rbind, Allcor_summary_smallerStartSize)
# All_wait_cor_summary_smallerStartSize<-do.call(rbind,All_wait_cor_summary_smallerStartSize)
# Allsummary_smallerStartSize<-do.call(rbind,Allsummary_smallerStartSize)
# Alldepth_wait_cor_summary_smallerStartSize<-do.call(rbind,Alldepth_wait_cor_summary_smallerStartSize)
# 
# fwrite(Allcor_summary_smallerStartSize, file = paste0(output_dir_data, "/Allcor_summary_smallerStartSize.csv"), row.names = FALSE)
# fwrite(All_wait_cor_summary_smallerStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_smallerStartSize.csv"), row.names = FALSE)
# fwrite(Allsummary_smallerStartSize, file = paste0(output_dir_data, "/Allsummary_smallerStartSize.csv"), row.names = FALSE)
# fwrite(Alldepth_wait_cor_summary_smallerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_smallerStartSize.csv"), row.names = FALSE)



#### Work On plot_corr_waiting_time_versus_start_size_ModifColor
wait_cor_summary_smallerStartSizeK1<-read_csv(paste0(output_dir_data, "/wait_cor_summary_smallerStartSizeK1.csv"), guess_max = 1E4)

col=cols_list[[1]]
Kval<-64
plot_corr_waiting_time_versus_start_size_ModifColor(wait_cor_summary_smallerStartSizeK1, 
                                                    output_filename = paste0(col, "_wait_correlations_K",Kval, "smallerStartSize" ), output_dir = output_dir_plots_WaitCorrelation, title=paste0("K=", Kval))


wait_cor_summary_smallerStartSizeK2<-read_csv(paste0(output_dir_data, "/wait_cor_summary_smallerStartSizeK2.csv"), guess_max = 1E4)

col=cols_list[[1]]
Kval<-512
plot_corr_waiting_time_versus_start_size_ModifColor(wait_cor_summary_smallerStartSizeK2, 
                                                    output_filename = paste0(col, "_wait_correlations_K",Kval, "smallerStartSize" ), output_dir = output_dir_plots_WaitCorrelation, title=paste0("K=", Kval))

wait_cor_summary_smallerStartSizeK3<-read_csv(paste0(output_dir_data, "/wait_cor_summary_smallerStartSizeK3.csv"), guess_max = 1E4)

col=cols_list[[1]]
Kval<-4096
plot_corr_waiting_time_versus_start_size_ModifColor(wait_cor_summary_smallerStartSizeK3, 
                                                    output_filename = paste0(col, "_wait_correlations_K",Kval, "smallerStartSize" ), output_dir = output_dir_plots_WaitCorrelation, title=paste0("K=", Kval))


# Create All_wait_cor_summary_smallerStartSize
All_wait_cor_summary_smallerStartSize<-vector(mode="list", length=3)

for(Kiter in c(1:3)){
  
  All_wait_cor_summary_smallerStartSize[[Kiter]]<-read_csv(paste0(output_dir_data, "/wait_cor_summary_smallerStartSizeK", Kiter, ".csv"), guess_max = 1E4)

}
All_wait_cor_summary_smallerStartSize<-do.call(rbind,All_wait_cor_summary_smallerStartSize)
fwrite(All_wait_cor_summary_smallerStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_smallerStartSize.csv"), row.names = FALSE)


########### Samller start size plot waiting time correlation ########### 

output_dir_plots_WaitCorrelation_SmallerStartSize<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_SmallerStartSize"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_SmallerStartSize), dir.create(output_dir_plots_WaitCorrelation_SmallerStartSize), FALSE)


#Let's do the plots of waiting time correlation, fixing mu and s but varying K
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)

cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
col=cols_list[[1]]

for(MuIndex in seq_along(MuValues)){
  for(SIndex in seq_along(SValues)){
    tmpSummary<-All_wait_cor_summary_smallerStartSize[which(All_wait_cor_summary_smallerStartSize$mu_driver_birth==MuValues[MuIndex] & All_wait_cor_summary_smallerStartSize$s_driver_birth==SValues[SIndex]),]
    
    plot_corr_waiting_time_versus_start_size_Title(tmpSummary, 
                                                   output_filename = paste0(col, "_wait_correlations_Mu",MuValues[MuIndex], "_S",SValues[SIndex], "SmallerStartSize"  ), output_dir = output_dir_plots_WaitCorrelation_SmallerStartSize,
                                                   title=paste0("Mu : ",MuValues[MuIndex], ", S :",SValues[SIndex]))
    
  }
}

for(SIndex in seq_along(SValues)){
  
  tmpSummary<-All_wait_cor_summary_smallerStartSize[which(All_wait_cor_summary_smallerStartSize$s_driver_birth==SValues[SIndex]),]
  
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
  
  png(paste0(output_dir_plots_WaitCorrelation_SmallerStartSize, paste0("/", col, "_wait_correlations_S",SValues[SIndex]  ), "SmallerStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
}

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary_smallerStartSize[which(All_wait_cor_summary_smallerStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
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
  
  png(paste0(output_dir_plots_WaitCorrelation_SmallerStartSize, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "SmallerStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}


#######  Medium Start size ##########
#Complete the hole in start_size
start_size_range <- seq(90000,100000, by=5000)
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value


for(Kiter in c(1:3)){
#for(Kiter in c(2:3)){
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- data[[Kiter]]
  
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  print("summary computed")
  
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
  
  fwrite(cor_summary, file = paste0(output_dir_data, "/cor_summary_mediumStartSizeK", Kiter,  ".csv"), row.names = FALSE)
  fwrite(wait_cor_summary, file = paste0(output_dir_data, "/wait_cor_summary_mediumStartSizeK", Kiter,".csv"), row.names = FALSE)
  fwrite(summary, file = paste0(output_dir_data, "/summary_mediumStartSizeK", Kiter,".csv"), row.names = FALSE)
  fwrite(depth_wait_cor_summary, file = paste0(output_dir_data, "/depth_wait_cor_summary_mediumStartSizeK", Kiter,".csv"), row.names = FALSE)
  
}

# Create All_wait_cor_summary_mediumStartSize
All_wait_cor_summary_mediumStartSize<-vector(mode="list", length=3)

for(Kiter in c(1:3)){
  
  All_wait_cor_summary_mediumStartSize[[Kiter]]<-read_csv(paste0(output_dir_data, "/wait_cor_summary_mediumStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  
}
All_wait_cor_summary_mediumStartSize<-do.call(rbind,All_wait_cor_summary_mediumStartSize)
fwrite(All_wait_cor_summary_mediumStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_mediumStartSize.csv"), row.names = FALSE)

#Let's put all All_wait_cor_summary together
All_wait_cor_summary_AllStartSize<-read_csv(paste0(output_dir_data, "/All_wait_cor_summary_smallerStartSize.csv"), guess_max = 1E4)
All_wait_cor_summary_AllStartSize<-rbind(All_wait_cor_summary_AllStartSize, read_csv(paste0(output_dir_data, "/All_wait_cor_summary_mediumStartSize.csv"), guess_max = 1E4))
All_wait_cor_summary_AllStartSize<-rbind(All_wait_cor_summary_AllStartSize, read_csv(paste0(output_dir_data, "/All_wait_cor_summary.csv"), guess_max = 1E4))

fwrite(All_wait_cor_summary_AllStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_AllStartSize.csv"), row.names = FALSE)


########### All start size plot waiting time correlation ########### 

output_dir_plots_WaitCorrelation_AllStartSize<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_AllStartSize"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_AllStartSize), dir.create(output_dir_plots_WaitCorrelation_AllStartSize), FALSE)


#Let's do the plots of waiting time correlation, fixing mu and s but varying K
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)

cols_list <-c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
#col=cols_list[[1]]
col=cols_list[[2]]

for(SIndex in seq_along(SValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$s_driver_birth==SValues[SIndex]),]
  
  gg<-ggplot(tmpSummary)+
    #geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(K), group = K ))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, color=as.factor(K), group = K ))+
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
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_S",SValues[SIndex]  ), "AllStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
}

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  gg<-ggplot(tmpSummary)+
    #geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(K), group = K ))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, color=as.factor(K), group = K ))+
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
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "AllStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue"))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green"))+
    scale_color_manual(name="", values=c("blue"="blue", "green"="green"), labels=c("DriverDiversity", "DriverEdgeDiversity"))+
    facet_grid( K ~ s_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))

  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "facetK_AllStartSize.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
  
}

for(SIndex in seq_along(SValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$s_driver_birth==SValues[SIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue"))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K, color="green" ))+
    scale_color_manual(name="", values=c("blue"="blue", "green"="green"), labels=c("DriverDiversity", "DriverEdgeDiversity"))+
    facet_grid(K ~ mu_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle(paste0("fitness effect = ",SValues[SIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_S",SValues[SIndex]  ), "facetK_AllStartSize.png"), width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
}


for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue"))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green"))+
    geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red"))+
    geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow"))+
    scale_color_manual(name="", values=c("blue"="blue", "green"="green", "red"="red", "yellow"="yellow"), labels=c("DriverDiversity", "DriverEdgeDiversity", "Cor_MeanBirthRate", "Cor_Drivers"))+
    facet_grid( K ~ s_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "facetK_AllStartSize_AllCor.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
  
}



for(KIndex in seq_along(KValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$K==KValues[KIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(mu_driver_birth), group = mu_driver_birth ))+
    facet_grid(. ~ s_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity vs waiting time"))+
    ggtitle(paste0("Deme size = ",KValues[KIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))+
    scale_color_discrete(name="Driver mutation rate")
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_K",KValues[KIndex]  ), "AllStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}

for(KIndex in seq_along(KValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$K==KValues[KIndex]),]
  
  gg<-ggplot(tmpSummary)+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(s_driver_birth), group = s_driver_birth ))+
    facet_grid(. ~ mu_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity vs waiting time"))+
    ggtitle(paste0("Deme size = ",KValues[KIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))+
    scale_color_discrete(name="Fitness effect")
  
  png(paste0(output_dir_plots_WaitCorrelation_AllStartSize, paste0("/", col, "_wait_correlations_K",KValues[KIndex]  ), "_facetFitness_AllStartSize.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}


#########  reduce max strat_size to 400000 ##########  

output_dir_plots_WaitCorrelation_StartSizeBis400000<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_StartSizeBis400000"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_StartSizeBis400000), dir.create(output_dir_plots_WaitCorrelation_StartSizeBis400000), FALSE)

for(SIndex in seq_along(SValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$s_driver_birth==SValues[SIndex]),]
  
  gg<-ggplot(filter(tmpSummary, 
                    start_size <= 400000))+
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
  
  png(paste0(output_dir_plots_WaitCorrelation_StartSizeBis400000, paste0("/", col, "_wait_correlations_S",SValues[SIndex]  ), "StartSizeBis400000.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
}

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  gg<-ggplot(filter(tmpSummary, 
                    start_size <= 400000))+
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
  
  png(paste0(output_dir_plots_WaitCorrelation_StartSizeBis400000, paste0("/", col, "_wait_correlations_M",MuValues[MuIndex]  ), "StartSizeBis400000.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}

for(KIndex in seq_along(KValues)){
  
  tmpSummary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$K==KValues[KIndex]),]
  
  gg<-ggplot(filter(tmpSummary, 
                    start_size <= 400000))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, color=as.factor(s_driver_birth), group = s_driver_birth ))+
    facet_grid(. ~ mu_driver_birth)+
    theme_bw()+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity vs waiting time"))+
    ggtitle(paste0("Deme size = ",KValues[KIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90))+
    scale_color_discrete(name="Fitness effect")
  
  png(paste0(output_dir_plots_WaitCorrelation_StartSizeBis400000, paste0("/", col, "_wait_correlations_K",KValues[KIndex]  ), "_facetFitness_StartSizeBis400000.png"), width = 500, height = 500, res = 100)
  print(gg)
  dev.off()
  
  
}


#plot in function of K
output_dir_plots_WaitCorrelation_K<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_K"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_K), dir.create(output_dir_plots_WaitCorrelation_K), FALSE)


col=cols_list[[1]]
for(KIndex in seq_along(KValues)){
  
  Kval<-KValues[KIndex]
  tmp_wait_cor_summary<-All_wait_cor_summary_AllStartSize[which(All_wait_cor_summary_AllStartSize$K==KValues[KIndex]),]
  
  plot_corr_waiting_time_versus_start_size_ModifColor(tmp_wait_cor_summary, 
                                                      output_filename = paste0(col, "_wait_correlations_K",Kval, "plot" ), output_dir = output_dir_plots_WaitCorrelation_K, title=paste0("K=", Kval))
  
  tmp_wait_cor_summary<-filter(tmp_wait_cor_summary, 
                               start_size <= 400000)
  plot_corr_waiting_time_versus_start_size_ModifColor(tmp_wait_cor_summary, 
                                                      output_filename = paste0(col, "_wait_correlations_K",Kval, "StartSizeBis400000_plot" ), output_dir = output_dir_plots_WaitCorrelation_K, title=paste0("K=", Kval))
  
  
}



##########  Metric plots ##########  

#data loading and saving
driver_genotype_properties<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  driver_genotype_properties[[Kiter]] <- read_csv(paste0(output_dir_data, "/driver_genotype_propertiesK", Kval, ".csv"), guess_max = 1E4)
}

Alldriver_genotype_properties<-do.call(rbind, driver_genotype_properties)

fwrite(Alldriver_genotype_properties, file = paste0(output_dir_data, "/Alldriver_genotype_properties.csv"), row.names = FALSE)

Alldriver_genotype_properties<-read_csv(paste0(output_dir_data, "/Alldriver_genotype_properties.csv"), guess_max = 1E4)

# DIversity computation

par_names <- colnames(Alldriver_genotype_properties)[1:18]
Alldriver_genotype_properties <- group_by_at(Alldriver_genotype_properties, par_names) %>% mutate(Diversity = inv_Simpson_index(Population/sum(Population)),
                                                Diversity2 = inv_Simpson_index2(Population, sum(Population)),
                                                Diversity_Descendants = inv_Simpson_index(Descendants/sum(Descendants))) %>%
  slice(1) %>% ungroup()

#

############## Augmented data
data_output<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  data_output[[Kiter]] <- read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), guess_max = 1E4)
  data_output[[Kiter]] <- data_output[[Kiter]]  %>% group_by(SimulationNb) %>%filter(Generation == max(Generation)) %>% ungroup()

  data_output[[Kiter]] <-as.data.frame(data_output[[Kiter]])
  
}
data_output<-do.call(rbind,data_output)

dataForMetricPlots<-full_join(data_output, Alldriver_genotype_properties)
# sum(is.na(dataForMetricPlots$Drivers))

# Add SimulationType
# dataForMetricPlots<-mutate(dataForMetricPlots,
#                            SimulationNb=0,
#                            SimulationType=0)
# for(i in seq_along(SimulationsDefinition[,1])){
#   
#   rowIndex<-which( (dataForMetricPlots$K== SimulationsDefinition$K[i]) &
#                      (dataForMetricPlots$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
#                      (dataForMetricPlots$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
#                      (dataForMetricPlots$seed==SimulationsDefinition$seed[i]))
#   dataForMetricPlots[rowIndex, 
#                 which(colnames(dataForMetricPlots)=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
#   dataForMetricPlots[rowIndex, 
#                 which(colnames(dataForMetricPlots)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
#   
# }

fwrite(dataForMetricPlots, file = paste0(output_dir_data, "/dataForMetricPlots.csv"), row.names = FALSE)

dataForMetricPlots<-read_csv(paste0(output_dir_data, "/dataForMetricPlots.csv"), guess_max = 1E4)



## output_dir metric 
output_dir_plots_Metric<-"plots/Batch_ForecastingPaper_bis/MetricPlots"
ifelse(!dir.exists(output_dir_plots_Metric), dir.create(output_dir_plots_Metric), FALSE)

## color
#color=interaction(K, s_driver_birth, mu_driver_birth))
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#with values=sample(color, 36))
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))

##labels
Lab<-vector(mode="character", length=length(unique(dataForMetricPlots$SimulationType)))
Param_range=unique(dataForMetricPlots$SimulationType)

  for(SimType in seq_along(Param_range)){
    
    MuVal<-dataForMetricPlots[which(dataForMetricPlots$SimulationType==Param_range[SimType])[1], which(colnames(dataForMetricPlots) %in% c("mu_driver_birth"))]
    SVal<-dataForMetricPlots[which(dataForMetricPlots$SimulationType==Param_range[SimType])[1], which(colnames(dataForMetricPlots) %in% c("s_driver_birth"))]
    KVal<-dataForMetricPlots[which(dataForMetricPlots$SimulationType==Param_range[SimType])[1], which(colnames(dataForMetricPlots) %in% c("K"))]
    
    Lab[SimType]<-paste0("K", KVal,"M", MuVal, "S", SVal )
    
  }



##Driver mutation diversity
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
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(SimulationType)) , dataForMetricPlots, alpha=0.5) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Parameters", values=(colfunc(36)), labels=Lab)
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

### facet by K
g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(SimulationType)) , dataForMetricPlots, alpha=0.5) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(.~K)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Parameters", values=(colfunc(36)), labels=Lab)
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend_facetK.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

## facet mu
g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(SimulationType)) , dataForMetricPlots, alpha=0.5) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(. ~ mu_driver_birth)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Parameters", values=(colfunc(36)), labels=Lab)
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend_facetMu.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

#facet by mu and K
g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(SimulationType)) , dataForMetricPlots, alpha=0.5) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(K ~ mu_driver_birth)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Parameters", values=(colfunc(36)), labels=Lab)
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend_facetMuAndK.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(s_driver_birth)) , dataForMetricPlots, alpha=0.3) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(K ~ mu_driver_birth)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend_facetMuAndK_colorS.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity_Descendants,color=as.factor(s_driver_birth)) , dataForMetricPlots, alpha=0.3) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(K ~ mu_driver_birth, scales="free_x")+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationDiversity_MetricPlot_Legend_facetMuAndK_colorS_FreeScale.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()





## Clonal diversity
g1 <- ggplot() +
  geom_point(aes(x=Drivers+1, y=Diversity,color=as.factor(s_driver_birth)) , dataForMetricPlots, alpha=0.3) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 40), name = "Mean driver mutations per cell") +
  facet_grid(K ~ mu_driver_birth)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/ClonalDiversity_MetricPlot_Legend_facetMuAndK_colorS.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


###  clonal versus mutation driver diversity
g1 <- ggplot() +
  geom_point(aes(x=Diversity, y=Diversity_Descendants,color=as.factor(s_driver_birth)) , dataForMetricPlots, alpha=0.3) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 70), name = "Clonal diversity") +
  facet_grid(K ~ mu_driver_birth)+
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
  geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationVsClonal_MetricPlot_facetMuAndK_colorS.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

g1 <- ggplot() +
  geom_point(aes(x=Diversity, y=Diversity_Descendants,color=as.factor(s_driver_birth)) , dataForMetricPlots, alpha=0.3) +
  scale_y_log10(limits = c(1, 70), name = "Driver mutation diversity") +
  #scale_x_continuous(limits = c(0, 13), trans = "sqrt", breaks = c(0.1, 1, 5, 10)) +
  scale_x_log10(limits = c(1, 70), name = "Clonal diversity") +
  facet_grid(K ~ mu_driver_birth)+
  geom_abline(slope=1, intercept = 0, linetype="dashed",color = "grey")+
  theme_bw()+
  # geom_line(aes(x = x, y = y1), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y2), line_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x, y = y3), line_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x2, y = y1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x3, y = y2), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u1, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = u2, y = v1), curve_df, lty = 2, color = "grey") +
  #geom_line(aes(x = x1, y = zmin), curve_df, lty = 2, color = "grey") +
  # geom_line(aes(x = x1, y = zmax), curve_df, lty = 2, color = "black") +
  # geom_line(aes(x = linex, y = liney), curve_df, lty = 2, color = "black") +
  # geom_line(aes(x = linex, y = bumps), curve_df, lty = 2, color = "grey") +
  # geom_vline(xintercept = 2, lty = 2, color = "grey") +
  # geom_line(aes(x = x4, y = y4), curve_df, lty = 2, color = "#9999FF") +
  # geom_line(aes(x = x5, y = y5), curve_df, lty = 2, color = "#99fF99") +
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "kidney")) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "lung"), shape = 15) + 
  # geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == 1, dataset == "breast"), shape = 17) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0)) + 
  #geom_line(aes(x = n, y = D, group = tumour), real_points) + 
  #geom_point(aes(x = n, y = D), filter(real_points, minimal == 0), color = "grey") + 
  #geom_text(aes(x = text_x, y = text_y, label = tumour), filter(real_points, minimal == 1), nudge_x=0.01, nudge_y=0.005) +
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))
# +  
#   theme(legend.position="none")
# pdf("/newMetricsInvestigation/Plot_DiversityDescendant/test.pdf", width = 3.5, height = 3.5)
png(paste0(output_dir_plots_Metric,"/DriverMutationVsClonal_MetricPlot_facetMuAndK_colorS_themebw.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()



######## Diversity in function of the number of cells ######## 
#select for each simulationType 20 simulations and plot the Diversity curves in function of the number of cells.
set.seed(42)

NbRandomSeeds =20

SelectedSimulations_df<-by(SimulationsDefinition, INDICES=SimulationsDefinition$SimulationType, FUN=sample_n, NbRandomSeeds)

SelectedSimulations_df<-as.data.frame(do.call(rbind,SelectedSimulations_df))

SelectedSimulations<-unique(SelectedSimulations_df$SimulationNb)

## Augmented data
data<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){
  
  Kval<-unique(tmp$Var1)[Kiter]
  data[[Kiter]] <- read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), guess_max = 1E4)
  data[[Kiter]]<- select(data[[Kiter]], c(colnames(data[[Kiter]])[1:18],"SimulationNb", "SimulationType", "DriverDiversity", "Generation","NumCells", "Drivers" , "DriverEdgeDiversity", "MeanBirthRate"))
  data[[Kiter]]<-as.data.frame(data[[Kiter]])
  }

data<-do.call(rbind, data)


Filtered_data<-data%>% filter(SimulationNb  %in% SelectedSimulations)
#FilteredDataForMetricPlots<- dataForMetricPlots %>% filter(SimulationNb  %in% SelectedSimulations)

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
output_dir_plots_Metric<-"plots/Batch_ForecastingPaper_bis/MetricPlots"
ifelse(!dir.exists(output_dir_plots_Metric), dir.create(output_dir_plots_Metric), FALSE)


g1<-ggplot(Filtered_data)+
  geom_line(aes(x=NumCells, y=DriverDiversity, color=as.factor(s_driver_birth), group= SimulationNb))+
  facet_grid(K ~ mu_driver_birth)+
  theme_bw()+
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))

png(paste0(output_dir_plots_Metric,"/DriverMutationDiversityVsNumCells.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


# Plot separately for each fitness effect
SValues<-unique(SimulationsDefinition$s_driver_birth)
for(SIndex in seq_along(SValues)){
  g1<-ggplot(subset(Filtered_data,
                    s_driver_birth==SValues[SIndex]
                    ))+
    geom_path(aes(x=NumCells, y=DriverDiversity, group= SimulationNb), color=colfunc(4)[SIndex])+
    facet_grid(K ~ mu_driver_birth)+
    theme_bw()
  
  png(paste0(output_dir_plots_Metric,"/DriverMutationDiversityVsNumCells_S",SValues[SIndex], ".png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
}

# remove the NA value to link the points
g1<-ggplot(Filtered_data[!is.na(Filtered_data$DriverDiversity),])+
  geom_line(aes(x=NumCells, y=DriverDiversity, color=as.factor(s_driver_birth), group= SimulationNb), alpha=0.5)+
  facet_grid(K ~ mu_driver_birth)+
  theme_bw()+
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))+
  guides(alpha=NULL)

png(paste0(output_dir_plots_Metric,"/DriverMutationDiversityVsNumCells_NaRemoved.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


SValues<-unique(SimulationsDefinition$s_driver_birth)
for(SIndex in seq_along(SValues)){
  g1<-ggplot(subset(Filtered_data[!is.na(Filtered_data$DriverDiversity),],
                    s_driver_birth==SValues[SIndex]
  ))+
    geom_path(aes(x=NumCells, y=DriverDiversity, group= SimulationNb), color=colfunc(4)[SIndex])+
    facet_grid(K ~ mu_driver_birth)+
    theme_bw()
  
  png(paste0(output_dir_plots_Metric,"/DriverMutationDiversityVsNumCells_S",SValues[SIndex], "_NaRemoved.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
}
  
SValues<-unique(SimulationsDefinition$s_driver_birth)
for(SIndex in seq_along(SValues)){
  g1<-ggplot(subset(Filtered_data[!is.na(Filtered_data$DriverDiversity),],
                    s_driver_birth==SValues[SIndex]
  ))+
    geom_path(aes(x=NumCells, y=DriverDiversity, group= SimulationNb))+
    facet_grid(K ~ mu_driver_birth)+
    theme_bw()
  
  png(paste0(output_dir_plots_Metric,"/DriverMutationDiversityVsNumCells_S",SValues[SIndex], "_NaRemoved_Black.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
}


# remove the NA value to link the points
g1<-ggplot(Filtered_data[!is.na(Filtered_data$DriverEdgeDiversity),])+
  geom_line(aes(x=NumCells, y=DriverEdgeDiversity, color=as.factor(s_driver_birth), group= SimulationNb), alpha=0.5)+
  facet_grid(K ~ mu_driver_birth)+
  theme_bw()+
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))+
  guides(alpha=NULL)

png(paste0(output_dir_plots_Metric,"/DriverEdgeDiversityVsNumCells_NaRemoved.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()

SValues<-unique(SimulationsDefinition$s_driver_birth)
for(SIndex in seq_along(SValues)){
  g1<-ggplot(subset(Filtered_data[!is.na(Filtered_data$DriverEdgeDiversity),],
                    s_driver_birth==SValues[SIndex]
  ))+
    geom_path(aes(x=NumCells, y=DriverEdgeDiversity, group= SimulationNb), color=colfunc(4)[SIndex])+
    facet_grid(K ~ mu_driver_birth)+
    theme_bw()
  
  png(paste0(output_dir_plots_Metric,"/DriverEdgeDiversityVsNumCells_S",SValues[SIndex], "_NaRemoved.png"), width = 1000, height = 1000, res = 100)
  
  print(g1)
  dev.off()
  
}


## Study of the fitness effect evolution with the num of cells

g1<-ggplot(Filtered_data[!is.na(Filtered_data$MeanBirthRate),])+
  geom_line(aes(x=NumCells, y=MeanBirthRate, color=as.factor(s_driver_birth), group= SimulationNb), alpha=0.5)+
  facet_grid(K ~ mu_driver_birth)+
  theme_bw()+
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))+
  guides(alpha=NULL)

png(paste0(output_dir_plots_Metric,"/MeanBirthRateVsNumCells_NaRemoved.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


######## median instead of mean fitness effect ######## 
output_dir_plots_FitnessEffect<-"plots/Batch_ForecastingPaper_bis/FitnessEffectStudy"
ifelse(!dir.exists(output_dir_plots_FitnessEffect), dir.create(output_dir_plots_FitnessEffect), FALSE)

driver_genotype_properties<-vector(mode="list", length=length(unique(tmp$Var1)))
for(Kiter in seq_along(unique(tmp$Var1))){

  Kval<-unique(tmp$Var1)[Kiter]
  driver_genotype_properties[[Kiter]] <- read_csv(paste0(output_dir_data, "/driver_genotype_propertiesK", Kval, ".csv"), guess_max = 1E4)
}
Alldriver_genotype_properties<-do.call(rbind, driver_genotype_properties)

Alldriver_genotype_properties<-mutate(Alldriver_genotype_properties, 
                                      SimulationType=0, 
                                      SimulationNb=0)

for(i in seq_along(SimulationsDefinition[, 1])){
  
  rowIndex<-which( (Alldriver_genotype_properties$K== SimulationsDefinition$K[i]) &
                     (Alldriver_genotype_properties$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                     (Alldriver_genotype_properties$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]) &
                     (Alldriver_genotype_properties$seed==SimulationsDefinition$seed[i])
                     )
  Alldriver_genotype_properties[rowIndex,
                                which(colnames(Alldriver_genotype_properties)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
  
  Alldriver_genotype_properties[rowIndex,
                                which(colnames(Alldriver_genotype_properties)=="SimulationNb")]<-SimulationsDefinition$SimulationNb[i]
  
  
}

fwrite(Alldriver_genotype_properties, file = paste0(output_dir_data, "/Alldriver_genotype_properties_Augmented.csv"), row.names = FALSE)


Alldriver_genotype_properties<-Alldriver_genotype_properties
MeanMedianBirthRate<-Alldriver_genotype_properties %>% group_by(SimulationNb) %>% mutate(MedianBirthRate = median(BirthRate), 
                                                                          MeanBirthRate = mean(BirthRate)) %>% ungroup()

MeanMedianBirthRatePerSimulation<-unique(MeanMedianBirthRate[, c(colnames( MeanMedianBirthRate)[1:18], 
                                                "SimulationType","SimulationNb","MedianBirthRate","MeanBirthRate")])


g1<-ggplot(MeanMedianBirthRatePerSimulation)+
  geom_point(aes(x=MeanBirthRate, y=MedianBirthRate, color=as.factor(s_driver_birth)), alpha=0.5)+
  facet_grid(K ~ mu_driver_birth)+
  theme_bw()+
  scale_color_manual(name = "Fitness effect", values=(colfunc(4)))+
  guides(alpha=NULL)+
  geom_abline(slope=1, intercept = 0, linetype="dashed",color = "grey")

png(paste0(output_dir_plots_FitnessEffect,"/MeanVsMedianFitnessEffect.png"), width = 1000, height = 1000, res = 100)

print(g1)
dev.off()


######## forecasting in function of biopsy samples ######## 

#start_size_range<-c(10000,25000, 40000, 50000, 70000, 90000, 100000, 125000,250000,375000,500000, 750000) 
start_size_range<-c(125000,250000,375000,500000, 750000)
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value

KValues<-unique(SimulationsDefinition$K)


for(Kiter in c(0,2)){
  
  
  
  #Kval<-KValues[Kiter]
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- read_csv(paste0(output_dir_data, "/dataAugmentedK", Kiter, ".csv"), guess_max = 1E4)
  
  if(Kiter==2){
    #those simulations don't achieve a final size larger thant 750000
    FilteredData<-FilteredData %>%filter( ! (K==4096 & mu_driver_birth==1e-6 & s_driver_birth==0.05 & seed %in% c(52, 1038, 7024, 38009)))
    
  }
  
  summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  print("summary computed")
  
  cols_list<-colnames(FilteredData)[grepl("DriverDiversity",colnames(FilteredData))]
  
  
  
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
  
  fwrite(depth_wait_cor_summary, file = paste0(output_dir_data, "/depth_wait_cor_summary_LargerStartSizeK", (Kiter+1),".csv"), row.names = FALSE)
  fwrite(cor_summary, file = paste0(output_dir_data, "/cor_summary_LargerStartSizeK", (Kiter+1),".csv"), row.names = FALSE)
  fwrite(wait_cor_summary, file = paste0(output_dir_data, "/wait_cor_summary_LargerStartSizeK", (Kiter+1),".csv"), row.names = FALSE)
  fwrite(summary, file = paste0(output_dir_data, "/summary_LargerStartSizeK", (Kiter+1),".csv"), row.names = FALSE)
   
  
  
  # Allcor_summary_smallerStartSize[[Kiter]]<-cor_summary
  # All_wait_cor_summary_smallerStartSize[[Kiter]]<-wait_cor_summary
  # Allsummary_smallerStartSize[[Kiter]]<-summary
  # Alldepth_wait_cor_summary_smallerStartSize[[Kiter]]<-Alldepth_wait_cor_summary
}


Alldepth_wait_cor_summary_SmallerStartSize<-vector(mode="list", length=3)
Alldepth_wait_cor_summary_MediumStartSize<-vector(mode="list", length=3)
Alldepth_wait_cor_summary_LargerStartSize<-vector(mode="list", length=3)
for(Kiter in c(1:3)){
  Alldepth_wait_cor_summary_SmallerStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/depth_wait_cor_summary_smallerStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  Alldepth_wait_cor_summary_MediumStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/depth_wait_cor_summary_mediumStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  Alldepth_wait_cor_summary_LargerStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/depth_wait_cor_summary_LargerStartSizeK", Kiter, ".csv"), guess_max = 1E4)

  }
Alldepth_wait_cor_summary_SmallerStartSize<-do.call(rbind, Alldepth_wait_cor_summary_SmallerStartSize)
Alldepth_wait_cor_summary_MediumStartSize<-do.call(rbind, Alldepth_wait_cor_summary_MediumStartSize)
Alldepth_wait_cor_summary_LargerStartSize<-do.call(rbind, Alldepth_wait_cor_summary_LargerStartSize)


fwrite(Alldepth_wait_cor_summary_SmallerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_SmallerStartSize.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary_MediumStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_MediumStartSize.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary_LargerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_LargerStartSize.csv"), row.names = FALSE)

#Let's put all Aldepthl_wait_cor_summary together
Alldepth_wait_cor_summary_AllStartSize<-rbind(Alldepth_wait_cor_summary_SmallerStartSize,Alldepth_wait_cor_summary_MediumStartSize)
Alldepth_wait_cor_summary_AllStartSize<-rbind(Alldepth_wait_cor_summary_AllStartSize,Alldepth_wait_cor_summary_LargerStartSize )
                                              
fwrite(Alldepth_wait_cor_summary_AllStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_AllStartSize.csv"), row.names = FALSE)



##### read


Alldepth_wait_cor_summary_AllStartSize<-read_csv(paste0(output_dir_data, "/Alldepth_wait_cor_summary_AllStartSize.csv"), guess_max = 1E4)

##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)
## color
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))  

## plots for forecasting in function of biopsy depth Y
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)
#for(MuIndex in seq_along(MuValues)){
for(MuIndex in c(2:length(MuValues))){
  
  tmpSummary<-Alldepth_wait_cor_summary_AllStartSize[which(Alldepth_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  
  SampleSize<-c("1","4", "1Big", "4Big")
  
  for(Samples in seq_along(SampleSize)){
    
    gg<-ggplot(tmpSummary)
    
    for(Val in c(0:10)){
      gg<-gg+geom_line(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), color=as.factor(Val)))
    }
    
    gg<-gg+
      scale_color_manual(name="Depth",
                         values=colfunc(length(c(0:10))),
                         labels=as.character(c(0:10)))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom", SampleSize[Samples] , "SamplesAtDepthX vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepthX_wait_correlations_M",MuValues[MuIndex]  ), "facetK_AllStartSize_AllCor.png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    
    
  }
}
  


  
## plots for forecasting in function of number of biopsies 

           
           
for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-Alldepth_wait_cor_summary_AllStartSize[which(Alldepth_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  
  
  
  for(Depth in c(0:10)){
    
    SampleSize<-c("1","4", "1Big", "4Big")
    gg<-ggplot(tmpSummary)
    
    for(Samples in seq_along(SampleSize)){
      
      gg<-gg+geom_line(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Depth), color=as.factor(Samples)))
    
    }
    
    gg<-gg+
      scale_color_manual(name="Biopsy size",
                         values=colfunc(length(SampleSize)),
                         labels=SampleSize)+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFromXSamplesAtDepth",Depth ,"vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFromXSamplesAtDepth",Depth ,"_wait_correlations_M",MuValues[MuIndex]  ), "facetK_AllStartSize_AllCor.png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    
  }
  
  
}


######## forecasting in function of random biopsy samples ######## 
min_count <- 20

#column of the seed 
seed_index <- which(colnames(data) == "seed")
# column of tall other paramaters
pars_without_seed <- (1:num_parameters)[which((1:num_parameters) != seed_index)]

start_size_range_list<-list(seq(10000,85000, by=15000),seq(90000,100000, by=5000),  c(125000,250000,375000,500000, 750000))
Size_start_size_range<-c("Smaller", "Medium", "Larger")

gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value

KValues<-unique(SimulationsDefinition$K)


#for(SSR in seq_along(start_size_range_list)){
  for(SSR in c(2:length(start_size_range_list))){
  
  start_size_range<-start_size_range_list[[SSR]]
  
  for(Kiter in c(0:2)){
  #for(Kiter in c(1,2)){
    
    
    print(paste("working on Kiter=",Kiter ))
    
    FilteredData<- read_csv(paste0(output_dir_data, "/dataAugmentedK", Kiter, ".csv"), guess_max = 1E4)
    
    
    if(Kiter==2){
      #those simulations don't achieve a final size larger thant 750000
      FilteredData<-FilteredData %>%filter( ! (K==4096 & mu_driver_birth==1e-6 & s_driver_birth==0.05 & seed %in% c(52, 1038, 7024, 38009)))
      
    }
    
    summary <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
    
    print("summary computed")
    
    cols_list<-paste0("DriverDiversityFrom", c("1","4", "1Big", "4Big"), "RandomSamples")
    
   
    #the following lines generates summary dataframe of correlations with "waiting_time"
    wait_cor_summary <- get_wait_cor_summary(summary, cols_list,
                                             num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
    
     
    fwrite(wait_cor_summary, file = paste0(output_dir_data, "/wait_cor_summary_RandomSample_",Size_start_size_range[SSR] , "StartSizeK", (Kiter+1),".csv"), row.names = FALSE)
    
    
  }
  
}


Alldepth_wait_cor_summary_RandomSample_SmallerStartSize<-vector(mode="list", length=3)
Alldepth_wait_cor_summary_RandomSample_MediumStartSize<-vector(mode="list", length=3)
Alldepth_wait_cor_summary_RandomSample_LargerStartSize<-vector(mode="list", length=3)
for(Kiter in c(1:3)){
  Alldepth_wait_cor_summary_RandomSample_SmallerStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/wait_cor_summary_RandomSample_SmallerStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  Alldepth_wait_cor_summary_RandomSample_MediumStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/wait_cor_summary_RandomSample_MediumStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  Alldepth_wait_cor_summary_RandomSample_LargerStartSize[[Kiter]]<-read_csv(paste0(output_dir_data,"/wait_cor_summary_RandomSample_LargerStartSizeK", Kiter, ".csv"), guess_max = 1E4)
  
}
Alldepth_wait_cor_summary_RandomSample_SmallerStartSize<-do.call(rbind, Alldepth_wait_cor_summary_RandomSample_SmallerStartSize)
Alldepth_wait_cor_summary_RandomSample_MediumStartSize<-do.call(rbind, Alldepth_wait_cor_summary_RandomSample_MediumStartSize)
Alldepth_wait_cor_summary_RandomSample_LargerStartSize<-do.call(rbind, Alldepth_wait_cor_summary_RandomSample_LargerStartSize)


fwrite(Alldepth_wait_cor_summary_RandomSample_SmallerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_RandomSample_SmallerStartSize.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary_RandomSample_MediumStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_RandomSample_MediumStartSize.csv"), row.names = FALSE)
fwrite(Alldepth_wait_cor_summary_RandomSample_LargerStartSize, file = paste0(output_dir_data, "/Alldepth_wait_cor_summary_RandomSample_LargerStartSize.csv"), row.names = FALSE)

#Let's put all Aldepthl_wait_cor_summary together
All_wait_cor_summary_RandomSample_AllStartSize<-rbind(Alldepth_wait_cor_summary_RandomSample_SmallerStartSize,Alldepth_wait_cor_summary_RandomSample_MediumStartSize)
All_wait_cor_summary_RandomSample_AllStartSize<-rbind(All_wait_cor_summary_RandomSample_AllStartSize,Alldepth_wait_cor_summary_RandomSample_LargerStartSize )

fwrite(All_wait_cor_summary_RandomSample_AllStartSize, file = paste0(output_dir_data, "/All_wait_cor_summary_RandomSample_AllStartSize.csv"), row.names = FALSE)


#plots

##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)

All_wait_cor_summary_RandomSample_AllStartSize<-read_csv(paste0(output_dir_data, "/All_wait_cor_summary_RandomSample_AllStartSize.csv"), guess_max = 1E4 )
Alldepth_wait_cor_summary_AllStartSize<-read_csv(paste0(output_dir_data, "/Alldepth_wait_cor_summary_AllStartSize.csv"), guess_max = 1E4)

MuValues<-unique(SimulationsDefinition$mu_driver_birth)
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))  

for(MuIndex in seq_along(MuValues)){
  
  tmpSummary<-Alldepth_wait_cor_summary_AllStartSize[which(Alldepth_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
  tmpSummary<-tmpSummary[-which(tmpSummary$start_size %in% c(500,1000,2000,5000, 20000, 80000)), ]
  
  tmpSummary_RandomSample<-All_wait_cor_summary_RandomSample_AllStartSize[which(All_wait_cor_summary_RandomSample_AllStartSize$mu_driver_birth==MuValues[MuIndex]), ]
  tmpSummary_RandomSample<-tmpSummary_RandomSample[-which(tmpSummary_RandomSample$K==64 & tmpSummary_RandomSample$start_size %in% c(25000, 55000, 70000, 85000)),]
  MergetmpSummary<-merge(tmpSummary, tmpSummary_RandomSample, by= intersect(colnames(tmpSummary), colnames(tmpSummary_RandomSample)))
  
  SampleSize<-c("1","4", "1Big", "4Big")
  
  for(Samples in seq_along(SampleSize)){
    
    # gg<-ggplot()+
    #   geom_line(data=tmpSummary_RandomSample, aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "RandomSamples")), linetype="dashed", color="black")
    
    gg<-ggplot(MergetmpSummary)+
      geom_line(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "RandomSamples")), color="black")
    
    
    for(Depth in c(0:5)){
 
      gg<-gg+geom_line(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Depth), color=as.factor(Depth)))
      
    }
    
  
  
    gg<-gg+
      scale_color_manual(name="Depth",
                         values=colfunc(length(c(0:5))),
                         labels=as.character(c(0:5)))+
      # scale_linetype_manual(name="",
      #                       values=c("twodash","solid"),
      #                       labels=c("Random", "Selected"))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom",SampleSize[Samples] , "Samplesvs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "Samples_wait_correlations_M",MuValues[MuIndex]  ), "facetK_AllStartSize_SelectedVsRandom.png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    
  }
  
  
}
           


## debugging
tmpSummary<-Alldepth_wait_cor_summary_AllStartSize[which(Alldepth_wait_cor_summary_AllStartSize$mu_driver_birth==MuValues[MuIndex]),]
tmpSummary_RandomSample<-All_wait_cor_summary_RandomSample_AllStartSize[which(All_wait_cor_summary_RandomSample_AllStartSize$mu_driver_birth==MuValues[MuIndex]), ]

test_1<-unique(tmpSummary[, which(colnames(tmpSummary) %in% c("K", "mu_driver_birth", "s_driver_birth", "start_size") )])
test_2<-unique(tmpSummary_RandomSample[, which(colnames(tmpSummary_RandomSample) %in% c("K", "mu_driver_birth", "s_driver_birth", "start_size") )])

Jointest<-inner_join(test_1, test_2)

length(which(! tmpSummary$start_size %in% c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000))))
length(which(! tmpSummary_RandomSample$start_size %in% c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000))))

length(which( tmpSummary$start_size %in% c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000))))

tmpSummary$start_size[which(! tmpSummary$start_size %in% c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000)))]

setdiff(unique(tmpSummary$start_size), c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000)))
setdiff(c(seq(10000,85000, by=15000), seq(90000,100000, by=5000), c(125000,250000,375000,500000, 750000)),unique(tmpSummary$start_size))
setdiff(unique(tmpSummary_RandomSample$start_size), unique(tmpSummary$start_size))


for(Kt in c(64, 512, 4096)){
  print(paste("K=", Kt))
  tt<-tmpSummary[which(tmpSummary$K==Kt), ]
  
  for(st in c(0.05, 0.1, 0.15, 0.2)){
    print(paste("s=", st))
    ttt<-tt[which(tt$s_driver_birth==st), ]
    print(paste("nb start_size :", length(unique(ttt$start_size))))
    
  }
  
}

unique(tmpSummary[which(tmpSummary$K==64), which(colnames(tmpSummary)=="start_size")])
unique(tmpSummary[which(tmpSummary$K==64), "start_size"])

colnames(cbind_tmpSummary)

### Scatter plots to figure out what is going wrong in the correlation plot
DataForScatterPlotK512<- read_csv(paste0(output_dir_data, "/dataAugmentedK1.csv"), guess_max = 1E4)
##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)


SampleSize<-c("1","4", "1Big", "4Big")

for(Samples in seq_along(SampleSize)){

  #for(Val in c(0:10)){
  Val=0
  gg<-ggplot()
  
    
    subsetDataForScatterPlotK512<-filter(DataForScatterPlotK512, 
                                         DataForScatterPlotK512$NumCells %in% c(150000,250000, 350000, 1e06), 
                                         DataForScatterPlotK512$K==512, 
                                         DataForScatterPlotK512$s_driver_birth==0.1, 
                                         DataForScatterPlotK512$mu_driver_birth==1e-06)
    
    subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb",  paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
    
    for(NumCellsAtMeasurement in c(150000,250000, 350000, 1e06)){
      gg<-gg+
        geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y=FinalTumourDiversity, color=as.factor(SimulationNb)), position="jitter")
    }
    
    
    subsetDataForScatterPlotK512<-tidyr::spread(subsetDataForScatterPlotK512, key="NumCells", value=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val))
    
    colnames(subsetDataForScatterPlotK512)<-c("SimulationNb", "DiversityTumourAtMeasurement", "FinalTumourDiversity")
    
    #gg<-ggplot(subsetDataForScatterPlotK512)
    gg<-gg+
      geom_point(aes(x=DiversityTumourAtMeasurement, y=FinalTumourDiversity, color=as.factor(SimulationNb)), position="jitter")
    
  
    subsetDataForScatterPlotK512<-filter(DataForScatterPlotK512, 
                                         DataForScatterPlotK512$NumCells %in% c(TumourSizeAtmeasur, 1e06), 
                                         DataForScatterPlotK512$K==512, 
                                         DataForScatterPlotK512$s_driver_birth==0.1, 
                                         DataForScatterPlotK512$mu_driver_birth==1e-06)
    
    subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb",  paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
                              
    subsetDataForScatterPlotK512<-tidyr::spread(subsetDataForScatterPlotK512, key="NumCells", value=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val))
    
    colnames(subsetDataForScatterPlotK512)<-c("SimulationNb", "DiversityTumourAtMeasurement", "FinalTumourDiversity")
    
    #gg<-ggplot(subsetDataForScatterPlotK512)
    gg<-gg+
      geom_point(aes(x=DiversityTumourAtMeasurement, y=FinalTumourDiversity, color=as.factor(SimulationNb)), position="jitter")+
      ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n Tumour size at measurement = 250000"))+
      theme_bw()+
      guides(color=FALSE)+
      theme(plot.title = element_text(hjust = 0.5))
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds.png")),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    

}

########### All summary ########
#summary contains the variable waiting_time
Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)


SampleSize<-c("1","4", "1Big", "4Big")

for(Samples in seq_along(SampleSize)){
  
  #for(Val in c(0:10)){
  Val=0

  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$NumCells >= 150000),
                                       (Allsummary$NumCells <= 350000), 
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
    gg<-ggplot(subsetDataForScatterPlotK512)+
      geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color="NumCells"), position="jitter")+
      ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n Versus waiting time"))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds.png")),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
  
}

# for comparison
Samples=4
subsetDataForScatterPlotK512<-filter(Allsummary, 
                                     (Allsummary$NumCells >= 150000),
                                     (Allsummary$NumCells <= 350000), 
                                     Allsummary$K==512, 
                                     Allsummary$s_driver_birth==0.05, 
                                     Allsummary$mu_driver_birth==1e-06)

subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]

gg<-ggplot(subsetDataForScatterPlotK512)+
  geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color="NumCells"), position="jitter")+
  ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n Versus waiting time"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds.png")),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

Samples=4
subsetDataForScatterPlotK512<-filter(Allsummary, 
                                     (Allsummary$NumCells >= 150000),
                                     (Allsummary$NumCells <= 350000), 
                                     Allsummary$K==512, 
                                     Allsummary$s_driver_birth==0.05, 
                                     Allsummary$mu_driver_birth==1e-06)

subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]

gg<-ggplot(subsetDataForScatterPlotK512)+
  geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color="NumCells"), alpha=0.5)+
  ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.05"))+
  theme_bw()+
  guides(alpha=FALSE)+
  scale_color_continuous(type = "viridis")+
  theme(plot.title = element_text(hjust = 0.5))

png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds.png")),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

  
SampleSize<-c("1","4", "1Big", "4Big")

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size >= 125000),
                                       (Allsummary$start_size <= 375000), 
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color=sprintf("factor(%s)",subsetDataForComparison_ScatterPlot$start_size)), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForComparison_ScatterPlot<-filter(Allsummary, 
                                              (Allsummary$start_size >= 125000),
                                              (Allsummary$start_size <= 375000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.05, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForComparison_ScatterPlot<-subsetDataForComparison_ScatterPlot[, which(colnames(subsetDataForComparison_ScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom4BigSamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForComparison_ScatterPlot)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom4BigSamplesAtDepth",Val), y="waiting_time", color=sprintf("factor(%s)",subsetDataForComparison_ScatterPlot$start_size)), alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom4BigSamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.05"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  

  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparison.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}


#

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size >= 125000),
                                       (Allsummary$start_size <= 375000), 
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","SimulationNb" ,"waiting_time", "start_size",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  #tmp<-by(subsetDataForScatterPlotK512, INDICES=subsetDataForScatterPlotK512$start_size, FUN=find_correlations, "waiting_time", paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), 20)
  
  colnames(subsetDataForScatterPlotK512)[which(colnames(subsetDataForScatterPlotK512)== paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val))]<-"var2"
  #CorSpearman<-by(subsetDataForScatterPlotK512, INDICES=subsetDataForScatterPlotK512$start_size, FUN=cor, subsetDataForScatterPlotK512$waiting_time, subsetDataForScatterPlotK512$var2, method = "spearman")
  
  CorSpearman_tmp<-group_by(subsetDataForScatterPlotK512, start_size) %>% mutate(CorSpearman=cor(waiting_time, var2,method = "spearman" ))
  CorSpearman_tmp<-unique( CorSpearman_tmp[,c("start_size", "CorSpearman")])
  
  if(Samples==1){
        #CorrelationReComputed<-do.call(rbind,tmp)
    CorSpearman_df<-as.data.frame(CorSpearman_tmp[, "CorSpearman" ])
  }else{
    #CorrelationReComputed<-cbind(CorrelationReComputed, do.call(rbind,tmp))
    CorSpearman_df<-cbind(CorSpearman_df, as.data.frame(unique( CorSpearman_tmp[,"CorSpearman"])))
  }
  
 
   
}
colnames(CorrelationReComputed)<-SampleSize
colnames(CorSpearman_df)<-SampleSize
rownames(CorSpearman_df)<-c("125000", "250000", "375000")







for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size >= 125000),
                                       (Allsummary$start_size <= 375000), 
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color="start_size"), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotK64<-filter(Allsummary, 
                                              (Allsummary$start_size >= 125000),
                                              (Allsummary$start_size <= 375000),
                                              Allsummary$K==64, 
                                              Allsummary$s_driver_birth==0.1, 
                                              Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK64<-subsetDataForScatterPlotK64[, which(colnames(subsetDataForScatterPlotK64) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotK64)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", color="start_size"), alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K64 mu=1e-06 s=0.1"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonWithSmallerK.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotK64<-filter(Allsummary, 
                                      (Allsummary$start_size == 250000),
                                      Allsummary$K==64, 
                                      Allsummary$s_driver_birth==0.1, 
                                      Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK64<-subsetDataForScatterPlotK64[, which(colnames(subsetDataForScatterPlotK64) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotK64)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K64 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonUniqueStartSize.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotComparison<-filter(Allsummary, 
                                      (Allsummary$start_size == 125000),
                                      Allsummary$K==512, 
                                      Allsummary$s_driver_birth==0.05, 
                                      Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotComparison<-subsetDataForScatterPlotComparison[, which(colnames(subsetDataForScatterPlotComparison) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotComparison)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=125000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonBestCorrelation.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}

######## confidence interval for correlation coefficient ########
#summary contains the variable waiting_time
Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)


cols_list<-colnames(Allsummary)[grepl("DriverDiversity",colnames(Allsummary))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10))

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]

cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")

#the following lines generates summary dataframe of correlations with "waiting_time"
min_count=20
wait_cor_summary <- get_wait_cor_summary_modified(Allsummary, cols_list, 
                                         num_parameters = num_parameters, min_count = min_count)


wait_cor_summary_exactNull <- get_wait_cor_summary_modified(Allsummary, cols_list, 
                                                            num_parameters = num_parameters, min_count = min_count)

wait_cor_summary_exactFALSE <- get_wait_cor_summary_modified(Allsummary, cols_list, 
                                                            num_parameters = num_parameters, min_count = min_count)

#wait_cor_summary_RVAideMemoire <- get_wait_cor_summary_modified(Allsummary, cols_list, 
                                                             #num_parameters = num_parameters, min_count = min_count)

#fwrite(wait_cor_summary_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summary_RVAideMemoire.csv"), row.names = FALSE)


wait_cor_summary_RVAideMemoire<-read_csv(paste0(output_dir_data, "/wait_cor_summary_RVAideMemoire.csv"), guess_max = 1E4)
# #debugging
# 
# subsetDataForScatterPlotComparison<-filter(Allsummary, 
#                                            (Allsummary$start_size == 125000),
#                                            Allsummary$K==512, 
#                                            Allsummary$s_driver_birth==0.05, 
#                                            Allsummary$mu_driver_birth==1e-06)
# test<-get_wait_cor_summary_modified(subsetDataForScatterPlotComparison, cols_list, 
#                               num_parameters = num_parameters, min_count = min_count)
# test<-get_wait_cor_summary_modified(Allsummary, cols_list, 
#                                     num_parameters = num_parameters, min_count = min_count)
# 
# 
# test<-subsetDataForScatterPlotComparison %>% 
#   mutate_(variance = interp(~var(var1), var1 = as.name(factor2))) %>% 
#   filter(variance > 0) %>% # to avoid warnings when all values of factor2 are identical
#   mutate_(count1 = interp(~length(var1), var1 = as.name(factor1)), 
#           count2 = interp(~length(var2), var2 = as.name(factor2))) %>% 
#   filter(count1 >= min_count, count2 >= min_count) %>% 
#   summarise_(temp_name = interp(~cor(var1, var2, method = "spearman"), var1 = as.name(factor1), var2 = as.name(factor2)),
#              temp_name_pval = interp(~cor_pval(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2)),
#              temp_name_cilow = interp(~corCI_low(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2)),
#              temp_name_cihigh = interp(~corCI_high(var1, var2), var1 = as.name(factor1), var2 = as.name(factor2)))%>% 
#   rename_(.dots = setNames(c("temp_name", "temp_name_pval","temp_name_cilow", "temp_name_cihigh") , paste0(result_name,  c("", "_pVal", "_CI_low", "_CI_high")) ))
# 
# 
# 
# 
# 
# 
# 
# result <- do.call(rbind,lapply(cols_list, function(x) {
#   cor.result<-cor.test(data$waiting_time,data[[x]])
#   # pvalue <- cor.result$p.value
#   # estimate <- cor.result$conf.int[1]
#   return(data.frame(pvalue = pvalue, CI_low = cor.result$conf.int[1], CI_high = cor.result$conf.int[2] ))
# })
# )
# 
# 


##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)
## color
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))  

## plots for forecasting in function of biopsy depth Y
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)
#for(MuIndex in seq_along(MuValues)){
for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary[which(wait_cor_summary$mu_driver_birth==MuValues[MuIndex]),]
  
  
  SampleSize<-c("1","4", "1Big", "4Big")
  
  for(Samples in seq_along(SampleSize)){
    
    gg<-ggplot(tmpSummary)
    
    Val=0
    gg<-gg+geom_point(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), color=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val,"_pVal")))

    
    gg<-gg+
      #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
      scale_colour_gradient(name="p value", low="green",high="red", breaks=seq(0, 1, 0.05), labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom", SampleSize[Samples] , "SamplesAtDepth0 vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth0_wait_correlations_M",MuValues[MuIndex]  ), "pValue.png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    
    
  }
}


for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary[which(wait_cor_summary$mu_driver_birth==MuValues[MuIndex]),]
  
  
  SampleSize<-c("1","4", "1Big", "4Big")
  
  for(Samples in seq_along(SampleSize)){
    
    Val=0
    tmpSampleSummary<-mutate(tmpSummary, 
                       pValSmaller0_05 = ifelse(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val,"_pVal")]] <= 0.05,tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val,"_pVal")]], as.logical(NA) ),
                       pValSmaller0_01 = ifelse(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val,"_pVal")]] <= 0.01,tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val,"_pVal")]], as.logical(NA) ))
  
    gg1<-ggplot(tmpSampleSummary)+
      geom_point(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), color="pValSmaller0_05"))
    
    
    gg1<-gg1+
      #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
      scale_colour_gradient(name="p value",
                            low="green",high="blue",
                            na.value = "red", 
                            breaks=seq(0, 0.05, 0.01))+
                            #labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom", SampleSize[Samples] , "SamplesAtDepth0 vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex], "\n Red points : p-value larger than 0.05"))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth0_wait_correlations_M",MuValues[MuIndex]  ), "pValueSmaller005.png"),  width = 1000, height = 1000, res = 100)
    print(gg1)
    dev.off()
    
    if(length(which(! is.na(tmpSampleSummary$pValSmaller0_01))) > 0){
      gg2<-ggplot(tmpSampleSummary)+
        geom_point(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), color="pValSmaller0_01"))
      
      
      gg2<-gg2+
        #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
        scale_colour_gradient(name="p value",
                              low="green",high="blue",
                              na.value = "red", 
                              breaks=seq(0, 0.01, 0.005))+
        #labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
        facet_grid( K ~ s_driver_birth)+
        theme_bw()+
        ylim(-1, 1)+
        geom_hline(aes(yintercept=0), linetype="dashed")+
        xlab("tumour size at measurement")+
        ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom", SampleSize[Samples] , "SamplesAtDepth0 vs waiting time"))+
        ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex],  "\n Red points : p-value larger than 0.01"))+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 90))
      
      
      png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth0_wait_correlations_M",MuValues[MuIndex]  ), "pValueSmaller001.png"),  width = 1000, height = 1000, res = 100)
      print(gg2)
      dev.off()
    }
    
    
    
  }
}


##output dir
output_dir_plots_WaitCorrelation_biopsy<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_biopsy"
ifelse(!dir.exists(output_dir_plots_WaitCorrelation_biopsy), dir.create(output_dir_plots_WaitCorrelation_biopsy), FALSE)
## color
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))  

## plots for forecasting in function of biopsy depth Y
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)
for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  
  SampleSize<-c("1","4", "1Big", "4Big")
  
  for(Samples in seq_along(SampleSize)){
    
    gg<-ggplot(tmpSummary)
    
    Val=0
    gg<-gg+
      geom_line(aes_string(x="start_size", y=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))+
      geom_errorbar(aes_string(x="start_size", 
                               ymin=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "_CI_low"), 
                               ymax=paste0("Cor_DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "_CI_high")),
                    width=.2,
                    position=position_dodge(0.05))
    
    
    gg<-gg+
      #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
      #scale_colour_gradient(name="p value", low="green",high="red", breaks=seq(0, 1, 0.05), labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", "DriverDiversityFrom", SampleSize[Samples] , "SamplesAtDepth0 vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth0_wait_correlations_M",MuValues[MuIndex]  ), "ErrorBar.png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
    
    
  }
}

## add jitter and  transparency to these plots
colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position=position_jitter(width = 0.1, height = 0, seed = 1), alpha=0.7)+
    geom_hex(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", fill = "..density.."), binwidth=c(0.01, 0.01) ) +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotComparison<-filter(Allsummary, 
                                             (Allsummary$start_size == 125000),
                                             Allsummary$K==512, 
                                             Allsummary$s_driver_birth==0.05, 
                                             Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotComparison<-subsetDataForScatterPlotComparison[, which(colnames(subsetDataForScatterPlotComparison) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotComparison)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position=position_jitter(width = 0.1, height = 0, seed = 1), alpha=0.7)+
    #geom_dotplot(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binaxis="y", binwidth=0.01, stackdir = "center", stackratio = 0.75)+
    geom_hex(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", fill = "..density.."), binwidth=c(0.01, 0.01) ) +
    # stat_density2d(h = 0.01, bins = 100,
    #                aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time",  fill = "..level.."), geom = "polygon") +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=125000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonBestCorrelation_GeomHexDensity.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}


for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position=position_jitter(width = 0.1, height = 0, seed = 1), alpha=0.7)+
    geom_hex(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binwidth=c(0.01, 0.01) ) +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotComparison<-filter(Allsummary, 
                                             (Allsummary$start_size == 125000),
                                             Allsummary$K==512, 
                                             Allsummary$s_driver_birth==0.05, 
                                             Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotComparison<-subsetDataForScatterPlotComparison[, which(colnames(subsetDataForScatterPlotComparison) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotComparison)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position=position_jitter(width = 0.1, height = 0, seed = 1), alpha=0.7)+
    #geom_dotplot(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binaxis="y", binwidth=0.01, stackdir = "center", stackratio = 0.75)+
    geom_hex(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binwidth=c(0.01, 0.01) ) +
    # stat_density2d(h = 0.01, bins = 100,
    #                aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time",  fill = "..level.."), geom = "polygon") +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=125000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonBestCorrelation_GeomHexCount.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}

for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    stat_density2d(h = 0.5, bins = 10,
                   aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", fill="..level.."), geom = "polygon") +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotComparison<-filter(Allsummary, 
                                             (Allsummary$start_size == 125000),
                                             Allsummary$K==512, 
                                             Allsummary$s_driver_birth==0.05, 
                                             Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotComparison<-subsetDataForScatterPlotComparison[, which(colnames(subsetDataForScatterPlotComparison) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotComparison)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position=position_jitter(width = 0.1, height = 0, seed = 1), alpha=0.7)+
    #geom_dotplot(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binaxis="y", binwidth=0.01, stackdir = "center", stackratio = 0.75)+
    #geom_hex(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), binwidth=c(0.01, 0.01) ) +
    stat_density2d(h = 0.5, bins = 10,
                   aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time", fill="..level.."), geom = "polygon") +
    colorscale+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=125000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonBestCorrelation_DensityMap.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}

#### p-values with the other diversity measures

for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  
  DiversityMeasure<-c("DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  
  for(Meas in seq_along(DiversityMeasure)){
    
    
    tmpSampleSummary<-mutate(tmpSummary, 
                             pValSmaller0_05 = ifelse(tmpSummary[[paste0("Cor_", DiversityMeasure[Meas] ,"_pVal")]] <= 0.05,tmpSummary[[paste0("Cor_", DiversityMeasure[Meas] ,"_pVal")]], as.logical(NA) ),
                             pValSmaller0_01 = ifelse(tmpSummary[[paste0("Cor_", DiversityMeasure[Meas] ,"_pVal")]] <= 0.01,tmpSummary[[paste0("Cor_", DiversityMeasure[Meas] ,"_pVal")]], as.logical(NA) ))
    
    gg1<-ggplot(tmpSampleSummary)+
      geom_point(aes_string(x="start_size", y=paste0("Cor_", DiversityMeasure[Meas]), color="pValSmaller0_05"))
    
    
    gg1<-gg1+
      #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
      scale_colour_gradient(name="p value",
                            low="green",high="blue",
                            na.value = "red", 
                            breaks=seq(0, 0.05, 0.01))+
      #labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw()+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n ", DiversityMeasure[Meas], " vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex], "\n Red points : p-value larger than 0.05"))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90))
    
    
    png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", DiversityMeasure[Meas], "_wait_correlations_M",MuValues[MuIndex]  ), "pValueSmaller005.png"),  width = 1000, height = 1000, res = 100)
    print(gg1)
    dev.off()
    
    if(length(which(! is.na(tmpSampleSummary$pValSmaller0_01))) > 0){
      gg2<-ggplot(tmpSampleSummary)+
        geom_point(aes_string(x="start_size", y=paste0("Cor_", DiversityMeasure[Meas]), color="pValSmaller0_01"))
      
      
      gg2<-gg2+
        #scale_colour_gradient2(name="p value", low="green", mid="blue", high="red", midpoint=0.05, breaks=seq(0, 1, 0.05))+
        scale_colour_gradient(name="p value",
                              low="green",high="blue",
                              na.value = "red", 
                              breaks=seq(0, 0.01, 0.005))+
        #labels=c("0", "0.05", "", "0.15", "", "0.25", "", "0.35", "", "0.45", "", "0.55", "", "0.65", "", "0.75", "", "0.85", "", "0.95", ""))+
        facet_grid( K ~ s_driver_birth)+
        theme_bw()+
        ylim(-1, 1)+
        geom_hline(aes(yintercept=0), linetype="dashed")+
        xlab("tumour size at measurement")+
        ylab(paste0("correlation coefficient:\n ", DiversityMeasure[Meas], " vs waiting time"))+
        ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex], "\n Red points : p-value larger than 0.01"))+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 90))
      
      
      png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", DiversityMeasure[Meas], "_wait_correlations_M",MuValues[MuIndex]  ), "pValueSmaller001.png"),  width = 1000, height = 1000, res = 100)
      print(gg2)
      dev.off()
    }
    
    
    
  }
}


for(Samples in seq_along(SampleSize)){
  
  Val=0
  
  subsetDataForScatterPlotK512<-filter(Allsummary, 
                                       (Allsummary$start_size == 250000),
                                       Allsummary$K==512, 
                                       Allsummary$s_driver_birth==0.1, 
                                       Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK512<-subsetDataForScatterPlotK512[, which(colnames(subsetDataForScatterPlotK512) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg1<-ggplot(subsetDataForScatterPlotK512)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K512 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subsetDataForScatterPlotK64<-filter(Allsummary, 
                                      (Allsummary$start_size == 250000),
                                      Allsummary$K==64, 
                                      Allsummary$s_driver_birth==0.1, 
                                      Allsummary$mu_driver_birth==1e-06)
  
  subsetDataForScatterPlotK64<-subsetDataForScatterPlotK64[, which(colnames(subsetDataForScatterPlotK64) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time",   paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val)))]
  
  gg2<-ggplot(subsetDataForScatterPlotK64)+
    geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    ggtitle(paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val, "\n K64 mu=1e-06 s=0.1 start_size=250000"))+
    theme_bw()+
    guides(alpha=FALSE)+
    scale_color_continuous(type = "viridis")+
    ylim(0, 0.7)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckCorrelation_",SampleSize[Samples] , "Samples_depth", Val, "allseeds_comparisonUniqueStartSize.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  
}


########  plots to study the positiv unsignificant correlation coefficients #### 


#K=64 mu=1e-06 s=0.15 start_size=750000
KOfInterest=64
start_sizeOfInterest=750000
s_driver_birthOfInterest =0.15
mu_driver_birthOfInterest = 1e-06

  subsetDataForScatterPlot<-filter(Allsummary, 
                                       (Allsummary$start_size == start_sizeOfInterest),
                                       Allsummary$K==KOfInterest, 
                                       Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                       Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
##############
  #K=64 mu=1e-06 s=0.2 start_size=500000
  KOfInterest=64
  start_sizeOfInterest=500000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  ##############
  #K=64 mu=1e-06 s=0.2 start_size=750000
  KOfInterest=64
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  ##750000 vs 500000
  #K=64 mu=1e-06 s=0.2 start_size=750000 vs 500000
  KOfInterest=64
  start_sizeOfInterest=500000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    ylim(0, 0.15)+
    xlim(1, 10)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #K=64 mu=1e-06 s=0.2 start_size=750000
  KOfInterest=64
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    ylim(0, 0.15)+
    xlim(1, 10)+
    theme(plot.title = element_text(hjust = 0.5))
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize750000vs500000.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  

  ##############
  #K=512 mu=1e-06 s=0.05 start_size=500000
  KOfInterest=512
  start_sizeOfInterest=500000
  s_driver_birthOfInterest =0.05
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  ##############
  #K=512 mu=1e-06 s=0.05 start_size=750000
  KOfInterest=512
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.05
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  ##############
  #K=512 mu=1e-06 s=0.05 start_size=750000 vs 500000
  KOfInterest=512
  start_sizeOfInterest=500000
  s_driver_birthOfInterest =0.05
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    ylim(0, 0.3)+
    xlim(1, 8)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #K=512 mu=1e-06 s=0.05 start_size=500000
  KOfInterest=512
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.05
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    ylim(0, 0.3)+
    xlim(1, 8)+
    theme(plot.title = element_text(hjust = 0.5))
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize750000vs500000.png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
  
  
  ##############
  #K=4096 mu=1e-06 s=0.2 start_size=750000
  KOfInterest=4096
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg2<-ggplot(subsetDataForScatterPlot)+
    #geom_point(aes_string(x=paste0("DriverDiversityFrom",SampleSize[Samples] , "SamplesAtDepth",Val), y="waiting_time"), alpha=0.5)+
    # geom_hex(aes(x=DriverDiversity, y=waiting_time, fill = ..density..), binwidth=c(0.01, 0.01) ) +
    geom_count(aes(x=DriverDiversity, y=waiting_time))+
    colorscale+
    ggtitle(paste0("DriverDiversity \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  ".png")),  width = 1000, height = 1000, res = 100)
  grid.arrange(gg1, gg2, nrow=2)
  dev.off()
 
  
  
  ###########
  #K=64 mu=1e-06 s=0.2 start_size=125000
  KOfInterest=64
  start_sizeOfInterest=125000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  #subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=MeanBirthRate), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity vs MeanBirthRate \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    #ylim(0, 0.1)+
    theme(plot.title = element_text(hjust = 0.5))
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  "DriverVsMeanBirthRate.png")),  width = 1000, height = 1000, res = 100)
  print(gg1)
  dev.off()
  
  cor(subsetDataForScatterPlot$DriverDiversity, subsetDataForScatterPlot$MeanBirthRate, method="spearman")
  cor.test(subsetDataForScatterPlot$DriverDiversity, subsetDataForScatterPlot$MeanBirthRate, method="spearman", alternative = "two.sided", exact=FALSE)

  ##############
  #K=64 mu=1e-06 s=0.2 start_size=750000
  KOfInterest=64
  start_sizeOfInterest=750000
  s_driver_birthOfInterest =0.2
  mu_driver_birthOfInterest = 1e-06
  
  subsetDataForScatterPlot<-filter(Allsummary, 
                                   (Allsummary$start_size == start_sizeOfInterest),
                                   Allsummary$K==KOfInterest, 
                                   Allsummary$s_driver_birth==s_driver_birthOfInterest, 
                                   Allsummary$mu_driver_birth==mu_driver_birthOfInterest)
  
  
  subsetDataForScatterPlot<-subsetDataForScatterPlot %>% filter(gap == min(gap, na.rm = TRUE))
  
  #subsetDataForScatterPlot<-subsetDataForScatterPlot[, which(colnames(subsetDataForScatterPlot) %in% c("NumCells","start_size", "SimulationNb" ,"waiting_time", "DriverDiversity"))]
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=MeanBirthRate), position="jitter", alpha=0.5)+
    ggtitle(paste0("DriverDiversity vs MeanBirthRate \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    theme(plot.title = element_text(hjust = 0.5))
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  "DriverVsMeanBirthRate.png")),  width = 1000, height = 1000, res = 100)
  print(gg1)
  dev.off()
  
  gg1<-ggplot(subsetDataForScatterPlot)+
    geom_point(aes(x=DriverDiversity, y=waiting_time), alpha=0.5)+
    ggtitle(paste0("DriverDiversity vs waiting_time \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
    theme_bw()+
    guides(alpha=FALSE)+
    ylim(0, 0.06)+
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio=1)
  
  png(paste0(output_dir_plots_WaitCorrelation_biopsy, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  "DriverVsMeanBirthRate_Ylim.png")),  width = 1000, height = 1000, res = 100)
  print(gg1)
  dev.off()
  
################Muller plots  ################
  #search for the simulation with large meanBirthRate and small driverDiversity
  as.data.frame(subsetDataForScatterPlot[which(subsetDataForScatterPlot$MeanBirthRate==max(subsetDataForScatterPlot$MeanBirthRate)),])
  
  seedOfInterest=as.numeric(subsetDataForScatterPlot[which(subsetDataForScatterPlot$MeanBirthRate==max(subsetDataForScatterPlot$MeanBirthRate)),
                                which(colnames(subsetDataForScatterPlot)=="seed")])
  
  pathToData<-"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_0/mu_driver_birth_0/s_driver_birth_3"
  
  for(seedIt in c(0:99)){
    tmpParam<-read_delim(paste0(pathToData, "/seed_", seedIt, "/parameters.dat"), "\t", escape_double = FALSE, trim_ws = TRUE)
    if(tmpParam$seed==seedOfInterest){
      print(paste0("associated seed ", seedIt))
      seedFolder=seedIt
      break
    }
    
  }
  

     plot_all_images(paste0(pathToData, "/seed_", seedFolder),
                  output_filename = paste0("MullerPlot_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest, "seedFolder", seedFolder ),
                  file_type = "png", 
                  output_dir = output_dir_plots_WaitCorrelation_biopsy,
                  trim = -1)
     
     Muller_df <- muller_df_from_file(paste0(pathToData,"/seed_",seedFolder,"/driver_phylo.dat"))
     #time compute the num of cells
     Muller_df<-group_by(Muller_df, Generation) %>% mutate(Time = sum(Population)) %>% ungroup()
     Muller_df$Generation<-NULL
     #colnames(Muller_df)[which(colnames(Muller_df)=="NumCells")]<-"Time"
     if(class(Muller_df) != "data.frame") return(NA)
     
     long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                       "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                       "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                       "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                       "#5E738F", "#D1A33D")
     dd <- 0:25
     dd.col <- long_palette
     names(dd.col) <- dd
     
     h1 <- ggmuller::Muller_plot(Muller_df, colour_by = "col_index", palette = dd.col, xlab ="Number of cells")
     h1<-h1+ylim(0,1)
     h2 <- Muller_pop_plot(Muller_df, colour_by = "col_index", palette = dd.col, xlab ="Number of cells")
     
     png(paste0(output_dir_plots_WaitCorrelation_biopsy, 
                paste0("/", "MullerPlot_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest, "seedFolder", seedFolder, "_NumCells.png")),  width = 1000, height = 1000, res = 100)
     grid.arrange(h1, h2, nrow=2)
     dev.off()
     
     
     
#search for the simulation with low meanBirthRate and large driverDiversity
     as.data.frame(subsetDataForScatterPlot[which(subsetDataForScatterPlot$DriverDiversity==max(subsetDataForScatterPlot$DriverDiversity)),])
     
     seedOfInterest=as.numeric(subsetDataForScatterPlot[which(subsetDataForScatterPlot$DriverDiversity==max(subsetDataForScatterPlot$DriverDiversity)),
                                                        which(colnames(subsetDataForScatterPlot)=="seed")])
     
     pathToData<-"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_0/mu_driver_birth_0/s_driver_birth_3"
     
     for(seedIt in c(0:99)){
       tmpParam<-read_delim(paste0(pathToData, "/seed_", seedIt, "/parameters.dat"), "\t", escape_double = FALSE, trim_ws = TRUE)
       if(tmpParam$seed==seedOfInterest){
         print(paste0("associated seed ", seedIt))
         seedFolder=seedIt
         break
       }
       
     }
     
     
     plot_all_images(paste0(pathToData, "/seed_", seedFolder),
                     output_filename = paste0("MullerPlot_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest, "seedFolder", seedFolder ),
                     file_type = "png", 
                     output_dir = output_dir_plots_WaitCorrelation_biopsy,
                     trim = -1)
     
     #search for the simulation with lowest meanBirthRate a
     
     seedOfInterest=as.numeric(subsetDataForScatterPlot[which(subsetDataForScatterPlot$MeanBirthRate==min(subsetDataForScatterPlot$MeanBirthRate)),
                                                        which(colnames(subsetDataForScatterPlot)=="seed")])
     
     pathToData<-"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_0/mu_driver_birth_0/s_driver_birth_3"
     
     for(seedIt in c(0:99)){
       tmpParam<-read_delim(paste0(pathToData, "/seed_", seedIt, "/parameters.dat"), "\t", escape_double = FALSE, trim_ws = TRUE)
       if(tmpParam$seed==seedOfInterest){
         print(paste0("associated seed ", seedIt))
         seedFolder=seedIt
         break
       }
       
     }
     
     
     plot_all_images(paste0(pathToData, "/seed_", seedFolder),
                     output_filename = paste0("MullerPlot_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest, "seedFolder", seedFolder ),
                     file_type = "png", 
                     output_dir = output_dir_plots_WaitCorrelation_biopsy,
                     trim = -1)
     
     Muller_df <- muller_df_from_file(paste0(pathToData,"/seed_",seedFolder,"/driver_phylo.dat"))
     #time compute the num of cells
     Muller_df<-group_by(Muller_df, Generation) %>% mutate(Time = sum(Population)) %>% ungroup()
     Muller_df$Generation<-NULL
     #colnames(Muller_df)[which(colnames(Muller_df)=="NumCells")]<-"Time"
     if(class(Muller_df) != "data.frame") return(NA)
     
     long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                       "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                       "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                       "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                       "#5E738F", "#D1A33D")
     dd <- 0:25
     dd.col <- long_palette
     names(dd.col) <- dd
     
     h1 <- ggmuller::Muller_plot(Muller_df, colour_by = "col_index", palette = dd.col, xlab ="Number of cells")
     h1<-h1+ylim(0,1)
     h2 <- Muller_pop_plot(Muller_df, colour_by = "col_index", palette = dd.col, xlab ="Number of cells")
     
     png(paste0(output_dir_plots_WaitCorrelation_biopsy, 
                paste0("/", "MullerPlot_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest, "seedFolder", seedFolder, "_NumCells.png")),  width = 1000, height = 1000, res = 100)
     grid.arrange(h1, h2, nrow=2)
     dev.off()
     
     
######### Correlation of Combined dataset ######### 
######### Combined by mu_driver_birth
Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)
     
cols_list<-colnames(Allsummary)[grepl("DriverDiversity",colnames(Allsummary))]
     
# ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
#                            paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
#                            paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
#                            paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10))
     
#cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]
     
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
     
#the following lines generates summary dataframe of correlations with "waiting_time"
min_count=20
     
wait_cor_summary_CombinedMutationRate<-get_wait_cor_summary_CombinedMutationRate(Allsummary, cols_list, 
                                               num_parameters = num_parameters, min_count = min_count)

fwrite(wait_cor_summary_CombinedMutationRate, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate.csv"))

#### add smaller start size range
start_size_range<-62500
min_count=20
gap_range <- (1:100)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
# data come from "dataAugmented"
#col_list are the same as the ones used for Allsummary
Allsummary_SS62500=vector(mode="list", length=3)
for(Kiter in c(1:3)){
  
  print(paste("working on Kiter=",Kiter ))
  
  FilteredData<- data[[Kiter]]
  
  Allsummary_SS62500[[Kiter]] <- get_summary(FilteredData, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  print("summary computed")
  
  
   
}

Allsummary_SS62500<-do.call(rbind, Allsummary_SS62500)
Allsummary_withSS62500<-rbind(Allsummary_SS62500,Allsummary )
fwrite(Allsummary_withSS62500, file = paste0(output_dir_data, "/Allsummary_withSS62500.csv"))



wait_cor_summary_CombinedMutationRate_withSS62500<-get_wait_cor_summary_CombinedMutationRate(Allsummary_withSS62500, cols_list, 
                                                                                        num_parameters = num_parameters, min_count = min_count)
fwrite(wait_cor_summary_CombinedMutationRate_withSS62500, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500.csv"))

 
wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved<-get_wait_cor_summary_CombinedMutationRate(Allsummary_withSS62500, cols_list, 
                                                                                             num_parameters = num_parameters, min_count = min_count)
fwrite(wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved.csv"))


wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved<-get_wait_cor_summary_CombinedFitnessEffect(Allsummary_withSS62500, cols_list, 
                                                                                                       num_parameters = num_parameters, min_count = min_count)
fwrite(wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved, file = paste0(output_dir_data, "/wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved.csv"))







# Creation of wait_cor summary for smaller start size with the pvalues
Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)

cols_list<-colnames(Allsummary)[grepl("DriverDiversity",colnames(Allsummary))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10))

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]

cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")


#the following lines generates summary dataframe of correlations with "waiting_time"
min_count=20

#for(Kiter in c(1:3)){
Kiter=2
  SummarySmall_Tmp<-read_csv(paste0(output_dir_data, "/summary_smallerStartSizeK",Kiter, ".csv" ))
  SummaryMedium_Tmp<-read_csv(paste0(output_dir_data, "/summary_mediumStartSizeK",Kiter, ".csv" ))
  
  wait_cor_summarySmall_RVAideMemoire<-get_wait_cor_summary_modified(SummarySmall_Tmp,cols_list, 
                                                                     num_parameters = num_parameters, min_count = min_count)
  wait_cor_summaryMedium_RVAideMemoire<-get_wait_cor_summary_modified(SummaryMedium_Tmp,cols_list, 
                                                                     num_parameters = num_parameters, min_count = min_count)
  
  fwrite(wait_cor_summarySmall_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summarysmallerStartSizeK",Kiter, "_RVAideMemoire.csv"), row.names = FALSE)
  
  fwrite(wait_cor_summaryMedium_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summarymediumStartSizeK",Kiter, "_RVAideMemoire.csv"), row.names = FALSE)
  
#}


Kiter=2
Allsummary_withSS62500<-read_csv(paste0(output_dir_data, "/Allsummary_withSS62500.csv"), guess_max = 1E4)
SummarySmall_Tmp<-read_csv(paste0(output_dir_data, "/summary_smallerStartSizeK",Kiter, ".csv" ))
SummaryMedium_Tmp<-read_csv(paste0(output_dir_data, "/summary_mediumStartSizeK",Kiter, ".csv" ))


cols_list<-colnames(Allsummary_withSS62500)[grepl("DriverDiversity",colnames(Allsummary_withSS62500))]

ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
                      paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10))

cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]

cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")


#the following lines generates summary dataframe of correlations with "waiting_time"
min_count=20



Svalues<-unique(Allsummary_withSS62500$s_driver_birth)
for(SIter in sort(Svalues) ){
  

  
    SummaryWithSmallerSS_reducedS<-rbind(subset(Allsummary_withSS62500,
                                                 Allsummary_withSS62500$s_driver_birth==SIter),
                                          subset(SummarySmall_Tmp,
                                                 SummarySmall_Tmp$s_driver_birth==SIter))
    SummaryWithSmallerSS_reducedS<-rbind(SummaryWithSmallerSS_reducedS,
                                          subset(SummaryMedium_Tmp,
                                                 SummaryMedium_Tmp$s_driver_birth==SIter))

  wait_cor_summary_WithSmallerSS_reducedS<-get_wait_cor_summary_CombinedMutationRate(SummaryWithSmallerSS_reducedS, cols_list, 
                                                                                      num_parameters = num_parameters, min_count = min_count)
  
  
  fwrite(wait_cor_summary_WithSmallerSS_reducedS, file = paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_S",  SIter, ".csv"))
  
  
  
}
for(SIter in sort(Svalues) ){

  if(SIter==min(Svalues)){
    paste("first loop")
    wait_cor_summary_WithSmallerSS_CombinedMutationRate<-read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_S",  SIter, ".csv"), guess_max = 1E4)

  }else{
    wait_cor_summary_WithSmallerSS_CombinedMutationRate<-rbind(wait_cor_summary_WithSmallerSS_CombinedMutationRate,
                                                                         read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_S",  SIter, ".csv"), guess_max = 1E4)
    )
  }

}

fwrite(wait_cor_summary_WithSmallerSS_CombinedMutationRate, file = paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_CombinedMutationRate_allSvalues.csv"))



Mu<-unique(Allsummary_withSS62500$mu_driver_birth)
for(MuIter in Mu ){
  SummaryWithSmallerSS_reducedMu<-rbind(subset(Allsummary_withSS62500,
                                               Allsummary_withSS62500$mu_driver_birth==MuIter),
                                        subset(SummarySmall_Tmp,
                                               SummarySmall_Tmp$mu_driver_birth==MuIter))
  SummaryWithSmallerSS_reducedMu<-rbind(SummaryWithSmallerSS_reducedMu,
                                        subset(SummaryMedium_Tmp,
                                               SummaryMedium_Tmp$mu_driver_birth==MuIter))
  wait_cor_summary_WithSmallerSS_reducedMu<-get_wait_cor_summary_CombinedFitnessEffect(SummaryWithSmallerSS_reducedMu, cols_list,
                                                                                      num_parameters = num_parameters, min_count = min_count)


  fwrite(wait_cor_summary_WithSmallerSS_reducedMu, file = paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_Mu",  MuIter, ".csv"))



}

for(MuIter in sort(Mu) ){

  if(MuIter==min(Mu)){
    paste("first loop")
    wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect<-read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_Mu",  MuIter, ".csv"), guess_max = 1E4)

  }else{
    wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect<-rbind(wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect,
                                                                         read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_Mu",  MuIter, ".csv"), guess_max = 1E4)
    )
  }

}

fwrite(wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect, file = paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect.csv"))




# SummaryWithSmallerSS<-rbind(Allsummary_withSS62500, SummarySmall_Tmp)
# SummaryWithSmallerSS<-rbind(SummaryWithSmallerSS, SummaryMedium_Tmp)
# 
# wait_cor_summary_CombinedMutationRate_WithSmallerSS<-get_wait_cor_summary_CombinedMutationRate(Allsummary_withSS62500, cols_list, 
#                                                                                              num_parameters = num_parameters, min_count = min_count)
# fwrite(wait_cor_summary_CombinedMutationRate_withSS62500, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500.csv"))
# 
# 
# wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved<-get_wait_cor_summary_CombinedMutationRate(Allsummary_withSS62500, cols_list, 
#                                                                                                        num_parameters = num_parameters, min_count = min_count)
# fwrite(wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved.csv"))
# 
# 
# wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved<-get_wait_cor_summary_CombinedFitnessEffect(Allsummary_withSS62500, cols_list, 
#                                                                                                          num_parameters = num_parameters, min_count = min_count)
# fwrite(wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved, file = paste0(output_dir_data, "/wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved.csv"))
# 
# 



####### New correlation functions #######

######## Recompute summary and correlations
Kval=512
Kiter=1
SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)
data <- read_csv(paste0(output_dir_data, "/datatestK", Kiter, ".csv"), guess_max = 1E4)

data <- add_columns(data, num_parameters = num_parameters)
data <- add_relative_time(data, start_size = 250000, num_parameters = num_parameters)

data<-mutate(data,
             SimulationNb =0,
             SimulationType=0)


Kiter=2
for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
  
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
fwrite(data, file = paste0(output_dir_data, "/dataAugmentedK", Kval, "_NewCorrelation.csv"), row.names = FALSE)
print(paste0("Wrote ", output_dir_data, "/dataAugmentedK", Kval, "_NewCorrelation.csv to file"))
# data<-read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, "_NewCorrelation.csv"), guess_max = 1E4)

start_size_range<-c(10000,30000,60000,  90000,  125000, 250000, 375000, 500000, 625000, 750000 )
gap_range<-0.1
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
min_count<-20

summary <- get_summary(data, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
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
  cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  
  
  # the following line generates summary dataframe of correlations with "outcome"
  cor_summary <- get_cor_summary(summary, cols_list, num_parameters = num_parameters, min_count = min_count)
  
  #the following lines generates summary dataframe of correlations with "waiting_time"
  wait_cor_summary <- get_wait_cor_summary_modified(summary, cols_list, 
                                           num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  # depth_wait_cor_summary <- get_wait_cor_summary_modified(summary, depth_cols_list, 
  #                                                num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
  # 
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
  

fwrite(cor_summary, file = paste0(output_dir_data, "/Allcor_summary_NewCorrelation_K512.csv"), row.names = FALSE)
fwrite(wait_cor_summary, file = paste0(output_dir_data, "/All_wait_cor_summary_NewCorrelation_K512.csv"), row.names = FALSE)
fwrite(summary, file = paste0(output_dir_data, "/Allsummary_NewCorrelation_K512.csv"), row.names = FALSE)





wait_cor_summary_CombinedFitnessEffectK512<-get_wait_cor_summary_CombinedFitnessEffect(summary, cols_list,
                                                                                     num_parameters = num_parameters, min_count = min_count)

fwrite(wait_cor_summary_CombinedFitnessEffectK512, file = paste0(output_dir_data, "/wait_cor_summary_CombinedFitnessEffect_K512_NewCorrelation.csv"))

wait_cor_summary_CombinedMutationRateK512<-get_wait_cor_summary_CombinedMutationRate(summary, cols_list,
                                                                                   num_parameters = num_parameters, min_count = min_count)

fwrite(wait_cor_summary_CombinedMutationRateK512, file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_K512_NewCorrelation.csv"))






######  correlation with growth rate ######  

Allsummary_NewCorrelation_K512<-read_csv(file = paste0(output_dir_data, "/Allsummary_NewCorrelation_K512.csv"), guess_max = 1E4 )

final_size <- 1E6
min_count<-20
Allsummary_NewCorrelation_K512<-mutate(Allsummary_NewCorrelation_K512, 
                                       AverageGrowthRate = (final_size-start_size) / waiting_time )


cols_list<-colnames(Allsummary_NewCorrelation_K512)[grepl("DriverDiversity",colnames(Allsummary_NewCorrelation_K512))]

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
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")


#the following lines generates summary dataframe of correlations with "waiting_time"
wait_cor_summary <- get_AverageGrowthRate_cor_summary_modified(Allsummary_NewCorrelation_K512, cols_list, 
                                                  num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells

wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)

for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
  
  rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
                     (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                     (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
  wait_cor_summary[rowIndex, 
                   which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
  
  
}


fwrite(wait_cor_summary, file = paste0(output_dir_data, "/All_AverageGrowthRate_cor_summary_NewCorrelation_K512.csv"), row.names = FALSE)
fwrite(Allsummary_NewCorrelation_K512, file = paste0(output_dir_data, "/Allsummary_NewCorrelation_withAverageGrowthRate_K512.csv"), row.names = FALSE)







##### 
Allsummary_NewCorrelation_K512<-read_csv(paste0(output_dir_data, "/Allsummary_NewCorrelation_withAverageGrowthRate_K512.csv"), guess_max = 1E4)

MeanSelectiveSweep<-Allsummary_NewCorrelation_K512 %>% group_by(SimulationType) %>% mutate(GeneralMeanSweep = mean(mean_autocor, na.rm=TRUE),
                                                                                           GeneralMedianSweep = median(mean_autocor, na.rm=TRUE)
                                                                                           )
MeanSelectiveSweep<-unique(MeanSelectiveSweep[, colnames(MeanSelectiveSweep) %in% c("K", "s_driver_birth","mu_driver_birth", "GeneralMeanSweep", "GeneralMedianSweep" )])



######  correlation meanSweep value with growth rate ######  

Allsummary_NewCorrelation_K512<-read_csv(paste0(output_dir_data, "/Allsummary_NewCorrelation_withAverageGrowthRate_K512.csv"), guess_max = 1E4)

final_size <- 1E6
min_count<-20
Kiter<-2
cols_list<-"mean_autocor"


#the following lines generates summary dataframe of correlations with "waiting_time"
wait_cor_summary <- get_AverageGrowthRate_cor_summary_modified(Allsummary_NewCorrelation_K512, cols_list, 
                                                               num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells

wait_cor_summary<-mutate(wait_cor_summary,  SimulationType=0)

for(i in c((1200*(Kiter-1)+1): (1200*Kiter) ) ){
  
  rowIndex<-which( (wait_cor_summary$K== SimulationsDefinition$K[i]) &
                     (wait_cor_summary$s_driver_birth== SimulationsDefinition$s_driver_birth[i]) &
                     (wait_cor_summary$mu_driver_birth==SimulationsDefinition$mu_driver_birth[i]))
  wait_cor_summary[rowIndex, 
                   which(colnames(wait_cor_summary)=="SimulationType")]<-SimulationsDefinition$SimulationType[i]
  
  
}


fwrite(wait_cor_summary, file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_NewCorrelation_K512.csv"), row.names = FALSE)



wait_cor_summary_CombinedFitnessEffectK512<-get_AverageGrowthRate_cor_summary_CombinedFitnessEffect(Allsummary_NewCorrelation_K512, cols_list,
                                                                                       num_parameters = num_parameters, min_count = min_count)

fwrite(wait_cor_summary_CombinedFitnessEffectK512, file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_CombinedFitnessEffect_K512_NewCorrelation.csv"))

wait_cor_summary_CombinedMutationRateK512<-get_AverageGrowthRate_cor_summary_CombinedMutationRate(Allsummary_NewCorrelation_K512, cols_list,
                                                                                     num_parameters = num_parameters, min_count = min_count)

fwrite(wait_cor_summary_CombinedMutationRateK512, file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_CombinedMutationRate_K512_NewCorrelation.csv"))



test<-Allsummary_NewCorrelation_K512[which(Allsummary_NewCorrelation_K512$mu_driver_birth==1e-04 & Allsummary_NewCorrelation_K512$start_size ==750000 & Allsummary_NewCorrelation_K512$s_driver_birth==0.15),]

####### Random mu #######

start_size_range<-c(125000)
gap_range<-0.1
final_size <-c(250000, 375000, 500000, 625000, 750000,875000,  1E6) # waiting time is measured until tumour reaches this NumCells value
min_count<-20
data<-read_csv(paste0("/data/Batch_RandomMu_ForecastingPaper_Working", "/data_BatchRandomMuAugmentedK512.csv"), guess_max = 1E4)

Summary_allFinalSize<-vector(mode="list", length= length(final_size))

for(FinSizeIt in seq_along(final_size) ){
  
  Summary_allFinalSize[[FinSizeIt]] <- get_summary(data, start_size_range, gap_range, final_size[FinSizeIt], num_parameters = num_parameters) %>% mutate(FinalSize =final_size[FinSizeIt])
  
}
Summary_allFinalSize<-do.call(rbind, Summary_allFinalSize)

Summary_allFinalSize<-mutate(Summary_allFinalSize, 
                             AverageGrowthRate = (FinalSize-start_size) / waiting_time, 
                             InverseWaitingTime = 1/ waiting_time
                             )

fwrite(Summary_allFinalSize, file = paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000.csv"))

Summary_allFinalSize<-read_csv(paste0(output_dir_data, "/Summary_allFinalSize_fromStartSize125000.csv"), guess_max = 1E4)

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
cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers", "mean_autocor", "mu_driver_birth")


AverageGrowthRate_cor_summary_CombinedMutationRateK512<-get_AverageGrowthRate_cor_summary_CombinedMutationRate_VaryingFinalSize(Summary_allFinalSize, cols_list,
                                                                                                               num_parameters = num_parameters, min_count = min_count)
fwrite(AverageGrowthRate_cor_summary_CombinedMutationRateK512, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRate_K512_RandomMu.csv"))


InverseWaitingTime_cor_summary_CombinedMutationRateK512<-get_InverseWaitingTime_cor_summary_CombinedMutationRate(Summary_allFinalSize, cols_list,
                                                                                                                 num_parameters = num_parameters, min_count = min_count)
fwrite(InverseWaitingTime_cor_summary_CombinedMutationRateK512, file = paste0(output_dir_data, "/InverseWaitingTime_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRate_K512_RandomMu.csv"))


######## Recompute summary and correlations for all K ######## 
Kvalue=c(64, 4096)
# Kiter=c(0,3)
SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)
SimulationsDefinition<-as.data.frame(SimulationsDefinition)

for(Kit in seq_along(c(0,2)) ){
  
  Kiter = c(0,2)[Kit]
  Kval=Kvalue[Kit]
  data <- read_csv(paste0(output_dir_data, "/datatestK", Kiter, ".csv"), guess_max = 1E4)
  
  data <- add_columns(data, num_parameters = num_parameters)
  data <- add_relative_time(data, start_size = 250000, num_parameters = num_parameters)
  
  SimulationsDefinition_tmp<-SimulationsDefinition[c((1200*(Kiter)+1): (1200*(Kiter+1) ) ) , ]
  
  data<-merge(data,SimulationsDefinition_tmp, by=c("K", "mu_driver_birth", "s_driver_birth", "seed") )
  
  #as it takes a while, save the modification
  fwrite(data, file = paste0(output_dir_data, "/dataAugmentedK", Kval, "_NewCorrelation.csv"), row.names = FALSE)
  print(paste0("Wrote ", output_dir_data, "/dataAugmentedK", Kval, "_NewCorrelation.csv to file"))
  
  start_size_range<-c(10000,30000,60000,  90000,  125000, 250000, 375000, 500000, 625000, 750000 )
  gap_range<-0.1
  final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value
  min_count<-20
  
  if(Kiter==2){
    #those simulations don't achieve a final size larger than 6250000
    data<-data %>%filter( ! (SimulationNb %in% c(2499,2460)))
    #those simulations don't achieve a final size larger thant 7500000
    data<-data %>%filter( ! (SimulationNb %in% c(2414, 2404)))
    
  }
  
  summary <- get_summary(data, start_size_range, gap_range, final_size, num_parameters = num_parameters) 
  
  summary<-mutate(summary,
                  AverageGrowthRate = (final_size-start_size) / waiting_time, 
                  InverseWaitingTime = 1/ waiting_time)
  
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
  cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
  cols_list<-cols_list[! cols_list %in% "mu_driver_birth" ]
  
  
  #the following lines generates summary dataframe of correlations with "waiting_time"
  wait_cor_summary <- get_wait_cor_summary_modified(summary, cols_list, 
                                                    num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  AverageGrowthRate_cor_summary<-get_AverageGrowthRate_cor_summary_modified(summary, cols_list, 
                                                                            num_parameters = num_parameters, min_count = min_count) # summary dataframe of correlations with "waiting_time", including all cells
  
  AverageGrowthRate_cor_summary<-merge(AverageGrowthRate_cor_summary, unique(SimulationsDefinition_tmp[, colnames(SimulationsDefinition_tmp) %in% c("K", "mu_driver_birth", "s_driver_birth", "SimulationType")]), by=c("K", "mu_driver_birth", "s_driver_birth" ))
  
  wait_cor_summary<-merge(wait_cor_summary, unique(SimulationsDefinition_tmp[, colnames(SimulationsDefinition_tmp) %in% c("K", "mu_driver_birth", "s_driver_birth", "SimulationType")]), by=c("K", "mu_driver_birth", "s_driver_birth"))
  
  fwrite(wait_cor_summary, file = paste0(output_dir_data, "/waitingTime_cor_summary_NewCorrelation_K",Kval,".csv"), row.names = FALSE)
  fwrite(summary, file = paste0(output_dir_data, "/Summary_NewCorrelation_K",Kval,".csv"), row.names = FALSE)
  fwrite(AverageGrowthRate_cor_summary, file = paste0(output_dir_data, "/AverageGrowthRate_cor_summary_NewCorrelation_K",Kval,".csv"), row.names = FALSE)
  
}


# K=512
# 
# Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)
# 
# cols_list<-colnames(Allsummary)[grepl("DriverDiversity",colnames(Allsummary))]
# 
# ToRemoveFomColList<-c(paste0("DriverDiversityFrom1SamplesAtDepth", 5:10),
#                       paste0("DriverDiversityFrom4SamplesAtDepth", 5:10),
#                       paste0("DriverDiversityFrom1BigSamplesAtDepth", 5:10),
#                       paste0("DriverDiversityFrom4BigSamplesAtDepth", 5:10))
# 
# cols_list<-cols_list[! cols_list %in% ToRemoveFomColList ]
# 
# cols_list <-c(cols_list, "DriverEdgeDiversity", "MeanBirthRate", "Drivers")
# 
# remove(Allsummary)
# 
# #the following lines generates summary dataframe of correlations with "waiting_time"
# min_count=20
# 
# 
# Kiter=2
# ## first Large Start size
# Allsummary_withSS62500<-read_csv(paste0(output_dir_data, "/Allsummary_withSS62500.csv"), guess_max = 1E4)
# Allsummary_withSS62500<-subset(Allsummary_withSS62500, 
#                                Allsummary_withSS62500$K==512)
# 
# wait_cor_summaryLarger_RVAideMemoire<-get_wait_cor_summary_modified(Allsummary_withSS62500,cols_list, 
#                                                                     num_parameters = num_parameters, min_count = min_count)
# fwrite(wait_cor_summaryLarger_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summaryLargerStartSizeK",K, "_RVAideMemoire_NewCorrelation.csv"), row.names = FALSE)
# 
# remove(Allsummary_withSS62500)
# 
# # small SS
# SummarySmall_Tmp<-read_csv(paste0(output_dir_data, "/summary_smallerStartSizeK",Kiter, ".csv" ))
# wait_cor_summarySmall_RVAideMemoire<-get_wait_cor_summary_modified(SummarySmall_Tmp,cols_list, 
#                                                                    num_parameters = num_parameters, min_count = min_count)
# fwrite(wait_cor_summarySmall_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summarysmallerStartSizeK",K, "_RVAideMemoire_NewCorrelation.csv"), row.names = FALSE)
# 
# remove(SummarySmall_Tmp)
# 
# # Medium SS
# SummaryMedium_Tmp<-read_csv(paste0(output_dir_data, "/summary_mediumStartSizeK",Kiter, ".csv" ))
# 
# wait_cor_summaryMedium_RVAideMemoire<-get_wait_cor_summary_modified(SummaryMedium_Tmp,cols_list, 
#                                                                     num_parameters = num_parameters, min_count = min_count)
# 
# fwrite(wait_cor_summaryMedium_RVAideMemoire, file = paste0(output_dir_data, "/wait_cor_summarymediumStartSizeK",K, "_RVAideMemoire_NewCorrelation.csv"), row.names = FALSE)
# 
# remove(SummaryMedium_Tmp)

AllSimulationNb<-unique(summary[which(summary$start_size==10000), which(colnames(summary) =="SimulationNb")])
setdiff(AllSimulationNb,unique(summary[which(summary$start_size==750000), which(colnames(summary) =="SimulationNb")]))

test<-data[which(data$SimulationNb %in% c(2499)),] %>% group_by(SimulationNb)
max(test$NumCells)
test<-data[which(data$SimulationNb %in% c(2460)),] %>% group_by(SimulationNb)
max(test$NumCells)
test<-data[which(data$SimulationNb %in% c(2414)),] %>% group_by(SimulationNb)
max(test$NumCells)
test<-data[which(data$SimulationNb %in% c(2404)),] %>% group_by(SimulationNb)
max(test$NumCells)





 ######### Correlation between Clonal diversity  and mean birth rate ######### 


  
Kvalue=c(64,512, 4096)
#Kiter=c(0,3)
SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)
SimulationsDefinition<-as.data.frame(SimulationsDefinition)
min_count<-20

for(Kit in seq_along(c(0:2)) ){
  
  Kiter = c(0:2)[Kit]
  Kval=Kvalue[Kit]
  
  if(Kval==512){
    summary<-read_csv(file = paste0(output_dir_data, "/Allsummary_NewCorrelation_withAverageGrowthRate_K512.csv"), guess_max = 1E4)
 
  }else{
    summary<-read_csv(file = paste0(output_dir_data, "/Summary_NewCorrelation_K",Kval,".csv"), guess_max = 1E4)
    
  }
    
  ClonalDiversityVsMeanBirthRate_correlation<-get_MeanBirthRate_correlation(summary, "DriverDiversity", num_parameters, min_count)
  
  SimulationsDefinition_tmp<-SimulationsDefinition[c((1200*(Kiter)+1): (1200*(Kiter+1) ) ) , ]
  
  ClonalDiversityVsMeanBirthRate_correlation<-merge(ClonalDiversityVsMeanBirthRate_correlation, unique(SimulationsDefinition_tmp[, colnames(SimulationsDefinition_tmp) %in% c("K", "mu_driver_birth", "s_driver_birth", "SimulationType")]), by=c("K", "mu_driver_birth", "s_driver_birth" ))
  
  fwrite(ClonalDiversityVsMeanBirthRate_correlation, file = paste0(output_dir_data, "/ClonalDiversityVsMeanBirthRate_correlation_K",Kval,".csv"), row.names = FALSE)
  
}


