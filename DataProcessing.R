library(demonanalysis)
library(readr)
library(dplyr)
library(data.table)
library(moments)

input_dir <- "all_results/Batch_TTT_100Seeds" # folder containing results of the batch
input_dir_data <- "data/Batch_TTT_100Seeds" # folder containing data files
output_dir_data <- "forecasting_analysis"

n_cores <- NA

all_output_Rob <- function(input_dir, include_diversities = TRUE, df_type = "output", max_generation = FALSE, vaf_cut_off = NA, generation = NA, numcells = NA, n_cores = NA,ExitCode4=FALSE) {
  
  df_type_list <- c("output", "driver_genotype_properties", "genotype_properties", 
                    "allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts", "driver_phylo",
                    "diversities")
  stopifnot(df_type %in% df_type_list)
  
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_check <- function(x, res) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) print(paste0(full_dir, " ", msg), quote = FALSE)
    else print("Failed")
  }
  apply_combinations(final_values, each_check)
  
  print("Finished checking", quote = FALSE)
  
  num_parameters <- count_parameters(input_dir)
  each_df <- function(x) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    print(paste0(full_dir, " ", msg), quote = FALSE)
    
    #Modifications to deal with treatment => use combine_dfs with argument ExitCode4
    if(!identical(msg, character(0))){
      
      if(msg == "Exit code 0"){
        print("Starting code for Exit code 0")
        return(combine_dfs_Rob(full_dir, include_diversities, 
                           df_type, max_generation, vaf_cut_off, generation, numcells, num_parameters, ExitCode4=FALSE))
        
      }else if (msg == "Exit code 4"){
        # if some of the simulations got "Exit code 4" and that we want to include them into the analysis, i.e by using ExitCode4=TRUE,
        # then we need to call combine_dfs with argument ExitCode4=TRUE 
        print("Starting code for Exit code 4")
        if(ExitCode4){
          return(combine_dfs_Rob(full_dir, include_diversities, 
                             df_type, max_generation, vaf_cut_off, generation, numcells, num_parameters, ExitCode4=TRUE ))
        }
        
      }
    }
    return(data.frame())
  }
  
  print(paste("final_values = ", final_values))
  final_values[3] <- 2 ### REMOVE THIS!!!
  print(paste("final_values = ", final_values))
  
  if(is.na(n_cores)){
    intermediate <- apply_combinations_Rob(final_values, each_df)
    print("Intermediate result:")
    print(class(intermediate))
    print(length(intermediate))
    print(head(intermediate))
    res <- rbindlist(intermediate)
  } else {
    res <- rbindlist(apply_combinations_parallel(n_cores, final_values, each_df))
  }
  
  # report seed counts:
  print("Number of seeds:", quote = FALSE)
  print(count_seeds(res, num_parameters))
  
  return(res)
}

apply_combinations_Rob <- function(vec, fn, ...){
  print("Starting")
  vecs <- mapply(seq, 0, vec, SIMPLIFY = FALSE) # a list of n sequences, where n = length(vec)
  print(length(vecs))
  tmp <- do.call(expand.grid, vecs) # a data frame where each row is a permuation of values from the n sequences
  print(dim(tmp))
  print(head(tmp))
  res <- apply(tmp, 1, fn, ...) # the result of applying fn to each row of tmp
  print("Done")
  return(res)
}

combine_dfs_Rob <- function(full_dir, include_diversities = TRUE, df_type = "output", max_generation = FALSE, vaf_cut_off = NA, generation = NA, numcells = NA, num_parameters = NA, ExitCode4=FALSE) {
  
  print("Starting combine_dfs_Rob")
  
  if(substr(full_dir, nchar(full_dir), nchar(full_dir)) == "/") full_dir <- substr(full_dir, 1, nchar(full_dir) - 1)
  
  file_pars <- paste0(full_dir, "/parameters.dat")
  file_out <- paste0(full_dir, "/output.dat")
  file_div <- paste0(full_dir, "/output_diversities.dat")
  file_driver_phylo <- paste0(full_dir, "/driver_phylo.dat")
  file_allele_counts <- paste0(full_dir, "/output_allele_counts.dat")
  
  df_out <- fread(file_out)
  
  if(ExitCode4){
    
    #first check that the new data frame as the column "treated", else this will not work
    if(! "Treated" %in% colnames(df_out)){
      warning(paste0("file ",file_out, " does not contain the variable Treated, but ExitCode4=TRUE, which is incompatible !" ))
    }
    df_out<-subset(df_out, df_out$Treated==0) #remove the last line corresponding to the update once no cells are remaining
  }
  
  df_pars <- fread(file_pars)
  
  print("End of first part")
  
  if (df_type == "output"){
    # procedure for 'traditional' df_out (output.dat)
    if(include_diversities){
      
      print("include_diversities is TRUE")
      
      df_div <- fread(file_div)
      
      if(ExitCode4){
        
        #first check that the data frame as the column "treated", else this will not work
        if(! "Treated" %in% colnames(df_div)){
          warning(paste0("file ",file_div, " does not contain the variable Treated, but ExitCode4=TRUE, which is incompatible !" ))
        }
        
        df_div<-subset(df_div, df_div$Treated==0) #remove the last line corresponding to the update once no cells are remaining
      }
    } 
    
    df_driver_phylo <- fread(file_driver_phylo)
    
    print("df_div:")
    print(df_div[1:12, 1:5])
    print("End of second part")
    
    if(ExitCode4){
      
      #first check that the data frame as the column "treated", else this will not work
      if(! "Treated" %in% colnames(df_driver_phylo)){
        warning(paste0("file ",file_driver_phylo, " does not contain the variable Treated, but ExitCode4=TRUE, which is incompatible !" ))
      }
      
      df_driver_phylo<-subset(df_driver_phylo, df_driver_phylo$Treated==0)#remove the last lines corresponding to the update once no cells are remaining
    }
    
    
    df_driver_phylo <- filter(df_driver_phylo, CellsPerSample == -1, NumSamples == 1, SampleDepth == -1)
    df_driver_phylo <- df_driver_phylo[!duplicated(df_driver_phylo), ]
    
    print("About to do get_population_df")
    
    pop_df <- get_population_df_Rob(df_driver_phylo) #should run with the version of get_population_df dealing with treatment.
    
    print("Done get_population_df")
    
    if(include_diversities){
      
      print("include_diversities is TRUE again")
      print("df_div:")
      print(df_div[1:12, 1:5])
      print("df_out:")
      print(df_out[1:12, 1:5])
      
      temp <- merge(df_out, df_div, all = TRUE)
      
      if("Treated" %in% colnames(temp)) {
        temp<-temp[order(Generation, Treated), ]#ensure that the output file is correctly ordered
      } else {
        temp<-temp[order(Generation), ]#ensure that the output file is correctly ordered
      }
      print("temp after merging and reordering:")
      print(temp[1:12, 1:5])
      
    }else{ 
      temp <- df_out
    }
    
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
    
    # adds maxgen and gen_adj columns
    if(is.na(num_parameters)) num_parameters <- count_parameters(full_dir)
    temp <- add_columns(temp, num_parameters)  #should run with the version of add_columns dealing with treatment.
    
    # add sweep_seq columns (specific for output.dat?)
    sweep_seq <- sweep_sequence(pop_df, lag_type = "proportions", breaks = 10)
    temp <- mutate(temp, mean_autocor = mean(sweep_seq), 
                   log_mean_autocor = log(mean(sweep_seq)), 
                   sqrt_mean_autocor = sqrt(mean(sweep_seq)), 
                   skewness = skewness(sweep_seq))
  } else if (df_type %in% c("allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts", "diversities")){
    temp <- fread(paste0(full_dir, "/output_", df_type, ".dat"))
    
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else if (df_type %in% c("driver_genotype_properties", "genotype_properties")){
    temp <- fread(paste0(full_dir, "/output_", df_type, ".dat"))
    # data contains columns Descendants (or AlleleCount in older versions) and pop_size
    colnames(temp)[colnames(temp) == "AlleleCount"] <- "Descendants"
    calc_VAF <- function(data){
      alpha <- 1
      coverage <- data$pop_size * alpha
      VAF <- data$Descendants / coverage
      return(VAF)
    }
    temp$pop_size <- (df_out %>% filter(Generation == max(Generation)) %>% select(NumCells) %>% slice(1))$NumCells
    temp$VAF <- calc_VAF(temp)
    if(!is.na(vaf_cut_off)) {
      temp <- temp %>% filter(VAF >= vaf_cut_off | Population > 0)
    }
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else if (df_type %in% c("driver_phylo")){
    
    df_driver_phylo <- fread(file_driver_phylo)
    
    if(ExitCode4){
      
      #first check that the data frame as the column "treated", else this will not work
      if(! "Treated" %in% colnames(df_driver_phylo)){
        warning(paste0("file ",file_driver_phylo, " does not contain the variable Treated, but ExitCode4=TRUE, which is incompatible !" ))
      }
      
      df_driver_phylo<-subset(df_driver_phylo, df_driver_phylo$Treated==0)#remove the last lines corresponding to the update once no cells are remaining
      
    }
    
    df_driver_phylo <- filter(df_driver_phylo, CellsPerSample == -1, NumSamples == 1, SampleDepth == -1)
    df_driver_phylo <- df_driver_phylo[!duplicated(df_driver_phylo), ]
    temp <- df_driver_phylo
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else {
    stop("no valid df_type argument was passed")
  }
  
  # filter if requested:
  temp <- filter_by_generation_or_numcells(temp, full_dir, generation, numcells)
  # only use last generation
  if(max_generation){
    temp <- temp %>% filter(Generation == max(Generation))
  }
  
  print(paste0("Result of combine_dfs has dimensions ", dim(temp)[1], " x ", dim(temp)[2]), quote = FALSE)
  
  return(temp)
}

get_population_df_Rob <- function(df) {
  
  original_colname <- "Generation"
  # rename Time column (original name will be restored later):
  if("Time" %in% colnames(df) && !("Generation" %in% colnames(df))) {
    colnames(df)[colnames(df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  
  # check column names:
  if(!("Generation" %in% colnames(df)) | !("Identity" %in% colnames(df)) | !("Population" %in% colnames(df))) 
    stop("colnames(df) must contain Generation (or Time), Identity and Population")
  
  . <- NULL # avoid check() note
  Population <- NULL # avoid check() note
  
  max_gen <- max(df$Generation) # final generation
  max_gen_ids <- filter_(df, ~Generation == max_gen)$Identity # vector containing all identities at final generation
  df <- filter_(df, ~Identity %in% max_gen_ids) # filter df to include only identities present at final generation
  n <- length(unique(df$Identity)) # number of unique identities in df after filtering
  
  print("Mark 1")
  print(c(length(df$Generation), n))
  
  InteractionGenerationNumCells <- unique(df[,c('Generation','NumCells')])
  
  # data frame containing all combinations of generations, NumCells and identities
  # Should also deal with the case where max(res$Generation) is achieved for 2 values of NumCells (ex 999991 and 10e6)
  master <- data.frame(Generation = rep(InteractionGenerationNumCells$Generation, 
                                        each = n),
                       NumCells =rep(InteractionGenerationNumCells$NumCells, 
                                     each = n),
                       Identity = unique(df$Identity)) 
  
  print(dim(master))
  print(dim(df))
  
  res <- left_join(master, df, by = c("Generation","NumCells",  "Identity")) %>% 
    mutate(Population = ifelse(Population %in% NA, 0, Population))
  
  print("Mark 2")
  
  cols <- colnames(df)[!(colnames(df) %in% c("Generation", "NumCells","Identity", "Population"))]
  
  # Deal with some problems occuring when the dimension of the columns are not compatibles.
  #warnings are given if the first case (taht should always happen) don't occur.
  for (col in cols){
    
    # the first case should always happen
    if(  length(res[, col]) %% length(res[res$Generation == max(res$Generation),col]) == 0){
      res[, col] <- res[res$Generation == max(res$Generation), 
                        col]
    }else if(length(unique(res[res$Generation == max(res$Generation),col]))==1){
      res[, col] <- unique(res[res$Generation == max(res$Generation), 
                               col])
      warning(paste0("unique value for ", col, " but pb of modulo which values : ",  length(res[, col]) %% length(res[res$Generation == max(res$Generation),col]), "instead of 0"))
    }else{
      warning(paste0("pb of dimension for column ", col, " modulo values ", length(res[, col]) %% length(res[res$Generation == max(res$Generation),col]) ))
      break
    }
    
  } 
  
  print("Mark 3")
  
  # restore original time column name:
  colnames(res)[colnames(res) == "Generation"] <- original_colname
  
  return(res)
}

print("About to read")

data_old <- read_csv(paste0(input_dir_data, "/dataWithExitCode4.csv"), guess_max = 1E4)

print(paste("data_old:", dim(data_old)))
print(data_old[1:12, ])

remove(data_old)

data <- all_output_Rob(input_dir, include_diversities = TRUE, max_generation = FALSE, 
                   n_cores = n_cores, ExitCode4=TRUE) # combined data for a batch of simulations, excluding diversity columns

print("Created data")
fwrite(data, file = paste0(output_dir_data, "/data.csv"), row.names = FALSE)
print("Wrote data")

data<-data %>% 
  mutate(JustAfterTTT_tmp=(c(0, diff(data$Treated))), 
         JustBeforeTTT_tmp=c(diff(data$Treated), 0)) %>% 
  mutate(JustAfterTTT= ifelse(JustAfterTTT_tmp==1, 1, 0), 
         JustBeforeTTT=ifelse(JustBeforeTTT_tmp==1, 1, 0))

print("Added columns")

#creation of simulationDefinition
SimulationsDefinition<-as.data.frame(data)
SimulationsDefinition<-unique(SimulationsDefinition[, which(colnames(SimulationsDefinition) %in% c("K", "mu_driver_birth", "s_driver_birth", "seed"))])
SimulationsDefinition<-SimulationsDefinition[with(SimulationsDefinition, order(mu_driver_birth,s_driver_birth, seed)), ]
SimulationsDefinition<-SimulationsDefinition %>% mutate(SimulationNumber = seq_along(SimulationsDefinition[, 1]))

print("About to merge")

data <- merge(data, SimulationsDefinition, by= c("K", "mu_driver_birth", "s_driver_birth", "seed"))

print("About to write modified data")
fwrite(data, file = paste0(output_dir_data, "/data_modified.csv"), row.names = FALSE)
print("All finished")


