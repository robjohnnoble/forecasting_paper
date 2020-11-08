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
#library('plot.matrix')
require(lattice)
library(viridis)
##########################  function to add aes and aes_string ########################## 
`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

readSeedFromParamFiles<-function(input_dir){
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  vecs <- mapply(seq, 0, final_values, SIMPLIFY = FALSE) # a list of n sequences, where n = length(vec)
  PathValues <- do.call(expand.grid, vecs) # a data frame where each row is a permuation of values from the n sequences
  
  finalParamFile<-c()
  for(i in c(1:dim(PathValues)[1])){
    #finalParamFile<-rbind(finalParamFile,each_msg(PathValues[i,]) )
    full_dir <- make_dir(input_dir, pars, PathValues[i,])
    tmp<-fread(paste0(full_dir, "parameters.dat"))
    tmp<-mutate(tmp, seedFolder =PathValues[i,dim(PathValues)[2]])
    finalParamFile<-rbind(finalParamFile,tmp)
  }
  
  return(finalParamFile)
}

get_Muller_df_modified<-function (edges, pop_df, cutoff = 0, start_positions = 0.5, threshold = NA,
                                  add_zeroes = NA, smooth_start_points = NA){
  Population <- NULL
  Generation <- NULL
  Identity <- NULL
  original_colname <- "Generation"
  if ("Time" %in% colnames(pop_df) && !("Generation" %in% colnames(pop_df))) {
    colnames(pop_df)[colnames(pop_df) == "Time"] <- "Generation"
    original_colname <- "Time"
  }
  if (dim(pop_df)[1] != length(unique(pop_df$Identity)) * length(unique(pop_df$Generation))) {
    added_rows <- expand.grid(Identity = unique(pop_df$Identity),
                              Generation = unique(pop_df$Generation))
    added_props <- group_by(pop_df, Identity) %>% slice(1) %>%
      ungroup() %>% select(-one_of("Generation", "Population"))
    added_rows <- merge(added_rows, added_props, all = TRUE)
    pop_df <- merge(added_rows, pop_df, all = TRUE)
    pop_df[is.na(pop_df$Population), "Population"] <- 0
    pop_df <- arrange_(pop_df, ~Generation)
    warning("missing population sizes replaced by zeroes")
  }
  if (!missing(add_zeroes)) {
    warning("argument add_zeroes is deprecated (it is now always TRUE).",
            call. = FALSE)
  }
  if (!missing(smooth_start_points)) {
    warning("argument smooth_start_points is deprecated (it is now always TRUE).",
            call. = FALSE)
  }
  if (!missing(threshold)) {
    warning("argument threshold is deprecated (use cutoff instead, noting that genotypes whose abundance never exceeds the cutoff value are removed, \n            whereas previously genotypes whose abundance never exceeded *twice* the threshold value were removed).",
            call. = FALSE)
    if (missing(cutoff))
      cutoff <- threshold * 2
  }
  if (!("Generation" %in% colnames(pop_df)) | !("Identity" %in%
                                                colnames(pop_df)) | !("Generation" %in% colnames(pop_df)))
    stop("colnames(pop_df) must contain Generation (or Time), Identity and Population")
  if (!is.na(edges)[1]) {
    set1 <- unique(pop_df$Identity)
    set2 <- unique(edges$Identity)
    set3 <- unique(edges$Parent)
    if (length(setdiff(set1, set2)) != 1)
      stop("Identity values in edges must match Identity values in pop_df, excluding the original genotype (which has no parent)")
    if (length(setdiff(set3, set2)) != 1)
      stop("Parent values in edges must also appear as Identity values in edges, excluding the original genotype (which has no parent)")
  }
  if (!is.na(edges)[1]) {
    if ("phylo" %in% class(edges)) {
      collapse.singles(edges)
      edges <- edges$edge
    }
    edges <- na.omit(edges)
    colnames(edges) <- c("Parent", "Identity")
    if (is.factor(edges$Parent))
      edges$Parent <- levels(edges$Parent)[edges$Parent]
    if (is.factor(edges$Identity))
      edges$Identity <- levels(edges$Identity)[edges$Identity]
  }
  pop_df <- add_start_points(pop_df, start_positions)
  pop_df <- arrange_(pop_df, ~-Population)
  pop_df <- arrange_(pop_df, ~Generation)
  lookup <- group_by_(pop_df, ~Identity) %>% filter_(~Population >
                                                       0 | Generation == max(Generation)) %>% slice(1) %>% arrange_(~Generation) %>%
    ungroup()
  lookup <- mutate(lookup, Age = 1:dim(lookup)[1]) %>% select_(~-c(Generation,
                                                                   Population))
  if (is.factor(lookup$Identity))
    lookup$Identity <- levels(lookup$Identity)[lookup$Identity]
  lookup <- select_(lookup, ~c(Identity, Age))
  pop_df <- pop_df %>% group_by_(~Generation) %>% mutate(Frequency = (Population/sum(Population))/2) %>%
    ungroup()
  pop_df$Population <- pop_df$Population/2
  pop_df$Frequency[is.nan(pop_df$Frequency)] <- 0
  Muller_df <- rbind(pop_df, pop_df)
  Muller_df <- arrange_(Muller_df, ~Generation)
  if (!is.na(edges)[1]) {
    edges <- filter_(edges, ~Identity %in% lookup$Identity)
    edges <- left_join(edges, lookup, by = "Identity")
    edges <- select_(edges, ~-Identity)
    colnames(edges) <- c("Parent", "Identity")
    edges <- arrange_(edges, ~Identity)
    colnames(lookup)[1] <- "Parent"
    edges <- left_join(edges, lookup, by = "Parent")
    edges$Parent <- edges$Age
    edges <- select_(edges, ~-Age)
    path <- path_vector(edges)
    path <- rev(path)
    path <- left_join(data.frame(Age = path), lookup, by = "Age")$Parent
  }
  else path <- c(unique(pop_df$Identity)[1], unique(pop_df$Identity)[1])
  
  # filter first for rare genotypes
  if (cutoff > 0) {
    Muller_df <- Muller_df %>% group_by_(~Identity) %>% filter_(~max(Frequency) >=
                                                                  cutoff/2)
    Muller_df <- Muller_df %>% group_by_(~Generation) %>%
      mutate(Frequency = Population/sum(Population)) %>%
      ungroup()
  }
  
  Muller_df <- reorder_by_vector(Muller_df, path)
  
  Muller_df$Group_id <- factor(Muller_df$Group_id, levels = rev(unlist(as.data.frame(Muller_df %>%
                                                                                       filter_(~Generation == max(Generation)) %>% select_(~Group_id)),
                                                                       use.names = FALSE)))
  colnames(Muller_df)[colnames(Muller_df) == "Generation"] <- original_colname
  return(as.data.frame(Muller_df))
}

#' Read a "phylo" dataframe and process it for plotting
#' 
#' @param file file name including path
#' 
#' @return a dataframe formatted for plotting
#' 
#' @export
#' @importFrom readr read_delim
#' @import dplyr
#' @import ggmuller
#' 
#' @examples
#' Muller_df <- muller_df_from_file(system.file("extdata", 
#' "driver_phylo.dat", package = "demonanalysis", mustWork = TRUE))
muller_df_from_file_cutoff <- function(file, cutoff=0) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    return(NA)
  }
  phylo <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
  phylo <- filter(phylo, CellsPerSample == -1)
  edges <- get_edges(phylo)
  if(dim(edges)[1] == 0) edges <- NA
  pop_df <- get_population_df(phylo)
  pop_df <- pop_df %>% mutate(col_index = pop_df$Identity)
  pop_df$col_index[pop_df$col_index > 0] <- pop_df$col_index[pop_df$col_index > 0] %% 25 + 1
  pop_df$col_index <- as.character(pop_df$col_index)
  return(get_Muller_df_modified(edges, pop_df, cutoff=cutoff))
}

#' Create a set of grid and Muller plots
#' 
#' @param path folder containing the input files
#' @param trim how many rows and columns to remove from grids; if trim < 0 (default) then all rows and columns containing NA are removed
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' 
#' @return either an image file or a plot displyed on screen
#' 
#' @export
#' @import ggmuller
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' 
#' @examples
#' plot_all_images(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
plot_all_images_cutoff <- function(path, output_filename = NA, file_type = "png", output_dir = NA, trim = -1, include_genotype_plots = TRUE, cutoff=0) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  Muller_df <- muller_df_from_file_cutoff(paste0(path, "driver_phylo.dat"), cutoff = cutoff)
  if(class(Muller_df) != "data.frame") return(NA)
  
  long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                    "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                    "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                    "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                    "#5E738F", "#D1A33D")
  dd <- 0:25
  dd.col <- long_palette
  names(dd.col) <- dd
  
  b_grid <- image_df_from_grid_file(paste0(path, "output_birthratesgrid.dat"), trim)
  
  if(include_genotype_plots) {
    h1 <- Muller_plot(Muller_df, colour_by = "col_index", palette = dd.col)
    h2 <- Muller_pop_plot(Muller_df, colour_by = "col_index", palette = dd.col)
  }
  min_birth_rate <- min(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  max_birth_rate <- max(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  h3 <- Muller_plot(Muller_df, colour_by = "BirthRate", add_legend = FALSE) + 
    scale_fill_distiller(palette = "RdBu", direction = -1, 
                         limits = c(min_birth_rate, max_birth_rate)) + 
    theme(line = element_blank(), rect = element_blank()) + 
    scale_x_continuous(breaks = c(0, round(max(Muller_df$Generation)))) + 
    scale_y_continuous(breaks = NULL)
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_driversgrid.dat"), trim)
    image_df[which(image_df$z > 0), "z"] <- as.character(image_df[which(image_df$z > 0), "z"] %% 25 + 1)
    g1 <- grid_plot(image_df, palette = dd.col, discrete = TRUE)
  }
  
  g2 <- grid_plot(b_grid, add_legend = TRUE, legend_title = "Mean cell\nproliferation rate   ") + 
    scale_fill_distiller(name = "Mean cell\nproliferation rate   ", palette ="RdBu", 
                         direction = -1, na.value="white", 
                         limits = c(min_birth_rate, max_birth_rate))
  
  image_df <- image_df_from_grid_file(paste0(path, "output_passengersgrid.dat"), trim)
  g3 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Mean passenger\nmutations per cell")
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_popgrid.dat"), trim)
    image_df[which(image_df$z == 0), "z"] <- NA
    g4 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Tumour cells\nper gland")
  }
  
  if(!is.na(output_filename)) print(paste0("Created all plots for file ", output_filename), quote = FALSE)
  
  if(include_genotype_plots) {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
      else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 10)
    }
    lay <- rbind(c(1,1,2),
                 c(3,3,3),
                 c(4,4,5),
                 c(NA,6,7))
    print(grid.arrange(h1, g1, h2, h3, g2, g3, g4, layout_matrix = lay, heights = c(1, 1, 0.75, 0.75)))
  } else {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 180, res = 100)
      else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 1.8)
    }
    lay <- rbind(c(1,2,3))
    print(grid.arrange(h3, g2, g3, layout_matrix = lay))
  }
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}

plot_all_images_cutoff_modified <- function(path, output_filename = NA, file_type = "png", output_dir = NA, trim = -1, include_genotype_plots = TRUE, cutoff=0) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  Muller_df <- muller_df_from_file_cutoff(paste0(path, "driver_phylo.dat"), cutoff = cutoff)
  if(class(Muller_df) != "data.frame") return(NA)
  
  long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                    "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                    "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                    "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                    "#5E738F", "#D1A33D")
  dd <- 0:25
  dd.col <- long_palette
  names(dd.col) <- dd
  
  b_grid <- image_df_from_grid_file(paste0(path, "output_birthratesgrid.dat"), trim)
  
  if(include_genotype_plots) {
    h1 <- Muller_plot(Muller_df, colour_by = "col_index", palette = dd.col)
    
    png(paste0(output_dir,output_filename,"_h1.png"), width = 1000, height = 1000, res = 100)
    print(h1)
    dev.off()
    
    h2 <- Muller_pop_plot(Muller_df, colour_by = "col_index", palette = dd.col)
    
    png(paste0(output_dir,output_filename,"_h2.png"), width = 1000, height = 1000, res = 100)
    print(h2)
    dev.off()
  }
  min_birth_rate <- min(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  max_birth_rate <- max(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  h3 <- Muller_plot(Muller_df, colour_by = "BirthRate", add_legend = FALSE) + 
    scale_fill_distiller(palette = "RdBu", direction = -1, 
                         limits = c(min_birth_rate, max_birth_rate)) + 
    theme(line = element_blank(), rect = element_blank()) + 
    scale_x_continuous(breaks = c(0, round(max(Muller_df$Generation)))) + 
    scale_y_continuous(breaks = NULL)
  
  png(paste0(output_dir,output_filename,"_h3.png"), width = 1000, height = 1000, res = 100)
  print(h3)
  dev.off()
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_driversgrid.dat"), trim)
    image_df[which(image_df$z > 0), "z"] <- as.character(image_df[which(image_df$z > 0), "z"] %% 25 + 1)
    g1 <- grid_plot(image_df, palette = dd.col, discrete = TRUE)
    
    png(paste0(output_dir,output_filename,"_g1.png"), width = 1000, height = 1000, res = 100)
    print(g1)
    dev.off()
  }
  
  g2 <- grid_plot(b_grid, add_legend = TRUE, legend_title = "Mean cell\nproliferation rate   ") + 
    scale_fill_distiller(name = "Mean cell\nproliferation rate   ", palette ="RdBu", 
                         direction = -1, na.value="white", 
                         limits = c(min_birth_rate, max_birth_rate))
  png(paste0(output_dir,output_filename,"_g2.png"), width = 1000, height = 1000, res = 100)
  print(g2)
  dev.off()
  
  image_df <- image_df_from_grid_file(paste0(path, "output_passengersgrid.dat"), trim)
  g3 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Mean passenger\nmutations per cell")
  
  png(paste0(output_dir,output_filename,"_g3.png"), width = 1000, height = 1000, res = 100)
  print(g3)
  dev.off()
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_popgrid.dat"), trim)
    image_df[which(image_df$z == 0), "z"] <- NA
    g4 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Tumour cells\nper gland")
    
    png(paste0(output_dir,output_filename,"_g4.png"), width = 1000, height = 1000, res = 100)
    print(g4)
    dev.off()
  }
  
  
  if(!is.na(output_filename)) print(paste0("Created all plots for file ", output_filename), quote = FALSE)
  
  if(include_genotype_plots) {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
      else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 10)
    }
    lay <- rbind(c(1,1,2),
                 c(3,3,3),
                 c(4,4,5),
                 c(NA,6,7))
    print(grid.arrange(h1, g1, h2, h3, g2, g3, g4, layout_matrix = lay, heights = c(1, 1, 0.75, 0.75)))
  } else {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 180, res = 100)
      else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 1.8)
    }
    lay <- rbind(c(1,2,3))
    print(grid.arrange(h3, g2, g3, layout_matrix = lay))
  }
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}

########################## Setting the folders ########################## 

subfolder_name <- "Batch_ForecastingPaper_bis"

n_cores <- 4


input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files

##output dir
output_dir_plots_paper<-"plots/Batch_ForecastingPaper_bis/WaitCorrelation_plotsPapers"
ifelse(!dir.exists(output_dir_plots_paper), dir.create(output_dir_plots_paper), FALSE)

####### Variables ####### 

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

################# Data loading ################# 

SimulationsDefinition<-read_csv(paste0(output_dir_data, "/SimulationsDefinition.csv"), guess_max = 1E4)
SimulationsDefinition<-as.data.frame(SimulationsDefinition)

#associated variables
MuValues<-unique(SimulationsDefinition$mu_driver_birth)
SValues<-unique(SimulationsDefinition$s_driver_birth)
KValues<-unique(SimulationsDefinition$K)


Allsummary<-read_csv("/data/Batch_ForecastingPaper_bis/Allsummary.csv", guess_max = 1E4)

Allsummary_withSS62500<-read_csv(paste0(output_dir_data, "/Allsummary_withSS62500.csv"), guess_max = 1E4)


All_wait_cor_summary_RandomSample_AllStartSize<-read_csv(paste0(output_dir_data, "/All_wait_cor_summary_RandomSample_AllStartSize.csv"), guess_max = 1E4 )
Alldepth_wait_cor_summary_AllStartSize<-read_csv(paste0(output_dir_data, "/Alldepth_wait_cor_summary_AllStartSize.csv"), guess_max = 1E4)

#also contains pvalue and Confidence Interval
wait_cor_summary_RVAideMemoire<-read_csv(paste0(output_dir_data, "/wait_cor_summary_RVAideMemoire.csv"), guess_max = 1E4)
# contain p-values
wait_cor_summary_CombinedMutationRate<-read_csv(paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate.csv"), guess_max = 1E4)

wait_cor_summary_CombinedMutationRate_withSS62500<-read_csv(paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500.csv"), guess_max = 1E4)
wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved<-read_csv(paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved.csv"), guess_max = 1E4)
wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved<-read_csv(paste0(output_dir_data, "/wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved.csv"), guess_max = 1E4)


wait_cor_summarysmallerStartSize<-read_csv(file = paste0(output_dir_data, "/wait_cor_summarysmallerStartSizeK2_RVAideMemoire.csv"), guess_max = 1E4)
wait_cor_summarymediumStartSize<-read_csv(file = paste0(output_dir_data, "/wait_cor_summarymediumStartSizeK2_RVAideMemoire.csv"), guess_max = 1E4)


wait_cor_summary_WithSmallerSS_CombinedMutationRate<-read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_CombinedMutationRate_allSvalues.csv"), guess_max = 1E4)

wait_cor_summary_WithSmallerSS_CombinedFitnessEffect<-read_csv(paste0(output_dir_data, "/wait_cor_summary_WithSmallerSS_reducedMu_CombinedFitnessEffect.csv"), guess_max = 1E4)


dataNewCorrelation<-read_csv(file = paste0(output_dir_data, "/dataAugmentedK512_NewCorrelation.csv"), guess_max = 1E4)
Allcor_summary_NewCorrelation<-read_csv(file = paste0(output_dir_data, "/Allcor_summary_NewCorrelation_K512.csv"), guess_max = 1E4)
All_wait_cor_summary_NewCorrelation<-read_csv(file = paste0(output_dir_data, "/All_wait_cor_summary_NewCorrelation_K512.csv"), guess_max = 1E4)
Allsummary_NewCorrelation<-read_csv(file = paste0(output_dir_data, "/Allsummary_NewCorrelation_K512.csv"),  guess_max = 1E4)
wait_cor_summary_CombinedMutationRateK512<-read_csv(file = paste0(output_dir_data, "/wait_cor_summary_CombinedMutationRate_K512_NewCorrelation.csv"), guess_max = 1E4)
wait_cor_summary_CombinedFitnessEffectK512<-read_csv(file = paste0(output_dir_data, "/wait_cor_summary_CombinedFitnessEffect_K512_NewCorrelation.csv"), guess_max = 1E4)



# correlation with Average growth rate K=512
All_AverageGrowthRate_summary_NewCorrelation_K512<-read_csv(paste0(output_dir_data, "/All_AverageGrowthRate_cor_summary_NewCorrelation_K512.csv"), guess_max = 1E4)

CorelationAverageGrowthRateWithMeanSweep_K512<-read_csv(file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_NewCorrelation_K512.csv"), guess_max = 1E4)

CorelationAverageGrowthRateWithMeanSweep_K512_CombinedFitnessEffect<-read_csv(file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_CombinedFitnessEffect_K512_NewCorrelation.csv"), guess_max = 1E4)
CorelationAverageGrowthRateWithMeanSweep_K512_CombinedMutationRate<-read_csv(file = paste0(output_dir_data, "/CorelationAverageGrowthRateWithMeanSweep_CombinedMutationRate_K512_NewCorrelation.csv"), guess_max = 1E4)


# correlation with Average growth rate K=64
All_AverageGrowthRate_summary_NewCorrelation_K64<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_NewCorrelation_K64.csv"), guess_max = 1E4)
# correlation with Average growth rate K=4096
All_AverageGrowthRate_summary_NewCorrelation_K4096<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_NewCorrelation_K4096.csv"), guess_max = 1E4)








## random mu
waitingTime_cor_summary_CombinedMutationRateK512<-read_csv(paste0("/data/Batch_RandomMu_ForecastingPaper_Working", "/waitingTime_cor_summary_CombinedMutationRate_K512_RandomMu.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_CombinedMutationRateK512<-read_csv(paste0("/data/Batch_RandomMu_ForecastingPaper_Working", "/AverageGrowthRate_cor_summary_CombinedMutationRate_K512_RandomMu.csv"), guess_max = 1E4)
InverseWaitingTime_cor_summary_CombinedMutationRateK512<-read_csv(paste0("/data/Batch_RandomMu_ForecastingPaper_Working", "/InverseWaitingTime_cor_summary_CombinedMutationRate_K512_RandomMu.csv"), guess_max = 1E4)

AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512<-read_csv(paste0(output_dir_data, "/AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRate_K512_RandomMu.csv"), guess_max = 1E4)





#random s
waitingTime_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0("/data/Batch_RandomS_ForecastingPaper_Working", "/waitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)
AverageGrowthRate_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0("/data/Batch_RandomS_ForecastingPaper_Working", "/AverageGrowthRate_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)
InverseWaitingTime_cor_summary_CombinedFitnessEffectK512<-read_csv(paste0("/data/Batch_RandomS_ForecastingPaper_Working", "/InverseWaitingTime_cor_summary_CombinedFitnessEffect_K512_RandomS.csv"), guess_max = 1E4)
####### Which Diversity measure would be the best ? ####### 


for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  
  for(DivMeas in seq_along(DiversityMeasure)){
    print(DiversityMeasure[DivMeas])
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                       tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
    
    
  }
  
  
  gg<-ggplot(tmpSummary)+
    #first all points which are not significant
    # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
    # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
    # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
    # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
    
    geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
    #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
    #then all points which are  significant
    geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
    #first all points which are not significant
    # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
    # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
    # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
    # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
    # #then all points which are  significant
    # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
    # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
    # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
    # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
    # scale_shape_identity()+
    #first the dashed line with all points
    geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
    #now only the points significant at 0.05 level
    geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
    #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
    scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
    #                     values = c("solid"=16),
    #                     labels=c("solid"="TRUE"))+
    scale_size_manual(2)+
    guides(linetype=FALSE)+
    facet_grid( K ~ s_driver_birth)+
    theme_bw(base_size = 16)+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))
  
  
  
  
  
  png(paste0(output_dir_plots_paper, paste0("/wait_correlations_M",MuValues[MuIndex]  ), "AllDiversityMeasures.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
  
}




####### Which Diversity measure would be feasible ? #######


for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  
  Depth_4Big<-c(0:4)
  
  for(DepthIt in c(0:4)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
    
    
  }
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )
  
  
  gg<-ggplot(tmpSummary)+
    #first all points which are not significant
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
    #then all points which are  significant
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
    geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
    #first the dashed line with all points
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
    geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
    #now only the points significant at 0.05 level
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
    geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
    #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
    scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
                       labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
    #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
    #                     values = c("solid"=16),
    #                     labels=c("solid"="TRUE"))+
    # scale_size_manual(2)+
    facet_grid( K ~ s_driver_birth)+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 5),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))+
    
    
    
    
    
    
    
    png(paste0(output_dir_plots_paper, paste0("/wait_correlations_M",MuValues[MuIndex]  ), "DriverDiversityFrom4BigSamples_DifferentDepth.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
  
  
}

# removing depth 0
for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  
  Depth_4Big<-c(1:4)
  Denomination<-c(letters[seq_along(Depth_4Big)], letters[(length(Depth_4Big)+1)])
  
  
  for(DepthIt in seq_along(Depth_4Big) ){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],  "_Significant",SignificanceLevel )
    
    
  }
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )
  
  
  gg<-ggplot(tmpSummary)
  
  for(DepthIt in seq_along(Depth_4Big) ){
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],", as.logical(NA))"), 
                               color=shQuote(Denomination[DepthIt])),
                 size=2)
    print("first geom point")
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse( Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[DepthIt])),
                 size=2)
    print("second geom point")
    
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt]),
                              color=shQuote(Denomination[DepthIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse( Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom4BigSamplesAtDepth",Depth_4Big[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[DepthIt])))
    print("second geom line")
    
  }
  
  
  gg<-gg+
    geom_point(aes(x=start_size,group = K, shape ="dashed")+
                 aes_string( y =paste0("ifelse(!Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel,
                                       ", Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA) )"), 
                             color=shQuote(Denomination[length(Denomination)])),
               size=2)
  
  gg<-gg+
    geom_point(aes(x=start_size,group = K, shape ="solid")+
                 aes_string( y =paste0("ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel,
                                       ", Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA) )"),
                             color=shQuote(Denomination[length(Denomination)])),
               size=2)
  
  
  gg<-gg+
    geom_line(aes(x= start_size, group = K, linetype="dashed")+
                aes_string( y ="Cor_DriverDiversityFrom4BigRandomSamples", 
                            color=shQuote(Denomination[length(Denomination)])))
  
  gg<-gg+
    geom_line(aes(x= start_size, group = K , linetype="solid")+
                aes_string( y =paste0("ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel,
                                      ", Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA) )"), 
                            color=shQuote(Denomination[length(Denomination)])))
  
  
  #add scales
  gg<-gg+
    scale_color_manual(name="Biopsy size", 
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="black"),
                       labels=c("a"="Depth 1", "b"="Depth 2", "c"="Depth 3", "d"="Depth 4", "e"="Random Samples"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( K ~ s_driver_birth)+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    xlab("tumour size at measurement")+
    ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
    ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))+
    
    
    png(paste0(output_dir_plots_paper, paste0("/wait_correlations_M",MuValues[MuIndex]  ), "DriverDiversityFrom4BigSamples_DifferentDepth_NoDepth0.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}






####### How small could be the biopsy sample? #######


for(MuIndex in 1){
  
  tmpSummary<-wait_cor_summary_RVAideMemoire[which(wait_cor_summary_RVAideMemoire$mu_driver_birth==MuValues[MuIndex]),]
  
  SampleSize<-c("1","4", "1Big", "4Big")
  Depth<-c(0, 1)
  Denomination<-c("a", "b", "c", "d")
  
  
  for(DepthIt in seq_along(Depth)){
    
    for(SizeIt in seq_along(SampleSize)){
      
      tmpSummary<-mutate(tmpSummary,
                         tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
      colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
    }
    
    gg<-ggplot(tmpSummary)
    
    for(SizeIt in seq_along(SampleSize)){
      
      gg<-gg+
        geom_point(aes(x=start_size,group = K, shape ="dashed")+
                     aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                           ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                                 color=shQuote(Denomination[SizeIt])),
                   size=2)
      print("first geom point")
      
      gg<-gg+
        geom_point(aes(x=start_size,group = K, shape ="solid")+
                     aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                           ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                                 color=shQuote(Denomination[SizeIt])),
                   size=2)
      print("second geom point")
      
      # gg<-gg+
      #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
      #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
      #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
      #                           color=shQuote(Denomination[SizeIt])))
      gg<-gg+
        geom_line(aes(x= start_size, group = K, linetype="dashed")+
                    aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                                color=shQuote(Denomination[SizeIt])))
      print("first geom line")
      gg<-gg+
        geom_line(aes(x= start_size, group = K , linetype="solid")+
                    aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                          ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                                color=shQuote(Denomination[SizeIt])))
      print("second geom line")
      
    }
    
    print("add scales")
    gg<-gg+
      scale_color_manual(name="Biopsy size", 
                         values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                         labels=c("a"="1 sample", "b"="1Big sample", "c"="4 samples","d"="4Big samples"))+
      scale_linetype_manual(name="",
                            values = c("dashed"="dashed", "solid"="solid"))+
      scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                         values = c("dashed"=17, "solid"=16),
                         labels=c("dashed"="FALSE", "solid"="TRUE"))+
      facet_grid( K ~ s_driver_birth)+
      theme_bw(base_size = 16)+
      guides(linetype=FALSE,
             color= guide_legend(title.theme = element_text(size = 12),
                                 label.theme = element_text(size=10)), 
             shape= guide_legend(title.theme = element_text(size = 12),
                                 label.theme = element_text(size=10)))+
      ylim(-1, 1)+
      geom_hline(aes(yintercept=0), linetype="dashed")+
      xlab("tumour size at measurement")+
      ylab(paste0("correlation coefficient:\n DriverDiversity at depth ",Depth[DepthIt], " vs waiting time"))+
      ggtitle(paste0("Driver mutation rate = ",MuValues[MuIndex]))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))+
      
      
      png(paste0(output_dir_plots_paper, paste0("/wait_correlations_M",MuValues[MuIndex]  ), "DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], ".png"),  width = 1000, height = 1000, res = 100)
    print(gg)
    dev.off()
  }
}  

#######  Combining mutation rate####### 
####### Which Diversity measure would be the best ? #######     



#tmpSummary<-wait_cor_summary_CombinedMutationRate
tmpSummary<-wait_cor_summary_CombinedMutationRate_withSS62500

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( K ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_CombinedMutationRate.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

####### Which Diversity measure would be feasible ? ####### 

tmpSummary<-wait_cor_summary_CombinedMutationRate_withSS62500

Depth_4Big<-c(0:4)

for(DepthIt in c(0:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
                     labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  # scale_size_manual(2)+
  facet_grid( K ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 5),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))







png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

####### How small could be the biopsy sample? #######


#tmpSummary<-wait_cor_summary_CombinedMutationRate_withSS62500
tmpSummary<-wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved
SampleSize<-c("1","4", "1Big", "4Big")
Depth<-c(0, 1)
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    if( (Depth[DepthIt]==0 & SizeIt !=length(SampleSize)) || (Depth[DepthIt]==1 & SizeIt ==1)  ){
      gg<-gg+
        geom_point(aes(x=start_size,group = K, shape ="dashed")+
                     aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                           ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                                 color=shQuote(Denomination[SizeIt])),
                   size=2)
      print("first geom point")
    }
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size", 
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 sample", "b"="1Big sample", "c"="4 samples","d"="4Big samples"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( K ~ s_driver_birth)+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n DriverDiversity at depth ",Depth[DepthIt], " vs waiting time"))+
    ggtitle(paste0("Combined driver mutation rates"))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], ".png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}


# ## Search for bug
# tmpSummary[which(tmpSummary$K==512 & tmpSummary$s_driver_birth==0.05), "Cor_DriverDiversityFrom1SamplesAtDepth1"]
# Allsummary_withSS62500[which(Allsummary_withSS62500$K==512 & Allsummary_withSS62500$s_driver_birth==0.05), c("DriverDiversityFrom1SamplesAtDepth1", "waiting_time")]
# 
# wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved[which(wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved$K==512 & wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved$s_driver_birth==0.05), "Cor_DriverDiversityFrom1SamplesAtDepth1"]


#######  Combining fitness effect ####### 
####### Which Diversity measure would be the best ? #######     


tmpSummary<-wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( K ~ mu_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined fitness effects "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnessEffect_AllDiversityMeasures_facetMu.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

####### Which Diversity measure would be feasible ? ####### 

tmpSummary<-wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved

Depth_4Big<-c(0:4)

for(DepthIt in c(0:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
                     labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  # scale_size_manual(2)+
  facet_grid( K ~ mu_driver_birth)+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 5),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
  ggtitle(paste0("Combined fitness effects"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))



png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnessEffect_DriverDiversityFrom4BigSamples_DifferentDepth_facetMu.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

####### How small could be the biopsy sample? #######


tmpSummary<-wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved
SampleSize<-c("1","4", "1Big", "4Big")
Depth<-c(0, 1)
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size", 
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 sample", "b"="1Big sample", "c"="4 samples","d"="4Big samples"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( K ~ mu_driver_birth)+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n DriverDiversity at depth ",Depth[DepthIt], " vs waiting time"))+
    ggtitle(paste0("Combined fitness effects"))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnesseEffect_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_facetMu.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}





####### Final Version     ####### 
#######  Fig1.b

tmpSummary<-wait_cor_summary_RVAideMemoire[which( (wait_cor_summary_RVAideMemoire$mu_driver_birth %in% c(1e-06, 1e-04)) &
                                                    (wait_cor_summary_RVAideMemoire$K==512) &
                                                    (wait_cor_summary_RVAideMemoire$s_driver_birth %in% c(0.05, 0.2)))
                                           ,]


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( mu_driver_birth ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Deme size = 512"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))





png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Figure 2  #####    
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


gg1<-ggplot(subsetDataForScatterPlot)+
  geom_point(aes(x=DriverDiversity, y=waiting_time), alpha=0.5)+
  ggtitle(paste0("DriverDiversity vs waiting_time \n K=",KOfInterest, ", mu=",mu_driver_birthOfInterest, " , s=",s_driver_birthOfInterest, "\n start_size=", start_sizeOfInterest))+
  theme_bw()+
  guides(alpha=FALSE)+
  ylim(0, 0.06)+
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio=1)

png(paste0(output_dir_plots_paper, paste0("/", "ScatterPlotToCheckPositiveCorrelation_K",KOfInterest, "s", s_driver_birthOfInterest, "startSize",start_sizeOfInterest,  "DriverVsMeanBirthRate_Ylim.png")),  width = 1000, height = 1000, res = 100)
print(gg1)
dev.off()


##### Figure 3  #####    

tmpSummary<-wait_cor_summary_RVAideMemoire[which( (wait_cor_summary_RVAideMemoire$mu_driver_birth %in% c(1e-06, 1e-04)) &
                                                    (wait_cor_summary_RVAideMemoire$K==512) &
                                                    (wait_cor_summary_RVAideMemoire$s_driver_birth %in% c(0.05, 0.2)))
                                           ,]


Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
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
  #geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  #geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  scale_color_manual(name="", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  # scale_size_manual(2)+
  facet_grid( mu_driver_birth ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 5),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
  ggtitle(paste0("Deme size = 512"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))+
  
  
  
  
  
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_K512_DriverDiversityFrom4BigSamples_DifferentDepth_FinalVersion1.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 4b  #####     

tmpSummary<-subset(wait_cor_summary_CombinedMutationRate_withSS62500, 
                   wait_cor_summary_CombinedMutationRate_withSS62500$K==512)

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates  \n Deme size=512 "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_K512.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()




gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates  \n Deme size=512 "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_K512_sqaure.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 0.1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates  \n Deme size=512 "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_K512_Ylim.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 4a  ##### 

tmpSummary<-subset(wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved, 
                   wait_cor_summary_CombinedFitnessEffect_withSS62500_NaRemoved$K==512)

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=DiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ mu_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Combined fitness effects  \n Deme size=512 "))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnessEffect_AllDiversityMeasures_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 5a  #####     

tmpSummary<-subset(wait_cor_summary_CombinedMutationRate_withSS62500, 
                   wait_cor_summary_CombinedMutationRate_withSS62500$K==512)

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
  #geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  #geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  # scale_size_manual(2)+
  facet_grid( . ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 5),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 0.1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n DriverDiversity from 4Big samples vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates \n Deme size=512"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)







png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_Ylim.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 5b       ##### 


tmpSummary<-subset(wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved, 
                   wait_cor_summary_CombinedMutationRate_withSS62500_NaRemoved$K==512)

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==0 & SizeIt !=length(SampleSize)) || (Depth[DepthIt]==1 & SizeIt ==1)  ){
    # gg<-gg+
    #   geom_point(aes(x=start_size,group = K, shape ="dashed")+
    #                aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                      ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
    #                            color=shQuote(Denomination[SizeIt])),
    #              size=2)
    # print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 sample", "b"="1Big sample", "c"="4 samples","d"="4Big samples"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( . ~ s_driver_birth)+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)), 
           shape= guide_legend(title.theme = element_text(size = 12),
                               label.theme = element_text(size=10)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n DriverDiversity at depth ",Depth[DepthIt], " vs waiting time"))+
    ggtitle(paste0("Combined driver mutation rates \ Deme size = 512"))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1))
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

gg<-gg+
  scale_color_manual(name="Biopsy size",
                     values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                     labels=c("a"="1 sample", "b"="1Big sample", "c"="4 samples","d"="4Big samples"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  facet_grid( . ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  ylim(-1, 0.1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n DriverDiversity at depth ",Depth[DepthIt], " vs waiting time"))+
  ggtitle(paste0("Combined driver mutation rates \ Deme size = 512"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512Ylim_Square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


######### Add all start size      ######### 

Allwait_cor_summary_RVAideMemoire<-wait_cor_summary_RVAideMemoire[, !grepl("CI", colnames(wait_cor_summary_RVAideMemoire))]
Allwait_cor_summary_RVAideMemoire<-rbind(Allwait_cor_summary_RVAideMemoire, wait_cor_summarysmallerStartSize)
Allwait_cor_summary_RVAideMemoire<-rbind(Allwait_cor_summary_RVAideMemoire, wait_cor_summarymediumStartSize)

####### Final Final Version     ####### 
#######  Fig1.b      ####### 

tmpSummary<-Allwait_cor_summary_RVAideMemoire[which( (Allwait_cor_summary_RVAideMemoire$mu_driver_birth %in% c(1e-06, 1e-04)) &
                                                       (Allwait_cor_summary_RVAideMemoire$K==512) &
                                                       (Allwait_cor_summary_RVAideMemoire$s_driver_birth %in% c(0.05, 0.2)))
                                              ,]


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( mu_driver_birth ~ s_driver_birth)+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000, 60000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle(paste0("Deme size = 512"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures_AllStartSize.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


mu_driver_birth.lab <- paste0(expression(mu), " = ",unique(tmpSummary$mu_driver_birth))
#mu_driver_birth.lab <- paste0(expression(mu, " = ", unique(tmpSummary$mu_driver_birth)))
names(mu_driver_birth.lab) <- unique(tmpSummary$mu_driver_birth)
s_driver_birth.lab <- paste0("s = ",unique(tmpSummary$s_driver_birth))
names(s_driver_birth.lab) <- unique(tmpSummary$s_driver_birth)

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  # facet_grid( mu_driver_birth ~ s_driver_birth, labeller = label_parsed(mu_driver_birth = paste0(expression(mu), " = ",unique(tmpSummary$mu_driver_birth)),
  #                                                                   s_driver_birth=paste0("s = ",unique(tmpSummary$s_driver_birth))))+
  facet_grid( mu_driver_birth ~ s_driver_birth, labeller = labeller(mu_driver_birth = mu_driver_birth.lab,
                                                                    s_driver_birth=s_driver_birth.lab))+
  # facet_grid( mu_driver_birth ~ s_driver_birth, labeller = label_parsed(mu_driver_birth = mu_driver_birth.lab,
  #                                                                   s_driver_birth=s_driver_birth.lab))+
  # 
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000, 60000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)








png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures_AllStartSize_NewLabelNoTitle.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( mu_driver_birth ~ s_driver_birth, labeller =label_bquote(rows = mu == .(mu_driver_birth), 
                                                                       cols= s== .(s_driver_birth)))+
  
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  # scale_x_continuous(name="tumour size at measurement", breaks =c(10000, 60000,125000, 250000, 375000, 500000, 625000,  750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)




png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures_AllStartSize_muLabelNoTitle.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()



##### Figure 3b  #####     


tmpSummary<-subset(wait_cor_summary_WithSmallerSS_CombinedMutationRate, 
                   wait_cor_summary_WithSmallerSS_CombinedMutationRate$K==512)

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

s_driver_birth.lab <- paste0("s = ",unique(tmpSummary$s_driver_birth))
names(s_driver_birth.lab) <- unique(tmpSummary$s_driver_birth)


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=labeller(s_driver_birth=s_driver_birth.lab))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_K512_square_WithSmallerSS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Figure 3a  ##### 

tmpSummary<-subset(wait_cor_summary_WithSmallerSS_CombinedFitnessEffect, 
                   wait_cor_summary_WithSmallerSS_CombinedFitnessEffect$K==512)

for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnessEffect_AllDiversityMeasures_square_WithSmallerSS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Figure 3d  #####     

tmpSummary<-subset(wait_cor_summary_WithSmallerSS_CombinedMutationRate, 
                   wait_cor_summary_WithSmallerSS_CombinedMutationRate$K==512)

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
  #geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  #geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_WithSmallerSS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Figure 3c       ##### 


tmpSummary<-subset(wait_cor_summary_WithSmallerSS_CombinedMutationRate, 
                   wait_cor_summary_WithSmallerSS_CombinedMutationRate$K==512)

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==0 & SizeIt !=length(SampleSize)) || (Depth[DepthIt]==1 & SizeIt ==1)  ){
    # gg<-gg+
    #   geom_point(aes(x=start_size,group = K, shape ="dashed")+
    #                aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                      ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
    #                            color=shQuote(Denomination[SizeIt])),
    #              size=2)
    # print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_WithSmallerSS.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}




####### Supp Figures       ####### 

gg<-ggplot(subset(Allsummary_withSS62500, 
                  Allsummary_withSS62500$start_size ==500000))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/DiversityVsSurvivalTime_SS500000.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_withSS62500, 
                  Allsummary_withSS62500$start_size ==500000))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10DiversityVsSurvivalTime_SS500000.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10DiversityVsSurvivalTime_SS500000K512s005.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 ) ))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10DiversityVsSurvivalTime_SS500000K512AllS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigRandomSamples, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigRandomSamplesVsSurvivalTime_SS500000.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsSurvivalTime_SS500000K512s005.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Survival time")+
  facet_grid(.~mu_driver_birth)+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsSurvivalTime_SS500000K512s005_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=DriverDiversity, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Clonal Diversity [log10]")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsLog10ClonalDiversity_SS500000K512s005.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_withSS62500, 
                  (Allsummary_withSS62500$start_size ==500000 &
                     Allsummary_withSS62500$K==512 & 
                     Allsummary_withSS62500$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=DriverDiversity, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Clonal Diversity [log10]")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  facet_grid(.~mu_driver_birth)+
  scale_x_log10()+
  scale_y_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsLog10ClonalDiversity_SS500000K512s005_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


# #search for simulatioin with ==2
# test<-subset(Allsummary_withSS62500, 
#        (Allsummary_withSS62500$start_size ==500000 &
#           Allsummary_withSS62500$K==512 & 
#           Allsummary_withSS62500$s_driver_birth==0.05) )
# colnames(test)
# 
# as.data.frame(unique(test[which(test$DriverDiversityFrom1BigSamplesAtDepth1==2), c(c(1:16), 18) ]))
# SimulationsWithDiversity2<-unique(test[which(test$DriverDiversityFrom1BigSamplesAtDepth1==2), c(1:18) ])
# fwrite(SimulationsWithDiversity2, file = paste0(output_dir_data, "/SimulationsWith_DriverDiversityFrom1BigSamplesAtDepth1_2.csv"))
# 



##### Figure 1a ######

# data<-subset(Allsummary_withSS62500, 
#               (Allsummary_withSS62500$K==512 &
#                  Allsummary_withSS62500$s_driver_birth==0.05 &
#                  Allsummary_withSS62500$mu_driver_birth==1e-06))
# 
# dataSS500000<-subset(Allsummary_withSS62500, 
#                      (Allsummary_withSS62500$start_size ==500000 &
#                         Allsummary_withSS62500$K==512 &
#                         Allsummary_withSS62500$s_driver_birth==0.05 &
#                         Allsummary_withSS62500$mu_driver_birth==1e-06))
# fwrite(dataSS500000, file = paste0(output_dir_data, "/data_SS500000.csv"))
# 
# parameters_dataSS500000<-unique(dataSS500000[, c(1:18)])
# fwrite(parameters_dataSS500000, file = paste0(output_dir_data, "/parameters_dataSS500000.csv"))
# 
# gg<-ggplot(data)+
#   geom_line(aes(x=new_time, y=NumCells, group=seed, colour = MeanBirthRate )) + 
#   scale_colour_distiller(palette = "RdYlBu", name = "Mean cell division rate") +
#   scale_x_continuous(name = "time") +
#   scale_y_continuous(name = "tumour size (number of cells)") +
#   theme_classic()
# 
# png(paste0(output_dir_plots_paper, "/TunourTrajectoryK512mu1e06s005SS500000.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()

data<-vector(mode="list", length=length(c(0,1,2)))
#for(Kiter in seq_along(unique(c(0,1,2)))){

#Kval<-c(0,1,2)[Kiter]
Kval<-1
data_tmp<- as.data.frame(read_csv(paste0(output_dir_data, "/dataAugmentedK", Kval, ".csv"), guess_max = 1E4))
data<-subset(data_tmp,
             (data_tmp$K==512 &
                data_tmp$s_driver_birth==0.05 &
                data_tmp$mu_driver_birth==1e-06))
remove(data_tmp)
#}

#data<-do.call(rbind, data)


gg<-ggplot(na.omit(data))+
  geom_line(aes(x=new_time, y=NumCells, group=seed, colour = MeanBirthRate )) + 
  scale_colour_distiller(palette = "RdYlBu", name = "Mean cell division rate") +
  scale_x_continuous(name = "time") +
  scale_y_continuous(name = "tumour size (number of cells)") +
  theme_classic()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/TumourTrajectoryK512mu1e06s005.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

######## New correlation version  (depending on waiting time and not on gen_adj)    ######## 

########## Fig1.a      ########## 
data<-subset(dataNewCorrelation,
             (dataNewCorrelation$K==512 &
                dataNewCorrelation$s_driver_birth==0.05 &
                dataNewCorrelation$mu_driver_birth==1e-06))

remove(dataNewCorrelation)

gg<-ggplot(na.omit(data))+
  geom_line(aes(x=new_time, y=NumCells, group=seed, colour = MeanBirthRate )) + 
  scale_colour_distiller(palette = "RdYlBu", name = "Mean cell division rate") +
  scale_x_continuous(name = "time") +
  scale_y_continuous(name = "tumour size (number of cells)") +
  theme_classic()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/TumourTrajectoryK512mu1e06s005_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(na.omit(data))+
  geom_line(aes(x=new_time, y=NumCells, group=seed, colour = 100*rank_MeanBirthRate0 )) + 
  scale_colour_distiller(palette = "RdYlBu", name = "Mean cell division rate (rank)") +
  scale_x_continuous(name = "time") +
  scale_y_continuous(name = "tumour size (number of cells)") +
  theme_classic(base_size=16)+
  theme(aspect.ratio=1)+
  guides(color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))

png(paste0(output_dir_plots_paper, "/TumourTrajectoryK512mu1e06s005_NewCorrelation_ColorRankMeanBirthRate0.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


########## Fig1.b      ########## 

tmpSummary<-All_wait_cor_summary_NewCorrelation[which( (All_wait_cor_summary_NewCorrelation$mu_driver_birth %in% c(1e-06, 1e-04)) &
                                                         (All_wait_cor_summary_NewCorrelation$K==512) &
                                                         (All_wait_cor_summary_NewCorrelation$s_driver_birth %in% c(0.05, 0.2)))
                                                ,]


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( mu_driver_birth ~ s_driver_birth, labeller =label_bquote(rows = mu == .(mu_driver_birth), 
                                                                       cols= s== .(s_driver_birth)))+
  
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  # scale_x_continuous(name="tumour size at measurement", breaks =c(10000, 60000,125000, 250000, 375000, 500000, 625000,  750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)




png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures_AllStartSize_muLabelNoTitle_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

########## Fig3.a      ########## 

tmpSummary<-wait_cor_summary_CombinedFitnessEffectK512


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedFitnessEffect_AllDiversityMeasures_square_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3b  #####     


tmpSummary<-wait_cor_summary_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

s_driver_birth.lab <- paste0("s = ",unique(tmpSummary$s_driver_birth))
names(s_driver_birth.lab) <- unique(tmpSummary$s_driver_birth)


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #   geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200", "#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  facet_grid( . ~ s_driver_birth, labeller=labeller(s_driver_birth=s_driver_birth.lab))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_AllDiversityMeasures_K512_square_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3c       ##### 


tmpSummary<-wait_cor_summary_CombinedMutationRateK512

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==0 & SizeIt !=length(SampleSize)) || (Depth[DepthIt]==1 & SizeIt ==1)  ){
    # gg<-gg+
    #   geom_point(aes(x=start_size,group = K, shape ="dashed")+
    #                aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                      ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
    #                            color=shQuote(Denomination[SizeIt])),
    #              size=2)
    # print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}


##### Figure 3d  #####     

tmpSummary<-wait_cor_summary_CombinedMutationRateK512

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
  #geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  #first the dashed line with all points
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth0, group = K, color="a", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth1, group = K , color="b", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth2, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth3, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigSamplesAtDepth4, group = K , color="e", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_DriverDiversityFrom4BigRandomSamples, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  #geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Diversity vs waiting time"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/wait_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


##### Supp figure #######

gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==500000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10DiversityVsSurvivalTime_SS500000K512s005_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==250000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversity, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10DiversityVsSurvivalTime_SS250000K512s005_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

test<-subset(Allsummary_NewCorrelation, 
             (Allsummary_NewCorrelation$start_size ==250000 &
                Allsummary_NewCorrelation$K==512 & 
                Allsummary_NewCorrelation$s_driver_birth==0.05) )
test<-test[, colnames(test) %in% c("DriverDiversity", "waiting_time" , "mu_driver_birth")]
fwrite(test, file = paste0(output_dir_data, "/Data_Log10DiversityVsSurvivalTime_SS250000K512s005_NewCorrelation.csv"))




gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==500000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Survival time")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsSurvivalTime_SS500000K512s005_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==500000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Survival time")+
  facet_grid(.~mu_driver_birth)+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsSurvivalTime_SS500000K512s005_NewCorrelation_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==250000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=waiting_time, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Survival time")+
  facet_grid(.~mu_driver_birth)+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsSurvivalTime_SS250000K512s005_NewCorrelation_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==500000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=DriverDiversity, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Clonal Diversity [log10]")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsLog10ClonalDiversity_SS500000K512s005_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==500000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=DriverDiversity, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Clonal Diversity [log10]")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  facet_grid(.~mu_driver_birth)+
  scale_x_log10()+
  scale_y_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsLog10ClonalDiversity_SS500000K512s005_NewCorrelation_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(subset(Allsummary_NewCorrelation, 
                  (Allsummary_NewCorrelation$start_size ==250000 &
                     Allsummary_NewCorrelation$K==512 & 
                     Allsummary_NewCorrelation$s_driver_birth==0.05) ))+
  geom_point(aes(x=DriverDiversityFrom1BigSamplesAtDepth1, y=DriverDiversity, color=as.factor(mu_driver_birth) ), alpha=0.5, size=1)+
  xlab("Clonal Diversity, \n 1 biopsy of 1000 cells, depth 1 [log10]")+
  ylab("Clonal Diversity [log10]")+
  scale_color_manual(name="Mutation rate",values=c("#00B5EC", "#009F71", "#FAC200"))+
  theme_bw()+
  facet_grid(.~mu_driver_birth)+
  scale_x_log10()+
  scale_y_log10()+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/Log10Diversity1BigSamplesAtDepth1VsLog10ClonalDiversity_SS250000K512s005_NewCorrelation_facet.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##########   Correlation average growth rate    ########## 
########## Fig1.b      ########## 

tmpSummary<-All_AverageGrowthRate_summary_NewCorrelation_K512[which( (All_AverageGrowthRate_summary_NewCorrelation_K512$mu_driver_birth %in% c(1e-06, 1e-04)) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$K==512) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$s_driver_birth %in% c(0.05, 0.2)))
                                                              ,]


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

gg<-ggplot(tmpSummary)+
  
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c"), linetype="solid")+
  
  # #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), 
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K, 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), 
                 group = K , 
                 fill="c", 
                 stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  #then all points which are  significant
  geom_point(aes(x= start_size, 
                 y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K , 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), 
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), 
                 group = K , 
                 fill="c", 
                 stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
  
  scale_size_manual(point_size)+
  guides(linetype=FALSE,
         color=FALSE,
         fill = guide_legend(override.aes = list(stroke = 0))) + 
  
  continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
                   palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                   breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
  
  facet_grid( mu_driver_birth ~ s_driver_birth, labeller =label_bquote(rows = mu == .(mu_driver_birth), 
                                                                       cols= s== .(s_driver_birth)))+
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n predictor variables vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)




png(paste0(output_dir_plots_paper, "/wait_correlations_K512_AllDiversityMeasures_AllStartSize_muLabelNoTitle_CorrelationAverageGrowthRate_Revised.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


####### Fig2    ####### 


data<-subset(dataNewCorrelation,
             (dataNewCorrelation$K==512 &
                dataNewCorrelation$s_driver_birth==0.05 &
                dataNewCorrelation$mu_driver_birth==1e-04 ))

remove(dataNewCorrelation)

# take the value around 750000 for which both the MeanBirthRate and the DriverDiversity are given 
data<-subset(data, data$NumCells <= 750000)
data<-data[! is.na(data$MeanBirthRate), ]
data<-data[! is.na(data$DriverDiversity), ]
data<-data %>% group_by(SimulationNb) %>% mutate(maxNumCell =max(NumCells)) %>% mutate(isMax = (NumCells==maxNumCell)) %>% ungroup()

gg<-ggplot(subset(data, 
                  data$isMax==TRUE))+
  geom_point(aes(x=DriverDiversity, y=MeanBirthRate))+
  xlab("Clonal diversity")+
  ylab("Mean cell division rate")+
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/NegCorrelationClonalDiverMeanCellDivRate.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

FileWithAllParam<-readSeedFromParamFiles("/all_results/Batch_ForecastingPaper_bis")


ParamPanelb<-data[which(data$MeanBirthRate==max(data$MeanBirthRate)), c(1:18)]
ParamPanelb<-merge(ParamPanelb, FileWithAllParam)
PathPanelb<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_2/s_driver_birth_0/seed_", ParamPanelb$seedFolder, "/")
Check<-fread(paste0(PathPanelb, "parameters.dat"))                         
#plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE)
plot_all_images_cutoff(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)
plot_all_images_cutoff(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff01", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.1)
plot_all_images_cutoff_modified(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff01", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.1)



ParamPanelc<-data[which(data$DriverDiversity==max(data$DriverDiversity)), c(1:18)]
ParamPanelc<-merge(ParamPanelc, FileWithAllParam)
PathPanelc<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_2/s_driver_birth_0/seed_", ParamPanelc$seedFolder, "/")
Check<-fread(paste0(PathPanelc, "parameters.dat"))                         
#plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE)
plot_all_images_cutoff(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)
plot_all_images_cutoff_modified(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)


##### Plot heatmap sweep values       ##### 

Allsummary[which(Allsummary$K==4096 & Allsummary$mu_driver_birth==0.0001 & Allsummary$s_driver_birth==0.2), which(colnames(Allsummary)=="SimulationType")]<-35

MeanSelectiveSweep<-Allsummary %>% group_by(SimulationType) %>% mutate(GeneralMeanSweep = mean(mean_autocor, na.rm=TRUE),
                                                                       GeneralMedianSweep = median(mean_autocor, na.rm=TRUE)
)
MeanSelectiveSweep<-unique(MeanSelectiveSweep[, colnames(MeanSelectiveSweep) %in% c("K", "s_driver_birth","mu_driver_birth", "GeneralMeanSweep", "GeneralMedianSweep" )])

remove(Allsummary)

for(Kit in unique(MeanSelectiveSweep$K)){
  
  MeanData<-subset(MeanSelectiveSweep, 
                   MeanSelectiveSweep$K==Kit)
  MeanData$GeneralMedianSweep<-NULL
  MeanData<-as.data.frame(MeanData)
  
  gg<-ggplot(MeanData, aes(x=as.factor(mu_driver_birth), y=as.factor(s_driver_birth), fill = GeneralMeanSweep)) +
    geom_tile() + 
    xlab("Driver mutation rate")+
    ylab("Mean driver mutation fitness effect")+
    ggtitle(paste0("Deme size = ",Kit))+
    theme_bw(base_size = 16)+
    theme(aspect.ratio=1,
          plot.title = element_text(hjust = 0.5))+
    guides(fill= guide_legend(title.theme = element_text(size = 12),
                              label.theme = element_text(size=10)))+
    #scale_fill_gradient(name = "Mean clonal turnover")+
    scale_fill_viridis(name = "Mean clonal turnover",  discrete = FALSE, limits=c(0.02, 0.26))
  #begin = 0.04, end =0.25,
  
  png(paste0(output_dir_plots_paper, "/HeatmapMeanSweepValues_K", Kit, ".png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
}

for(Kit in unique(MeanSelectiveSweep$K)){
  
  MedianData<-subset(MeanSelectiveSweep, 
                     MeanSelectiveSweep$K==Kit)
  MedianData$GeneralMeanSweep<-NULL
  MedianData<-as.data.frame(MedianData)
  
  gg<-ggplot(MedianData, aes(x=as.factor(mu_driver_birth), y=as.factor(s_driver_birth), fill = GeneralMedianSweep)) +
    geom_tile() + 
    xlab("Driver mutation rate")+
    ylab("Median driver mutation fitness effect")+
    ggtitle(paste0("Deme size = ",Kit))+
    theme_bw(base_size = 16)+
    theme(aspect.ratio=1,
          plot.title = element_text(hjust = 0.5))+
    guides(fill= guide_legend(title.theme = element_text(size = 12),
                              label.theme = element_text(size=10)))+
    #scale_fill_gradient(name = "Median clonal turnover")+
    scale_fill_viridis(name = "Median clonal turnover",  discrete = FALSE, limits=c(0.02, 0.26))
  #begin = 0.04, end =0.25,
  
  png(paste0(output_dir_plots_paper, "/HeatmapMedianSweepValues_K", Kit, ".png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
  
}


#as.data.frame(test[1, c(c(1:18), which(colnames(test) %in% c("SimulationNb", "MeanBirthRate" ) ))])
#   unique(test[which(test$SimulationNb==2022), which(colnames(test)==c( "MeanBirthRate" ))])
#   
#   dataTest<-subset(data,
#                (data$K==512 &
#                   data$s_driver_birth==0.05 &
#                   data$mu_driver_birth==1e-04 ))
#   test<-dataTest[which(dataTest$seed==2000), ]
# test[which(is.na(test$MeanBirthRate)),which(colnames(test)=="MeanBirthRate")]
# which(test$NumCells==375000)
# test[c(653:656),which(colnames(test)%in% c("MeanBirthRate" , "NumCells"))]
# test[which(test$NumCells==750000),which(colnames(test)%in% c("MeanBirthRate" , "NumCells"))]
# test[c(810:813),which(colnames(test)%in% c("MeanBirthRate" , "NumCells"))]

########## Correlation average growth rate with Mean Sweep value        ########## 

tmpSummary<-CorelationAverageGrowthRateWithMeanSweep_K512

tmpSummary<-mutate(tmpSummary,
                   Cor_mean_autocor_pVal_Significant =(tmpSummary[["Cor_mean_autocor_pVal"]] <= SignificanceLevel),
                   Cor_mean_autocor_pVal_Unsignificant =(tmpSummary[["Cor_mean_autocor_pVal"]] > SignificanceLevel))

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( Cor_mean_autocor_pVal_Unsignificant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K, color="blue", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  scale_color_manual(name="", values=c("#00B5EC"), labels="Mean clonal turnover")+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( mu_driver_birth ~ s_driver_birth,
              labeller =label_bquote(rows = mu == .(mu_driver_birth),
                                     cols= s== .(s_driver_birth)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n  Mean Clonal turnover vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/CorrelationAverageGrowthRateWithMeanSweep_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


# Combined Fitness Effect
tmpSummary<-CorelationAverageGrowthRateWithMeanSweep_K512_CombinedFitnessEffect

tmpSummary<-mutate(tmpSummary,
                   Cor_mean_autocor_pVal_Significant =(tmpSummary[["Cor_mean_autocor_pVal"]] <= SignificanceLevel),
                   Cor_mean_autocor_pVal_Unsignificant =(tmpSummary[["Cor_mean_autocor_pVal"]] > SignificanceLevel))

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( Cor_mean_autocor_pVal_Unsignificant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  #first the dashed line with all points
  #geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K, color="blue", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  scale_color_manual(name="", values=c("#00B5EC"), labels="Mean clonal turnover")+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,
              labeller =label_bquote(cols= mu== .(mu_driver_birth)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n  Mean Clonal turnover vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/CorrelationAverageGrowthRateWithMeanSweep_K512_square_CombinedFitnessEffect.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

# Combined Mutation rate
tmpSummary<-CorelationAverageGrowthRateWithMeanSweep_K512_CombinedMutationRate

tmpSummary<-mutate(tmpSummary,
                   Cor_mean_autocor_pVal_Significant =(tmpSummary[["Cor_mean_autocor_pVal"]] <= SignificanceLevel),
                   Cor_mean_autocor_pVal_Unsignificant =(tmpSummary[["Cor_mean_autocor_pVal"]] > SignificanceLevel))

gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  #geom_point(aes(x= start_size, y =ifelse( Cor_mean_autocor_pVal_Unsignificant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  #first the dashed line with all points
  #geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K, color="blue", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(Cor_mean_autocor_pVal_Significant, Cor_mean_autocor, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  scale_color_manual(name="", values=c("#00B5EC"), labels="Mean clonal turnover")+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ s_driver_birth,
              labeller =label_bquote(cols= s== .(s_driver_birth)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n  Mean Clonal turnover vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/CorrelationAverageGrowthRateWithMeanSweep_K512_square_CombinedMutationRate.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3  with inverse waiting time       ##### 

########## Fig3.a      ########## 

tmpSummary<-InverseWaitingTime_cor_summary_CombinedFitnessEffectK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}



gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variable vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedFitnessEffect_AllDiversityMeasuresWithSweep_square_randomS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


#remove sewwp and mean number of driver
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  # 
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variable vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedFitnessEffect_AllDiversityMeasuresNOSweepNOMeanNbDriver_square_randomS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3b (Rob: Maybe an old version? New version is on line 6965)  #####     


tmpSummary<-InverseWaitingTime_cor_summary_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("blue"="#00B5EC", "green"="#009F71", "red"="#FAC200", "yellow"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_AllDiversityMeasures_K512_square_NewCorrelation.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()




# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
# scale_shape_identity()+
#first the dashed line with all points
geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DiversityMeasuresWithSweep_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

# without Sweep values and mean number of drivers
gg<-ggplot(tmpSummary)+
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
  #geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  #geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  #geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DiversityMeasuresNOSweepNOMeanDriverNb_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

############ Figure 3c      ############ 
tmpSummary<-InverseWaitingTime_cor_summary_CombinedMutationRateK512

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==1 & SizeIt ==1)  ){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Clonal diversity vs inverse waiting time"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

##### Figure 3d  #####     

tmpSummary<-InverseWaitingTime_cor_summary_CombinedMutationRateK512

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Clonal diversity vs inverse waiting time"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3  with inverse waiting time   and without s=0.15    ##### 
##### Figure 3b  #####     


tmpSummary<-subset(InverseWaitingTime_cor_summary_CombinedMutationRateK512, 
                   ! InverseWaitingTime_cor_summary_CombinedMutationRateK512$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("blue"="#00B5EC", "green"="#009F71", "red"="#FAC200", "yellow"="#F35E00"), labels=LabelsDiversityMeasure)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_AllDiversityMeasures_K512_square_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()




# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DiversityMeasuresWithSweep_K512_square_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

# without Sweep values and mean number of drivers
gg<-ggplot(tmpSummary)+
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
  #geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  #geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  #geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs inverse waiting time"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DiversityMeasuresNOSweepNOMeanDriverNb_K512_square_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

############ Figure 3c      ############ 
tmpSummary<-subset(InverseWaitingTime_cor_summary_CombinedMutationRateK512, 
                   ! InverseWaitingTime_cor_summary_CombinedMutationRateK512$s_driver_birth==0.15)

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==1 & SizeIt ==1)  ){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Clonal diversity vs inverse waiting time"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMuNOS015.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

##### Figure 3d  #####     
tmpSummary<-subset(InverseWaitingTime_cor_summary_CombinedMutationRateK512, 
                   ! InverseWaitingTime_cor_summary_CombinedMutationRateK512$s_driver_birth==0.15)

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Clonal diversity vs inverse waiting time"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/InverseWaitTime_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMuNOS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
##### Figure 3  with Average growth rate      ##### 

########## Fig3.a      ########## 

tmpSummary<-AverageGrowthRate_cor_summary_CombinedFitnessEffectK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}



gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #scale_x_continuous(name="tumour size at measurement", breaks =c(62500, 125000, 250000, 375000, 500000, 750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Predictor variable vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)





png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedFitnessEffect_AllDiversityMeasuresWithSweep_square_randomS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

point_size = 2

#remove sewwp and mean number of driver
gg<-ggplot(tmpSummary)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c"), linetype="solid")+
  
  # #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), 
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K, 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  #then all points which are  significant
  geom_point(aes(x= start_size, 
                 y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K , 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), 
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), 
                 group = K , 
                 fill="c", 
                 stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
  
  scale_size_manual(point_size)+
  guides(linetype=FALSE,
         color=FALSE,
         fill = guide_legend(override.aes = list(stroke = 0))) + 
  
  continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
                   palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                   breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
  
  theme_bw(base_size = 12)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n predictor variables vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)+
  facet_grid( . ~ mu_driver_birth,labeller =label_bquote(cols = mu == .(mu_driver_birth)) )






png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedFitnessEffect_AllDiversityMeasuresNOSweepNOMeanNbDriver_square_randomS_revised.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### Figure 3b  #####     


tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DiversityMeasuresWithSweep_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

gg <- three_predictors_plot(tmpSummary) + 
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))

# without Sweep values and mean driver number
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  # now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  # 
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DiversityMeasuresNOSweepNOMeanDriverNb_K512_square_revised.png"),  width = 1500, height = 1000, res = 100)
print(gg)
dev.off()


############ Figure 3c      ############ 
tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==1 & SizeIt ==1)  ){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

##### Figure 3d  #####     

tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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

gg <- biopsy_size_plot(tmpSummary) + facet_grid( . ~ s_driver_birth,labeller =label_bquote( cols= s== .(s_driver_birth)))

png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu_revised.png"),  width = 1500, height = 1000, res = 100)
print(gg)
dev.off()


##### Figure 3  with Average growth rate without panel s=0.15     ##### 
##### Figure 3b  #####     


tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   ! AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.15)


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DiversityMeasuresWithSweep_K512_square_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()



# without Sweep values and mean driver number
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= start_size, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  # now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= start_size, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  # 
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DiversityMeasuresNOSweepNOMeanDriverNb_K512_square_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


############ Figure 3c      ############ 
tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==1 & SizeIt ==1)  ){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
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
    ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu_NoS015.png"),  width = 1000, height = 1000, res = 100)
  print(gg)
  dev.off()
}

##### Figure 3d  #####     

tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512

Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu_NoS015.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()


#### Forecasting with fixeed start size and varying final size ########


tmpSummary<-AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
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
  # #geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099"), labels=LabelsDiversityMeasureWithSweep)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  scale_x_continuous(name="Final tumour size", breaks =c(250000, 375000, 500000, 625000,  750000, 875000, 1E6))+
  ylab(paste0("correlation coefficient:\n Predictor variables at start_size=125000 vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_DiversityMeasuresWithSweep_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

tmpSummary<-AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweepAndMuDriverBirth)){
  print(DiversityMeasureWithSweepAndMuDriverBirth[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweepAndMuDriverBirth[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweepAndMuDriverBirth[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweepAndMuDriverBirth[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweepAndMuDriverBirth[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # #geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(mu_driver_birth_Significant0.05, Cor_mu_driver_birth, as.logical(NA)), group = K , color="f", shape ="solid"), size=2)+
  
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_mu_driver_birth, group = K , color="f", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(mu_driver_birth_Significant0.05, Cor_mu_driver_birth, as.logical(NA)), group = K , color="f", linetype="solid"))+
  
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00","e"="#990099", "f"="#B22222"), labels=LabelsDiversityMeasureWithSweepAndMuDriverBirth)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  scale_x_continuous(name="Final tumour size", breaks =c(250000, 375000, 500000, 625000,  750000, 875000, 1E6))+
  ylab(paste0("correlation coefficient:\n Predictor variables at start_size=125000 vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)), 
         shape= guide_legend(title.theme = element_text(size = 12),
                             label.theme = element_text(size=10)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1),
        aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000_CombinedMutationRate_DiversityMeasuresWithSweepAndMutationRate_K512_square.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

##### to see how parameter values affect the strength of correlations      ##### 

All_AverageGrowthRate_summary_NewCorrelation_K64_SS250000<-subset(All_AverageGrowthRate_summary_NewCorrelation_K64, 
                                                                  All_AverageGrowthRate_summary_NewCorrelation_K64$start_size==250000
)

All_AverageGrowthRate_summary_NewCorrelation_K512_SS250000<-subset(All_AverageGrowthRate_summary_NewCorrelation_K512, 
                                                                   All_AverageGrowthRate_summary_NewCorrelation_K512$start_size==250000)

All_AverageGrowthRate_summary_NewCorrelation_K4096_SS250000<-subset(All_AverageGrowthRate_summary_NewCorrelation_K4096, 
                                                                    All_AverageGrowthRate_summary_NewCorrelation_K4096$start_size==250000)


remove(All_AverageGrowthRate_summary_NewCorrelation_K64)
remove(All_AverageGrowthRate_summary_NewCorrelation_K512)
remove(All_AverageGrowthRate_summary_NewCorrelation_K4096)

All_AverageGrowthRate_summary_NewCorrelation_K64_SS250000$mean_AverageGrowthRate<-NULL
All_AverageGrowthRate_summary_NewCorrelation_K4096_SS250000$mean_AverageGrowthRate<-NULL


df<-rbind(All_AverageGrowthRate_summary_NewCorrelation_K64_SS250000, All_AverageGrowthRate_summary_NewCorrelation_K512_SS250000)
df<-rbind(df, All_AverageGrowthRate_summary_NewCorrelation_K4096_SS250000)

sum_df <- df%>%
  group_by(K, mu_driver_birth, s_driver_birth) %>%
  summarise(MedianSweepValue = median(mean_autocor),
            MedianCorDriverDiversity = median(Cor_DriverDiversity), 
            FirstQuantileSweepValue =quantile(mean_autocor,probs =0.25),
            ThirdQuantileSweepValue =quantile(mean_autocor,probs =0.75),
            FirstQuantileCorDriverDiversity =quantile(Cor_DriverDiversity,probs =0.25),
            ThirdQuantileCorDriverDiversity =quantile(Cor_DriverDiversity,probs =0.75))

# gg<-ggplot(data = sum_df, aes(x = MedianSweepValue, y = MedianCorDriverDiversity)) +
#   geom_point()+
#   geom_errorbar(aes(ymin=FirstQuantileCorDriverDiversity, ymax=ThirdQuantileCorDriverDiversity), width=.2,position=position_dodge(.9))
#        
# png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity)) +
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

#####

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(K) )) +
  geom_point()+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByK_facetS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

#####
gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(K))) +
  geom_point()+
  facet_grid( . ~ mu_driver_birth, labeller=label_bquote(cols= mu== .(mu_driver_birth)))


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByK_facetMu.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(mu_driver_birth) )) +
  geom_point()+
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByMu_facetS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####
gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(mu_driver_birth))) +
  geom_point()+
  facet_grid( . ~ K, labeller=label_bquote(cols= K== .(K)))


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByMu_facetK.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####
gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(K))) +
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)

png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByK.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(s_driver_birth))) +
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByS.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, color= as.factor(mu_driver_birth))) +
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_colorByMu.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()
#####


gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(s_driver_birth), 
                          fill= as.factor(s_driver_birth), 
                          shape = as.factor(mu_driver_birth),
                          size = as.factor(log10(K)))) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Mean driver fitness effect",values=c("#00B5EC","#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Mean driver fitness effect",values=c("#00B5EC","#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Driver mutation rate", values = c(21, 22, 24) )+
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_AllParamShowed.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()



gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(s_driver_birth), 
                          fill= as.factor(s_driver_birth), 
                          shape = as.factor(mu_driver_birth),
                          size = as.factor(log10(K)), alpha=0.1)) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Mean driver fitness effect",values=c("#00B5EC","#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Mean driver fitness effect",values=c("#00B5EC","#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Driver mutation rate", values = c(21, 22, 24) )+
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_AllParamShowed_Alpha.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth)
                          #,size = as.factor(log10(K))
                          )) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) + 
  #scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  geom_point(size = 4)+
  xlab("mean clonal turnover index")+
  ylab("correlation coefficient:\n clonal diversity vs average growth rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1) +
  facet_grid( . ~ K, labeller=label_bquote(cols= K== .(K)))


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_AllParamShowed_ColorMu_FacetK.png"),  width = 1500, height =1000, res = 100)
print(gg)
dev.off()

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth),
                          size = as.factor(log10(K)), alpha=0.3)) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  geom_point()+
  xlab("Mean Sweep value")+
  ylab("correlation coefficient:\n Clonal diversity vs average growth rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MedianMeanSweepValueVsCorDriverDiversity_AllParamShowed_ColorMu_Alpha.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()


##### Last update of the figure      ##### 
########## Fig1.b      ########## 

tmpSummary<-All_AverageGrowthRate_summary_NewCorrelation_K512[which( (All_AverageGrowthRate_summary_NewCorrelation_K512$mu_driver_birth %in% c(1e-05)) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$K==512) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$s_driver_birth %in% c(0.1)))
                                                              ,]


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

point_size = 4

gg<-ggplot(tmpSummary)+
  # dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="a"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="b"), linetype="solid")+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="c"), linetype="solid")+
  
  # #first all points which are not significant
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), 
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K, 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  #then all points which are  significant
  geom_point(aes(x= start_size, 
                 y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),
                 group = K, 
                 fill="a", 
                 stroke=DriverDiversity_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), 
                 group = K , 
                 fill="b", 
                 stroke=DriverEdgeDiversity_Significant0.05 * (thick_width - thin_width) + thin_width), 
             shape=21, size=point_size)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), 
                 group = K , 
                 fill="c", 
                 stroke=MeanBirthRate_Significant0.05 * (thick_width - thin_width) + thin_width),
             shape=21, size=point_size)+
  
  scale_fill_manual(name="Predictor variable", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="solid", "solid"="solid"))+
  scale_color_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"))+
  
  scale_size_manual(point_size)+
  guides(linetype=FALSE,
         color=FALSE,
         fill = guide_legend(override.aes = list(stroke = 0))) + 
  
  continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel), 
                   palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                   breaks = c(thin_width, thick_width), labels = c("no", "yes"))+
  
  theme_bw(base_size = 16)+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n predictor variables vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)




# png(paste0(output_dir_plots_paper, "/wait_correlations_K512mu1e-5s01_AllDiversityMeasures_AllStartSize_muLabelNoTitle_CorrelationAverageGrowthRate.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()

png(paste0(output_dir_plots_paper, "/wait_correlations_K512mu1e-5s01_AllDiversityMeasures_AllStartSize_muLabelNoTitle_CorrelationAverageGrowthRate_revised.png"),  width = 800, height = 800, res = 100)
print(gg)
dev.off()


###### Fig 2 ######
data<-subset(dataNewCorrelation,
             (dataNewCorrelation$K==512 &
                dataNewCorrelation$s_driver_birth==0.05 &
                dataNewCorrelation$mu_driver_birth==1e-04 ))

remove(dataNewCorrelation)

# take the value around 750000 for which both the MeanBirthRate and the DriverDiversity are given 
data<-subset(data, data$NumCells <= 750000)
data<-data[! is.na(data$MeanBirthRate), ]
data<-data[! is.na(data$DriverDiversity), ]
data<-data %>% group_by(SimulationNb) %>% mutate(maxNumCell =max(NumCells)) %>% mutate(isMax = (NumCells==maxNumCell)) %>% ungroup()

gg<-ggplot(subset(data, 
                  data$isMax==TRUE))+
  geom_point(aes(x=DriverDiversity, y=MeanBirthRate))+
  xlab("Clonal diversity")+
  ylab("Mean cell division rate")+
  theme_bw(base_size = 24)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/NegCorrelationClonalDiverMeanCellDivRate_LargerSize.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()


#test
datatest<-subset(dataNewCorrelation,
                 (dataNewCorrelation$K==512 &
                    dataNewCorrelation$s_driver_birth==0.05 &
                    dataNewCorrelation$mu_driver_birth==1e-04 ))

remove(dataNewCorrelation)
data<-datatest[! is.na(datatest$MeanBirthRate), ]
data<-data[! is.na(data$DriverDiversity), ]
data<-data %>% group_by(SimulationNb) %>% mutate(Diff750000 =abs(750000-NumCells)) %>% mutate(minDiff = min(Diff750000)) %>% mutate(isMinDiff = (Diff750000==minDiff)) %>% ungroup()

as.data.frame(data[which(data$SimulationNb==2022), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])
as.data.frame(data[which(data$SimulationNb==2098), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])

gg<-ggplot(subset(data, 
                  data$isMinDiff==TRUE))+
  geom_point(aes(x=DriverDiversity, y=MeanBirthRate))+
  xlab("Clonal diversity")+
  ylab("Mean cell division rate")+
  theme_bw(base_size = 24)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/NegCorrelationClonalDiverMeanCellDivRate_LargerSize_UsingMinDiff.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()

FileWithAllParam<-readSeedFromParamFiles("/all_results/Batch_ForecastingPaper_bis")

dataMinDiff<-subset(data, 
                    data$isMinDiff==TRUE)
ParamPanelb<-dataMinDiff[which(dataMinDiff$MeanBirthRate==max(dataMinDiff$MeanBirthRate)), c(1:18)]
ParamPanelb<-merge(ParamPanelb, FileWithAllParam)
PathPanelb<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_2/s_driver_birth_0/seed_", ParamPanelb$seedFolder, "/")
Check<-fread(paste0(PathPanelb, "parameters.dat"))                         
#plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE)
plot_all_images_cutoff(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)
plot_all_images_cutoff(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff01", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.1)
plot_all_images_cutoff_modified(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff01", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.1)



ParamPanelc<-dataMinDiff[which(dataMinDiff$DriverDiversity==max(dataMinDiff$DriverDiversity)), c(1:18)]
ParamPanelc<-merge(ParamPanelc, FileWithAllParam)
PathPanelc<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_2/s_driver_birth_0/seed_", ParamPanelc$seedFolder, "/")
Check<-fread(paste0(PathPanelc, "parameters.dat"))                         
#plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE)
plot_all_images_cutoff(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)
plot_all_images_cutoff_modified(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff001", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0.01)



##### Figure 3  with Average growth rate only panel s=0.1     ##### 
##### Figure 3b  #####     


tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.1)


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

#only clonal diversity and clonal diversity at edge and mean cell division rate
gg <- three_predictors_plot(tmpSummary)

png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DiversityMeasuresNOSweepNOMeanDriverNb_K512_square_Onlys01_revised.png"),  width = 800, height = 800, res = 100)
print(gg)
dev.off()

#### fig 3c       #### 
tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.1)

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-1
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    #if( (Depth[DepthIt]==1 & SizeIt ==1)  ){
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="dashed")+
                   aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("first geom point")
    #}
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( . ~ s_driver_birth,
                labeller =label_bquote( cols= s== .(s_driver_birth)))+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)), 
           shape= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=12)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu_OnlyS01.png"),  width = 700, height = 700, res = 100)
  print(gg)
  dev.off()
}

##### Figure 3d  #####     

tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.1)


Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour boundary", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), 
        aspect.ratio=1)

gg <- biopsy_size_plot(tmpSummary)

png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu_onlyS01_revised.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()


#### fig 3c with depth 4   #### 
# tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
#                    AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.15)

tmpSummary<-AverageGrowthRate_cor_summary_CombinedMutationRateK512 

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-4
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    if( !(SampleSize[SizeIt] ==4)  ){
      gg<-gg+
        geom_point(aes(x=start_size,group = K, shape ="dashed")+
                     aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                           ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                                 color=shQuote(Denomination[SizeIt])),
                   size=2)
      print("first geom point")
    }
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1,000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( . ~ s_driver_birth,
                labeller =label_bquote( cols= s== .(s_driver_birth)))+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
           shape= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=12)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  # png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu_NoS015.png"),  width = 1000, height = 1000, res = 100)
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu_allS.png"),  width = 1000, height = 1000, res = 100)
  
  print(gg)
  dev.off()
}

#### fig 3c with depth 4    #### 
tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.1)

SampleSize<-c("1","4", "1Big", "4Big")
Depth<-4
Denomination<-c("a", "b", "c", "d")


for(DepthIt in seq_along(Depth)){
  
  for(SizeIt in seq_along(SampleSize)){
    
    tmpSummary<-mutate(tmpSummary,
                       tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt] ,"_pVal")]] <= SignificanceLevel))
    colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel )
  }
  
  gg<-ggplot(tmpSummary)
  
  for(SizeIt in seq_along(SampleSize)){
    
    if((SampleSize[SizeIt] ==1)  ){
      gg<-gg+
        geom_point(aes(x=start_size,group = K, shape ="dashed")+
                     aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                           ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                                 color=shQuote(Denomination[SizeIt])),
                   size=2)
      print("first geom point")
    }
    
    
    gg<-gg+
      geom_point(aes(x=start_size,group = K, shape ="solid")+
                   aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                         ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
                               color=shQuote(Denomination[SizeIt])),
                 size=2)
    print("second geom point")
    
    # gg<-gg+
    #   geom_line(aes(x= start_size, group = K, linetype="dashed")+
    #               aes_string( y =paste0("ifelse(! ", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
    #                                     ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"),
    #                           color=shQuote(Denomination[SizeIt])))
    gg<-gg+
      geom_line(aes(x= start_size, group = K, linetype="dashed")+
                  aes_string( y =paste0("Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt]),
                              color=shQuote(Denomination[SizeIt])))
    print("first geom line")
    gg<-gg+
      geom_line(aes(x= start_size, group = K , linetype="solid")+
                  aes_string( y =paste0("ifelse(", "Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],  "_Significant",SignificanceLevel,
                                        ", Cor_DriverDiversityFrom",SampleSize[SizeIt], "SamplesAtDepth",Depth[DepthIt],", as.logical(NA))"), 
                              color=shQuote(Denomination[SizeIt])))
    print("second geom line")
    
  }
  
  print("add scales")
  gg<-gg+
    scale_color_manual(name="Biopsy size",
                       values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00"),
                       labels=c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1, 000 cells"))+
    scale_linetype_manual(name="",
                          values = c("dashed"="dashed", "solid"="solid"))+
    scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                       values = c("dashed"=17, "solid"=16),
                       labels=c("dashed"="FALSE", "solid"="TRUE"))+
    facet_grid( . ~ s_driver_birth,
                labeller =label_bquote( cols= s== .(s_driver_birth)))+
    theme_bw(base_size = 16)+
    guides(linetype=FALSE,
           color= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=14)),
           shape= guide_legend(title.theme = element_text(size = 14),
                               label.theme = element_text(size=12)))+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
    ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
          aspect.ratio=1)
  
  
  png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFromAllSamplesSize_Depth", Depth[DepthIt], "_K512_NewCorrelation_RandomMu_OnlyS01.png"),  width = 700, height = 700, res = 100)
  print(gg)
  dev.off()
}


###### Fig 2  with mu=1e-06, s=0.05######

data_subset<-subset(dataNewCorrelation,
                    (dataNewCorrelation$K==512 &
                       dataNewCorrelation$s_driver_birth==0.05 &
                       dataNewCorrelation$mu_driver_birth==1e-06 ))

remove(dataNewCorrelation)

data<-data_subset[! is.na(data_subset$MeanBirthRate), ]
data<-data[! is.na(data$DriverDiversity), ]
data<-data %>% group_by(SimulationNb) %>% mutate(Diff750000 =abs(750000-NumCells)) %>% mutate(minDiff = min(Diff750000)) %>% mutate(isMinDiff = (Diff750000==minDiff)) %>% ungroup()

# as.data.frame(data[which(data$SimulationNb==2022), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])
# as.data.frame(data[which(data$SimulationNb==2098), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])

gg<-ggplot(subset(data, 
                  data$isMinDiff==TRUE))+
  geom_point(aes(x=DriverDiversity, y=MeanBirthRate))+
  xlab("Clonal diversity")+
  ylab("Mean cell division rate")+
  theme_bw(base_size = 24)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/NegCorrelationClonalDiverMeanCellDivRate_LargerSize_UsingMinDiff_mu1e-6s005.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()

FileWithAllParam<-readSeedFromParamFiles("/all_results/Batch_ForecastingPaper_bis")

dataMinDiff<-subset(data, 
                    data$isMinDiff==TRUE)
ParamPanelb<-dataMinDiff[which(dataMinDiff$MeanBirthRate==max(dataMinDiff$MeanBirthRate)), c(1:18)]
ParamPanelb<-merge(ParamPanelb, FileWithAllParam)
PathPanelb<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_0/s_driver_birth_0/seed_", ParamPanelb$seedFolder, "/")
#/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_0/s_driver_birth_0/seed_52/
Check<-fread(paste0(PathPanelb, "parameters.dat"))                         
plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff1e-3", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-3)
plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff1e-1", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-1)
plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_NOCutoff", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0)
plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_NOCutoff_genotypePlots", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE, cutoff = 0)




ParamPanelc<-dataMinDiff[which(dataMinDiff$DriverDiversity==max(dataMinDiff$DriverDiversity)), c(1:18)]
ParamPanelc<-merge(ParamPanelc, FileWithAllParam)
PathPanelc<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_0/s_driver_birth_0/seed_", ParamPanelc$seedFolder, "/")
#/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_0/s_driver_birth_0/seed_22/
Check<-fread(paste0(PathPanelc, "parameters.dat"))                         
plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-3", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-3)
plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_NoCutoff", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0)
plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-1", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-1)
plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_NoCutoff_genotypePlots", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE, cutoff = 0)


#### Forecasting with fixeed start size and varying final size ########


tmpSummary<-AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512


for(DivMeas in seq_along(DiversityMeasureWithSweep)){
  print(DiversityMeasureWithSweep[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasureWithSweep[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasureWithSweep[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}

gg <- three_predictors_plot(tmpSummary, x_axis = "FinalSize") + 
  facet_grid( . ~ s_driver_birth, labeller=label_bquote(cols= s== .(s_driver_birth)))

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # #geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  # 
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512_square_revised.png"),  width = 1500, height = 1000, res = 100)
print(gg)
dev.off()


tmpSummary<-subset(tmpSummary, 
                   tmpSummary$s_driver_birth==0.1)

gg <- three_predictors_plot(tmpSummary, x_axis = "FinalSize")

# with Sweep values
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # #geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  #geom_point(aes(x= FinalSize, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape ="solid"), size=2)+
  geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape ="solid"), size=2)+
  # geom_point(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", shape ="solid"), size=2)+
  # 
  #first all points which are not significant
  # geom_point(aes(x= FinalSize, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =17))+
  # geom_point(aes(x= FinalSize, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="a", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", shape =16))+
  # geom_point(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= FinalSize, y =Cor_DriverDiversity, group = K, color="a", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_DriverEdgeDiversity, group = K , color="b", linetype="dashed"))+
  geom_line(aes(x= FinalSize, y =Cor_MeanBirthRate, group = K , color="c", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_Drivers, group = K , color="d", linetype="dashed"))+
  # geom_line(aes(x= FinalSize, y =Cor_mean_autocor, group = K , color="e", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= FinalSize, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="a", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= FinalSize, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="c", linetype="solid"))+
  # geom_line(aes(x= FinalSize, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="d", linetype="solid"))+
  # geom_line(aes(x= FinalSize, y =ifelse(mean_autocor_Significant0.05, Cor_mean_autocor, as.logical(NA)), group = K , color="e", linetype="solid"))+
  # 
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200"), labels=LabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
  ylab(paste0("correlation coefficient:\n Predictor variables vs average growth rate"))+
  ggtitle("")+
  guides(linetype=FALSE,
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512S01_square_revised.png"),  width = 800, height = 800, res = 100)
print(gg)
dev.off()

png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512S01_square_Size800.png"),  width = 800, height = 800, res = 100)
print(gg)
dev.off()
png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512S01_square_Size900.png"),  width = 900, height = 900, res = 100)
print(gg)
dev.off()
png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512S01_square_Size1000.png"),  width = 1000, height = 1000, res = 100)
print(gg)
dev.off()

png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlationsWithAllFinalSize_fromStartSize125000__CombinedMutationRate_3DiversityMeasures_K512S01_square_Size850.png"),  width = 850, height = 850, res = 100)
print(gg)
dev.off()


##### correlation meanBrithRate     and clonal diversity  ##### 
Kvalue=c(64,512, 4096)
for(Kit in seq_along(c(0:2)) ){
  
  Kval=Kvalue[Kit]
  ClonalDiversityVsMeanBirthRate_correlation_tmp<-read_csv( file = paste0(output_dir_data, "/ClonalDiversityVsMeanBirthRate_correlation_K",Kval,".csv"),  guess_max = 1E4)
  
  ClonalDiversityVsMeanBirthRate_correlation_tmp<-subset(ClonalDiversityVsMeanBirthRate_correlation_tmp, 
                                                         ClonalDiversityVsMeanBirthRate_correlation_tmp$start_size==250000)
  
  if(Kit==1){
    df<-ClonalDiversityVsMeanBirthRate_correlation_tmp
  }else{
    df<-rbind(df, ClonalDiversityVsMeanBirthRate_correlation_tmp)
  }
  
}

# df <- df %>%
#   group_by(K, mu_driver_birth, s_driver_birth) %>%
#   summarise(MedianSweepValue = median(mean_autocor)) %>% ungroup()

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth),
                          size = as.factor(log10(K)), alpha=0.3)) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  geom_point()+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("Mean clonal turnover index")+
  ylab("correlation coefficient:\n Clonal diversity vs Mean cell division rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanSweepValueVsCOrClonalDiversityMeanBirthRate_AllParamShowed_ColorMu_Alpha.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth),
                          size = as.factor(log10(K)))) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  geom_point()+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("Mean clonal turnover index")+
  ylab("correlation coefficient:\n Clonal diversity vs Mean cell division rate")+
  theme_bw(base_size = 20)+
  guides(fill=FALSE,  
         color= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=7)), 
         shape= guide_legend(title.theme = element_text(size = 18),
                             label.theme = element_text(size=16),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 18),
                            label.theme = element_text(size=16)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanSweepValueVsCOrClonalDiversityMeanBirthRate_AllParamShowed_ColorMu.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()

###### Fig 2  with mu=1e-06, s=0.05######

data_subset<-subset(dataNewCorrelation,
                    (dataNewCorrelation$K==512 &
                       dataNewCorrelation$s_driver_birth==0.2 &
                       dataNewCorrelation$mu_driver_birth==1e-05 ))

remove(dataNewCorrelation)

data<-data_subset[! is.na(data_subset$MeanBirthRate), ]
data<-data[! is.na(data$DriverDiversity), ]
data<-data %>% group_by(SimulationNb) %>% mutate(Diff750000 =abs(750000-NumCells)) %>% mutate(minDiff = min(Diff750000)) %>% mutate(isMinDiff = (Diff750000==minDiff)) %>% ungroup()

# as.data.frame(data[which(data$SimulationNb==2022), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])
# as.data.frame(data[which(data$SimulationNb==2098), which(colnames(data) %in% c("NumCells", "Diff750000", "minDiff", "isMinDiff"))])

gg<-ggplot(subset(data, 
                  data$isMinDiff==TRUE))+
  geom_point(aes(x=DriverDiversity, y=MeanBirthRate))+
  xlab("Clonal diversity")+
  ylab("Mean cell division rate")+
  theme_bw(base_size = 24)+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/NegCorrelationClonalDiverMeanCellDivRate_LargerSize_UsingMinDiff_mu1e-5s02.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()


FileWithAllParam<-readSeedFromParamFiles("/all_results/Batch_ForecastingPaper_bis")

dataMinDiff<-subset(data, 
                    data$isMinDiff==TRUE)
ParamPanelb<-dataMinDiff[which(dataMinDiff$MeanBirthRate==max(dataMinDiff$MeanBirthRate)), c(1:18)]
ParamPanelb<-merge(ParamPanelb, FileWithAllParam)
PathPanelb<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_", ParamPanelb$seedFolder, "/")
#PathPanelb<-"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_15/"
Check<-fread(paste0(PathPanelb, "parameters.dat"))                         

plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_NOCutoff_genotypePlots_mu1e-5s02", file_type = "png", 
                output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE, cutoff = 0,min_birth_rate = 1, max_birth_rate = 10)

# plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff1e-3", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-3)
# plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_Cutoff1e-1", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-1)
# plot_all_images(PathPanelb, output_filename = "MullerPlots_largeMeanBirthRate_smallDriverDiversity_NOCutoff", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0)
# 



ParamPanelc<-dataMinDiff[which(dataMinDiff$DriverDiversity==max(dataMinDiff$DriverDiversity)), c(1:18)]
ParamPanelc<-merge(ParamPanelc, FileWithAllParam)
PathPanelc<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_", ParamPanelc$seedFolder, "/")
"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_28/"
Check<-fread(paste0(PathPanelc, "parameters.dat"))                         
# plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-3", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-3)
# plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_NoCutoff", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0)
# plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-1", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-1)
plot_all_images(PathPanelc, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_NoCutoff_genotypePlots_mu1e-5s02", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE, cutoff = 0)

#ParamPaneld<-dataMinDiff[which(dataMinDiff$MeanBirthRate==min(dataMinDiff$MeanBirthRate)), c(1:18)]

ParamPaneld<-dataMinDiff[which(dataMinDiff$DriverDiversity==sort(dataMinDiff$DriverDiversity,partial=(length(dataMinDiff$DriverDiversity)-1) )[(length(dataMinDiff$DriverDiversity)-1)]), c(1:18)]
ParamPaneld<-merge(ParamPaneld, FileWithAllParam)
PathPaneld<-paste0("/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_", ParamPaneld$seedFolder, "/")
PathPaneld<-"/all_results/Batch_ForecastingPaper_bis/log2_deme_carrying_capacity_1/mu_driver_birth_1/s_driver_birth_3/seed_14/"
Check<-fread(paste0(PathPaneld, "parameters.dat"))                         
# plot_all_images(PathPaneld, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-3", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-3)
# plot_all_images(PathPaneld, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_NoCutoff", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 0)
# plot_all_images(PathPaneld, output_filename = "MullerPlots_largeDriverDiversity_smallMeanBirthRate_Cutoff1e-1", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = FALSE, cutoff = 1e-1)
plot_all_images(PathPaneld, output_filename = "MullerPlots_SecondlargestDriverDiversity_smallMeanBirthRate_NoCutoff_genotypePlots_mu1e-5s02", file_type = "png", output_dir = output_dir_plots_paper, trim = -1, include_genotype_plots = TRUE, cutoff = 0,min_birth_rate = 1, max_birth_rate = 10)


# compute correlation between DriverDiversity and MeanBirthRate

CorDriverDiversityMeanBirthRate<-cor(dataMinDiff$DriverDiversity, dataMinDiff$MeanBirthRate, method = "spearman", use="na.or.complete")
pvalCorDriverDiversityMeanBirthRate <- cor_pval(dataMinDiff$DriverDiversity, dataMinDiff$MeanBirthRate)
######## Fig 1b new legend      ######## 

tmpSummary<-All_AverageGrowthRate_summary_NewCorrelation_K512[which( (All_AverageGrowthRate_summary_NewCorrelation_K512$mu_driver_birth %in% c(1e-05)) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$K==512) &
                                                                       (All_AverageGrowthRate_summary_NewCorrelation_K512$s_driver_birth %in% c(0.1)))
                                                              ,]


for(DivMeas in seq_along(DiversityMeasure)){
  print(DiversityMeasure[DivMeas])
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] <= SignificanceLevel),
                     tmpVar2 =(tmpSummary[[paste0("Cor_", DiversityMeasure[DivMeas], "_pVal")]] > SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0(DiversityMeasure[DivMeas], "_Significant",SignificanceLevel )
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar2")]<-paste0(DiversityMeasure[DivMeas], "_Unsignificant",SignificanceLevel )
  
  
}
gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y = Cor_DriverDiversity, group = K, color="blue", shape ="dashed"))+
  # geom_point(aes(x= start_size, y = Cor_DriverEdgeDiversity, group = K , color="green", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", shape ="dashed"))+
  # geom_point(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", shape ="dashed"))+
  
  geom_point(aes(x= start_size, y =ifelse( DriverDiversity_Unsignificant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Unsignificant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Unsignificant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="dashed"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse( Drivers_Unsignificant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="dashed"), size=2)+
  #then all points which are  significant
  geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape ="solid"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape ="solid"), size=2)+
  #geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape ="solid"), size=2)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse( ! DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse(! DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =17))+
  # geom_point(aes(x= start_size, y =ifelse( ! Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =17))+
  # #then all points which are  significant
  # geom_point(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)), group = K, color="blue", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", shape =16))+
  # geom_point(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", shape =16))+
  # scale_shape_identity()+
  #first the dashed line with all points
  geom_line(aes(x= start_size, y =Cor_DriverDiversity, group = K, color="blue", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_DriverEdgeDiversity, group = K , color="green", linetype="dashed"))+
  geom_line(aes(x= start_size, y =Cor_MeanBirthRate, group = K , color="red", linetype="dashed"))+
  #geom_line(aes(x= start_size, y =Cor_Drivers, group = K , color="yellow", linetype="dashed"))+
  #now only the points significant at 0.05 level
  geom_line(aes(x= start_size, y =ifelse(DriverDiversity_Significant0.05, Cor_DriverDiversity, as.logical(NA)),  group = K, color="blue", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(DriverEdgeDiversity_Significant0.05, Cor_DriverEdgeDiversity, as.logical(NA)), group = K , color="green", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(MeanBirthRate_Significant0.05, Cor_MeanBirthRate, as.logical(NA)), group = K , color="red", linetype="solid"))+
  #geom_line(aes(x= start_size, y =ifelse(Drivers_Significant0.05, Cor_Drivers, as.logical(NA)), group = K , color="yellow", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  scale_color_manual(name="", values=c("#00B5EC", "#009F71", "#FAC200"), labels= NewLabelsDiversityMeasure_NoMeanDriverNb)+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
  scale_shape_manual(name=paste0("Significant at \n alpha level " , SignificanceLevel),
                     values = c("dashed"=17, "solid"=16),
                     labels=c("dashed"="FALSE", "solid"="TRUE"))+
  # scale_shape_manual(name=paste0("Significant at alpha level " , SignificanceLevel),
  #                     values = c("solid"=16),
  #                     labels=c("solid"="TRUE"))+
  scale_size_manual(2)+
  guides(linetype=FALSE)+
  theme_bw(base_size = 16)+
  facet_grid( mu_driver_birth ~ s_driver_birth, labeller =label_bquote(rows = mu == .(mu_driver_birth), 
                                                                       cols= s== .(s_driver_birth)))+
  
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  # scale_x_continuous(name="tumour size at measurement", breaks =c(10000, 60000,125000, 250000, 375000, 500000, 625000,  750000))+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n predictor variable vs average growth rate"))+
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




# png(paste0(output_dir_plots_paper, "/wait_correlations_K512mu1e-5s01_AllDiversityMeasures_AllStartSize_muLabelNoTitle_CorrelationAverageGrowthRate.png"),  width = 1000, height = 1000, res = 100)
# print(gg)
# dev.off()

png(paste0(output_dir_plots_paper, "/wait_correlations_K512mu1e-5s01_AllDiversityMeasures_AllStartSize_muLabelNoTitle_CorrelationAverageGrowthRate_NewLegend.png"),  width = 800, height = 800, res = 100)
print(gg)
dev.off()

#### fig 3d with depth 4  new legend  #### 
tmpSummary<-subset(AverageGrowthRate_cor_summary_CombinedMutationRateK512, 
                   AverageGrowthRate_cor_summary_CombinedMutationRateK512$s_driver_birth==0.1)


Depth_4Big<-c(1:4)

for(DepthIt in c(1:4)){
  
  tmpSummary<-mutate(tmpSummary,
                     tmpVar =(tmpSummary[[paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt ,"_pVal")]] <= SignificanceLevel))
  colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth",DepthIt,  "_Significant",SignificanceLevel )
  
  
}

tmpSummary<-mutate(tmpSummary,
                   tmpVar =(tmpSummary$Cor_DriverDiversityFrom4BigRandomSamples_pVal <= SignificanceLevel))
colnames(tmpSummary)[which(colnames(tmpSummary)=="tmpVar")]<-paste0("Cor_DriverDiversityFrom4BigRandomSamples_Significant",SignificanceLevel )


gg<-ggplot(tmpSummary)+
  #first all points which are not significant
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth0_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth0, as.logical(NA)), group = K, color="a", shape ="dashed"), size=2)+
  geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", shape ="dashed"), size=2)+
  # geom_point(aes(x= start_size, y =ifelse(! Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", shape ="dashed"), size=2)+
  # #then all points which are  significant
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
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth1_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth1, as.logical(NA)), group = K , color="b", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth2_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth2, as.logical(NA)), group = K , color="c", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth3_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth3, as.logical(NA)), group = K , color="d", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigSamplesAtDepth4_Significant0.05, Cor_DriverDiversityFrom4BigSamplesAtDepth4, as.logical(NA)), group = K , color="e", linetype="solid"))+
  geom_line(aes(x= start_size, y =ifelse(Cor_DriverDiversityFrom4BigRandomSamples_Significant0.05, Cor_DriverDiversityFrom4BigRandomSamples, as.logical(NA)), group = K , color="f", linetype="solid"))+
  #scale_color_manual(name="", values=colfunc(length(DiversityMeasure)), labels=DiversityMeasure)+
  # scale_color_manual(name="", values=c("a"="#00B5EC", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="#990099", "f"="black"),
  #                    labels=c("a"="Depth 0", "b"="Depth 1", "c"="Depth 2", "d"="Depth 3", "e"="Depth 4", "f"="Random Samples"))+
  # #labels=c(paste("Depth", c(0:4)), "Random Samples"))+
  scale_color_manual(name="Distance from \n tumour edge", values=c("b"="#00B5EC", "c"="#009F71", "d"="#FAC200", "e"="#F35E00", "f"="black"),
                     labels=c("b"="10%", "c"="20%", "d"="30%", "e"="40%", "f"="Random Sample"))+
  scale_linetype_manual(name="",
                        values = c("dashed"="dashed", "solid"="solid"))+
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
         color= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=14)), 
         shape= guide_legend(title.theme = element_text(size = 14),
                             label.theme = element_text(size=12)))+
  ylim(-1, 1)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_x_continuous(name="tumour size at measurement", breaks =c(10000,125000, 250000, 375000, 500000, 625000,  750000))+
  ylab(paste0("correlation coefficient:\n Clonal diversity vs average growth rate"))+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), 
        aspect.ratio=1)



png(paste0(output_dir_plots_paper, "/AverageGrowthRate_correlations_CombinedMutationRate_DriverDiversityFrom4BigSamples_DifferentDepth_Square_K512_RandomMu_onlyS01_NewLegend.png"),  width = 700, height = 700, res = 100)
print(gg)
dev.off()


####### Fig 2a after revision ##### 

Kvalue=c(64,512, 4096)
for(Kit in seq_along(c(0:2)) ){
  
  Kval=Kvalue[Kit]
  ClonalDiversityVsMeanBirthRate_correlation_tmp<-read_csv( file = paste0(output_dir_data, "/ClonalDiversityVsMeanBirthRate_correlation_K",Kval,".csv"),  guess_max = 1E4)
  
  ClonalDiversityVsMeanBirthRate_correlation_tmp<-subset(ClonalDiversityVsMeanBirthRate_correlation_tmp, 
                                                         ClonalDiversityVsMeanBirthRate_correlation_tmp$start_size==250000)
  
  if(Kit==1){
    df<-ClonalDiversityVsMeanBirthRate_correlation_tmp
  }else{
    df<-rbind(df, ClonalDiversityVsMeanBirthRate_correlation_tmp)
  }
  
}

# df <- df %>%
#   group_by(K, mu_driver_birth, s_driver_birth) %>%
#   summarise(MedianSweepValue = median(mean_autocor)) %>% ungroup()
df<-subset(df, ! df$s_driver_birth==0.15)

gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth))) +
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  geom_point(size=3)+
  facet_grid(.~K)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("mean clonal turnover index")+
  ylab("correlation coefficient:\n Clonal diversity vs mean cell division rate")+
  theme_bw(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  guides(fill=FALSE, 
         alpha=FALSE, 
         color= guide_legend(title.theme = element_text(size = 16),
                             label.theme = element_text(size=14),
                             override.aes = list(size=5)), 
         shape= guide_legend(title.theme = element_text(size = 16),
                             label.theme = element_text(size=14),
                             override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanSweepValueVsCOrClonalDiversityMeanBirthRate_AllParamShowed_ColorMu_FacetK.png"),  width = 1500, height =1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth),
                          size = as.factor(log10(K)))) +
  scale_size_manual(name= "Deme size", values= c(2,4,6), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  geom_point(colour="black", shape=21, stroke = 1)+
  facet_grid(.~s_driver_birth)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("Mean clonal turnover index")+
  ylab("correlation coefficient:\n Clonal diversity vs Mean cell division rate")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  guides(color=guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14),
                            override.aes = list(size=5)), 
         alpha=FALSE, 
         fill= guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14),
                            override.aes = list(size=5)), 
         # shape= guide_legend(title.theme = element_text(size = 16),
         #                     label.theme = element_text(size=14),
         #                     override.aes = list(size=5)),
         size= guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14)
         ))+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanSweepValueVsCOrClonalDiversityMeanBirthRate_AllParamShowed_ColorMu_FacetS.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()


gg<-ggplot(data = df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth),
                          shape = as.factor(log10(K)))) +
  scale_shape_manual(name= "Deme size", values= c(21, 22, 23), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  #geom_point(colour="black", shape=21, stroke = 1)+
  geom_point(size=3, colour="black", stroke = 0.5)+
  facet_grid(.~s_driver_birth)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("Mean clonal turnover index")+
  ylab("correlation coefficient:\n Clonal diversity vs Mean cell division rate")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  guides(color=guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14),
                            override.aes = list(size=5)), 
         alpha=FALSE, 
         fill= guide_legend(title.theme = element_text(size = 16),
                            label.theme = element_text(size=14),
                            override.aes = list(size=5)), 
         shape= guide_legend(title.theme = element_text(size = 16),
                             label.theme = element_text(size=14),
                             override.aes = list(size=5))
         # size= guide_legend(title.theme = element_text(size = 16),
         #                    label.theme = element_text(size=14))
  )+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanSweepValueVsCOrClonalDiversityMeanBirthRate_AllParamShowed_ColorMu_FacetS_DemeSizeAsShape.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()

### Check the influence of max s_driver_birth
gg<-ggplot(data = df, aes(x = mean_MeanBirthRate, y = Cor_DriverDiversity, 
                          color=as.factor(mu_driver_birth), 
                          fill= as.factor(mu_driver_birth), 
                          shape = as.factor(s_driver_birth),
                          size = as.factor(log10(K)))) +
  scale_size_manual(name= "Deme size", values= c(3,4,5), label=c("64", "512","4096"))+
  scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
  scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
  #geom_point(colour="black", shape=21, stroke = 1)+
  geom_point()+
  #facet_grid(.~s_driver_birth)+
  geom_hline(yintercept = 0, linetype ="dashed")+
  xlab("Mean cell dividion rate (computed over each set of simulations)")+
  ylab("correlation coefficient:\n Clonal diversity vs Mean cell division rate")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  guides(
    alpha=FALSE, 
    shape= guide_legend(title.theme = element_text(size = 16),
                        label.theme = element_text(size=14),
                        override.aes = list(size=5)),
    fill= guide_legend(title.theme = element_text(size = 16),
                       label.theme = element_text(size=14),
                       override.aes = list(size=5)),
    color= guide_legend(title.theme = element_text(size = 16),
                        label.theme = element_text(size=14),
                        override.aes = list(size=5))
    # size= guide_legend(title.theme = element_text(size = 16),
    #                    label.theme = element_text(size=14))
  )+
  theme(aspect.ratio=1)


png(paste0(output_dir_plots_paper, "/MeanMeanBirthRateVsCorClonalDiversityGrowthRate.png"),  width = 1000, height =1000, res = 100)
print(gg)
dev.off()



