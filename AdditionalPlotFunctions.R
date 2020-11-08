library(demonanalysis)
library(readr)
library(dplyr)
library(ggplot2)

correlations_plot <- function(plot_df, group = "predictors", x_axis = "start_size", depth = 1, facet = "none", filename = NULL, 
                              width = 800, height = 800, res = 100, SignificanceLevel = 0.05, ylab = "average growth rate") {
  
  thin_width = 0.1
  thick_width = 2
  point_size = 4
  
  if(group == "predictors") {
    yvec <- c("Cor_DriverDiversity", "Cor_DriverEdgeDiversity", "Cor_MeanBirthRate")
  } else if(group == "biopsydepth") {
    yvec <- c(paste0("Cor_DriverDiversityFrom4BigSamplesAtDepth", 1:4), "Cor_DriverDiversityFrom4BigRandomSamples")
  } else if(group == "biopsysize") {
    yvec <- paste0("Cor_DriverDiversityFrom", c("1", "4", "1Big", "4Big"), "SamplesAtDepth", depth)
  }
  lvec <- vector()
  svec <- paste0(yvec, "_Significant0.05")
  
  for(i in 1:length(yvec)){
    plot_df <- mutate(plot_df, tmpVar = (plot_df[[paste0(yvec[i], "_pVal")]] <= 0.05))
    colnames(plot_df)[which(colnames(plot_df) == "tmpVar")] <- svec[i]
  }
  
  gg <- ggplot(plot_df)
  
  for(i in 1:length(yvec)) {
    col <- paste0("'", letters[i], "'")
    len <- dim(unique(plot_df[, svec[i]]))[1]
    if(is.null(len)) {
      lvec[i] <- length(unique(plot_df[, svec[i]]))
    } else lvec[i] <- len
    
    gg <- gg + geom_line(aes_string(x = x_axis, y = yvec[i], group = "K", color = col), linetype = "solid")
    
    if(lvec[i] == 1) {
      gg <- gg + geom_point(aes_string(x = x_axis,
                                       y = yvec[i],
                                       group = "K",
                                       fill = col,
                                       stroke = paste0(svec[i], " * (thick_width - thin_width) + thin_width")),
                            shape=21, size=point_size)
    } else {
      gg <- gg + geom_point(aes_string(x = x_axis,
                                       y = paste0("ifelse(!", svec[i], ", ", yvec[i], ", as.logical(NA))"),
                                       group = "K",
                                       fill = col,
                                       stroke = paste0(svec[i], " * (thick_width - thin_width) + thin_width")),
                            shape=21, size=point_size) +
        geom_point(aes_string(x = x_axis,
                              y = paste0("ifelse(", svec[i], ", ", yvec[i], ", as.logical(NA))"),
                              group = "K",
                              fill = col,
                              stroke = paste0(svec[i], " * (thick_width - thin_width) + thin_width")),
                   shape=21, size=point_size)
    }
  }
  
  if(prod(lvec) == 1) {
    gg <- gg + continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel),
                                palette = function(x){scales::rescale(x, to=thick_width / 2, from =thick_width / 2)},
                                breaks = thick_width / 2, labels = c("yes"))
  } else {
    gg <- gg + continuous_scale("stroke", scale_name = "", name = paste0("Significant at \n alpha level " , SignificanceLevel),
                                palette = function(x){scales::rescale(x, to=c(thin_width, thick_width), from =c(thin_width, thick_width))},
                                breaks = c(thin_width, thick_width), labels = c("no", "yes"))
  }
  
  if(group == "predictors") {
    col_labs <- c("a"="mediumpurple", "b"="#009F71", "c"="#FAC200")
    labels_vec <- c("Clonal diversity", "Clonal diversity at boundary", "Mean cell division rate")
    gg <- gg + 
      scale_fill_manual(name="Predictor variable", values = col_labs, labels = labels_vec)+
      scale_color_manual(name=paste0("Significant at \n alpha level ", SignificanceLevel), values = col_labs)+
      ylab(paste0("correlation coefficient:\n predictor variables vs ", ylab))
  } else if(group == "biopsydepth") {
    col_labs <- c("a"="mediumpurple", "b"="#009F71", "c"="#FAC200", "d"="#F35E00", "e"="gray30")
    labels_vec <- c("a"="10%", "b"="20%", "c"="30%", "d"="40%", "e"="Random sample")
    gg <- gg +
      scale_fill_manual(name="Distance from \n tumour boundary", values=col_labs, labels = labels_vec)+
      scale_color_manual(name=paste0("Significant at \n alpha level ", SignificanceLevel), values = col_labs)+
      ylab(paste0("correlation coefficient:\n clonal diversity vs ", ylab))
  } else if(group == "biopsysize") {
    col_labs <- c("a"="mediumpurple", "b"="#009F71", "c"="#FAC200", "d"="#F35E00")
    labels_vec <- c("a"="1 biopsy of 100 cells", "b"="1 biopsy of 1,000 cells", "c"="4 biopsies of 100 cells","d"="4 biopsies of 1,000 cells")
    gg <- gg +
      scale_fill_manual(name="Biopsy size", values=col_labs, labels = labels_vec)+
      scale_color_manual(name=paste0("Significant at \n alpha level ", SignificanceLevel), values = col_labs)+
      ylab(paste0("correlation coefficient:\n clonal diversity vs ", ylab))
  }
  
  if(x_axis == "start_size") {
    gg <- gg + scale_x_continuous(name="tumour size at measurement", breaks = c(10000, 125000*1:6))
  } else if(x_axis == "FinalSize") {
    gg <- gg + scale_x_continuous(name="endpoint size", breaks = 125000*2:8)
  } else if(x_axis == "s_driver_birth") {
    gg <- gg + scale_x_continuous(name="driver mutation effect", breaks = 0.05*(0:4))
  }
  
  gg <- gg + scale_size_manual(point_size)+
    scale_linetype_manual(name="", values = c("dashed"="solid", "solid"="solid"))+
    guides(linetype=FALSE,
           color=FALSE,
           fill = guide_legend(override.aes = list(stroke = 0))) + 
    theme_bw(base_size = 16)+
    ylim(-1, 1)+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    ggtitle("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1),
          aspect.ratio=1)
  
  if(facet == "mu") {
    gg <- gg + facet_grid(. ~ mu_driver_birth, labeller = label_bquote(cols = mu == .(mu_driver_birth)))
  } else if(facet == "s") {
    gg <- gg + facet_grid(. ~ s_driver_birth, labeller = label_bquote(cols = s == .(s_driver_birth)))
  } else if(facet == "both") {
    gg <- gg + facet_grid(mu_driver_birth ~ s_driver_birth, labeller = label_bquote(rows = mu == .(mu_driver_birth), cols = s == .(s_driver_birth)))
  } else if(facet == "K") {
    gg <- gg + facet_grid(. ~ K, labeller = label_bquote(cols = K == .(K)))
  }
  
  if(!is.null(filename)) {
    png(filename, width = width, height = height, res = res)
    print(gg)
    dev.off()
    return("Saved plot")
  } else {
    return(gg)
  }
}

clonal_turnover_plot <- function(plot_df, correlate = "average growth rate", filename = NULL, width = 800, height = 800, res = 100) {
  gg<-ggplot(data = plot_df, aes(x = mean_autocor, y = Cor_DriverDiversity, 
                            color=as.factor(mu_driver_birth), 
                            fill= as.factor(mu_driver_birth), 
                            shape = as.factor(s_driver_birth)
  )) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + 
    scale_fill_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
    scale_color_manual(name="Driver mutation rate",values=c("#009F71", "#FAC200", "#F35E00"))+
    scale_shape_manual(name="Mean driver fitness effect", values = c(21, 22, 23,  24) )+
    geom_point(size = 4)+
    xlab("mean clonal turnover index")+
    ylab(paste0("correlation coefficient:\n clonal diversity vs ", correlate))+
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
  
  if(!is.null(filename)) {
    png(filename, width = width, height = height, res = res)
    print(gg)
    dev.off()
    return("Saved plot")
  } else {
    return(gg)
  }
}

plot_all_images_test <- function(path, output_filename = NA, file_type = "png", output_dir = NA, trim = -1, include_genotype_plots = TRUE, 
                            cutoff = 0, min_birth_rate = NA, max_birth_rate = NA) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  Muller_df <- muller_df_from_file(paste0(path, "driver_phylo.dat"), cutoff = cutoff)
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
  if(is.na(min_birth_rate)) min_birth_rate <- min(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  if(is.na(max_birth_rate)) max_birth_rate <- max(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  h3 <- Muller_plot(Muller_df, colour_by = "BirthRate", add_legend = TRUE) + 
    scale_fill_distiller(palette = "RdBu", direction = -1, 
                         limits = c(min_birth_rate, max_birth_rate)) #+ 
  #theme(line = element_blank(), rect = element_blank()) + 
  #scale_x_continuous(breaks = c(0, round(max(Muller_df$Generation)))) + 
  #scale_y_continuous(breaks = NULL)
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_driversgrid.dat"), trim)
    image_df[which(image_df$z > 0), "z"] <- as.character(image_df[which(image_df$z > 0), "z"] %% 25 + 1)
    g1 <- grid_plot_test(image_df, palette = dd.col, discrete = TRUE)
  }
  
  g2 <- grid_plot_test(b_grid, add_legend = FALSE, legend_title = "Mean cell\ndivision rate   ") + 
    scale_fill_distiller(name = "Mean cell\ndivision rate   ", palette ="RdBu", 
                         direction = -1, na.value="white", 
                         limits = c(min_birth_rate, max_birth_rate))
  
  image_df <- image_df_from_grid_file(paste0(path, "output_passengersgrid.dat"), trim)
  g3 <- grid_plot_test(image_df, add_legend = TRUE, legend_title = "Mean passenger\nmutations per cell")
  
  if(include_genotype_plots) {
    image_df <- image_df_from_grid_file(paste0(path, "output_popgrid.dat"), trim)
    image_df[which(image_df$z == 0), "z"] <- NA
    g4 <- grid_plot_test(image_df, add_legend = TRUE, legend_title = "Tumour cells\nper gland")
  }
  
  if(!is.na(output_filename)) print(paste0("Created all plots for file ", output_filename), quote = FALSE)
  
  if(include_genotype_plots) {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
      else cairo_pdf(paste0(output_dir,output_filename,".pdf"), width = 7, height = 7, fallback_resolution = 800)
    }
    lay <- rbind(c(1,1,2),
                 c(3,3,3),
                 c(4,4,5),
                 c(NA,6,7))
    print(grid.arrange(h1, g1, h2, h3, g2, g3, g4, layout_matrix = lay, heights = c(1, 1, 1, 0.75)))
  } else {
    if(!is.na(output_filename) & !is.na(output_dir)) {
      if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 180, res = 100)
      else pdf(paste0(output_dir,output_filename,".pdf"), width = 7, height = 1.8)
    }
    lay <- rbind(c(1,2,3))
    print(grid.arrange(h3, g2, g3, layout_matrix = lay))
  }
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}

grid_plot_test <- function(image_df, palette = NA, discrete = FALSE, add_legend = FALSE, legend_title = "") {
  if(length(image_df) == 1) return(ggplot(data.frame()) + theme_classic())
  h2 <- ggplot(image_df, aes(x, y, fill = z)) + 
    geom_tile() +
    theme(legend.position = ifelse(add_legend, "right", "none")) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  if(discrete) {
    if(!is.na(palette[1])) h2 <- h2 + scale_fill_manual(name = legend_title, values = palette, na.value="white") +
        scale_color_manual(values = palette, na.value="white")
  }
  else {
    if(is.na(palette[1])) {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette ="RdBu", direction = -1, na.value="white") + 
        scale_color_distiller(palette ="RdBu", na.value="white")
    }
    else {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette = palette, direction = -1, na.value="white") + 
        scale_color_distiller(palette = palette, na.value="white")
    }
  }
  return(h2)
}

All_AverageGrowthRate_summary_NewCorrelation_K64 <- read_csv("AverageGrowthRate_cor_summary_NewCorrelation_K64.csv", guess_max = 1E4)
All_AverageGrowthRate_summary_NewCorrelation_K64$mean_AverageGrowthRate<-NULL
All_AverageGrowthRate_summary_NewCorrelation_K512 <- read_csv("All_AverageGrowthRate_cor_summary_NewCorrelation_K512.csv", guess_max = 1E4)
All_AverageGrowthRate_summary_NewCorrelation_K4096 <- read_csv("AverageGrowthRate_cor_summary_NewCorrelation_K4096.csv", guess_max = 1E4)
All_AverageGrowthRate_summary_NewCorrelation_K4096$mean_AverageGrowthRate<-NULL
ClonalDiversityVsMeanBirthRate_correlation_K64 <- read_csv( file = "ClonalDiversityVsMeanBirthRate_correlation_K64.csv")
ClonalDiversityVsMeanBirthRate_correlation_K512 <- read_csv( file = "ClonalDiversityVsMeanBirthRate_correlation_K512.csv")
ClonalDiversityVsMeanBirthRate_correlation_K4096 <- read_csv( file = "ClonalDiversityVsMeanBirthRate_correlation_K4096.csv")
AverageGrowthRate_cor_summary_CombinedMutationRateK512 <- read_csv("AverageGrowthRate_cor_summary_CombinedMutationRate_K512_RandomMu.csv", guess_max = 1E4)
AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512 <- read_csv("AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRate_K512_RandomMu.csv", guess_max = 1E4)
AverageGrowthRate_cor_summary_CombinedFitnessEffectK512 <- read_csv("AverageGrowthRate_cor_summary_CombinedFitnessEffect_K512_RandomS.csv", guess_max = 1E4)
df_summary <- read_csv("df_summary.csv", guess_max = 1E4)

df1 <- filter(All_AverageGrowthRate_summary_NewCorrelation_K512, mu_driver_birth == 1e-05, K == 512, s_driver_birth == 0.1)
df2 <- filter(rbind(rbind(All_AverageGrowthRate_summary_NewCorrelation_K64, All_AverageGrowthRate_summary_NewCorrelation_K512), All_AverageGrowthRate_summary_NewCorrelation_K4096), start_size == 250000)
df3a <- filter(AverageGrowthRate_cor_summary_CombinedMutationRateK512, s_driver_birth==0.1)
df3b <- filter(AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512, s_driver_birth==0.1)
dfS1 <- filter(All_AverageGrowthRate_summary_NewCorrelation_K512, mu_driver_birth %in% c(1e-6, 1e-4), K == 512, s_driver_birth %in% c(0.05, 0.2))
dfS3 <- filter(rbind(rbind(ClonalDiversityVsMeanBirthRate_correlation_K64, ClonalDiversityVsMeanBirthRate_correlation_K512), ClonalDiversityVsMeanBirthRate_correlation_K4096), start_size == 250000)
dfS4a <- AverageGrowthRate_cor_summary_CombinedMutationRateK512
dfS4b <- AverageGrowthRate_cor_summary_allFinalSize_fromStartSize125000_CombinedMutationRateK512
dfS5 <- AverageGrowthRate_cor_summary_CombinedFitnessEffectK512

final_size <- 1e6
df_summary<-mutate(df_summary,
                   AverageGrowthRate = (final_size-start_size) / waiting_time, 
                   ExactAverageGrowthRate = (final_size-NumCells) / waiting_time)
cols_list<-colnames(df_summary)[grepl("DriverDiversity",colnames(df_summary))]
wait_cor_summary <-get_wait_cor_summary(df_summary, c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "mu_driver_birth"), num_parameters = 18, min_count = 5, 
                                        VariablesToNotGroupBy = c("mu_driver_birth", "mu_passenger"))
wait_cor_summary_kendall <-get_wait_cor_summary(df_summary, c("DriverDiversity", "DriverEdgeDiversity", "MeanBirthRate", "mu_driver_birth"), num_parameters = 18, min_count = 5, 
                                                VariablesToNotGroupBy = c("mu_driver_birth", "mu_passenger"), method = "kendall")

# Fig 1b:
correlations_plot(df1, filename = "Fig1b.png", width = 800 * 8, height = 800 * 8, res = 800)

# Fig 2a:
clonal_turnover_plot(df2, filename = "Fig2a.png", width = 1500 * 8, height = 1000 * 8, res = 800)

p1 <- "/Users/rnoble/Documents/MontpellierDocuments/Data/ForecastingJune2020/seed_15"
plot_all_images_test(p1, output_filename = "MullerPlots1", file_type = "pdf", 
                output_dir = getwd(), trim = -1, include_genotype_plots = TRUE, cutoff = 0,min_birth_rate = 1, max_birth_rate = 10)
p2 <- "/Users/rnoble/Documents/MontpellierDocuments/Data/ForecastingJune2020/seed_14"
plot_all_images_test(p2, output_filename = "MullerPlots2", file_type = "pdf", 
                output_dir = getwd(), trim = -1, include_genotype_plots = TRUE, cutoff = 0,min_birth_rate = 1, max_birth_rate = 10)

# Fig 3:
correlations_plot(df3a, filename = "Fig3a.png", width = 800 * 8, height = 800 * 8, res = 800)
correlations_plot(df3b, x_axis = "FinalSize", filename = "Fig3b.png", width = 800 * 8, height = 800 * 8, res = 800)
correlations_plot(df3a, group = "biopsysize", depth = 4, filename = "Fig3c.png", width = 800 * 8, height = 800 * 8, res = 800)
correlations_plot(df3a, group = "biopsydepth", depth = 4, filename = "Fig3d.png", width = 750 * 8, height = 750 * 8, res = 800)

# Fig 4:
correlations_plot(filter(wait_cor_summary, start_size == -2), group = "predictors", x_axis = "s_driver_birth", depth = 1, facet = "K", filename = "Fig4a.png", 
                  width = 1200 * 8, height = 800 * 8, res = 800, SignificanceLevel = 0.05, ylab = "\nprogression-free survival")
correlations_plot(filter(wait_cor_summary_kendall, start_size == -2), group = "predictors", x_axis = "s_driver_birth", depth = 1, facet = "K", filename = "Fig4a_kendall.png", 
                  width = 1200 * 8, height = 800 * 8, res = 800, SignificanceLevel = 0.05, ylab = "\nprogression-free survival")

# Fig S1:
correlations_plot(dfS1, filename = "FigS1.png", width = 1000 * 8, height = 1000 * 8, facet = "both", res = 800)

# Fig S3:
clonal_turnover_plot(dfS3, correlate = "mean cell division rate", filename = "FigS3.png", width = 1500 * 8, height = 1000 * 8, res = 800)

# Fig S4:
correlations_plot(dfS4a, filename = "FigS4a.png", width = 1500 * 8, height = 1000 * 8, res = 800, facet = "s")
correlations_plot(dfS4b, x_axis = "FinalSize", filename = "FigS4b.png", width = 1500 * 8, height = 1000 * 8, res = 800, facet = "s")
correlations_plot(dfS4a, group = "biopsysize", filename = "FigS4c.png", width = 1500 * 8, height = 1000 * 8, res = 800, facet = "s")
correlations_plot(dfS4a, group = "biopsydepth", depth = 4, filename = "FigS4d.png", width = 1500 * 8, height = 1000 * 8, res = 800, facet = "s")

# Fig S5:
correlations_plot(dfS5, filename = "FigS5.png", width = 1500 * 3/4 * 8, height = 1000 * 8, res = 800, facet = "mu")








