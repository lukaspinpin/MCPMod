rm(list = ls())
usethis::edit_r_profile()
library(detectseparation)
library(clustermq)
library(DoseFinding)
library(ggplot2)
library(logistf)
library(dplyr)
library(broom)
library(tidyr)

## Set LSF as scheduler and path to template
options(clustermq.scheduler = "lsf")

source("MCP_Rand.R") ## load auxiliary code
export <- list(logit=logit, inv_logit=inv_logit, PBR = PBR, generate_data=generate_data, 
               reg_mod=reg_mod, MCP_RD = MCP_RD, get_randomised_data=get_randomised_data, 
               MCP_Rand=MCP_Rand, MCTtest=MCTtest, detect_separation=detect_separation)

## setup scenarios (toy example)
doses <- c(0, 10, 25, 100)
covariate <- c(TRUE)
reg <- c("glm", "logistf")
test_type <- c("population", "randomisation_mean", "randomisation_residuals")

# Create the full factorial grid first
params_changed <- expand.grid(test_type=test_type, reg=reg, covariate=covariate)

# Filter for MCP benchmark configurations only
filtered_params <- params_changed[!(params_changed$test_type == "population"  & params_changed$reg !=  "glm"),]
filtered_params <- filtered_params[!(filtered_params$test_type == "population"  & filtered_params$covariate ==  FALSE),]

# Convert remaining combinations to a data frame
filtered_params <- as.data.frame(filtered_params)

# Define the custom function
proportion_less_than_005 <- function(x) {
  # Remove NAs from the input vector
  x <- na.omit(x)
  if (length(x) == 0) return(NA) # Handle case where no data points left after omitting NAs
  sum(x < 0.05) / length(x)
}

# Define the custom function
proportion_less_than_01 <- function(x) {
  # Remove NAs from the input vector
  x <- na.omit(x)
  if (length(x) == 0) return(NA) # Handle case where no data points left after omitting NAs
  sum(x < 0.1) / length(x)
}

MCR_Rand_dat <- function(seed=1, randMethod= "PBR", doses=doses, n=n,  B=7, p1=0.2, pK=0.2, nrand=nrand, timetrend = FALSE, filtered_params=filtered_params){
  n_methods <- dim(filtered_params)[1]
  p_vals <- rep(0, n_methods)
  run_time <- rep(0, n_methods)  # Initialize run_time to store the runtime of each method
  no_respon = FALSE; fail_MCP = FALSE; completesep = FALSE
  mods <- DoseFinding::Mods(emax = c(10, 50), 
                            sigEmax = rbind(c(5, 3), c(25, 3)),
                            #betaMod = c(0.2475, 2.025),
                            betaMod = c(0.1529, 0.5809),
                            doses = doses)
  plotMods(mods, trafo = inv_logit)
  x1 <- rnorm(sum(n), 0, 1)
  dat <- generate_data(seed=seed, doses=doses, mods=mods, randMethod=randMethod, B=B, n=n, x1=x1, p1=p1, pK=pK, timetrend = timetrend )
  result <- glm(y ~ 1, data = dat, family = binomial, method = "detect_separation")
  completesep <- result$outcome
  
  if(sum(dat$y)==0){
   p_vals = rep(NA, n_methods)
   run_time = p_vals
  } else {
    for(i in 1:n_methods){
      if(sum(dat[dat$dose == 0,]$y)==0){completesep=TRUE}
      start_time <- Sys.time()  # Capture start time
      
      p_vals[i] <- MCP_Rand(dat=dat, doses=doses,  n=n, 
                            test_type=filtered_params$test_type[i], reg=filtered_params$reg[i], 
                            randMethod="PBR", B=B, x1=x1, 
                            covariate=filtered_params$covariate[i], nrand=nrand)
      
      end_time <- Sys.time()  # Capture end time
      
      run_time[i] <- end_time - start_time  # Calculate runtime
      if(filtered_params$test_type[i]=="population" & is.na(p_vals[i]) ){fail_MCP = TRUE}
    }
  }
  
  list(p_vals = p_vals, run_time = run_time, no_respon = no_respon,
       fail_MCP = fail_MCP, completesep = completesep)  # Return both p-values and runtimes
}

simPlots <- function(results, paramscom, sim_matrix_p_vals=0){
  # Convert covariate to factor for better visualization
  results$reg <- factor(results$reg, labels = c("glm", "logistf"))
  
  # Add an additional column to combine "test_type" and "reg" for purely visual purposes
  results$combined <- interaction(results$test_type, results$reg)
  # Defining the colors
  colors <- c(rgb(215,25,28, maxColorValue = 255),
              rgb(253, 174, 97, maxColorValue = 255),
              rgb(255, 255, 51, maxColorValue = 255),
              rgb(171, 221, 164, maxColorValue = 255),
              rgb(43, 131, 186, maxColorValue = 255)) 
  # Define the legend labels
  methods <- c(
    "Population-based Test",
    "Randomisation Test (MLE)",
    "Randomisation Test (MLE; residual-based)",
    "Randomisation Test (penalized MLE)",
    "Randomisation Test (pen. MLE; residual-based)"
  )
  legend_labels <- methods
  
  
  # Create the bar plot using ggplot
  t1 <- ggplot(results, aes(x = combined, y = TypeIerror_Power*100, fill = combined)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) + # Add borders to bars
    scale_fill_manual(values = colors, labels = legend_labels) +
    labs(title = "Type I Error Rate for Different Tests",
         x = "",
         y = "Type I Error Rate in %",
         fill = "Legend") + 
    ylim(0,12) +
    theme_minimal(base_size = 15) + # Applying a minimal theme with larger base font size
    theme(
      axis.text.x = element_blank(),           # Remove the text under each bar
      axis.ticks.x = element_blank(),          # Remove the ticks under each bar
      legend.title = element_blank(),          # Remove legend title
      legend.text.align = 0,                   # Align legend text
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # Center the title, increase size, and make it bold
      legend.position = "bottom",               # Position the legend to the right
      legend.spacing.y = unit(0.5, 'cm'),      # Increase space between legend elements
      legend.text = element_text(size = 12)# Increase the legend text size
    ) +
    geom_text(aes(label = round(TypeIerror_Power*100, 3)), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 4)
  # Save the plot to a PDF file
  ggsave(paste0("t1_", paramscom, ".pdf"), plot = t1, width = 14, height = 7)
  
  
  # Create the bar plot using ggplot
  ptt <- ggplot(results, aes(x = combined, y = RunTime, fill = combined)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = colors, labels = legend_labels) +
    labs(title = "RunTime for Different Tests",
         x = "",
         y = "RunTime (Seconds)",
         fill = "Test Type") +
    theme_minimal(base_size = 15) + # Applying a minimal theme with larger base font size
    theme(
      axis.text.x = element_blank(),           # Remove the text under each bar
      axis.ticks.x = element_blank(),          # Remove the ticks under each bar
      legend.title = element_blank(),          # Remove legend title
      legend.text.align = 0,                   # Align legend text
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # Center the title, increase size, and make it bold
      legend.position = "bottom",               # Position the legend to the right
      legend.spacing.y = unit(0.5, 'cm'),      # Increase space between legend elements
      legend.text = element_text(size = 12)    # Increase the legend text size
    ) +
    geom_text(aes(label = round(RunTime, 3)), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 4)
  
  # Save the plot to a PDF file
  ggsave(paste0("ptt_", paramscom, ".pdf"), plot = ptt, width = 11, height = 5)
  
  # Convert matrix to long format
  # df <- sim_matrix_p_vals[-1,]
  df <- sim_matrix_p_vals
  df <- df %>%
    # mutate(method = methods[-1]) %>%
    mutate(method = methods) %>%
    pivot_longer(cols = -method, names_to = "variable", values_to = "p_value")
  # Factorize the method variable with specified levels' order
  df <- df %>% 
    mutate(method = factor(method, levels = methods))
  # mutate(method = factor(method, levels = methods[-1]))
  
  # Exclude NAs and create the plot
  plot_pval <- ggplot(df %>% filter(!is.na(p_value)), aes(x = p_value, fill = method)) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity", boundary = 0, breaks = seq(0, 1, length.out = 51)) +
    #scale_fill_manual(values = colors[-1]) +
    scale_fill_manual(values = colors) +
    labs(title = "Histogram p-values by test", x = "P-value", y = "Frequency", fill = "Test") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    facet_wrap(~ method, scales = "free_y", ncol = 2)+ guides(fill="none")
  
  
  ggsave(paste0("plot_pval_", paramscom, ".pdf"), plot = plot_pval, width = 11, height = 8)
}

# Define the vectors
n <- c(7, 14, 14, 14)
time_trend <- c(TRUE, FALSE)
#time_trend <- c(FALSE)
rand_method <- c("CR", "RA", "PBR")
#rand_method <- c( "PBR")
sample_size <- list(n, 2*n, 10*n)   
#sample_size <- list(n)  
pK <- seq(0.2, 0.9, 0.1)
#pK <- c(0.2, 0.8)

#grid <- expand.grid(
 #  pK = pK,
  #  rand_method = rand_method,
   # sample_size = I(sample_size),
   # time_trend = time_trend,
   # KEEP.OUT.ATTRS = FALSE,
   # stringsAsFactors = FALSE)

# Generate the combination grid
grid <- rbind(
      expand.grid(
      pK = pK,
      rand_method = rand_method[2:3],
      sample_size = I(sample_size[1]),
      time_trend = time_trend,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ), 
    expand.grid(
      pK = c(0.61, seq(0.2, 0.7, 0.1)),
      rand_method = rand_method[2:3],
      sample_size = I(sample_size[2]),
      time_trend = time_trend,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    expand.grid(
      pK = c(0.364,seq(0.2, 0.45, 0.05)),
      rand_method = rand_method,
      sample_size = I(sample_size[3]),
      time_trend = time_trend,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
)


# Loop over all combinations
for (i in 1:nrow(grid)) {
  current_combination <- grid[i,]
  # Extracting current values (you can use these values in your tasks)
  current_time_trend <- current_combination$time_trend
  current_rand_method <- current_combination$rand_method
  current_sample_size <- current_combination$sample_size[[1]]  # note the extraction of list element
  current_pK <- current_combination$pK
  
  print(current_combination)
  
  # Create the full factorial grid first
  param_constatnt <- list(nrand=1000, doses=doses, n=current_combination$sample_size[[1]],   
                          B=7,  timetrend = current_combination$time_trend, pK= current_combination$pK,
                          randMethod= current_combination$rand_method, filtered_params=filtered_params)
  
      # Number of Simulations
      nsim = 10000
      seeds <- data.frame(seed=seq(1,nsim))
      
      sim_out <- Q_rows(fun = MCR_Rand_dat, 
                        df = seeds,
                        const = param_constatnt,
                        seed = 23456,
                        n_jobs = nrow(seeds),
                        template = list(
                          walltime = 300,
                          job_name = "MCP_Rand",
                          log_file = "logs.txt"),
                        pkgs = c("data.table", "DoseFinding", "logistf"),
                        export = export)
      
      paramscom <- paste("nsim", nsim, "nrand", param_constatnt$nrand, "pK", round(current_combination$pK*100),
                         "rand_method", current_combination$rand_method, "N", sum(current_combination$sample_size[[1]]), 
                         'TimeTrend', current_time_trend, sep="_")
      filename <- paste("simulation", paramscom, "FULL", sep="_")
      saveRDS(sim_out, file= filename)
      
      percent_completesep = sum(do.call(rbind, lapply(sim_out, function(x) x$completesep)))/nsim
      percent_fail_MCP = sum(do.call(rbind, lapply(sim_out, function(x) x$fail_MCP)))/nsim
      percent_no_respon = sum(do.call(rbind, lapply(sim_out, function(x) x$no_respon)))/nsim
      print(c(percent_completesep, percent_fail_MCP, percent_no_respon))
      
      p_vals_list <- lapply(sim_out, function(x) x$p_vals)
      run_time_list <- lapply(sim_out, function(x) x$run_time)
      
      # Convert the list into a matrix
      sim_matrix_p_vals <- as.data.frame(t(do.call(rbind, p_vals_list)))
      run_times <- as.data.frame(t(do.call(rbind, run_time_list)))
      
      row.names(filtered_params) <- NULL
      # Convert the matrix to a data frame
      results <- cbind(filtered_params, RunTime = rowMeans(run_times), TypeIerror_Power=apply(sim_matrix_p_vals, 1, proportion_less_than_01))
      print(results)
      
      filename <- paste("simulation", paramscom, "RESULTS", sep="_")
      saveRDS(results, file= filename)
  
      #plot
      simPlots(results, paramscom, sim_matrix_p_vals=sim_matrix_p_vals)
}

###### Power Plot Figure 2 in Paper ######
## Other Plots can be created by changing imported file names 

pK_values <- as.integer(seq(0.2, 0.9, 0.1)*100)
#pK_values <- sort(as.integer(c(0.61, seq(0.2, 0.7, 0.1))*100))
#pK_values <- sort(as.integer(c(0.364,seq(0.2, 0.45, 0.05))*100))

# Initialize a list to store the last columns of the files
last_columns <- list()

# Loop through each pK value and read the corresponding RDS file
for (pK in pK_values) {
  # Construct the filename
  filename <- sprintf("/home/pinlu1/Documents/simulation_nsim_10000_nrand_1000_pK_%d_rand_method_PBR_N_49_TimeTrend_FALSE_RESULTS", pK)
  
  # Print the filename (optional, for debugging)
  cat("Reading file:", filename, "\n")
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the last column and save it in the list
  last_columns[[as.character(pK)]] <- data[[ncol(data)]]
}

# Optionally, combine all the last columns into a single data frame
# You can name the columns using the corresponding pK values
combined_last_columns <- do.call(cbind, last_columns)
colnames(combined_last_columns) <- names(last_columns)

# Specify the output file and open a PNG graphics device
pdf(file = "power_nsim_10000_nrand_1000_rand_method_PBR_N_49_TimeTrend_FALSE_RESULTS.pdf", width = 10, height = 7)

# Plot the data
##colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")
colors <- c(rgb(215,25,28, maxColorValue = 255),
            rgb(253, 174, 97, maxColorValue = 255),
            rgb(255, 255, 51, maxColorValue = 255),
            rgb(171, 221, 164, maxColorValue = 255),
            rgb(43, 131, 186, maxColorValue = 255))
methods <- c(
  "Population-based Test",
  "Randomisation Test (MLE)",
  "Randomisation Test (MLE; residual-based)",
  "Randomisation Test (penalized MLE)",
  "Randomisation Test (pen. MLE; residual-based)"
)
lty <- c(3,1,1,2,2)
matplot(pK_values, t(combined_last_columns), type = 'l', lty = lty, col = colors, lwd = 3,
        xlab = "Success Rate on Arm K-1 (%)",
        ylab = "Power",
        main = "Power Comparison of Different Tests",
        cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.4)
abline(h=seq(0,1,by=0.05), lty=3, col="lightgrey")
abline(v=seq(20,90,by=10), lty=3, col="lightgrey")
matlines(pK_values, t(combined_last_columns), type = 'l', lty = lty, col = colors, lwd = 2,)

# Add dots to indicate actual data points
for (i in 1:nrow(combined_last_columns)) {
  points(pK_values, combined_last_columns[i, ], col = colors[i], pch = 16, cex = 1)
}

abline(h = 0.8, col = "black", lty = 2)

# Add legend in the upper-left corner with more space between elements
legend("bottomright", legend = methods, col = colors, lty = lty, cex = 0.8, bty = "n", inset = 0.04, y.intersp = 1.5, lwd=2)

dev.off()

print(round(combined_last_columns*100,2))


#### Print percentage of complete separation and faild MCP-Mod
# for original data scenario but other can analysed by changing file name
filename <- sprintf("/home/pinlu1/Documents/simulation_nsim_10000_nrand_1000_pK_%d_rand_method_PBR_N_49_TimeTrend_FALSE_FULL", 20)

# Print the filename (optional, for debugging)
cat("Reading file:", filename, "\n")

# Read the RDS file
data <- readRDS(filename)

percent_completesep = sum(do.call(rbind, lapply(data, function(x) x$completesep)))/10000
percent_fail_MCP = sum(do.call(rbind, lapply(data, function(x) x$fail_MCP)))/10000
c(percent_completesep, percent_fail_MCP)
