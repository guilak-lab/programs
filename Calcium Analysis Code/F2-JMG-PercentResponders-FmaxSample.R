# RScript1
# Use CSV file from Fiji
# Note: (rows,columns) or (1,2), row "1" = frame, column "2" = cell detected by ROI mask; values in matrix correspond to Mean Gray Values (MGV = F i.e. fluorescence) from Fiji
# Use file.choose() to select your directory
input_directory <- dirname(file.choose(new=TRUE))
# Set the selected directory as the working directory
setwd(input_directory)

# Get a list of file names in the directory
file_list <- list.files(input_directory, full.names = TRUE)

# Create an empty list to store result vector tables generated in for loop
results = list()

# Loop through each file and process it
for (file in file_list) {
  # Read in single cvs file and process file as needed (i.e. exclude first column since it's frame numbers and not MVG values)
  csv_file <- read.csv(file)
  F_data <-  csv_file[ , 2:dim(csv_file)[2]]
  
  # Calculate the Fo (average F of the baseline frames) and corresponding standard deviation for each cell (i.e. column)
  # Extract the first 30 rows of each column in "F_data"; assign new matrix to "baseline_frs" 
  baseline_frs <-  F_data[1:30, ]
  # Average "baseline_frs" to get Fo; assign new vector to "Fo"
  Fo <-  apply(baseline_frs, 2, mean)
  # Calculate the standard deviation of Fo; assign to "Fo_std"
  Fo_std <-  apply(baseline_frs, 2, sd)
  
  # Generate normalizing matrix based on total frames (rows) and the "Fo" vector (columns) 
  # Replicate "Fo" vector for each row in "F_data" (e.g. 140 frames) to generate a baseline matrix; assign new matrix to "Fo_matrix"
  Fo_matrix <-  matrix(rep(Fo, dim(F_data)[1]), nrow = dim(F_data)[1], byrow = TRUE)
  
  # Calculate the relative fluorescence (note: F-Fo = dF, relative fluorescence = dF/Fo ); assign matrix to "rF"
  rF <-  (F_data - Fo_matrix)/Fo_matrix

# ### Save rF matrix (normalized matrix) and Fo matrix (baseline matrix) as CSV files
  write.csv(rF, file = gsub('.csv', '_normalized.csv', file))
  write.csv(Fo_matrix, file = gsub('.csv', '_baseline.csv', file))
  
  # Find the maximum F and rF for each cell (i.e. max response from each cell); assign to rFmax_cell"
  rFmax_cell <-  apply(rF, 2, max)
  
  # Establish threshold for each cell (i.e. Fo of cell + 3 x  Fo's standard deviation); assign to "threshold"
  threshold <- Fo + 3 * Fo_std
  
  # Define the condition to be met by cell to be considered a responder (i.e. rFmax value for a cell must be greater than its respective threshold); note "!is.na" excludes NA values
  responder <- !is.na(rFmax_cell) & !is.na(threshold) & rFmax_cell > threshold
  
  #Create new matrices with cells that meet the responder condition (responders) or not (nonresponders)
  responder_matrix_norm <- subset(rF, select = responder)
  responder_matrix_raw <- subset(F_data, select = responder)
  nonresponder_matrix_norm <- subset(rF, select = !responder)
  nonresponder_matrix_raw <-subset(F_data, select = !responder)
  
# # ### Save responder and nonresponder matrices as CSV files
#   write.csv(responder_matrix_norm, file = gsub('.csv', '_responders-norm.csv', file))
#   write.csv(responder_matrix_raw, file = gsub('.csv', '_responders-raw.csv', file))
#   write.csv(nonresponder_matrix_norm, file = gsub('.csv', '_nonresponders-norm.csv', file))
#   write.csv(nonresponder_matrix_raw, file = gsub('.csv', '_nonresponders-raw.csv', file))
  
  #Find the Fmax of each cell
  Fmax_cell <- apply(F_data, 2, max)
  
  # Generate a data frame containing: cell's rFmax, threshold, Fo, Fo's standard deviation, and responder status 
  response_summary <- data.frame(Fmax_cell = Fmax_cell, rFmax_cell = rFmax_cell, threshold = threshold, baseline_avg = Fo, baseline_std = Fo_std, responder = responder)
  
  # ### Save summary of response as CSV file
  write.csv(response_summary, file = gsub('.csv', '_response_table.csv', file))
  
  # Calculate the the maximum F and percent responders for the sample
  # Calculate the average Fmax for all cells (i.e. overall response of field of cells or sample) using rFmax_cell vector; assign to "rFmax_sample"
  vector <- c(Fmax_cell)
  Fmax_sample <- mean(vector[is.finite(vector)],na.rm = TRUE)
  
  # Count the number of responders (i.e. TRUE values) in "responders" logical variable; assign to "number_of_responders"
  number_of_responders <- sum(responder, na.rm = TRUE)
  
  # Count the number of total cells in F_data 
  total_cells <-  ncol(F_data)
  
  # Calculate the percent of responding cells
  percent_responders <- (number_of_responders /total_cells) * 100
  
  # Generate data frame for sample containing: total cells, number of responders, percent responders, and the relative maximum response (rFmax_sample) of the sample
  result <- data.frame(total_cells =  total_cells, responding_cells = number_of_responders, percent_of_responders = percent_responders,Fmax_sample = Fmax_sample)  
  # Combine the new data frame with the existing results in results table
  results [[file]] <- result
}

# Use file.choose() to select your directory
output_directory <- dirname(file.choose(new=TRUE))
# Set the selected directory as the working directory
setwd(output_directory)
#Combines all the results tables into one table 
results_table <-do.call(rbind, results)
### Save table as a CSV file
write.csv(results_table, file = "Results-Table-3std.csv")


