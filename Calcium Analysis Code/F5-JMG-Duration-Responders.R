# RScript4
# Use min-max normalized responder matrix from RScript3
# Use file.choose() to select your directory
input_directory <- dirname(file.choose(new = TRUE))
# Set the selected directory as the working directory
setwd(input_directory)

# Get a list of file names in the directory
file_list <- list.files(input_directory, full.names = TRUE)

# Create an empty list to store the resulting vector tables generated in the for loop
results <- list()

# Loop through each file and process it
for (file in file_list) {
  # Read in min-max normalized CSV file and process file as needed
  csv_file <- read.csv(file)
  responders <- csv_file[, 2:dim(csv_file)[2]]  
  
  # Find the peak values of cells
  peak_values <- apply(responders, 2, max)
  
  # Get baseline parameters from resp matrix
  baseline_frs <- responders[1:30, ]
  baseline_avg <- apply(baseline_frs, 2, mean)
  baseline_std <- apply(baseline_frs, 2, sd)
  
  # Establish threshold for each cell (i.e. Fo of cell + 3 x Fo's standard deviation); assign to "threshold"
  threshold <- baseline_avg + 3 * baseline_std
  
  # Define the condition to be met by cell to be considered a true peak
  true_peak <- !is.na(peak_values) & !is.na(threshold) & peak_values > threshold
  
  true_peak_matrix <- responders[, true_peak]
  
  half_peak_durations <- numeric(ncol(true_peak_matrix))  # Create empty vector to store durations
  
  for (i in 1:ncol(true_peak_matrix)) {
    # Find peak value from true peak matrix
    true_peak_max <- max(true_peak_matrix[, i])
    
    # Calculate half of the peak value
    half_peak_value <- true_peak_max / 2
    
    # Find the indices where the values cross the half-peak value
    crossing_indices <- which(diff(sign(true_peak_matrix[, i] - half_peak_value)) != 0)
    
    # Calculate the duration between the first and last crossing point
    if (length(crossing_indices) >= 2) {
      half_peak_duration <- crossing_indices[length(crossing_indices)] - crossing_indices[1]
    } else {
      half_peak_duration <- NA
    }
    
    # Store the half-peak duration for the current column
    half_peak_durations[i] <- half_peak_duration
  }
  
  #Write half-peak durations to a CSV file
  write.csv(half_peak_durations, file = gsub('.csv', '_cell_resp-duration.csv', file), row.names = TRUE)
  
  # Calculate average cell response duration for the current file
  avg_response_duration <- mean(half_peak_durations, na.rm = TRUE)
  
  # Store the result in the results list
  results[[file]] <- avg_response_duration
}

# Combine all the results into one table 
response_duration <- do.call(rbind, results)

# Save result as CSV
write.csv(response_duration, file = "Responders-Duration.csv")

