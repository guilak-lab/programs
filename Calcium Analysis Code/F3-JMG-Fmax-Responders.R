# RScript2
# Use raw responder matrix from RScript1
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
  # Read in raw responder matrix and process file as needed
  csv_file <- read.csv(file)
  responders <- csv_file[, 2:dim(csv_file)[2]]
  
  # Find the Fmax of each cell in the responders maxtrix 
  Fmax_cell <- apply(responders, 2, max)

  #Calculate average Fmax for responders
  Fmax_responders <-mean(Fmax_cell, na.rm = TRUE) 
  
  # Store the result in the results list
  results[[file]] <- Fmax_responders
}
# Combine all the results into one table 
response_duration <- do.call(rbind, results)

# Save result as CSV
write.csv(response_duration, file = "Responders-Fmax.csv")

