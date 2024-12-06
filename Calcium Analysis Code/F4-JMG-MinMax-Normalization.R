# RScript3
# Use raw responder matrix from RScript1 
# !!!!!!MAKE SURE YOU HAVE MORE THAN 2 COLUMNS!!!!!!
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

  # Function to normalize matrix by minimal and maximal values for columns
  normalize_minmax <- function(responders) {
    num_cols <- ncol(responders)    # Number of columns in the matrix
  
    # Check if the matrix has any columns
    if (num_cols > 0) {
      # Loop through each column and normalize the values
      for (col in 1:num_cols) {
        col_min <- min(responders[, col])    # Find the minimum value in the column
        col_max <- max(responders[, col])    # Find the maximum value in the column
        responders[, col] <- (responders[, col] - col_min) / (col_max - col_min)    # Normalize the column
       }
     }
    return(responders) # Return the normalized matrix
   }

  # Call the function to normalize the columns of the matrix
  normalized_responders <- normalize_minmax(responders)
  
  # Save the normalized matrix
  write.csv(normalized_responders, file = gsub('.csv', '_minmax-norm.csv', file))
}