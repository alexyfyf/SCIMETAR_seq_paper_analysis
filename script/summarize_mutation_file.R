# Function to summarize a single TSV file
summarize_mutation_file <- function(file_path) {
  # cat("Processing file:", file_path, "\n")
  
  # Read the TSV file
  data <- read.table(file_path, 
                     sep = "\t", 
                     header = FALSE, 
                     stringsAsFactors = FALSE,
                     col.names = c("read_id", "read_type", "mutation"))
  
  # Group by read_id and mutation to see read type combinations per read
  read_combinations <- data %>%
    group_by(read_id, mutation) %>%
    summarise(read_types = paste(sort(unique(read_type)), collapse = "_"), 
              .groups = "drop")
  
  # Count how many reads have each combination of read types for each mutation
  summary <- read_combinations %>%
    group_by(mutation, read_types) %>%
    summarise(read_count = n(), .groups = "drop") %>%
    arrange(mutation, desc(read_count))
  
  return(summary)
}