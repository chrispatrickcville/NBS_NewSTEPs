### Use to consolidate NewSTEPs ancillary reports. ###

# Make sure to enter environmental variables in NewSTEPs_variable_setting
# before running code below

###########

# Read all files indicating difference in analytes reported out from
# those expected to be called out

diff_data <- read_data(folder=output_path, patt="Difference_in_results_called", 
                       separator=",", type="other")

# Remove "X" column
diff_data$X <- NULL

# Remove any duplicated records
diff_data <- unique(diff_data)

# Fix column names
colnames(diff_data) <- gsub("\\.", " ", colnames(diff_data))
colnames(diff_data)[colnames(diff_data) == "Date of Call s "] <- "Date of Call(s)"

# Write csv
write.csv(diff_data, row.names=FALSE, paste0(output_path, "Aggregated_Difference_in_results_called.csv"))

###########

# Read all files with samples that were not indicated as reported out in crit data

missing_samp_crit <- read_data(folder=output_path, patt="Samples with critical results in sample data", 
                               separator=",", type="other")

# Remove "X" column
missing_samp_crit$X <- NULL

# Remove any duplicated records
missing_samp_crit <- unique(missing_samp_crit)

# Fix column names
colnames(missing_samp_crit)[colnames(missing_samp_crit) == "combined"] <- "Critical/Critical Group Analytes"

# Write csv
write.csv(missing_samp_crit, 
          row.names=FALSE, paste0(output_path, 
                                  "Aggregated_Samples with criticals not reported as called out.csv"))

###########

# Read all files with samples that were reported out as critical but not indicated as such in sample data

missing_samp_QI <- read_data(folder=output_path, patt="Samples with critical results in critical result", 
                             separator=",", type="other")

# Remove "X" column
missing_samp_QI$X <- NULL

# Remove any duplicated records
missing_samp_QI <- unique(missing_samp_QI)

# Write csv
write.csv(missing_samp_QI, 
          row.names=FALSE, paste0(output_path, 
                                  "Aggregated_Samples with results called out but not indicated as critical in data.csv"))