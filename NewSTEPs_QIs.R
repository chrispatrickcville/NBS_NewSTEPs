####### User Does Not Change This Code ######

### Run after setting variables in "NewSTEPs_variable_setting.R" ###

# Load packages and functions
source(paste0(ns_wd, "NewSTEPs_load_packages_and_functions.R"))

# Remove crit_data if it exists
if (exists("crit_data")) {
  rm(crit_data)
}

# Convert start and end date for use in file titles
start_date_file <- gsub("/", "-", start_date)
end_date_file <- gsub("/", "-", end_date)

# Read import file
Monthly_Template <- suppressWarnings(read.csv(Monthly_Template_path))

# Create data frame to store the sample IDs for any "Unknown" samples for any measure
unknowns <- data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("Source", "Measure", "Count", "SampleIDs"))), 
                       stringsAsFactors = FALSE)

# Check for correct columns
col_check(sample_data_path, "sample")

# Read sample data (if does not already exist) 
if (!exists("QI_data")) {
  QI_data <- read_data(folder=sample_data_path, patt=NULL, separator="|", type="sample", 
                       "BIRTHDATE", "COLLECTIONDATE", "RECEIVEDATE", "RELEASEDATE")
}

# Check that range of RECEIVEDATE dates overlaps the requested start and end date
date_comp_check(QI_data, "sample")

# Filter samples by desired dates for RECEIVEDATE and RECEIVED.DATE
QI_filter <- subset(QI_data, QI_data$RECEIVEDATE >= start_date & QI_data$RECEIVEDATE <= end_date)

##### Out of range data checking and preparation #####

# Call log data should be checked for all periods after and including June 2016
if (start_date >= as.Date("2016-06-01")) {
  
  # Check for correct columns 
  col_check(call_log_path, "call_log")
  
  # Read data, reformat Date.of.Call, and join all analytes
  # from call into single cell
  crit_data <- read_data(folder=call_log_path, patt=NULL, separator=",", type="call_log", 
                         "Date.of.Call")
  
}

# Approved log data should be checked for all periods before and including June 2016
if (start_date <= as.Date("2016-06-01")) {
  
  # Check for correct columns 
  col_check(approved_path, "approve_log")
  
  # Read data, reformat Approved.Date, rename columns, filter for critical results,
  # and restructure so grouped by EXTERNAL_ID and Date.of.Call with all analytes
  # in a single cell
  approved_data <- read_data(folder=approved_path, patt=NULL, separator=",", type="approve_log", 
                             "Approved.Date")
  
  # If crit_data exists, bind approved_data to it; otherwise, rename approved_data
  # as crit_data
  if (exists("crit_data")) {
    crit_data <- rbind(crit_data, approved_data)
  } else {
    crit_data <- approved_data
  }
  
}

# Get receipt date for out-of-range result data by joining on QI_data
crit_check <- inner_join(crit_data, QI_data[, c("SAMPLEID", "RECEIVEDATE")],
                         by = c("EXTERNAL_ID" = "SAMPLEID"))

# Check that range of RECEIVEDATE dates overlaps the requested start and end date
date_comp_check(crit_check, "critical result reporting")

# Get birthdate and receipt date for out-of-range result data. Inner
# join on filtered sample data will mean the out-of-range data will
# be filtered for the same range as the sample data
crit_data_filter <- inner_join(crit_data, QI_filter[, c("SAMPLEID", "BIRTHDATE", "RECEIVEDATE", "RECALL_FLAG")],
                               by = c("EXTERNAL_ID" = "SAMPLEID"))

# Remove records with RECALL_FLAG = Y (this may change later)
crit_data_filter <- crit_data_filter[crit_data_filter$RECALL_FLAG == "N", ]

# Reorganize columns and drop RECALL_FLAG column
crit_data_filter <- crit_data_filter[, c("EXTERNAL_ID", "BIRTHDATE", "RECEIVEDATE", 
                                         "Date.of.Call", "combined")]

# Find time differences in birth to report and receipt to report
crit_data_filter$diff_days_birth2report <- crit_data_filter$Date.of.Call - crit_data_filter$BIRTHDATE
crit_data_filter$diff_days_receipt2report <- crit_data_filter$Date.of.Call - crit_data_filter$RECEIVEDATE

# Report out any birth2report or receipt2report values that are negative 
bad_samples1 <- paste0(crit_data_filter$EXTERNAL_ID[!is.na(crit_data_filter$diff_days_birth2report) &
                                                      (crit_data_filter$diff_days_birth2report < 0 | 
                                                         crit_data_filter$diff_days_receipt2report < 0)], collapse = ", ")
if (bad_samples1 != "") {
  message <- paste0("\nAt least one report date for the critical result reporting data is incorrect\n(captured as being either before the birth date or before the receipt date,\nor both). See the following SAMPLEIDs/EXTERNAL_ID(s) in the critical result\nreporting data:\n\n",
                    bad_samples1, "\n\nPlease address this before running the report.")
  cat(message)
}

# Report out any birth2report or receipt2report values that are after the current date (also may want to add in some kind of max)
bad_samples2 <- paste0(crit_data_filter$EXTERNAL_ID[!is.na(crit_data_filter$Date.of.Call) &
                                                      crit_data_filter$Date.of.Call > Sys.Date()], collapse = ", ")
if (bad_samples2 != "") {
  message <- paste0("\nAt least one report date for the critical result reporting data is incorrect\n(captured as being after today's date). See the following SAMPLEIDs/EXTERNAL_ID(s)\nin the critical result reporting data:\n\n",
                    bad_samples2, "\n\nPlease address this before running the report.")
  cat(message)
}

# Stop reporting if either bad_samples1 or bad_samples2 are not equal to ""
if (bad_samples1 != "" | bad_samples2 != "") {
  stopQuietly()
}

#####################################################################################################
## Check to ensure that out-of-range samples in call log match out-of-range samples in sample data ##
#####################################################################################################

# Aggregate results in critical data by SAMPLEID
crit_filt_check <- crit_data_filter[!is.na(crit_data_filter$Date.of.Call), 
                                    c("EXTERNAL_ID", "combined")]
crit_filt_check <- aggregate(combined ~ EXTERNAL_ID, data = crit_filt_check, toString)
crit_filt_check$combined <- gsub(" ", "", crit_filt_check$combined)

# Subset QI_filter data for samples with Abnormal or Critical in any analyte columns
QI_filter_crits <- QI_filter[apply(QI_filter[, analytes_check], 1, 
                                   function(x) any(grepl("Critical|Abnormal", x))) & 
                               QI_filter$RECALL_FLAG != "Y", ]

# Get a list of all critical results for each sample
QI_filter_crits$combined <- as.character(apply(QI_filter_crits, 1, getCols, colValue="combined", 
                                               colsToCheck=analytes_check, strCheck="ritical"))

# Create analyte sets for disorders in which multiple abnormals result in a critical.
# NOTE: These are the ways these should be listed in the sample data, not necessarily
# the call log data.
MCAD <- c("C6", "C8", "C10", "C8_C10")
LCHAD <- c("C16_OH", "C16", "C18_1_OH")
VLCAD <- c("C14_1", "C14", "C16")
PROP <- c("C3", "C3_C2")
MCD <- c("C3", "C5_OH")
MUT_CB1AB <- c("C3", "C4_DC", "C3_C2")
BKT <- c("C5_1", "C5_OH")

# Create a list of all abnormal combinations
a_s <- list(MCAD, LCHAD, VLCAD, PROP, MCD, MUT_CB1AB, BKT)

# For each set of analytes where 2 or more abnormals make a critical,
# add these results to the list of criticals for each sample
for (a in a_s) {
  QI_filter_crits$combined <- as.character(apply(QI_filter_crits, 1, getCols, colValue="combined", 
                                                 colsToCheck=a, strCheck="Critical|Abnormal", cutoff=2))
}

# Remove any records where combined is not ""
QI_filter_crits <- QI_filter_crits[QI_filter_crits$combined != "", ]

# Remove any punctuation that appears in the combined field, to make
# it easier to match on crit data
QI_filter_crits$combined <- gsub("[_/:-]", "", QI_filter_crits$combined)

# Match QI_filter_crits with crit_data_filter on ID, find which rows are not in call log
not_in_crit_data <- anti_join(QI_filter_crits, crit_data_filter, by = c("SAMPLEID" = "EXTERNAL_ID"))
not_in_QI_data <- anti_join(crit_data_filter, QI_filter_crits, by = c("EXTERNAL_ID" = "SAMPLEID"))

# Join not_in_QI_data on whatever records are in the sample data 
not_in_QI_data <- left_join(not_in_QI_data[, c("EXTERNAL_ID", "Date.of.Call", "combined")], 
                            QI_data, by = c("EXTERNAL_ID" = "SAMPLEID"))

# Order dfs by RECEIVEDATE
not_in_crit_data <- not_in_crit_data[order(not_in_crit_data$RECEIVEDATE), ]
not_in_QI_data <- not_in_QI_data[order(not_in_QI_data$RECEIVEDATE), ]

# Rename columns in not_in_QI_data to indicate what comes from call log
colnames(not_in_QI_data)[colnames(not_in_QI_data) %in% c("Date.of.Call", "combined")] <-
  c("CALL_LOG_Date.of.Call", "CALL_LOG_Reported.Results")

# Write results to csv if they exist
if (nrow(not_in_crit_data) > 0) {
  write.csv(not_in_crit_data, paste0(output_path, start_date_file, "_", end_date_file, "_", 
                                     "Samples with critical results in sample data but not found in critical result reporting data.csv"))
}

if (nrow(not_in_QI_data) > 0) {
  write.csv(not_in_QI_data, paste0(output_path, start_date_file, "_", end_date_file, "_", 
                                   "Samples with critical results in critical result reporting data but not found in sample data.csv"))
}

# Output warning message if not_in_crit_data has more than 0 rows
if (nrow(not_in_crit_data) > 0) {
  e_verb <- ifelse(nrow(not_in_crit_data) == 1, "is", "are")
  e_adj <- ifelse(nrow(not_in_crit_data) == 1, "one", "several")
  e_plur <- ifelse(nrow(not_in_crit_data) == 1, "", "s")
  e_pro_cap <- ifelse(nrow(not_in_crit_data) == 1, "This", "These")
  e_pro <- ifelse(nrow(not_in_crit_data) == 1, "this", "these")
  e_hv <- ifelse(nrow(not_in_crit_data) == 1, "has", "have")
  e_messages <- paste(not_in_crit_data$SAMPLEID, collapse=", ")
  cat(sprintf("\nWARNING: There %s %s critical result%s in the sample data that %s not\nin the critical result reporting data. %s record%s %s been saved to the\n'Samples with critical results in sample data but not found in critical\nresult reporting data' file, and the NewSTEPs report will be generated\nwithout %s sample%s. Problem sample ID%s %s:\n%s\n", 
              e_verb, e_adj, e_plur, e_verb, e_pro_cap, e_plur, e_hv, e_pro, e_plur, e_plur, e_verb, e_messages))
}

# Output warning message if not_in_QI_data has more than 0 rows
if (nrow(not_in_QI_data) > 0) {
  e_verb <- ifelse(nrow(not_in_QI_data) == 1, "is", "are")
  e_adj <- ifelse(nrow(not_in_QI_data) == 1, "one", "several")
  e_plur <- ifelse(nrow(not_in_QI_data) == 1, "", "s")
  e_pro_cap <- ifelse(nrow(not_in_QI_data) == 1, "This", "These")
  e_pro <- ifelse(nrow(not_in_QI_data) == 1, "this", "these")
  e_hv <- ifelse(nrow(not_in_QI_data) == 1, "has", "have")
  e_messages <- paste(not_in_QI_data$EXTERNAL_ID, collapse=", ")
  cat(sprintf("\nWARNING: There %s %s critical result%s in the critical result\nreporting data that %s not in the sample data. %s record%s %s\nbeen saved to the 'Samples with critical results in critical result\nreporting data but not found in sample data' file. The NewSTEPs\nreport will INCLUDE %s sample%s. Problem External_ID%s %s:\n%s\n", 
              e_verb, e_adj, e_plur, e_verb, e_pro_cap, e_plur, e_hv, e_pro, e_plur, e_plur, e_verb, e_messages))
}

# Replace values in crit data to match sample data
crit_filt_check$combined <- gsub("XLE", "MSUD", crit_filt_check$combined)
crit_filt_check$combined <- gsub("MET", "HCU", crit_filt_check$combined)
crit_filt_check$combined <- gsub("TYRI", "TYR", crit_filt_check$combined)
crit_filt_check$combined <- gsub("SCID", "TREC", crit_filt_check$combined)
crit_filt_check$combined <- gsub("CFTR", "CF", crit_filt_check$combined)

# Rename columns and join data - don't need to do a full join, because already
# have the results of the SAMPLEID differences
colnames(crit_filt_check) <- c("SAMPLEID", "crit_combined")
crits <- inner_join(QI_filter_crits, crit_filt_check, by="SAMPLEID")

# Find items in sample data not in crit data
crits$misscrits <- mapply(
  function(x, y) {
    test <- setdiff(unique(unlist(strsplit(x, ","))), 
                    unique(unlist(strsplit(y, ","))))
    if (length(test) == 0) {
      test = ""
    } else {
      test = paste0(test, collapse=",")
    }
    return(test)
  },
  x = crits$combined,
  y = crits$crit_combined
)

# Find items in crit data not in QI data
crits$missQI <- mapply(
  function(x, y) {
    test <- setdiff(unique(unlist(strsplit(x, ","))), 
                    unique(unlist(strsplit(y, ","))))
    if (length(test) == 0) {
      test = ""
    } else {
      test = paste0(test, collapse=",")
    }
    return(test)
  },
  x = crits$crit_combined,
  y = crits$combined
)

# Reduce dataframe to rows with values for either missing QI or missing crit
crits <- crits[crits$misscrit != "" | crits$missQI != "",]

# Report out results if crits has more than one row
if (nrow(crits) > 0) {
  
  # Get dates for all calls for each sample
  crit_calls <- crit_data_filter[, c("EXTERNAL_ID", "Date.of.Call")]
  crit_calls <- aggregate(Date.of.Call ~ EXTERNAL_ID, data = crit_calls, toString)
  crit_calls$Date.of.Call <- gsub(" ", "", crit_calls$Date.of.Call)
  names(crit_calls)[names(crit_calls) == 'Date.of.Call'] <- 'Date of Call(s)'
  
  # Join crits on crit_data_filter to atttach call date information
  crits <- left_join(crits, crit_calls[, c("EXTERNAL_ID", "Date of Call(s)")], by=
                       c("SAMPLEID" = "EXTERNAL_ID"))
  
  # Warning message if any values appear in misscrits
  if (sum(crits$misscrits != "") > 0 ) {
    number <- ifelse(sum(crits$misscrits != "") == 1, "sample has", "samples have")
    cat(sprintf("\nWARNING: %s %s a subset of critical results that were not\nindicated as having been called out. See the 'Results not listed as reported\nout in CRIT data' column in output file ‘Difference_in_results_called.csv’\nfor details.\n",
                sum(crits$misscrits != ""), number))
  }
  
  # Warning message if any values appear in missQI
  if (sum(crits$missQI != "") > 0 ) {
    number <- ifelse(sum(crits$missQI != "") == 1, "sample has", "samples have")
    cat(sprintf("\nWARNING: %s %s a subset of critical results that were called\nout but were not indicated as critical/critical group in the sample data source.\nSee the ‘Results not listed as being critical but were reported out' column in\noutput file ‘Difference_in_results_called.csv’ for details.\n",
                sum(crits$missQI != ""), number))     
  }
  
  # Rename colums
  names(crits)[names(crits) %in% c("combined", "crit_combined", "misscrits",
                                   "missQI")] <- c("Critical/Critical Group results from SAMPLE data",
                                                   "Critical/Critical group results listed as reported in CRIT data",
                                                   "Results not listed as reported out in CRIT data",
                                                   "Results not listed as being critical but were reported out")
  
  # Write crits to csv
  write.csv(crits, paste0(output_path, start_date_file, "_", end_date_file, "_", "Difference_in_results_called.csv"))
  
}

####################################################################################
##  Check for analytes in critical call log but not captured in variable setting  ##
####################################################################################

# Get full list of unique analytes from out-of-range data
analytes <- unique(unlist(strsplit(crit_data_filter$combined, ",")))

# Create list of all variables that have "Time_Critical" or "Non_Time_Critical"
tn_list <- ls(pattern = "Time_Critical")
n_list <- tn_list[grepl("Non", tn_list)]
t_list <- tn_list[!grepl("Non", tn_list)]

# Get full list of analytes in the variable-setting R file. Note that searching for
# the 'Time_Critical' pattern will also catch all of the Non_Time_Critical variables.
analytes_var <- unlist(lapply(tn_list, function(x) eval(parse(text = x))))

# Get difference in set between analytes and Time_Critical/Non_Time_Critical sets
missing_analytes <- setdiff(analytes, analytes_var)

# Stop code if missing_analytes has any items
if (length(missing_analytes) != 0) {
  e_begin <- ifelse(length(missing_analytes) == 1, "One analyte", "Several analytes")
  e_verb <- ifelse(length(missing_analytes) == 1, "is", "are")
  e_adj <- ifelse(length(missing_analytes) == 1, "this", "these")
  e_messages <- paste(missing_analytes, collapse=", ")
  {stop(sprintf("%s in the critical result reporting data %s not listed in either the\n'Time_Critical' or 'Non_Time_Critical' variables:\n\n%s\n\nPlease add %s in 'NewSTEPs_variable_setting.R' before running reports.", 
                e_begin, e_verb, e_messages, e_adj)) }
}

############################################################

# Add a column for each of the disorder categories to the call 
# log data, which will indicate whether that row features a disorder 
# for that category

# Perform test of columns on entire set of disorder categories from tn_list
crit_data_filter <- cbind(crit_data_filter , lapply(tn_list, 
                                                    function(x) cbind(checkCol(crit_data_filter, x))))
colnames(crit_data_filter)[grepl("structure", colnames(crit_data_filter))] <- tn_list

############################################################

# Add empty row to import file
temprow_month <- matrix(c(rep.int("", length(Monthly_Template))), nrow=1, ncol=length(Monthly_Template))
temprow_month <- as.data.frame(temprow_month)
colnames(temprow_month) <- colnames(Monthly_Template)
Monthly_Template <- rbind(Monthly_Template,temprow_month)

### Monthly Report

# State
Monthly_Template$state <- "Virginia"

# Year (4-digit)
year <- format(as.Date(start_date, format="%d/%m/%Y"),"%Y")
Monthly_Template$year <- as.numeric(year[1])

# Month (as number 1-12)
month <- format(as.Date(start_date, format="%d/%m/%Y"),"%m")
Monthly_Template$month <- as.numeric(month[1])

# DBS Samples
Monthly_Template$dbsSamples <- nrow(QI_filter)

### Unsat statistics

unsat_data <- QI_filter[!is.na(QI_filter$UNSATCODE),]

# Change column names of unsat code dataframes
colnames(improper_collection_codes) <- c("codes")
colnames(improper_transport_codes) <- c("codes")

# Match codes with samples
improper_coll <- unsat_data[unsat_data$UNSATDESC %in% improper_collection_codes$codes,]
Monthly_Template$improperCollectionCount <- nrow(improper_coll)

improper_trans <- unsat_data[unsat_data$UNSATDESC %in% improper_transport_codes$codes,]
Monthly_Template$improperTransportCount <- nrow(improper_trans)

# Test that the sum of the improper columns matches the count of non-null UNSATDESC
# values that are not equal to "Insufficient Information"
sumCheck("improper", QI_filter[!is.na(QI_filter$UNSATDESC) & QI_filter$UNSATDESC != "" &
                                 QI_filter$UNSATDESC != "Insufficient Information" &
                                 QI_filter$UNSATDESC != "Parental Refusal", ])

################################################
#####     TIMING OF SAMPLE COLLECTION      #####
#####  (COLLECTION TIME MINUS BIRTH TIME)  #####
################################################

initialc <- QI_filter

# Format times
birthhr <- sprintf("%04d",initialc$BIRTHTIME)
collecthr <- sprintf("%04d",initialc$COLLECTIONTIME)

# Paste dates and times
initialc$final_birth<- with(initialc, paste0(BIRTHDATE, birthhr))
initialc$final_collect <- with(initialc, paste0(COLLECTIONDATE, collecthr))

# Change to date/time objects
initialc$final_birth <- as.POSIXct(strptime(initialc$final_birth, format = "%Y-%m-%d%H%M"))
initialc$final_collect <- as.POSIXct(strptime(initialc$final_collect, format = "%Y-%m-%d%H%M"))

initialc$time_diff <- as.numeric(difftime(initialc$final_collect, initialc$final_birth, units="hours"))
initialc$time_diff_days <- initialc$COLLECTIONDATE - initialc$BIRTHDATE

##### TIMING OF SAMPLE COLLECTION: INITIAL SAMPLES

# Filter by initial sample recall flag
initial <- initialc[initialc$RECALL_FLAG == "N",]

# Samples collected under 12 hours
less_12 <- initial[!is.na(initial$time_diff) & initial$time_diff < 12 & initial$time_diff >= 0,]
Monthly_Template$initialDbsCollectionCount.12 <- nrow(less_12)

# Samples collected between 12 and 24 hours of birth
twelve_twenty4 <- initial[!is.na(initial$time_diff) & initial$time_diff >= 12 & initial$time_diff <= 24,]
Monthly_Template$initialDbsCollectionCount12.24 <- nrow(twelve_twenty4)

# Samples collected between 1 and 2 days after birth
one_two_days <- initial[!is.na(initial$time_diff) & initial$time_diff > 24 & initial$time_diff <= 48,]
Monthly_Template$initialDbsCollectionCount1.2 <- nrow(one_two_days)

# Samples collected between 2 and 3 days after birth
two_three_days <- initial[!is.na(initial$time_diff) & initial$time_diff > 48 & initial$time_diff <= 72,]
Monthly_Template$initialDbsCollectionCount2.3 <- nrow(two_three_days)

# Samples collected greater than 3 days after birth
three_days <- initial[(!is.na(initial$time_diff) & initial$time_diff > 72) |
                        (is.na(initial$time_diff) & !is.na(initial$time_diff_days)) & initial$time_diff_days >= 4, ]
Monthly_Template$initialDbsCollectionCount.3 <- nrow(three_days)

# Unknown collection time (includes cases where time difference is negative or where the day count
# is less than 4 and either the collection time or the birth time are unknown)
Monthly_Template$initialDbsCollectionCountTimeUnknown <- sum(is.na(initial$time_diff_days)) +
  nrow(initial[!is.na(initial$time_diff) & initial$time_diff < 0, ]) + 
  nrow(initial[is.na(initial$time_diff) & !is.na(initial$time_diff_days) & initial$time_diff_days < 4, ])

# Validation test for sums to equal expected value
sumCheck("initialDbsCollectionCount", initial)

# Get sample IDs that resulted in 'unknown' values
un_samp <- c(initial$SAMPLEID[is.na(initial$time_diff_days)],
             initial$SAMPLEID[!is.na(initial$time_diff) & initial$time_diff < 0],
             initial$SAMPLEID[is.na(initial$time_diff) & !is.na(initial$time_diff_days) & 
                                initial$time_diff_days < 4])

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "initialDbsCollectionCount", 
                                    Monthly_Template$initialDbsCollectionCountTimeUnknown, un_samp)

##### TIMING OF SAMPLE COLLECTION: SUBSEQUENT SAMPLES

# Subsequent collection
subsequent <- initialc[initialc$RECALL_FLAG == "Y",]

# Samples collected under 7 days
less_7 <- subsequent[!is.na(subsequent$time_diff_days) & subsequent$time_diff_days < 7 & 
                       subsequent$time_diff_days >= 0,]
Monthly_Template$subsequentDbsCollectionCount.7 <- nrow(less_7)

# Samples collected between 7 and 10 days of birth
seven_ten <- subsequent[!is.na(subsequent$time_diff_days) & subsequent$time_diff_days >= 7 & 
                          subsequent$time_diff_days <= 10,]
Monthly_Template$subsequentDbsCollectionCount7.10 <- nrow(seven_ten)

# Samples collected between 10 and 14 days after birth # current import file
###### -----> 11 to 14 days on newest resource
elev_four <- subsequent[!is.na(subsequent$time_diff_days) & 
                          subsequent$time_diff_days >= 11 & subsequent$time_diff_days <= 14,]
Monthly_Template$subsequentDbsCollectionCount10.14 <- nrow(elev_four)

# Samples collected greater than 14 days after birth
###### -----> greater than 15 on newest resource
fifteen_days <- subsequent[!is.na(subsequent$time_diff_days) & subsequent$time_diff_days >= 15,]
Monthly_Template$subsequentDbsCollectionCount.14 <- nrow(fifteen_days)

# Unknown subsequent collection time (includes cases where time difference in days is negative)
Monthly_Template$subsequentDbsCollectionCountTimeUnknown <- sum(is.na(subsequent$time_diff_days)) +
  nrow(subsequent[!is.na(subsequent$time_diff_days) & subsequent$time_diff_days < 0, ])

# Validation test for sums to equal expected value
sumCheck("subsequentDbsCollectionCount", subsequent)

# Get sample IDs that resulted in 'unknown' values
un_samp <- c(subsequent$SAMPLEID[is.na(subsequent$time_diff_days)],
             subsequent$SAMPLEID[!is.na(subsequent$time_diff_days) & subsequent$time_diff_days < 0])

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "subsequentDbsCollectionCount", 
                                    Monthly_Template$subsequentDbsCollectionCountTimeUnknown, un_samp)

##################################################
#####       TIMING OF SAMPLE RECEIPT         #####
#####  (RECEIPT DATE MINUS COLLECTION DATE)  #####
##################################################

# Receipt Counts
receipt <- QI_filter

# Find differences between receipt date and collection date
receipt$time_diff_days <- receipt$RECEIVEDATE - receipt$COLLECTIONDATE

##### TIMING OF SAMPLE RECEIPT: INITIAL SAMPLES

# Filter by initial sample recall flag
receipt_init <- receipt[receipt$RECALL_FLAG == "N",]

# Initial receipt count on same day 
zero <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 0,]
Monthly_Template$initialDbsReceiptCountDay0 <- nrow(zero)

# Initial receipt day after collection
one <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 1,]
Monthly_Template$initialDbsReceiptCountDay1 <- nrow(one)

# Initial receipt day 2 after collection
two <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 2,]
Monthly_Template$initialDbsReceiptCountDay2 <- nrow(two)

# Initial receipt day 3 after collection
three <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 3,]
Monthly_Template$initialDbsReceiptCountDay3 <- nrow(three)

# Initial receipt day 4 after collection
four <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 4,]
Monthly_Template$initialDbsReceiptCountDay4 <- nrow(four)

# Initial receipt day 5 after collection
five <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 5,]
Monthly_Template$initialDbsReceiptCountDay5 <- nrow(five)

# Initial receipt day 6 after collection
six <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days == 6,]
Monthly_Template$initialDbsReceiptCountDay6 <- nrow(six)

# Initial receipt day 7 or greater after collection
seven_plus <- receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days >= 7,]
Monthly_Template$initialDbsReceiptCountDay7AndGreater <- nrow(seven_plus)

# Unknown inital receipt time (includes cases where time difference is negative)
Monthly_Template$initialDbsReceiptCountTimeUnknown <- sum(is.na(receipt_init$time_diff_days)) +
  nrow(receipt_init[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days < 0, ])

# Validation test for sums to equal expected value
sumCheck("initialDbsReceiptCount", receipt_init)

# Get sample IDs that resulted in 'unknown' values
un_samp <- c(receipt_init$SAMPLEID[is.na(receipt_init$time_diff_days)],
             receipt_init$SAMPLEID[!is.na(receipt_init$time_diff_days) & receipt_init$time_diff_days < 0])

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "initialDbsReceiptCount", 
                                    Monthly_Template$initialDbsReceiptCountTimeUnknown, un_samp)

##### TIMING OF SAMPLE RECEIPT: SUBSEQUENT SAMPLES

receipt_sub <- receipt[receipt$RECALL_FLAG == "Y",]

# Subsequent receipt same day as collection
zero <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 0,]
Monthly_Template$subsequentDbsReceiptCountDay0 <- nrow(zero)

# Subsequent receipt day after collection
one <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 1,]
Monthly_Template$subsequentDbsReceiptCountDay1 <- nrow(one)

# Subsequent receipt day 2 after collection
two <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 2,]
Monthly_Template$subsequentDbsReceiptCountDay2 <- nrow(two)

# Subsequent receipt day 3 after collection
three <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 3,]
Monthly_Template$subsequentDbsReceiptCountDay3 <- nrow(three)

# Subsequent receipt day 4 after collection
four <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 4,]
Monthly_Template$subsequentDbsReceiptCountDay4 <- nrow(four)

# Subsequent receipt day 5 after collection
five <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 5,]
Monthly_Template$subsequentDbsReceiptCountDay5 <- nrow(five)

# Subsequent receipt day 6 after collection
six <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days == 6,]
Monthly_Template$subsequentDbsReceiptCountDay6 <- nrow(six)

# Subsequent receipt day 7 or greater after collection
seven_plus <- receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days >= 7,]
Monthly_Template$subsequentDbsReceiptCountDay7AndGreater <- nrow(seven_plus)

# Unknown subsequent receipt time (includes cases where time difference is negative)
Monthly_Template$subsequentDbsReceiptCountTimeUnknown <- sum(is.na(receipt_sub$time_diff_days)) +
  nrow(receipt_sub[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days < 0, ])

# Validation test for sums to equal expected value
sumCheck("subsequentDbsReceiptCount", receipt_sub)

# Get sample IDs that resulted in 'unknown' values
un_samp <- c(receipt_sub$SAMPLEID[is.na(receipt_sub$time_diff_days)],
             receipt_sub$SAMPLEID[!is.na(receipt_sub$time_diff_days) & receipt_sub$time_diff_days < 0])

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "subsequentDbsReceiptCount", 
                                    Monthly_Template$subsequentDbsReceiptCountTimeUnknown, un_samp)

###############################################################
#####     TIMING OF REPORTING TIME-CRITICAL DISORDERS     #####
#####          (REPORT TIME MINUS RECEIPT TIME)           #####
###############################################################

# Receipt to report time critical same day
day_0 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 0, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay0 <- ifelse(nrow(day_0) == 0, 0, sum(day_0))

# Receipt to report time critical day after
day_1 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 1, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay1 <- ifelse(nrow(day_1) == 0, 0, sum(day_1))

# Receipt to report time critical day 2 
day_2 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 2, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay2 <- ifelse(nrow(day_2) == 0, 0, sum(day_2))

# Receipt to report time critical day 3 
day_3 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 3, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay3 <- ifelse(nrow(day_3) == 0, 0, sum(day_3))

# Receipt to report time critical day 4 
day_4 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 4, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay4 <- ifelse(nrow(day_4) == 0, 0, sum(day_4))

# Receipt to report time critical day 5 
day_5 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 5, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay5 <- ifelse(nrow(day_5) == 0, 0, sum(day_5))

# Receipt to report time critical day 6 
day_6 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 6, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay6 <- ifelse(nrow(day_6) == 0, 0, sum(day_6))

# Receipt to report time critical day 7 or greater 
day_7_plus <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                                 crit_data_filter$diff_days_receipt2report >= 7, t_list]
Monthly_Template$receiptToReportTimeCriticalCountDay7AndGreater <- ifelse(nrow(day_7_plus) == 0, 0, sum(day_7_plus))

# Receipt to report time critical unknown (includes cases where time difference is negative)
unknown <- crit_data_filter[is.na(crit_data_filter$diff_days_receipt2report) | 
                              crit_data_filter$diff_days_receipt2report < 0, t_list]
Monthly_Template$receiptToReportTimeCriticalCountTimeUnknown <- ifelse(nrow(unknown) == 0, 0, sum(unknown))

# Validation test for sums to equal expected value
sumCheckOutOfRange("receiptToReportTimeCriticalCount", crit_data_filter, time_critical=TRUE)

# If there are any unknowns in this set, get the EXTERNAL_IDs
if (Monthly_Template$receiptToReportTimeCriticalCountTimeUnknown > 0) {
  uns <- crit_data_filter[is.na(crit_data_filter$diff_days_receipt2report) | 
                            crit_data_filter$diff_days_receipt2report < 0, 
                          c("EXTERNAL_ID", t_list)]
  uns$sum <- rowSums(uns[t_list])
  un_samp <- paste(as.character(uns$EXTERNAL_ID[uns$sum > 0]), collapse=";")
} else {
  un_samp <- NA
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Critical reporting data", "receiptToReportTimeCriticalCount", 
                                    Monthly_Template$receiptToReportTimeCriticalCountTimeUnknown, un_samp)

################################################################
#####     TIMING OF REPORTING NON-TIME-CRITICAL DISODERS   #####
#####          (REPORT TIME MINUS RECEIPT TIME)            #####
################################################################

# Receipt to report same day
day_0 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 0, n_list]
Monthly_Template$receiptToReportPositiveCountDay0 <- ifelse(nrow(day_0) == 0, 0, sum(day_0))

# Receipt to report next day
day_1 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 1, n_list]
Monthly_Template$receiptToReportPositiveCountDay1 <- ifelse(nrow(day_1) == 0, 0, sum(day_1))

# Receipt to report day 2
day_2 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 2, n_list]
Monthly_Template$receiptToReportPositiveCountDay2 <- ifelse(nrow(day_2) == 0, 0, sum(day_2))

# Receipt to report day 3
day_3 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 3, n_list]
Monthly_Template$receiptToReportPositiveCountDay3 <- ifelse(nrow(day_3) == 0, 0, sum(day_3))

# Receipt to report day 4
day_4 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 4, n_list]
Monthly_Template$receiptToReportPositiveCountDay4 <- ifelse(nrow(day_4) == 0, 0, sum(day_4))

# Receipt to report day 5
day_5 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 5, n_list]
Monthly_Template$receiptToReportPositiveCountDay5 <- ifelse(nrow(day_5) == 0, 0, sum(day_5))

# Receipt to report day 6
day_6 <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                            crit_data_filter$diff_days_receipt2report == 6, n_list]
Monthly_Template$receiptToReportPositiveCountDay6 <- ifelse(nrow(day_6) == 0, 0, sum(day_6))

# Receipt to report day 7 or greater
day_7_plus <- crit_data_filter[!is.na(crit_data_filter$diff_days_receipt2report) & 
                                 crit_data_filter$diff_days_receipt2report >= 7, n_list]
Monthly_Template$receiptToReportPositiveCountDay7AndGreater <- ifelse(nrow(day_7_plus) == 0, 0, sum(day_7_plus))

# Receipt to report unknown
unknown <- crit_data_filter[is.na(crit_data_filter$diff_days_receipt2report) | 
                              crit_data_filter$diff_days_receipt2report < 0, n_list]
Monthly_Template$receiptToReportPositiveCountTimeUnknown <- ifelse(nrow(unknown) == 0, 0, sum(unknown))

# Validation test for sums to equal expected value
sumCheckOutOfRange("receiptToReportPositiveCount", crit_data_filter, time_critical=FALSE)

# If there are any unknowns in this set, get the EXTERNAL_IDs
if (Monthly_Template$receiptToReportPositiveCountTimeUnknown > 0) {
  uns <- crit_data_filter[is.na(crit_data_filter$diff_days_receipt2report) | 
                            crit_data_filter$diff_days_receipt2report < 0, 
                          c("EXTERNAL_ID", n_list)]
  uns$sum <- rowSums(uns[n_list])
  un_samp <- paste(as.character(uns$EXTERNAL_ID[uns$sum > 0]), collapse=";")
} else {
  un_samp <- NA
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Critical reporting data", "receiptToReportPositiveCount", 
                                    Monthly_Template$receiptToReportPositiveCountTimeUnknown, un_samp)

##################################################################
#####        TIMING OF REPORTING FOR ALL RESULTS             #####
#####         (RELEASE DATE MINUS RECEIVE DATE)              #####
#####        *** NOTE: This does NOT use the ***             #####
#####   *** date of call information for out-of-range ***    #####
#####     *** results - just the RELEASEDATE from ***        #####
#####                *** the sample data ***                 #####
##################################################################

# Create copy of QI_filter
report_filter <- QI_filter

# Subtract the RELEASE DATE from the RECEIVE DATE (in days)
report_filter$time_diff_days <- report_filter$RELEASEDATE - report_filter$RECEIVEDATE

##### TIMING OF REPORTING: INITIAL SAMPLES

# Subset report_filter for records with RECALL_FLAG = "N" (initial samples)
all_report <- report_filter[report_filter$RECALL_FLAG == "N",]

# Receipt to report same day
day_0 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 0,]
Monthly_Template$receiptToReportFirstCountDay0 <- nrow(day_0)

# Receipt to report day after
day_1 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 1,]
Monthly_Template$receiptToReportFirstCountDay1 <- nrow(day_1)

# Receipt to report day 2
day_2 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 2,]
Monthly_Template$receiptToReportFirstCountDay2 <- nrow(day_2)

# Receipt to report day 3
day_3 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 3,]
Monthly_Template$receiptToReportFirstCountDay3 <- nrow(day_3)

# Receipt to report day 4
day_4 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 4,]
Monthly_Template$receiptToReportFirstCountDay4 <- nrow(day_4)

# Receipt to report day 5
day_5 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 5,]
Monthly_Template$receiptToReportFirstCountDay5 <- nrow(day_5)

# Receipt to report day 6
day_6 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days == 6,]
Monthly_Template$receiptToReportFirstCountDay6 <- nrow(day_6)

# Receipt to report day 7 or greater
day_7 <- all_report[!is.na(all_report$time_diff_days) & 
                      all_report$time_diff_days >= 7,]
Monthly_Template$receiptToReportFirstCountDay7AndGreater <- nrow(day_7)

# Receipt to report day intial samples unknown (including cases where time difference is negative)
unknown <- all_report[is.na(all_report$time_diff_days) |
                        all_report$time_diff_days < 0, ]
Monthly_Template$receiptToReportFirstCountTimeUnknown <- nrow(unknown)

# Validation test for sums to equal expected value
sumCheck("receiptToReportFirstCount", all_report)

# Get sample IDs that resulted in 'unknown' values
un_samp <- all_report$SAMPLEID[is.na(all_report$time_diff_days) |
                                 all_report$time_diff_days < 0 ]

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "receiptToReportFirstCountDay", 
                                    Monthly_Template$receiptToReportFirstCountTimeUnknown, un_samp)

##### TIMING OF REPORTING: SUBSEQUENT SAMPLES

# Subset report_filter for records with RECALL_FLAG = "Y" (subsequent samples)
all_report <- report_filter[report_filter$RECALL_FLAG == "Y",]

# Receipt to report same day
day_0_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 0,]
Monthly_Template$receiptToReportSubsequentCountDay0 <- nrow(day_0_s)

# Receipt to report day after
day_1_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 1,]
Monthly_Template$receiptToReportSubsequentCountDay1 <- nrow(day_1_s)

# Receipt to report day 2
day_2_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 2,]
Monthly_Template$receiptToReportSubsequentCountDay2 <- nrow(day_2_s)

# Receipt to report day 3
day_3_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 3,]
Monthly_Template$receiptToReportSubsequentCountDay3 <- nrow(day_3_s)

# Receipt to report day 4
day_4_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 4,]
Monthly_Template$receiptToReportSubsequentCountDay4 <- nrow(day_4_s)

# Receipt to report day 5
day_5_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 5,]
Monthly_Template$receiptToReportSubsequentCountDay5 <- nrow(day_5_s)

# Receipt to report day 6
day_6_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days == 6,]
Monthly_Template$receiptToReportSubsequentCountDay6 <- nrow(day_6_s)

# Receipt to report day 7 or greater
day_7_s <- all_report[!is.na(all_report$time_diff_days) & 
                        all_report$time_diff_days >= 7,]
Monthly_Template$receiptToReportSubsequentCountDay7AndGreater <- nrow(day_7_s)

# Receipt to report unknown (includes cases where time difference is negative)
unknown_s <- all_report[is.na(all_report$time_diff_days) |
                          all_report$time_diff_days < 0,]
Monthly_Template$receiptToReportSubsequentCountTimeUnknown <- nrow(unknown_s)

# Validation test for sums to equal expected value
sumCheck("receiptToReportSubsequentCount", all_report)

# Get sample IDs that resulted in 'unknown' values
un_samp <- all_report$SAMPLEID[is.na(all_report$time_diff_days) |
                                 all_report$time_diff_days < 0]

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "receiptToReportSubsequentCount", 
                                    Monthly_Template$receiptToReportSubsequentCountTimeUnknown, un_samp)

##################################################################
#####          TIMING FROM BIRTH TO REPORTING                #####
#####            FOR TIME-CRITICAL DISORDERS                 #####
##################################################################

# Less than or equal to 2 days after birth (but not less than 0)
day_2 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) & 
                            crit_data_filter$diff_days_birth2report >= 0 &
                            crit_data_filter$diff_days_birth2report <= 2, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay2AndLess <- ifelse(nrow(day_2) == 0, 0, sum(day_2))

# Day 3 after birth
day_3 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 3, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay3 <- ifelse(nrow(day_3) == 0, 0, sum(day_3))

# Day 4 after birth
day_4 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 4, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay4 <- ifelse(nrow(day_4) == 0, 0, sum(day_4))

# Day 5 after birth
day_5 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 5, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay5 <- ifelse(nrow(day_5) == 0, 0, sum(day_5))

# Day 6 after birth
day_6 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 6, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay6 <- ifelse(nrow(day_6) == 0, 0, sum(day_6))

# Day 7 after birth
day_7 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 7, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay7 <- ifelse(nrow(day_7) == 0, 0, sum(day_7))

# Day 8 after birth
day_8 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 8, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay8 <- ifelse(nrow(day_8) == 0, 0, sum(day_8))

# Day 9 after birth
day_9 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 9, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay9 <- ifelse(nrow(day_9) == 0, 0, sum(day_9))

# Day 10 after birth or greater
day_10_plus <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                                  crit_data_filter$diff_days_birth2report >= 10, t_list]
Monthly_Template$birthToReportTimeCriticalCountDay10AndGreater <- ifelse(nrow(day_10_plus) == 0, 0, sum(day_10_plus))

# Birth to report time unknown (including cases where time difference is negative)
unknown <- crit_data_filter[is.na(crit_data_filter$diff_days_birth2report) |
                              crit_data_filter$diff_days_birth2report < 0, t_list]
Monthly_Template$birthToReportTimeCriticalCountUnknown <- ifelse(nrow(unknown) == 0, 0, sum(unknown))

# Validation test for sums to equal expected value
sumCheckOutOfRange("birthToReportTimeCriticalCount", crit_data_filter, time_critical=TRUE)

# If there are any unknowns in this set, get the EXTERNAL_IDs
if (Monthly_Template$birthToReportTimeCriticalCountUnknown > 0) {
  uns <- crit_data_filter[is.na(crit_data_filter$diff_days_birth2report) | 
                            crit_data_filter$diff_days_birth2report < 0, 
                          c("EXTERNAL_ID", t_list)]
  uns$sum <- rowSums(uns[t_list])
  un_samp <- paste(as.character(uns$EXTERNAL_ID[uns$sum > 0]), collapse=";")
} else {
  un_samp <- NA
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Critical reporting data", "birthToReportTimeCriticalCount", 
                                    Monthly_Template$birthToReportTimeCriticalCountUnknown, un_samp)

##################################################################
#####          TIMING FROM BIRTH TO REPORTING                #####
#####          FOR NON-TIME-CRITICAL DISORDERS               #####
##################################################################

# Less than or equal to 2 days after birth (but not less than 0)
day_2 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) & 
                            crit_data_filter$diff_days_birth2report >= 0 &
                            crit_data_filter$diff_days_birth2report <= 2, n_list]
Monthly_Template$birthToReportPositiveCountDay2AndLess <- ifelse(nrow(day_2) == 0, 0, sum(day_2))

# 3 days after birth
day_3 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 3, n_list]
Monthly_Template$birthToReportPositiveCountDay3 <- ifelse(nrow(day_3) == 0, 0, sum(day_3))

# 4 days after birth
day_4 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 4, n_list]
Monthly_Template$birthToReportPositiveCountDay4 <- ifelse(nrow(day_4) == 0, 0, sum(day_4))

# 5 days after birth
day_5 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 5, n_list]
Monthly_Template$birthToReportPositiveCountDay5 <- ifelse(nrow(day_5) == 0, 0, sum(day_5))

# 6 days after birth
day_6 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 6, n_list]
Monthly_Template$birthToReportPositiveCountDay6 <- ifelse(nrow(day_6) == 0, 0, sum(day_6))

# 7 days after birth
day_7 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 7, n_list]
Monthly_Template$birthToReportPositiveCountDay7 <- ifelse(nrow(day_7) == 0, 0, sum(day_7))

# 8 days after birth
day_8 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 8, n_list]
Monthly_Template$birthToReportPositiveCountDay8 <- ifelse(nrow(day_8) == 0, 0, sum(day_8))

# 9 days after birth
day_9 <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                            crit_data_filter$diff_days_birth2report == 9, n_list]
Monthly_Template$birthToReportPositiveCountDay9 <- ifelse(nrow(day_9) == 0, 0, sum(day_9))

# 10 days after birth or greater
day_10_plus <- crit_data_filter[!is.na(crit_data_filter$diff_days_birth2report) &
                                  crit_data_filter$diff_days_birth2report >= 10, n_list]
Monthly_Template$birthToReportPositiveCountDay10AndGreater <- ifelse(nrow(day_10_plus) == 0, 0, sum(day_10_plus))

# unknown time from birth to report (including cases where time difference is negative)
unknown <- crit_data_filter[is.na(crit_data_filter$diff_days_birth2report) |
                              crit_data_filter$diff_days_birth2report < 0, n_list]
Monthly_Template$birthToReportPositiveCountUnknown <- ifelse(nrow(unknown) == 0, 0, sum(unknown))

# Validation test for sums to equal expected value
sumCheckOutOfRange("birthToReportPositiveCount", crit_data_filter, time_critical=FALSE)

# If there are any unknowns in this set, get the EXTERNAL_IDs
if (Monthly_Template$birthToReportPositiveCountUnknown > 0) {
  uns <- crit_data_filter[is.na(crit_data_filter$diff_days_birth2report) | 
                            crit_data_filter$diff_days_birth2report < 0, 
                          c("EXTERNAL_ID", n_list)]
  uns$sum <- rowSums(uns[n_list])
  un_samp <- paste(as.character(uns$EXTERNAL_ID[uns$sum > 0]), collapse=";")
} else {
  un_samp <- NA
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Critical reporting data", "birthToReportPositiveCount", 
                                    Monthly_Template$birthToReportPositiveCountUnknown, un_samp)

##################################################################
#####          TIMING FROM BIRTH TO REPORTING                #####
#####                 FOR ALL SAMPLES                        #####
#####        *** NOTE: This does NOT use the ***             #####
#####   *** date of call information for out-of-range ***    #####
#####     *** results - just the RELEASEDATE from ***        #####
#####                *** the sample data ***                 #####
##################################################################

# Get difference in days between RELEASEDATE and BIRTHDATE
report_filter$time_diff_days_birth2report <- report_filter$RELEASEDATE - report_filter$BIRTHDATE

##### TIMING OF REPORTING: INITIAL SAMPLES
all_report_birth <- report_filter[report_filter$RECALL_FLAG == "N",]

# Time from birth to report less than or equal to 2 days
day_2_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report >= 0 &
                              all_report_birth$time_diff_days_birth2report <= 2,]
Monthly_Template$birthToReportFirstScreenCountDay2AndLess <- nrow(day_2_f)

# Time from birth to report day 3
day_3_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 3,]
Monthly_Template$birthToReportFirstScreenCountDay3 <- nrow(day_3_f)

# Time from birth to report day 4
day_4_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 4,]
Monthly_Template$birthToReportFirstScreenCountDay4 <- nrow(day_4_f)

# Time from birth to report day 5
day_5_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 5,]
Monthly_Template$birthToReportFirstScreenCountDay5 <- nrow(day_5_f)

# Time from birth to report day 6
day_6_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 6,]
Monthly_Template$birthToReportFirstScreenCountDay6 <- nrow(day_6_f)

# Time from birth to report day 7
day_7_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 7,]
Monthly_Template$birthToReportFirstScreenCountDay7 <- nrow(day_7_f)

# Time from birth to report day 8
day_8_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 8,]
Monthly_Template$birthToReportFirstScreenCountDay8 <- nrow(day_8_f)

# Time from birth to report day 9
day_9_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 9,]
Monthly_Template$birthToReportFirstScreenCountDay9 <- nrow(day_9_f)

# Time from birth to report day 10 or greater 
day_10_f <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                               all_report_birth$time_diff_days_birth2report >= 10,]
Monthly_Template$birthToReportFirstScreenCountDay10AndGreater <- nrow(day_10_f)

# Time from birth to report unknown
unknown_f <- all_report_birth[is.na(all_report_birth$time_diff_days_birth2report) |
                                all_report_birth$time_diff_days_birth2report < 0, ]
Monthly_Template$birthToReportFirstScreenCountUnknown <- nrow(unknown_f)

# Validation test for sums to equal expected value
sumCheck("birthToReportFirstScreenCount", all_report_birth)

# Get sample IDs that resulted in 'unknown' values
un_samp <- all_report_birth$SAMPLEID[is.na(all_report_birth$time_diff_days_birth2report) |
                                       all_report_birth$time_diff_days_birth2report < 0]

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "birthToReportFirstScreenCount", 
                                    Monthly_Template$birthToReportFirstScreenCountUnknown, un_samp)

##### TIMING OF REPORTING: SUBSEQUENT SAMPLES

all_report_birth <- report_filter[report_filter$RECALL_FLAG == "Y",]

# Time from birth to report less than or equal to 2 days
day_2_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) &
                              all_report_birth$time_diff_days_birth2report >= 0 &
                              all_report_birth$time_diff_days_birth2report <= 2,]
Monthly_Template$birthToReportSubsequentScreenCountDay2AndLess <- nrow(day_2_s)

# Time from birth to report day 3
day_3_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 3,]
Monthly_Template$birthToReportSubsequentScreenCountDay3 <- nrow(day_3_s)

# Time from birth to report day 4
day_4_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 4,]
Monthly_Template$birthToReportSubsequentScreenCountDay4 <- nrow(day_4_s)

# Time from birth to report day 5
day_5_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 5,]
Monthly_Template$birthToReportSubsequentScreenCountDay5 <- nrow(day_5_s)

# Time from birth to report day 6
day_6_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 6,]
Monthly_Template$birthToReportSubsequentScreenCountDay6 <- nrow(day_6_s)

# Time from birth to report day 7
day_7_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 7,]
Monthly_Template$birthToReportSubsequentScreenCountDay7 <- nrow(day_7_s)

# Time from birth to report day 8
day_8_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 8,]
Monthly_Template$birthToReportSubsequentScreenCountDay8 <- nrow(day_8_s)

# Time from birth to report day 9
day_9_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                              all_report_birth$time_diff_days_birth2report == 9,]
Monthly_Template$birthToReportSubsequentScreenCountDay9 <- nrow(day_9_s)

# Time from birth to report day 10 or greater 
day_10_s <- all_report_birth[!is.na(all_report_birth$time_diff_days_birth2report) & 
                               all_report_birth$time_diff_days_birth2report >= 10,]
Monthly_Template$birthToReportSubsequentScreenCountDay10AndGreater <- nrow(day_10_s)

# Time from birth to report unknown
unknown_s <- all_report_birth[is.na(all_report_birth$time_diff_days_birth2report) |
                                all_report_birth$time_diff_days_birth2report < 0, ]
Monthly_Template$birthToReportSubsequentScreenCountUnknown <- nrow(unknown_s)

# Validation test for sums to equal expected value
sumCheck("birthToReportSubsequentScreenCount", all_report_birth)

# Get sample IDs that resulted in 'unknown' values
un_samp <- all_report_birth$SAMPLEID[is.na(all_report_birth$time_diff_days_birth2report) |
                                       all_report_birth$time_diff_days_birth2report < 0]

# Change un_samp to NA if there are no IDs
if (length(un_samp) == 0) {
  un_samp <- NA
} else {
  un_samp = paste(as.character(un_samp), collapse=";")
}

# Add these to the dataframe of unknown values
unknowns[nrow(unknowns) + 1, ] <- c("Sample data", "birthToReportSubsequentScreenCount", 
                                    Monthly_Template$birthToReportSubsequentScreenCountUnknown, un_samp)

############################################
############################################

## Rename columns to match template
for (i in 1:ncol(Monthly_Template)) {
  col = colnames(Monthly_Template[i])
  # if column name has a period ...
  if (grepl("\\.", col)) {
    # and if column name has digits on either side of the period ...
    if (grepl("\\d\\.\\d", col)) {
      # replace the period with a dash
      colnames(Monthly_Template)[i] <- gsub("\\.", "-", col)
      # else if the column name has .3 or .14 ...
    } else if (grepl("\\.3", col) | grepl("\\.14", col)) {
      # replace the period with greater than
      colnames(Monthly_Template)[i] <- gsub("\\.", ">", col)
      # else replace the period with less than
    } else {
      colnames(Monthly_Template)[i] <- gsub("\\.", "<", col)
    }
  }
}

# Final file name and location
final_file_name <- paste(start_date_file, end_date_file, "NewSTEPs.csv", sep = "_")

# Final file location
final <- paste0(output_path, final_file_name)

# Create vertical version of file for data review purposes
review <- data.frame(t(Monthly_Template))
colnames(review) <- "Values"
write.csv(review, file = paste0(output_path, start_date_file, "_", end_date_file, "_NewSTEPs_FOR_REVIEW_ONLY.csv"))

# Write file to desired location
write.csv(Monthly_Template, file = final, row.names = FALSE)

# Write unknowns file to desired location
write.csv(unknowns, file = paste0(output_path, start_date_file, "_", end_date_file, "_NewSTEPs_Unknowns.csv"),
          quote=TRUE)

# Indicate reports are complete
cat("\nReports should be complete and saved to your output folder.")