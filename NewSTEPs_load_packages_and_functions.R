# Install libraries
libs <- c('data.table',
          'dplyr')

for (l in libs) {
  if(!is.element(l, .packages(all.available = TRUE)) ) {
    install.packages(l)
  }
  suppressPackageStartupMessages(library(l, character.only=TRUE))
}

# Set file separator
slash <- ifelse(comp_type == 'PC', '\\', '/')

# Read in example file for use in removing rows that appear from call logs
examps <- data.frame(suppressWarnings(read.csv(examples, stringsAsFactors = FALSE)))

# Strip odd characters from column names
colnames(examps) <- gsub("ï..", "", colnames(examps))

stopQuietly <- function(...) {
  
  # Stops a source file quietly (without printing an error message), used in cases
  # where we have multiple files that need to stop running, but only have one of them
  # throw an error.
  
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
  
} 

# Reformat start date and end date as dates
if (exists("start_date")) {
  start_date <- as.Date(start_date, "%m/%d/%Y")
  end_date <- as.Date(end_date, "%m/%d/%Y")
  
  if (is.na(start_date)) {
    message = "\nYour start date is invalid. Please enter a valid start date."
    cat(message)
    stopQuietly()
  }
  
  if (is.na(end_date)) {
    message = "\nYour end date is invalid. Please enter a valid end date."
    cat(message)
    stopQuietly()
  }
  
}

removeExamples <- function(examples, df) {
  
  # Removes from call log data any rows that appear
  # in the example file
  
  # Get subset of columns from df and examples that are likely to be formatted the same way. 
  df_sub = subset(df, select = c(EXTERNAL_ID, BIRTH.TIME, Critical.1:Critical.5))
  temp_sub = subset(examples, select = c(EXTERNAL_ID, BIRTH.TIME, Critical.1:Critical.5))
  
  # Bind both dataframes together
  both = rbind(df_sub, temp_sub)
  
  # Get which indices for duplicated rows in 'both' (reading from end of the dataframe to 
  # the beginning). This will give the indices in df that need to be removed (i.e., are
  # the same as the observations in 'examples')
  indices = which(duplicated(both, fromLast = TRUE))
  
  # Remove duplicated indices from df if indices has any values
  if (length(indices) > 0) {
    df = df[-indices, ]
  }
  
  # Return the filtered df
  return(df)
  
}

get_file_list <- function(folder) {
  
  # Returns list of .txt or .csv files in a folder
  
  files <- list.files(folder, pattern = "csv|txt")
  temp <- paste0(folder, slash, files)
  
  return(temp)
  
}

substrRight <- function(x, n){
  
  # Returns a string of specified length counting backwards
  
  return(substr(x, nchar(x)-n+1, nchar(x)))
  
}

date_reformat <- function(df, ...) {
  
  # Reformats set of columns as dates in a dataframe, first checks to see how year is 
  # formatted
  
  date_cols = list(...)
  
  for (col in unlist(date_cols)) {
    
    # check to see if last three characters in string contain a "/"
    # (thus indicating the year is in 2-digit rather than 4-digit format)
    
    if (grepl("/", substrRight(df[col, 1], 3))) {
      df[, col] = as.Date(df[, col], "%m/%d/%y")
    } else {
      df[, col] = as.Date(df[, col], "%m/%d/%Y")
    }
  }
  
  return(df)
}

date_repair <- function(df, ...) {
  
  # Repairs dates that have been read incorrectly (e.g., if year is listed as '0016' 
  # instead of '2016'). Dates for checking should have already been formatted as dates.
  
  date_cols = list(...)
  
  # Pull all records where the dates do not appear to have been formatted correctly
  for (col in unlist(date_cols)) {
    
    rows = which(!is.na(df[, col]) & df[, col] < as.Date("1970-01-01")) 
    
    for (row in rows) {
      
      if (as.numeric(substr(df[row, col], 3, 4)) > as.numeric(format(Sys.Date(), "%y"))) {
        cent = "19"
      } else {
        cent = "20"
      }
      
      repaired = gsub("00", cent, df[row, col])
      df[row, col] = as.Date(repaired, format = "%Y-%m-%d")
      
    }
    
  }
  
  return(df)
  
}

# Create list of analytes to check for in sample data
analytes_check <- c("T4", "TSH", "PKU", "MSUD", "HCU", "GALT", "GALH", "TGAL", "BIOT",
                    "HGB", "CAH", "C8", "C0", "C6", "C8_C10", "C10", "C14", "C14_1", "C16",
                    "C16_OH", "C18_1_OH", "C2", "C3", "C3_C2", "C4_DC", "C5", "C5_1", "C5_DC",
                    "C5_OH", "TYR", "CIT", "PHE_TYR", "IRT", "CF", "SUAC", "TREC")

col_check <- function(folder, type) {
  
  # Checks that a given list of files have the correct columns, returns
  # list of 'bad files' along with column names that are not found in the
  # files indicated.
  
  # Get file list
  temp = get_file_list(folder)
  
  # Loop through files to ensure they all have the correct columns
  bad_files = c()
  
  if (type == "sample") {
    cols = c("SAMPLEID", "BIRTHDATE", "BIRTHTIME", "COLLECTIONDATE", "COLLECTIONTIME", "RECEIVEDATE", 
             "RECEIVETIME", "RELEASEDATE", "UNSATCODE", "UNSATDESC", "RECALL_FLAG", "CATEGORY", 
             analytes_check)
    separator = "|"
  } else if (type == "call_log") {
    cols = c("EXTERNAL_ID", "Critical.1", "Critical.2", "Critical.3", "Critical.4", 
             "Critical.5", "Date.of.Call")
    separator = ","
  } else if (type == "approve_log") {
    cols = c("Control.ID", "Test.Name", "Approved.Date", "Task.Appr.Status")
    separator = ","
  }
  
  for (f in temp) {
    # Read in single row from each file to get the column names
    temp_file = suppressWarnings(read.table(f, nrows = 1, header = TRUE, sep=separator, 
                                            fill = TRUE))
    # Strip odd characters from column names
    colnames(temp_file) <- gsub("ï..", "", colnames(temp_file))
    # Find set difference between expected columns and columns in temp_file
    bf_temp = setdiff(cols, colnames(temp_file))
    if (length(bf_temp) != 0) {
      bad_files = c(bad_files, paste0(f, ":\n   missing columns - ", paste(bf_temp, collapse = ", "), "\n"))
    } 
    
  }
  
  if (length(bad_files) != 0) {
    e_begin <- ifelse(length(bad_files) == 1, "One file", "Several files")
    e_verb <- ifelse(length(bad_files) == 1, "does", "do")
    e_art <- ifelse(length(bad_files) == 1, "this", "these")
    e_plur <- ifelse(length(bad_files) == 1, "", "s")
    message = paste0(e_begin, " in the ", type, "_data_path location ", e_verb, 
                     " not have the required column\nheadings. Please check the column headings for ",
                     e_art, " file", e_plur, " before running reports:\n", paste0(bad_files, collapse=""))
    
    cat(message)
    stopQuietly()
  }
  
}

get_dates <- function(df) {
  
  # Gets minimum and maximum date for the column data will be 
  # filtered by. Used as input for date_compare.
  
  date_col = "RECEIVEDATE"
  
  # Get the earliest date from date_col
  min_date = min(df[, date_col], na.rm = TRUE)
  
  # Get the latest date from date_col
  max_date = max(df[, date_col], na.rm = TRUE)
  
  return(list(min_date, max_date))
  
}

date_compare <- function(data_type, start_or_end, compare_obj) {
  
  # Checks data to see if start date or end date desired by user
  # is outside of bounds of data source
  
  # compare_obj is a list that contains minimum and maximum date
  # from the data source (created by date_check function)
  
  comp_s = as.Date(compare_obj[1][[1]]) # minimum date
  comp_e = as.Date(compare_obj[2][[1]]) # maximum date
  
  if (start_or_end == "start") {
    adj <- "earliest"
    adv <- "earlier"
    date <- start_date
    compare_date <- comp_s
  } else {
    adj <- "latest"
    adv <- "later"
    date <- end_date
    compare_date <- comp_e
  }
  
  # Stop report completely if requested start date is after the end date in the data
  # (regardless of whether user has indicated a check for start or end)
  if (start_date > comp_e) {
    message <- paste0("\nYour desired start_date, ", start_date, ", is later than any dates for RECEIVEDATE\nin your ",
                      data_type, " data source (the latest of which is ", comp_e, 
                      ").\nBecause there is no overlap between your requested dates and the dates in\nyour data, the report generation will stop.")
    cat(message)
    stopQuietly()
  }
  
  # Stop report completely if requested end date is before the start date in the data
  # (regardless of whether user has indicated a check for start or end)
  if (end_date < comp_s) {
    message <- paste0("\nYour desired end_date, ", end_date, ", is earlier than any dates for RECEIVEDATE\nin your ",
                      data_type, " data source (the earliest of which is ", comp_s, 
                      ").\nBecause there is no overlap between your requested dates and the dates in\nyour data, the report generation will stop.")
    cat(message)
    stopQuietly()
  }
  
  # Allow user to choose to stop report if start date is after minimum date in the data or end date
  # is before maximum date in the data
  if ( (start_or_end == "start" & comp_s > start_date) || (start_or_end == "end" & comp_e < end_date) ) {
    message <- paste0("\nYour desired ", start_or_end, "_date, ", date, ", is ", adv, 
                      " than any dates for RECEIVEDATE\nin your ",
                      data_type, " data source (the ", adj , " of which is ", as.Date(compare_date),
                      ").\n\nDo you wish to proceed? (enter 'Y' or 'y' for yes)\n")
    cat(message)                
    ans <- readline(prompt = "> ")
    if (tolower(ans) != 'y') {
      cat("\nStopping per your request.")
      stopQuietly()
    } else {
      message <- paste0("\nContinuing with report generation with your selection for ", start_or_end, "_date.\nThis will take a few moments.\n")
      cat(message)
    }
  }
}

date_comp_check <- function(df, data_type) {
  
  # Compares start date and end date entered by user with
  # minimum and maximum dates of filter column.
  
  # Get minimum and maximum dates for filter column in df
  date_test <- get_dates(df)
  
  # Check that earliest date for filt_col in data is earlier than requested start_date
  date_compare(data_type, "start", date_test)
  
  # Check that latest date for filt_col in data is later than requested end_date
  date_compare(data_type, "end", date_test)
  
}

trim <- function (x) {
  
  # Remove spaces at the end or the beginning of any elements
  gsub("^\\s+|\\s+$", "", x)
  
}

read_data <- function(folder, patt=NULL, separator, type, ...) {
  
  # Returns dataframe of data. Optional arguments are columns to be reformatted as dates
  # (for use with csv and txt files).
  
  # Make list of columns to be reformatted as dates
  date_cols = list(...)
  
  # Get file list
  temp <- get_file_list(folder)
  
  # If pattern is identified, get files that contain pattern
  if (!is.null(patt)) {
    temp <- temp[grepl(patt, temp)]
  }
  
  # read in data and bind together
  lst <- lapply(temp, function(x) suppressWarnings(read.csv(x, stringsAsFactors = FALSE, 
                                                            header=TRUE, sep=separator, fileEncoding="latin1")))
  initial_dd <- rbindlist(lst, fill=TRUE)
  
  # read in data second time in order to get character vector of SAMPLEID
  sub_lst <- lapply(temp, function(x) suppressWarnings(read.csv(x, stringsAsFactors = FALSE, header=TRUE, 
                                                                sep=separator, fileEncoding="latin1",
                                                                colClasses="character")))
  
  sub_dd <- rbindlist(sub_lst, fill=TRUE)
  
  # replace initial_dd SAMPLEID with sub_dd version (to keep leading zeros)
  if ("SAMPLEID" %in% names(initial_dd)) {
    initial_dd$SAMPLEID <- as.character(sub_dd$SAMPLEID)
  }
  
  # change initial_dd to a dataframe
  initial_dd <- data.frame(initial_dd)
  
  # Strip odd characters from column names
  colnames(initial_dd) <- gsub("ï..", "", colnames(initial_dd))
  
  # Remove any records from temp_file that appear in 'examples' and strip spaces from critical results
  if (type == "call_log") {
    
    initial_dd <- removeExamples(examps, initial_dd)
    
    for (c in colnames(initial_dd)[grepl("Critical", colnames(initial_dd))]) {
      initial_dd[, c] <- trim(initial_dd[, c])
    }
    
  }
  
  # Reformat and repair any specified columns as dates
  if (length(date_cols) > 0) {
    initial_dd <- date_reformat(initial_dd, date_cols)
    initial_dd <- date_repair(initial_dd, date_cols)
  }
  
  # Reformat SUBMITTERID as character
  if("SUBMITTERID" %in% names(initial_dd)) {
    initial_dd$SUBMITTERID <- as.character(initial_dd$SUBMITTERID)
  }
  
  # Replace 9999 values in transit time column with NA
  if("TRANSITTIME" %in% names(initial_dd)) {
    initial_dd$TRANSIT_TIME[initial_dd$TRANSIT_TIME == 9999] <- NA
  }
  
  # If dataframe has CATEGORY column, remove any records that have category listed as "Proficiency", 
  # "Treatment", or "Treatment - PKU"
  remove_cats <- c("Proficiency","Treatment","Treatment - PKU")
  if (!is.null(initial_dd$CATEGORY)) {initial_dd <- initial_dd[!(initial_dd$CATEGORY %in% remove_cats),]}
  
  # For approve_log data:
  if (type == "approve_log") {
    
    # Filter for just results containing 'Critical'
    initial_dd <- initial_dd[grepl("critical", tolower(initial_dd$Task.Appr.Status)), ] 
    
    # Select and rename columns of interest
    initial_dd <- initial_dd %>%
      mutate(EXTERNAL_ID = Control.ID,
             combined = toupper(Test.Name),
             Date.of.Call = Approved.Date) %>%
      select(EXTERNAL_ID, combined, Date.of.Call)
    
    # Restructure data to group by EXTERNAL_ID and Date.of.Call
    initial_dd <- aggregate(combined ~ EXTERNAL_ID + Date.of.Call, data = initial_dd, toString)
    initial_dd$combined <- gsub(" ", "", initial_dd$combined)
  }
  
  # For call_log data:
  if (type == "call_log") {
    
    # Combine all analytes called for a sample into a single cell, separated by commas
    initial_dd$combined <- apply(subset(initial_dd, select=Critical.1:Critical.5), 1, function (x) 
      gsub(" ", "," ,sub("\\s+$", "", gsub("NA", "", paste(toupper(x), collapse=" ")))))
    
    # Select columns of interest
    initial_dd <- initial_dd %>%
      select(EXTERNAL_ID, combined, Date.of.Call)
    
  }
  
  # Remove any '-', '/', or ':' characters from combined field from critical data 
  # (to limit the number of different representations of analytes)
  if("combined" %in% names(initial_dd)) {
    initial_dd$combined <- gsub("[/:-]", "", initial_dd$combined)
  }
  
  # Reformat EXTERNAL_ID as character
  if("EXTERNAL_ID" %in% names(initial_dd)) {
    initial_dd$EXTERNAL_ID <- as.character(initial_dd$EXTERNAL_ID)
  }
  
  # Remove LINKID if it appears in the data - the same sample can have more
  # than one LINKID because of the transition from V9 to V10
  if("LINKID" %in% names(initial_dd)) {
    initial_dd <- subset(initial_dd, select=-LINKID)
  }
  
  # Remove any duplicate rows
  initial_dd <- unique(initial_dd)
  
  return(initial_dd)
  
}

checkString <- function(cell_to_check, check_set) {
  
  # Given a cell to check, returns 1 if any substring in that cell matches an element
  # in a set of values. Assumes that the check_cell has values separated by
  # commas (e.g., T4,TSH,CAH). Used to see whether a sample in the call log 
  # has any analyte listed for a particular category of disorders.
  
  val = ifelse(any(unlist(strsplit(cell_to_check, split=",")) %in% check_set), 1, 0)
  return(val)
}

checkCol <- function(df, colString) {
  
  # Given a dataframe and a column expressed as a string, perform checkString on a
  # set of values of interest and add a new column with the results of the check.
  # Used to add the results from checkString across an entire column.
  
  vals <- apply(df, 1, function (x) checkString(x, eval(parse(text = colString))))
  return(vals)
  
}

sumCheck = function(cols, df) {
  
  # Checks that the sum of calculated values across a set of columns in
  # the Monthly Template matches the expected value (nrows of a filtered
  # dataframe)
  
  expect_sum <- nrow(df)
  test_sum <- rowSums(Monthly_Template[grepl(cols, colnames(Monthly_Template))])
  
  # If the sums do not match, create an error message, print it out, and stop the running
  # of the report.
  if (expect_sum != test_sum) {
    message = paste0("\nThe sum of the amounts for the '", cols, 
                     "' columns is ", test_sum, ", which does not\nmatch the expected value, ", 
                     expect_sum, ". Please review the 'NewSTEPs_QIs.R' file\nby searching for '", cols, 
                     "' to determine the source of the problem.")
    cat(message)
    stopQuietly()
    
  }
  
}

sumCheckOutOfRange = function(cols, df, time_critical) {
  
  # Checks that the sum of calculated values across a set of columns for
  # out of range results (either time-critical or non-time-critical) in
  # the Monthly Template matches the expected value (sum of all time-criticals
  # or non-time-criticals)
  
  # Accepted values for time_critical_indic are TRUE (for time-critical out-of-range
  # results) or FALSE (for non-time-critical out-of-range results)
  
  # Test for time_critical being entered correctly
  if (time_critical == TRUE) {
    col_check = t_list
  } else if (time_critical == FALSE) {
    col_check = n_list
  } else {
    cat("\nPlease enter either TRUE or FALSE for the third argument in the 'sumCheckOutOfRange'\nfunction; TRUE if you wish to filter for time-critical results and FALSE if you\nwish to filter for non-time-critical results.")
    stopQuietly()
  }
  
  df <- data.frame(df)
  
  expect_sum <- sum(df[, col_check])
  test_sum <- rowSums(Monthly_Template[grepl(cols, colnames(Monthly_Template))])
  
  # If the sums do not match, create an error message, print it out, and stop the running
  # of the report.
  if (expect_sum != test_sum) {
    message = paste0("The sum of the amounts for the '", cols, 
                     "' columns is ", test_sum, ", which does not\nmatch the expected value, ", 
                     expect_sum, ". Please review the 'NewSTEPs_QIs.R' file\nby searching for '", cols, 
                     "' to determine the source of the problem.")
    cat(message)
    stopQuietly()
    
  }
  
}

getCols <- function(df_row, colValue, colsToCheck, strCheck, cutoff=1) {
  
  # Returns a list of all analyte columns that feature a certain value,
  # when the number of columns that contain that value is greater than
  # a cutoff value
  
  # df_row -      a row of sample data that features columns of analyte results
  # colValue -    which column in the df_row that has a stored list of analytes
  #               that have already been checked
  # colsToCheck - which set of columns to check for result values 
  # strCheck    - string to check for in analyte columns (e.g., 'Abnormal',
  #               'Critical', or both - 'Abnormal|Critical')
  # cutoff      - count of values equivalent to a Critical. For example,
  #               this should be set to 2 when checking for combinations of
  #               Abnormal or Critical results, but set to 1 when just checking
  #               for Critical results (since a single critical is treated as 
  #               a critical)
  
  # if length of vector containing string is greater than cutoff:
  if (sum(grepl(strCheck, df_row[colsToCheck])) >= cutoff) {
    
    # Creates string containing all column names that contain the strCheck
    cols <- paste0(names(unlist(lapply(df_row[colsToCheck], function(x) which(grepl(strCheck, x))))),
                   collapse=",")
    
  } else {
    
    cols <- ""
    
  }
  
  # If value of colValue is not NULL, NA, or "", and if cols != "", 
  # paste the current results to the previous value; otherwise cols = previous value
  if (!is.null(df_row[colValue]) & !is.na(df_row[colValue]) & as.character(df_row[colValue]) != "") {
    
    if (cols != "") {
      cols <- paste0(as.character(df_row[colValue]), ",", cols)
    } else {
      cols <- as.character(df_row[colValue])
    }
    
  }
  
  # Remove duplicates from cols
  if (cols != "") {
    col_list <- unique(unlist(strsplit(cols, ",")))
    cols <- paste(col_list, collapse=",")
    
  }
  
  return(cols)
  
}