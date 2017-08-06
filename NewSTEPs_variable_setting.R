# Change to match desired date range and locations on your computer

# Set computer type (enter 'PC' or 'MAC')
comp_type <- 'MAC'


# Read in directory where you have your NewSTEPs R files stored
# (remember to use "\\" for Windows, "/" for Mac)
ns_wd <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/r_files/NewSTEPs/"


# Read in 360 template
Monthly_Template_path <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/r_files/NewSTEPs/NewSTEPS360_monthly_import.csv"


# Read in example records (these are the 3 records that may appear in call log data)
examples <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/r_files/NewSTEPs/call_log_examples.csv"


# Path for call log data
call_log_path <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/NewSTEPS_data"


# Path for approved dates (for determining date of call for samples for period before 
# call logs were used, October 2016)
approved_path <- ""


# Read in data to be used for quality indicator calculations (change location
# to your own location)
sample_data_path <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/Sample"


# Output path (where to save the NewSTEPs report and the validation reports)
output_path <- "/Users/chrispatrick/Documents/Classes/Fall 2016/DS 6001/Newborn Screening/FINAL DELIVERABLES (and other code)/Output/NewSTEPs_final_reports/"


# Add or remove Unsat codes
improper_collection_codes <- as.data.frame(c("Improperly Collected", "Scratched or Abraded", "Wet", 
                                             "Oversaturated", "Contaminated", "No Blood", 
                                             "Infant > 6 months old", "Outdated filter paper card",
                                             "Insufficient Quantity", "Interfering Substances Present"))

improper_transport_codes <- as.data.frame(c("Old Sample > 10 days in transit", "Other"))

# List of time-critical vs. non time-critical disorders by category.
# NOTE: If an analyte needs to be added, change it to ALL CAPITALS and
# remove any punctuation (e.g., "phe/tyr" should be listed as "PHETYR").
Time_Critical_Metabolic     <- c("CIT", "MSUD", "XLE", "C14", "C141", "C16", "C16OH",
                                 "C181OH", "C8", "C6", "C8C10", "C10", "C3", "C5OH", 
                                 "C4DC", "C3C2", "C51", "C5", "C5DC")
Time_Critical_Galactosemia  <- c("TGAL", "GALT", "GALH")
Time_Critical_CAH           <- c("CAH")
Non_Time_Critical_Metabolic <- c("PKU", "PHETYR", "HCU", "MET", "SUAC", "TYR", "TYRI", "C0")
Non_Time_Critical_Hemo      <- c("HGB")
Non_Time_Critical_CF        <- c("CF", "CFTR", "IRT")
Non_Time_Critical_CH        <- c("T4", "TSH")
Non_Time_Critical_BIOT      <- c("BIOT")
Non_Time_Critical_SCID      <- c("SCID", "TREC")

######### Run Report ###############

# Choose start and end dates
start_date <- "07/01/2016"
end_date <- "07/31/2016"

source(paste0(ns_wd, "NewSTEPs_QIs.R"))
