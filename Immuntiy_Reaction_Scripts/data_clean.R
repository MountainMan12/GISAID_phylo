data <- read.delim("gisaid_hcov-clinical_filt.tsv", header = TRUE, sep = '\t')

data$Patient_age[data$Patient_age == "60's"] <- "61"

data$Patient_age[data$Patient_age == "40's"] <- "41"

data$Patient_age[data$Patient_age == "20 - 30"] <- "25"

data$Patient_age[data$Patient_age == "61, 11 months"] <- "62"

data$Patient_age[data$Patient_age == "53 years"] <- "53"

data$Patient_age[data$Patient_age == ">60"] <- "61"

data$Patient_age[data$Patient_age == "4 months"] <- "0.4"

data$Patient_age[data$Patient_age == "51 years"] <- "51"

data$Patient_age[data$Patient_age == "68 years"] <- "68"

data$Patient_age[data$Patient_age == "7 months"] <- "0.7"

data$Patient_age[data$Patient_age == "10-20"] <- "15"

data$Patient_age[data$Patient_age == "18-49"] <- "35"

data$Patient_age[data$Patient_age == "50-64"] <- "57"

data$Patient_age[data$Patient_age == "20s"] <- "21"

data$Patient_age[data$Patient_age == "2020-1975"] <- "45"

data$Patient_age[data$Patient_age == "2020-1986"] <- "34"

data$Patient_age[data$Patient_age == "2020-1976"] <- "44"

data$Patient_age[data$Patient_age == "60s"] <- "61"

data$Patient_age[data$Patient_age == "70s"] <- "71"

data$Patient_age[data$Patient_age == "50s"] <- "51"

data$Patient_age[data$Patient_age == "40s"] <- "41"

data$Patient_age[data$Patient_age == "80s"] <- "81"

data$Patient_age[data$Patient_age == "36, 11 months"] <- "37"

data$Patient_age[data$Patient_age == "2 months"] <- "0.2"

data$Patient_age[data$Patient_age == "1 month"] <- "0.1"

data$Patient_age[data$Patient_age == "31 years 6 months"] <- "32"

data$Patient_age[data$Patient_age == "59 years 1 months"] <- "59"

data$Patient_age[data$Patient_age == "33 years 5 months"] <- "34"

data$Patient_age[data$Patient_age == "51 years 3 months"] <- "51"

data$Patient_age[data$Patient_age == "75, 4 months"] <- "75"

data$Patient_age[data$Patient_age == "50-54"] <- "52"

data$Patient_age[data$Patient_age == "55-59"] <- "57"

data$Patient_age[data$Patient_age == "65-69"] <- "67"

data$Patient_age[data$Patient_age == "25-29"] <- "27"

data$Patient_age[data$Patient_age == "30-34"] <- "32"

data$Patient_age[data$Patient_age == "25-29"] <- "27"

data$Patient_age[data$Patient_age == "45-49"] <- "47"

data$Patient_age[data$Patient_age == "75-79"] <- "77"

write.csv(data, "gisaid_hcov-clinical.csv")


