# Load the necessary library
library(dplyr)
library(irr)
library(readxl)

# Read in the data
irrData <- read_excel("inclusionInterraterAgreemeentData.xlsx")

# Recode "K" to "Y" and calculate the Cohen's kappa and percentage agreement
irrData <- irrData %>% mutate(across(c(ratingFinalRater2, ratingFinalRater1), ~recode(factor(.), "K" = "Y")))

# Cohen's kappa
kappa <- kappa2(irrData[,c("ratingFinalRater1", "ratingFinalRater2")])

# Percentage agreement
perc_agreement <- agree(irrData[,c("ratingFinalRater1", "ratingFinalRater2")])

# Results
list(kappa = kappa, perc_agreement = perc_agreement)
