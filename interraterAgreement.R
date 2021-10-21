#+ setup, include=FALSE
library(readxl)
library(tidyverse)
library(irr)

rm(list = ls())
firstCoder <- read_delim("firstCoder.csv", delim = ";")
secondCoder <- read_delim("secondCoder.csv", delim = ";")

ratersData <- list(firstCoder, secondCoder)

variableNames <- c("focal", "design", "predictedDirection", "randomization", "sampleRestricted", "possibleCOI", "published", "controlsActive", "Att_type_exp_imp", "inLabAdministration", "controlledGameplay",  "persuasive_mechanic", "Commercial",
                   "nTotal", "nMale", "nFemale", "meanAge", "mExp", "sdExp", "mCtrl", "sdCtrl", "nExp", "nCtrl", "Fstat", "df2", "interventionDuration(m)")

ratersData <- lapply(ratersData, function(x){
  x %>% select(paperID, studyID, effectID, starts_with(variableNames)) %>%
    mutate(paperID = as.numeric(paperID),
           studyID = as.numeric(studyID),
           effectID = as.numeric(effectID),
           mExp = as.numeric(as.character(mExp)))})

datCut <- as_tibble(left_join(ratersData[[1]], ratersData[[2]], by = c("paperID", "studyID", "effectID"), suffix = c("_first", "_second"), keep = TRUE))
datCut <- datCut[,order(colnames(datCut),decreasing=TRUE)] %>% select(paperID_second, paperID_first, studyID_second, studyID_first, effectID_second, effectID_first, everything()) %>% filter(!is.na(paperID_second))
#View(datCut)

# datCut %>% select(starts_with(variableNames, ignore.case = F)) %>% apply(., 2, table)

kappas <- list()
for(i in 1:length(variableNames[14:26])){
  kappas[[i]] <- datCut %>% select(starts_with(variableNames[14:26][i], ignore.case = F)) %>% kappa2(.)
}
names(kappas) <- variableNames[14:26]

agreement <- list()
for(i in 1:length(variableNames[1:13])){
  agreement[[i]] <- datCut %>% select(starts_with(variableNames[1:13][i], ignore.case = F)) %>% agree(.)
}
names(agreement) <- variableNames[1:13]


#+ include=TRUE
#'# Cohen's kappa for metric variables
kappas
#'# Inter-rater agreement for categorical variables
agreement

#'# Mean IRR
#'
#' Mean Cohen's kappa for metric variables
meanKappa <- NA
for(i in 1:length(kappas)){
  meanKappa[i] <- kappas[[i]]$value
}
mean(meanKappa, na.rm = T)

#' Mean % agreement for metric variables
meanAgreeement <- NA
for(i in 1:length(agreement)){
  meanAgreeement[i] <- agreement[[i]]$value
}
mean(meanAgreeement, na.rm = T)
