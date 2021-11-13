#' ---
#' title: ""
#' author: ""
#' date: ""
#' output:
#'    html_document:
#'       html: true  
#'       toc: false
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

#+ include = FALSE
pdf.options(encoding='ISOLatin2.enc')
options(scipen=0.3, knitr.kable.NA = '')

list.of.packages <- c("knitr", "kableExtra", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

test <- "one-tailed"
source("functions.R")
source("esConversion.R")

# Create ES table for supplementary material
esTable <- data %>% 
  mutate(Paper = text_spec(paste(Author, ", ", pubYear, sep = ""), "html", 
                           tooltip = paste(data$studyTitle, "\n", data$doi, "Journal's H5 index: ", journalH5, "\n", "Citations: ", citations)),
         yi = round(yi, 2),
         CI = text_spec(paste("(", round(yi - 1.96*sqrt(vi), 2), ", ", round(yi + 1.96*sqrt(vi), 2), ")", sep = ""), "html", 
                        tooltip = paste("SE = ", round(sqrt(vi), 2), sep = "")),
         percFemale = text_spec(round(percFemale*100, 0), "html", 
                                tooltip = paste(nMale, " male; ", nFemale, " female", sep = "")),
         gamegenre = text_spec(gamegenre, "html", 
                               tooltip = paste(gameName)),
         interventionDuration.m. = text_spec(interventionDuration.m., "html", 
                                             tooltip = paste(posttestDelay.day., " days posttest delay", sep = "")),
         meanAge = round(meanAge, 2),
         Att_type_exp_imp = case_when(Att_type_exp_imp == 1 ~ "Explicit",
                                      Att_type_exp_imp == 0 ~ "Implicit"),
         controlsActive = case_when(controlsActive == 1 ~ "Active",
                                     controlsActive == 0 ~ "Passive"),
         persuasive_mechanic = dplyr::recode(persuasive_mechanic, 
                                      'Economical model' = "Meaningful feedback", 
                                      'Procedural rhetorics' = "Meaningful feedback",
                                      'Reward system' = "Meaningful feedback",
                                      'Procedural rhetorics/Reward system' = "Meaningful feedback")) %>%
  select("ID" = label, Paper, "N" = nTotal, "g" = yi, "95% CI" = CI,
         "%female" = percFemale, "Mean age" = meanAge, "Implicit/Explicit" = Att_type_exp_imp, 
         "Game genre" = gamegenre, "Persuasive mechanics" =  persuasive_mechanic, "Attitudes theme" = Attitudes_theme, 
         "Intervention duration (min.)" = interventionDuration.m., "Control group" = controlsActive)

#+ include = TRUE
#'## Table of effect sizes
#'#### For the paper: [Video Games and Attitude Change: A Meta-analysis](https://osf.io/4aeqt/).
#'
#'
#' Please click and hover over the former columns to see the latter information.
#' 
#' Paper -> *Title of the paper, DOI, journal's H5 index, and the number of citations*; 
#' 95% CI -> *Standard error of the estimate*; 
#' %female -> *Natural frequencies for males/females*; 
#' Game genre -> *Game name*; Intervention duration -> *Posttest delay in days*
#'
#+ echo = FALSE
esTable %>% kable(escape = FALSE) %>% kable_paper("hover") %>% 
  kable_styling( latex_options = c("striped", "HOLD_position"), full_width = T, font_size = 12, position = "left")
