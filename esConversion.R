# Vymazat posledny stlpec a spodne riadky
# odstranit vsetky "~?" a ";"
#odstranit . v F riadok 19; skryte znaky na konci v mEXP riadok 6-9
#vymazat riadok 131
# remove "in press"
# NIKDY NEKOPIROVAT HODNOTY DO TABULKY
# risk of bias stlpce vyhodit z data
# premenovat SE_Exp a SE_ctrl na seExp

# Read in the data
# install required R libraries if not installed already
list.of.packages <- c("car", "tidyverse", "psych", "metafor", "esc", "lme4", "ggplot2", "knitr", "puniform", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp", "magrittr", "readxl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
select <- dplyr::select

dat <- read_excel("codingSheet.xlsx", sheet = "Data")
#str(dat)
#view(dat)
dat <- dat %>% modify_at(., .at = c("pubYear"), .f = ~as.numeric(as.character(.)))
#grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", dat$mExp)

# Some data wrangling to get the right type of data (formatting of the raw dataset in Excel introduces a lot of junk otherwise)
dat$pReported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$pReported)))

# Compute gender ratio (% of female)
dat$percFemale <- dat$nFemale/(dat$nFemale + dat$nMale)

# Compute SDs from SEs and df2 from nTotal if not reported
dat <- dat %>% mutate(nExp = ifelse(is.na(nExp) & !is.na(nTotal), nTotal/2, nExp),
                      nCtrl = ifelse(is.na(nCtrl) & !is.na(nTotal), nTotal/2, nCtrl),
                      sdExp = ifelse(is.na(sdExp) & !is.na(seExp), seExp*sqrt(nExp), sdExp),
                      sdCtrl = ifelse(is.na(sdCtrl) & !is.na(seCtrl), seCtrl*sqrt(nCtrl), sdCtrl),
                      df2 = ifelse(is.na(df2) & !is.na(nTotal), nTotal - 2, df2)) # Change if other than 2 group designs are present

dat <- escalc(measure = "SMD", m1i = mExp, m2i = mCtrl, sd1i = sdExp, sd2i = sdCtrl, n1i = nExp, n2i = nCtrl, data = dat)

dat <- dat %>% mutate(yi = abs(yi) * predictedDirection)
dat <- dat %>% rowwise %>% mutate(p = as.numeric(round(tTestSummary(mExp, mCtrl, sdExp, sdCtrl, nExp, nCtrl, withinSS = FALSE)["p-value"], 5))) %>% data.frame()

# Create result ID
dat$study <- paste(dat$paperID, "/", dat$studyID, sep = "") 
dat$result <- 1:nrow(dat)

# Initialize new variables
dat$yiConv <- NA
dat$viConv <- NA
dat$useCellN <- NA

# F-Test between with df1 == 1 ---------------------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = case_when(!is.na(F) & !is.na(df1) & !is.na(df2) & df1 == 1 ~ "F1"),
                      ni = ifelse(finalDesign == "F1", df2 + 2, NA),
                      p = ifelse(finalDesign == "F1", 1-pf(F, df1, df2), NA))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$nExp + dat$nCtrl) >= (dat$ni - 2) & (dat$nExp + dat$nCtrl) <= (dat$ni + 2), 1, 0),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$yiConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$viConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$yiConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = nExp, grp2n = nCtrl, es.type = "g")$es
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$viConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = nExp, grp2n = nCtrl, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "F1") %>% select(yiConv, viConv, yi, vi, dReported, p, df2, nExp, nCtrl, ni, useCellN)

# t-tests between ---------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(t) & !is.na(df2), "tBtw", finalDesign),
                      ni = ifelse(finalDesign == "tBtw", df2 + 2, ni),
                      p = ifelse(finalDesign == "tBtw", 2*pt(abs(t), df2, lower.tail = FALSE), p))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$nExp + dat$nCtrl) >= (dat$ni - 2) & (dat$nExp + dat$nCtrl) <= (dat$ni + 2), 1, useCellN),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$yiConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$viConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$yiConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = nExp, grp2n = nCtrl, es.type = "g")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$viConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = nExp, grp2n = nCtrl, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "tBtw") %>% select(yiConv, viConv, yi, vi, dReported, p, df2, nExp, nCtrl, ni, useCellN)

# Correlation -------------------------------------------------------------
# Specify the design, compute es, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(r) & !is.na(df2), "cor", finalDesign),
                      yiConv = ifelse(finalDesign == "cor", ((2*abs(r))/sqrt(1-r^2))*(1 - (3/(4*df2 - 1))), yiConv),
                      ni = ifelse(finalDesign == "cor", df2 + 2, ni),
                      rVar = ifelse(finalDesign == "cor", dat %$% escalc(measure = "COR", ri = r, ni = df2 + 2)$vi, NA),
                      viConv = ifelse(finalDesign == "cor", (1 - (3/(4*df2 - 1))) * (4 * rVar/(1 - r^2)^3), viConv),
                      p = ifelse(finalDesign == "cor", 2*pt(abs(r*sqrt(df2 / (1 - r^2))), df2, lower.tail = FALSE), p))

# Show the converted ESs
dat %>% filter(finalDesign == "cor") %>% select(yiConv, r, rVar, viConv, p, ni, dReported)

 
# # # PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
# dat <- dat %>% mutate(gConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * dReported, gConv),
#                       gVarConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * ((reportedOverallN)/(reportedOverallN/2 * reportedOverallN/2) + (dReported^2)/(2 * (reportedOverallN))), gVarConv)) 

# Variable computations
# # Multiply the ES by -1 if not in the opposite direction
dat <- dat %>% mutate(yi = ifelse(is.na(yi) & (!is.na(yiConv) & !is.na(predictedDirection)), predictedDirection * yiConv, yi),
                      vi = ifelse(is.na(vi) & (!is.na(viConv) & !is.na(predictedDirection)), viConv, vi),
                      label = paste(paperID, "/", studyID, "/", effectID, sep = ""),
                      ni = ifelse(is.na(ni), nExp + nCtrl, ni),
                      nTerm = 2/ni,
                      persuasiveMechanicImpl = ifelse(persuasive_mechanic == "Perspective-taking", 0, ifelse(is.na(persuasive_mechanic), NA, 1)),
                      persuasiveMechanicExpl = ifelse(persuasive_mechanic == "Stereotype rehearsal", 0, ifelse(is.na(persuasive_mechanic), NA, 1)),
                      gamegenreBinary = ifelse(gamegenre %in% c("action game", "action violent game"), "action games", ifelse(is.na(gamegenre), NA, "non-action games")),
                      inLabAdministrationBinary <- ifelse(inLabAdministration == 1, 1, ifelse(is.na(inLabAdministration), NA, 0)),
                      Att_type_exp_imp <- factor(as.numeric(Att_type_exp_imp), levels = c(0, 1)))

# GRIM & GRIMMER Test
dat <- dat %>% mutate(items = ifelse(is.na(items), 0, items))
grim(dat)
grimmer(dat)

# Create dats objects
rob2 <- read_delim("ROB2_IRPG_beta_v8.csv", delim = ";")
data <- left_join(dat, rob2[,-8], by = "paperID") %>% mutate(overallBias = factor(overallBias, levels = c("Low", "Some concerns", "High"))) %>% filter(!is.na(yi))
datExplicit <- data %>% filter(Att_type_exp_imp == 1)
datImplicit <- data %>% filter(Att_type_exp_imp == 0)

# Remove outliers (based on the results from the maDiag script)
# dat <- dat %>% filter(!result %in% c())

# dat %>% select(label, yi, vi, yiConv, dReported, predictedDirection, mExp, mCtrl, sdExp, sdCtrl, seExp, seCtrl, nExp, nCtrl, df2, F, t, r) %>% view()

# dat %>% filter(ni != nTotal) %>% select(ni, nTotal, df2, nExp, nCtrl)
# dat %>% filter(abs(dReported - gConv) > .2) %>% select(result, gConv, dReported)
# 
# dat$studentSample <- ifelse(dat$sampleType == "student", 1, ifelse(dat$sampleType == "general", 0, NA))

# Remove outliers (based on the results from the maDiag script)
# dat <- dat %>% filter(!result %in% c(194))

# Create dat objects
# Mood and overall effect data objects

