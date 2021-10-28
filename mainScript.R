#' ---
#' title: "The effect of videogames on attitude change - a meta-analysis"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

#' **This is the supplementary analytic output for the paper The effect of videogames on attitudes - a meta-analysis **
#' 
#' **It reports detailed results for all models reported in the paper. The analytic R script by which this html report was generated can be found on the project's OSF page at: https://osf.io/4aeqt/**
#' 
#' ------------------------------------
#' 
#' **Brief information about the methods used in the analysis:**
#' 
#' **RMA results with model-based SEs**
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#' **RVE SEs with Satterthwaite small-sample correction**
#' Estimate based on a multilevel RE model with constant sampling correlation model (CHE - correlated hierarchical effects - working model) (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/). 
#' Interpretation of naive-meta-analysis should be based on these estimates.
#'
#' **Prediction interval**
#' Shows the expected range of true effects in similar studies.
#' As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between PI LB and PI UB.
#' Note that these are non-adjusted estimates. An unbiased newly conducted study will more likely fall in an interval centered around bias-adjusted estimate with a wider CI width.
#'
#' **Heterogeneity**
#' Tau can be interpreted as the total amount of heterogeneity in the true effects. 
#' I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates. Estimates calculated by two approaches are reported.
#' This is followed by separate estimates of between- and within-cluster heterogeneity and estimated intra-class correlation of underlying true effects.
#' 
#' **Proportion of significant results**
#' What proportion of effects were statistically at the alpha level of .05.
#' 
#' **ES-precision correlation**
#' Kendalls's correlation between the ES and precision.
#' 
#' **4/3PSM**
#' Applies a permutation-based, step-function 4-parameter selection model (one-tailed p-value steps = c(.025, .5, 1)). 
#' Falls back to 3-parameter selection model if at least one of the three p-value intervals contains less than 5 p-values.
#' For this meta-analysis, we applied 3-parameter selection model by default as there were only 11 independent effects in the opposite direction overall (6%), causing the estimates to be unstable across iterations.
#' pvalue = p-value testing H0 that the effect is zero. ciLB and ciUB are lower and upper bound of the CI. k = number of studies. steps = 3 means that the 4PSM was applied, 2 means that the 3PSM was applied.
#' 
#' **PET-PEESE**
#' Estimated effect size of an infinitely precise study. Using 4/3PSM as the conditional estimator instead of PET (can be changed to PET). If the PET-PEESE estimate is in the opposite direction, the effect can be regarded nil. 
#' By default (can be changed to PET), the function employs a modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/). 
#' It also uses the same RVE sandwich-type based estimator in a CHE (correlated hierarchical effects) working model with the identical random effects structure as the primary (naive) meta-analytic model. 
#' 
#' We report results for both, PET and PEESE, with the first reported one being the primary (based on the conditional estimator).
#' 
#' **WAAP-WLS**
#' The combined WAAP-WLS estimator (weighted average of the adequately powered - weighted least squares) tries to identify studies that are adequately powered to detect the meta-analytic effect. 
#' If there is less than two such studies, the method falls back to the WLS estimator (Stanley & Doucouliagos, 2015). If there are at least two adequately powered studies, WAAP returns a WLS estimate based on effects from only those studies.
#' 
#' type = 1: WAAP estimate, 2: WLS estimate. kAdequate = number of adequately powered studies
#' 
#' **p-uniform**
#' P-uniform* is a selection model conceptually similar to p-curve. It makes use of the fact that p-values follow a uniform distribution at the true effect size while it includes also nonsignificant effect sizes.
#' Permutation-based version of p-uniform method, the so-called p-uniform* (van Aert, van Assen, 2021).
#' 
#' **p-curve**
#' Permutation-based p-curve method. Output should be self-explanatory. For more info see p-curve.com
#' 
#' **Power for detecting SESOI and bias-corrected parameter estimates**
#' Estimates of the statistical power for detecting a smallest effect sizes of interest equal to .20, .50, and .70 in SD units (Cohen's d). 
#' A sort of a thought experiment, we also assumed that population true values equal the bias-corrected estimates (4/3PSM or PET-PEESE) and computed power for those.
#' 
#' **Handling of dependencies in bias-correction methods**
#' To handle dependencies among the effects, the 4PSM, p-curve, p-uniform are implemented using a permutation-based procedure, randomly selecting only one focal effect (i.e., excluding those which were not coded as being focal) from a single study and iterating nIterations times.
#' Lastly, the procedure selects the result with the median value of the ES estimate (4PSM, p-uniform) or median z-score of the full p-curve (p-curve).

#+ setup, include = FALSE
# NOTE: Please note that to run the script, you need the development versions of metafor and dmetar packages from github.
knitr::opts_chunk$set(echo = FALSE, warning = FALSE)

rm(list = ls())

# Settings ----------------------------------------------------------------

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
rmCor <- 0.5

# Assumed constant sampling correlation
rho <- 0.5

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
side <- "right"

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 500 # Set to 5 just to make code checking/running fast.
nIterationsPcurve <- 300
nIterationVWsensitivity <- 300 # Number of iterations for the Vevea & Woods (2005) step function model sensitivity analysis 

# Number of chains and iterations for Robust Bayesian model-averaging approach
runRobMA <- TRUE
robmaChains <- 2
robmaSamples <- 1000

# Controls for PET-PEESE
condEst <- FALSE

# Controls for the multiple-parameter selection models 

# Whether to apply a 4- or 3-parameter selection model. If fallback == TRUE, the procedure falls back to the 3-parameter selection model. 
# This should be selected when too few effects in the opposite side make the estimate unstable.
fallback <- TRUE
# Even when fallback == FALSE, the 4-parameter selection model still falls back to 3 parameters for the given iteration if,
# (1) it fails to converge or (2) the number of p-values in each of the step intervals gets smaller than minPvalues.
minPvalues <- 4

# Steps and delta parameters for Vevea & Woods selection models 
# Can be adjusted if a different selection process is assumed. 
# Please note that steps vector represents one-tailed p-values.
stepsDelta <- data.frame(
  steps =     c(.0025, .005, .0125, .025, .05, .10, .25, .50, 1),
  moderateSelection = c(1, 0.99, 0.97, 0.95, 0.80, 0.60, 0.50, 0.50, 0.50),
  severeSelection =   c(1, 0.99, 0.97, 0.95, 0.65, 0.40, 0.25, 0.25, 0.25),
  extremeSelection =  c(1, 0.98, 0.95, 0.90, 0.50, 0.20, 0.10, 0.10, 0.10))

# Sourcing and dat -----------------------------------------------------------------
#+ include = FALSE
source("functions.R")
source("pcurvePlotOption.R")
source("esConversion.R")
statcheck <- read_csv("statcheck.csv")
funnel <- metafor::funnel
forest <- metafor::forest

# Descriptives ------------------------------------------------------------

#'# Descriptives
#'
#'## Publication year
c("from" = min(dat$pubYear, na.rm = T), "to" = max(dat$pubYear, na.rm = T))

#'## Sample sizes
#'
#'### N of effects
dat %>% filter(!is.na(yi)) %>% nrow()

#'### N of studies
length(unique(dat$study)) # overall
dat %>% filter(!is.na(yi)) %$% length(unique(.$study)) # for which ES dat were available

#'###  N of papers
length(unique(dat$paperID))
dat %>% filter(!is.na(yi)) %$% length(unique(.$paperID))

#'### Median N across all the ES eligible for meta-analysis
median(dat$ni, na.rm = T)

#'### Total meta-analytic N
out <- list(NA)
for(i in unique(dat$study)){
  out[i] <- dat %>% filter(study == i) %>% select(ni) %>% max()
}
sum(unlist(out), na.rm = T)

#'### Mean gender ratio (percent female)
out <- list(NA)
for(i in unique(dat$study)){
  out[i] <- dat %>% filter(study == i) %>% select(percFemale) %>% unlist() %>% median()
}
mean(unlist(out), na.rm = T)
sd(unlist(out), na.rm = T)

#'### Weighted mean age of included samples
weighted.mean(dat$meanAge, dat$ni, na.rm = T)

# Meta-analysis -----------------------------------------------------------
#'# Meta-analysis results
#'
#'## H1: Videogames induce a change in players’ explicit (H1a) and implicit attitudes (H1b).
#'
#' Number of iterations run equal to `r nIterationsPcurve` for p-curve and `r nIterations` for all other bias correction functions.
rmaOverall <- data %>% rmaCustom()
(resultsOverall <- data %>% maResults(., rmaObject = rmaOverall, bias = T))

#fit <- RoBMA(d = data$yi, se = sqrt(data$vi), study_names = data$result, seed = 1,
#             chains = 2, iter = 4000)
#summary(fit)

#'### Table 
maResultsTable(resultsOverall, bias = T)

#'### Forest plot
par(mfrow = c(1,2))
data %$% metafor::forest(yi, vi, at = c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2, 2.5),
         xlim = c(-3.5,4.5), alim = c(-3, 3), xlab = "Hedges' g for overall effect",        ### adjust horizontal plot region limits
         order = order(vi),        ### order by size of yi
         slab = NA, annotate = FALSE, ### remove study labels and annotations
         efac = 0,                  ### remove vertical bars at end of CIs
         pch = 19,                  ### changing point symbol to filled circle
         col = "gray40",            ### change color of points/CIs
         psize = 2,                 ### increase point size
         cex.lab=1, cex.axis=.7,   ### increase size of x-axis title/labels
         lty=c("solid","blank"))  ### remove horizontal line at top of plot
addpoly(rmaOverall[[1]], row = -2, mlab = "", cex = 1, annotate = FALSE)

#'### Contour-enhanced funnel plot
funnel(rmaOverall[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline = 0, pch = 20, yaxis = "sei", xlab = "Hedges' g for overall effect", xlim = c(-1.2, 2), steps = 7, digits = c(1,2))

#'### p-curve plot
par(mfrow = c(1,1))
quiet(pcurveMod(metaResultPcurve, effect.estimation = FALSE, plot = TRUE))
title("Overall effect", cex.main = 1)

# Implicit / Explicit attitudes
namesObjects <- c("Implicit", "Explicit")
levels(data$Att_type_exp_imp) <- namesObjects
dataObjects <- list("Implicit" = datImplicit, "Explicit" = datExplicit)
rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

# Results
resultsEffType <- list(NA)
metaResultsPcurve <- list(NA)
for(i in 1:length(rmaObjects)){
  resultsEffType[[i]] <- maResults(data = dataObjects[[i]], rmaObject = rmaObjects[[i]])
  metaResultsPcurve[[i]] <- metaResultPcurve
}

resultsEffType <- setNames(resultsEffType, nm = namesObjects)
metaResultsPcurve <- setNames(metaResultsPcurve, nm = namesObjects)

#'## H1a: Videogames induce a change in players’ explicit attitudes (H1a).
resultsEffType$Explicit

#'## H1b: Videogames induce a change in players’ implicit attitudes (H1b).
resultsEffType$Implicit

#'### Table 
do.call(rbind, resultsEffType %>% map(maResultsTable, bias = T)) %>% noquote() %>% t()

#'### Plots
#'
#' Contour-enhanced funnel plots
par(mfrow = c(1,2))
datExplicit %$% funnel(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g")
title("Explicit attitudes", cex.main = 1)

datImplicit %$% funnel(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g")
title("Implicit attitudes", cex.main = 1)

#'### Forest plots for the explicit and implicit attitudes effects.
datExplicit %$% forest(yi, vi, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
                                                  #xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
                                                  order = order(vi),        ### order by size of yi
                                                  slab = NA, annotate = FALSE, ### remove study labels and annotations
                                                  efac = 0,                  ### remove vertical bars at end of CIs
                                                  pch = 19,                  ### changing point symbol to filled circle
                                                  col = "gray40",            ### change color of points/CIs
                                                  psize = 1.5,                 ### increase point size
                                                  cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                                                  lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Explicit attitudes", cex.main = 1)
addpoly(rmaObjects[[2]][[1]], row = .2, mlab = "", cex = 2, annotate = F)


datImplicit %$% forest(yi, vi, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
                                                  #xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
                                                  order = order(vi),        ### order by size of yi
                                                  slab = NA, annotate = FALSE, ### remove study labels and annotations
                                                  efac = 0,                  ### remove vertical bars at end of CIs
                                                  pch = 19,                  ### changing point symbol to filled circle
                                                  col = "gray40",            ### change color of points/CIs
                                                  psize = 1.5,                 ### increase point size
                                                  cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                                                  lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Implicit attitudes", cex.main = 1)
addpoly(rmaObjects[[1]][[1]], row = .2, mlab = "", cex = 2, annotate = F)

#'### p-curve plots
quiet(pcurveMod(metaResultsPcurve$Explicit, effect.estimation = FALSE, plot = TRUE))
title("Explicit attitudes", cex.main = 1)

quiet(pcurveMod(metaResultsPcurve$Implicit, effect.estimation = FALSE, plot = TRUE))
title("Implicit attitudes", cex.main = 1)

#'#### Comparison of effect types
#'
#' Model without covariates
viMatrix <- data %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaObjectExpImp <- data %$% rma.mv(yi ~ 0 + factor(Att_type_exp_imp), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(RVEmodelExpImp <- data %$% list("k" = table(Att_type_exp_imp),
                                "test" = coef_test(rmaObjectExpImp, vcov = "CR2", test = "z", cluster = study), 
                                "CIs" = conf_int(rmaObjectExpImp, vcov = "CR2", test = "z", cluster = study),
                                "RVE Wald test" = Wald_test(rmaObjectExpImp, constraints = constrain_equal(1:2), vcov = "CR2")))

#' Model with covariates
#' 
#' Controlling for design-related factors that are prognostic w.r.t. the effect sizes (i.e., might vary across moderator categories), namely risk of bias, published status, and age of the sample. 
rmaObjectExpImpCov <- data %$% rma.mv(yi ~ 0 + factor(Att_type_exp_imp) + overallBias + published +  meanAge + inLabAdministration, V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(RVEmodelExpImpCov <- data %$% list("test" = coef_test(rmaObjectExpImpCov, vcov = "CR2", test = "z", cluster = study), 
                                   "CIs" = conf_int(rmaObjectExpImpCov, vcov = "CR2", test = "z", cluster = study),
                                   "RVE Wald test" = Wald_test(rmaObjectExpImpCov, constraints = constrain_equal(1:2), vcov = "CR2")))


#'## H2: Duration of intervention (=“interventionDuration” ) is positively related to the magnitude of explicit (H2a) and implicit (H2b) attitude change. 
#'
#' The below reported meta-regressions are all implemented as a multivariate RVE-based models using the CHE working model (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/).
#' Testing of contrasts is carried out using a robust Wald-type test testing the equality of estimates across levels of the moderator. 
#'
#'### Overall effect moderated by the duration of intervention
rmaModDurInt <- data %$% rma.mv(yi ~ 0 + I(interventionDuration.m./60), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModDurInt <- data %$% list("test" = coef_test(rmaModDurInt, vcov = "CR2", test = "z", cluster = study), 
                               "CIs" = conf_int(rmaModDurInt, vcov = "CR2", test = "z", cluster = study)))

#'### Stability of the overal effect in time
viMatrixPosttestDelay <- data %>% filter(posttestDelay.day. > 0) %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaModPosttestDelay <- data %>% filter(posttestDelay.day. > 0) %$% rma.mv(yi ~ 0 + posttestDelay.day., V = viMatrixPosttestDelay, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModDurInt <- data %>% filter(posttestDelay.day. > 0) %$% list("k" = rmaModPosttestDelay$k.all, "test" = coef_test(rmaModPosttestDelay, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModPosttestDelay, vcov = "CR2", test = "z", cluster = study)))

#'### H2a: Explicit attitudes moderated by the duration of intervention
viMatrixExplicit <- datExplicit %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaModDurIntExplicit <- datExplicit %$% rma.mv(yi ~ 0 + I(interventionDuration.m./60), V = viMatrixExplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModDurIntExplicit <- datExplicit %$% list("test" = coef_test(rmaModDurIntExplicit, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModDurIntExplicit, vcov = "CR2", test = "z", cluster = study)))

#'### H2b: Implicit attitudes moderated by the duration of intervention
viMatrixImplicit <- datImplicit %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaModDurIntImplicit <- datImplicit %$% rma.mv(yi ~ 0 + I(interventionDuration.m./60), V = viMatrixImplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModDurIntImplicit <- datImplicit %$% list("test" = coef_test(rmaModDurIntImplicit, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModDurIntImplicit, vcov = "CR2", test = "z", cluster = study)))

#'## H3: The magnitude of explicit attitude change is smaller in studies using control groups with topic-related activities (Controls_string=”Activity related”) than in studies using control groups with topic-unrelated activities (Controls string = “Activity unrelated” || “Game diff mechanic” || “Game sim mechanic” || “Nothing” || “Combination”).   
data$controlsActivityRelated <- ifelse(data$Controls_string == "Activity related", 1, ifelse(is.na(data$Controls_string), NA, 0))
rmaModActRelExplicit <- data %>% filter(Att_type_exp_imp == "Explicit") %$% rma.mv(yi ~ 0 + factor(controlsActivityRelated), V = viMatrixExplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModActRelExplicit <- data %>% filter(Att_type_exp_imp == "Explicit") %$% list("test" = coef_test(rmaModActRelExplicit, vcov = "CR2", test = "z", cluster = study),
                                                                         "CIs" = conf_int(rmaModActRelExplicit, vcov = "CR2", test = "z", cluster = study),
                                                                         "RVE Wald test" = Wald_test(rmaModActRelExplicit, constraints = constrain_equal(1:2), vcov = "CR2")))

#'## H4: Videogames change player's implicit attitudes only if compared to control groups with unrelated activities towards the measured topic (H8a) (Controls string = “Activity unrelated” || “Game diff mechanic” || “Game sim mechanic” || “Nothing” || “Combination”), but not if compared to control groups using related activities (H8b) (Controls_string= ”Activity related”) towards the measured phenomena.
#'
#' No such implicit studies
# rmaModActRelImplicit <- data %>% filter(Att_type_exp_imp == 0) %$% rma.mv(yi ~ 0 + factor(controlsActivityRelated), V = viMatrixImplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
# (rveModActRelImplicit <- data %>% filter(Att_type_exp_imp == 0) %$% list("test" = coef_test(rmaModActRelImplicit, vcov = "CR2", test = "z", cluster = study), 
#                                                                           "CIs" = conf_int(rmaModActRelImplicit, vcov = "CR2", test = "z", cluster = study),
#                                                                           "RVE Wald test" = Wald_test(rmaModActRelImplicit, constraints = constrain_equal(1:2), vcov = "CR2")))

#'## H5: Implicit attitudes are more affected by games using Stereotype Rehearsal and Meaningful feedback (persuasive mechanic = Stereotype Rehearsal || Reward System || procedural rhetorics || economical model) as persuasive mechanics than by other game mechanics (persuasive mechanic = perspective-taking).
rmaModPersMechImpl <- datImplicit %$% rma.mv(yi ~ 0 + factor(persuasiveMechanicImpl), V = viMatrixImplicit, method = "REML", random = ~ 1|study/result)
(rveModPersMechImpl <- datImplicit %$% list("test" = coef_test(rmaModPersMechImpl, vcov = "CR2", test = "Satterthwaite", cluster = study), 
                                                                       "CIs" = conf_int(rmaModPersMechImpl, vcov = "CR2", test = "Satterthwaite", cluster = study),
                                                                       "RVE Wald test" = Wald_test(rmaModPersMechImpl, constraints = constrain_equal(1:2), vcov = "CR2")))

#'## H6: Explicit attitudes are more affected by games using Perspective-taking and Meaningful feedback (persuasive mechanic = perspective-taking || procedural rhetorics || economical model || reward system as persuasive mechanics than by other game mechanics (persuasive mechanic = Stereotype Rehearsal)).
#'
rmaModPersMechExpl <- datExplicit %$% rma.mv(yi ~ 0 + factor(persuasiveMechanicExpl), V = viMatrixExplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModPersMechExpl <- datExplicit %$% list("test" = coef_test(rmaModPersMechExpl, vcov = "CR2", test = "Satterthwaite", cluster = study), 
                                                                       "CIs" = conf_int(rmaModPersMechExpl, vcov = "CR2", test = "Satterthwaite", cluster = study),
                                                                       "RVE Wald test" = Wald_test(rmaModPersMechExpl, constraints = constrain_equal(1:2), vcov = "CR2")))

#'#### H6: Exploratory analysis, comparing only Perspective-taking (= 1) vs Stereotype rehearsal (= 0).
viMatrixExpl2 <- datExplicit %>% filter(persuasive_mechanic %in% c("Perspective-taking", "Stereotype rehearsal")) %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaModPersMechExpl2 <- datExplicit %>% filter(persuasive_mechanic %in% c("Perspective-taking", "Stereotype rehearsal")) %$% rma.mv(yi ~ 0 + factor(persuasiveMechanicExpl), V = viMatrixExpl2, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModPersMechExpl2 <- datExplicit %>% filter(persuasive_mechanic %in% c("Perspective-taking", "Stereotype rehearsal")) %$% list("test" = coef_test(rmaModPersMechExpl2, vcov = "CR2", test = "z", cluster = study), 
                                            "CIs" = conf_int(rmaModPersMechExpl2, vcov = "CR2", test = "z", cluster = study),
                                            "RVE Wald test" = Wald_test(rmaModPersMechExpl2, constraints = constrain_equal(1:2), vcov = "CR2")))


#'## H7: Action games (gamegenre =”action game” || “action violent game”) have a larger effect on implicit attitudes than on explicit attitudes.
viMatrixAction <- data %>% filter(gamegenre %in% c("action game", "action violent game")) %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaModActionGames <- data %>% filter(gamegenre %in% c("action game", "action violent game")) %$% rma.mv(yi ~ 0 + factor(Att_type_exp_imp), V = viMatrixAction, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModActionGames <- data %>% filter(gamegenre %in% c("action game", "action violent game")) %$% list("test" = coef_test(rmaModActionGames, vcov = "CR2", test = "z", cluster = study), 
                                                                                                       "CIs" = conf_int(rmaModActionGames, vcov = "CR2", test = "z", cluster = study),
                                                                                                       "RVE Wald test" = Wald_test(rmaModActionGames, constraints = constrain_equal(1:2), vcov = "CR2")))

#'## H8: Non-action games (gamegenre=”simulation/strategy” || “adventure”) have a larger effect on explicit attitudes than action games. 
rmaModNonactExpl <- datExplicit %$% rma.mv(yi ~ 0 + factor(gamegenreBinary), V = viMatrixExplicit, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModNonactExpl <- datExplicit %$% list("test" = coef_test(rmaModNonactExpl, vcov = "CR2", test = "z", cluster = study), 
                                                                     "CIs" = conf_int(rmaModNonactExpl, vcov = "CR2", test = "z", cluster = study),
                                                                     "RVE Wald test" = Wald_test(rmaModNonactExpl, constraints = constrain_equal(1:2), vcov = "CR2")))

#'##	Exploratory analyses
#'
#'### Overall effect moderated by age
rmaModAge <- data %$% rma.mv(yi ~ meanAge, V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModAge <- data %$% list("test" = coef_test(rmaModAge, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModAge, vcov = "CR2", test = "z", cluster = study)))

#'### Overall effect moderated by gender ratio
rmaModGender <- data %$% rma.mv(yi ~ percFemale, V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModGender <- data %$% list("test" = coef_test(rmaModGender, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModGender, vcov = "CR2", test = "z", cluster = study)))

#'## Methodological moderators
#'
#'### In-lab administration
rmaModInLab <- data %$% rma.mv(yi ~ 0 + factor(inLabAdministrationBinary), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModInLab <- data %$% list("test" = coef_test(rmaModInLab, vcov = "CR2", test = "z", cluster = study), 
                              "CIs" = conf_int(rmaModInLab, vcov = "CR2", test = "z", cluster = study),
                              "RVE Wald test" = Wald_test(rmaModInLab, constraints = constrain_equal(1:2), vcov = "CR2")))

#'### Commercial status
rmaModCommercial <- data %$% rma.mv(yi ~ 0 + factor(Commercial), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModCommercial <- data %$% list("test" = coef_test(rmaModCommercial, vcov = "CR2", test = "z", cluster = study), 
                                  "CIs" = conf_int(rmaModCommercial, vcov = "CR2", test = "z", cluster = study),
                                  "RVE Wald test" = Wald_test(rmaModCommercial, constraints = constrain_equal(1:2), vcov = "CR2")))


#'### Range restriction
rmaModRestricted <- data %$% rma.mv(yi ~ 0 + factor(sampleRestricted), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModRestricted <- data %$% list("k" = table(data$sampleRestricted),
                                   "test" = coef_test(rmaModRestricted, vcov = "CR2", test = "z", cluster = study), 
                                   "CIs" = conf_int(rmaModRestricted, vcov = "CR2", test = "z", cluster = study),
                                   "RVE Wald test" = Wald_test(rmaModRestricted, constraints = constrain_equal(1:2), vcov = "CR2")))

#'#### F-test of equality of variances
#'
#' Mean vi for restricted designs
(viRestricted <- data %>% filter(sampleRestricted == 1) %$% mean(vi, na.rm = T))
dfRestricted <- data %>% filter(sampleRestricted == 1, !is.na(vi)) %>% nrow() - 1
#' Mean vi for non-restricted designs
(viNonRestricted <- data %>% filter(sampleRestricted == 0) %$% mean(vi, na.rm = T))
dfNonRestricted <- data %>% filter(sampleRestricted == 0, !is.na(vi)) %>% nrow() - 1
#' F-statistics
viNonRestricted/viRestricted
#' F-test p-value
2 * (1 - pf(viNonRestricted/viRestricted, df1 = dfNonRestricted, df2 = dfRestricted, lower.tail = FALSE))

#'### Year of Publication
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEpubYear <- summary(lmer(scale(sqrt(vi)) ~ scale(journalH5) + scale(pubYear) + (1|study), data = data, REML = T))$coefficients)
#' Comment: all the variables were centered for easier interpretation of model coefficients. See the negative beta for Publication Year. The more recent a publication, the lower the variance (better precision), controlling for H5.
#'

#'#### Scatterplot year <-> precision
#'
#' Size of the points indicate the H5 index (the bigger the higher) of the journal that the ES is published in.
yearPrecisionScatter <- data %>%  ggplot(aes(pubYear, sqrt(vi))) + 
  geom_point(aes(size = journalH5), alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  scale_x_continuous() +
  xlab("Year of publication") +
  ylab("Imprecision") +
  theme_bw() +
  theme(legend.position = "none")
yearPrecisionScatter

#'### Citations
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEcitations <- summary(lmer(scale(sqrt(vi)) ~ scale(pubYear) + scale(journalH5) + scale(citations) + (1|study), data = data, REML = T))$coefficients)

#'#### Scatterplot precision <-> citations
#'
#' The relationship between precision (sqrt of variance) and number of citations (log).
citationsPrecisionScatter <- data %>% ggplot(aes(log(citations + 1), sqrt(vi))) + 
  geom_point(alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  xlim(0, 5) +
  xlab("Citations (log)") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")
citationsPrecisionScatter

#'### H5 index
#'
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEh5 <- summary(lmer(scale(sqrt(vi)) ~ scale(journalH5) + (1|study), data = data, REML = T))$coefficients)

#'#### Scatterplot precision <-> journal H5
#'
#' The relationship between precision (sqrt of variance) and H5 index of the journal.
h5PrecisionScatter <- data %>% ggplot(aes(journalH5, sqrt(vi))) + 
  geom_point(alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  xlim(0, 130) +
  xlab("H5 index") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")
h5PrecisionScatter

#'### Decline effect
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study.
(declineEff <- summary(lmer(scale(abs(as.numeric(rmaOverall[[1]][1]) - yi)) ~ scale(sqrt(vi)) + scale(pubYear) + (1|study), data = data, REML = T))$coefficients)

#'### Citation bias
#' Do more highly-cited studies report larger effect sizes?
data <- data %>% mutate(rConv = yi/sqrt(yi^2 + 4))
(LMEcitationsYi <- summary(lmer(rConv ~ scale(pubYear) + scale(journalH5) + scale(citations) + (1|study), data = data, REML = T))$coefficients)
(LMEcitationsYiVi <- summary(lmer(rConv ~ scale(pubYear) + scale(journalH5) + scale(citations) + scale(vi) + (1|study), data = data, REML = T))$coefficients)

#'## Sensitivity analyses
#'
#'### Risk of bias assessment
rob_summary(rob2, tool = "ROB2", overall = T)

#' Moderator analysis for risk of bias
rmaModRoB <- data %$% rma.mv(yi ~ 0 + factor(overallBias), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModRoB <- data %$% list("k" = table(data$overallBias), 
                            "test" = coef_test(rmaModRoB, vcov = "CR2", test = "z", cluster = study), 
                            "CIs" = conf_int(rmaModRoB, vcov = "CR2", test = "z", cluster = study),
                            "RVE Wald test" = Wald_test(rmaModRoB, constraints = constrain_equal(1:3), vcov = "CR2")))

#'### Excluding effects based on inconsistent means or SDs
#' 
#' GRIM test, testing the inconsistencies in means
table(data$inconsistenciesCountGRIM, useNA = "always")
#' GRIMMER test, testing the inconsistencies in SD (and means)
table(data$inconsistenciesCountGRIMMER, useNA = "always")
rmaModIncon <- data %$% rma.mv(yi ~ 0 + factor(as.logical(inconsistenciesCountGRIMMER)), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModIncon <- data %$% list("Count of GRIM/GRIMMER inconsistencies" = table(as.logical(data$inconsistenciesCountGRIMMER)),
                             "test" = coef_test(rmaModIncon, vcov = "CR2", test = "z", cluster = study), 
                             "CIs" = conf_int(rmaModIncon, vcov = "CR2", test = "z", cluster = study),
                             "RVE Wald test" = Wald_test(rmaModIncon, constraints = constrain_equal(1:2), vcov = "CR2")))

#'### Excluding effects reported in studies with possible conflicts of interests
rmaModCoI <- data %$% rma.mv(yi ~ 0 + factor(as.logical(data$possibleCOI)), V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(rveModCoI <- data %$% list("Count of effects reported in studies with CoI" = table(data$possibleCOI),
                              "test" = coef_test(rmaModCoI, vcov = "CR2", test = "z", cluster = study), 
                              "CIs" = conf_int(rmaModCoI, vcov = "CR2", test = "z", cluster = study),
                              "RVE Wald test" = Wald_test(rmaModCoI, constraints = constrain_equal(1:2), vcov = "CR2")))

#'### Numerical inconsistencies in reported p-values
nrow(statcheck) # how many results were analyzed
length(unique(statcheck$Source)) # how many papers reported results in APA format
prop.table(table(statcheck$Error)) # how many statcheck errors
table(statcheck$DecisionError)[2]/table(statcheck$Error)[2] # how many statcheck errors affected the decision
statcheck %>% filter(Error == TRUE) %>% select(Source) %>% unique() %>% nrow()/length(unique(statcheck$Source)) # How many papers contained statcheck errors

# Record session info
sessionInfo()
