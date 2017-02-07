suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("caTools"))

## This file performs fitting of the Beat AML inhibition data to DSS curves
## It is based heavily on the file curve_Fitting.R provided by
## Swapnil.

## References to Ritz et al. below are to
## Dose-response analysis using R
## Ritz, Baty, Streibig, and Gerhard (2015) PLoS ONE
## PMID:26717316

## References to Yadav et al. below are to
## Quantitative scoring of differential drug sensitivity for individually optimized anticancer therapies.
## Yadav, Pemovska, Szwajda, Kulesskiy, Kontro, Karjalainen, Majumder, Malani, Murumägi, Knowles, Porkka, Heckman, Kallioniemi, Wennerberg, and Aittokallio (2014) Scientific Reports
## PMID: 24898935

## dose: drug dose concentration (in natural, not log10, units)
## viability: cell viability under that dose (in percent--i.e., 100 = 100%)
## NB: viability = 100 - inhibition
## We will fit a log-logistic function, which is scale invariant through the parameter b,
## as noted in Ritz et al.
## Duplicated response (i.e., viability) will be randomly jittered.  Hence, you may want
## to set.seed before calling this function.
calc.ic50 <- function(dose, viability, file) {
  
  ## Jitter any duplicated data to avoid singularities during fitting
  dupes <- duplicated(viability)
  if(any(dupes)) {
    viability[dupes] <- viability[dupes] + rnorm(n = length(which(dupes)), mean = 0, sd = 0.0001)
  }
  
  ## Define inhibition in terms of viability
  inhibition <- 100 - viability
  
  ## Combine the data and sort by dose
  mat_tbl <- data.frame(logconc = log10(dose), dose = dose, viability = viability, inhibition = inhibition)
  mat_tbl <- mat_tbl[order(mat_tbl$logconc), ]
  
  ## Calculate the IC50.
  
  ## Use drc to get the initial parameters.  
  ## FIMM function:
  ##   y = d + ( a - d ) / ( 1 + ( 10^( b * ( c - x ) ) ) )
  ## Asymptotics for positive b: 
  ##              x -> infty --> y -> a
  ##              x -> -infty --> y -> d
  ## In an inhibition setting, FIMM sets b = 1, d (min asymptote) = 0, and determines a and c via the linear regression self starter
  ## drc function:
  ##   y = c' + ( d' - c' ) / ( 1 + ( exp( b' * ( log(x) - log(e') ) ) )
  ## here c' is the lower asymptote, d' is the upper asymptote, and b' is the slope (see Ritz et al). 
  ## Note that the drc slope b' has the opposite sign as the FIMM slope b.
  ## Hence, to ensure that y is a decreasing function (e.g., suitable for viability response where viability
  ## decreases with concentration) and c' is the lower asymptote, we need b' to be positive.
  ## Asymptotics:
  ##    x -> infy, b' < 0 --> y -> d'
  ##    x -> infy, b' > 0 --> y -> c'
  ##    x ->    0, b' < 0 --> y -> c'
  ##    x ->    0, b' > 0 --> y -> d'
  ## Hence, for an increasing function (x increases -> y increases), we will have b' < 0
  ## so that d' is the upper asymptote.  For a decreasing function (x increases -> y decreases),
  ## we will have b' > 0, so that c' is the lower asymptote.
  ## And this is stated by Ritz et al. here:
  ## 'the slope paramter b [what I'm called b'] should be negative for an increasing dose-response relationship
  ## to ensure the interpretation that c is the lower asymptote and d the upper asymptote.'
  
  ## For inhibition data (which is what we fit with drc), the function y will be increasing.  
  ## Hence, b' should be negative.
  ## In any case, c' and d' are the min and max, respectively.
  ## Parameters are in the order names = c("b'", "c'", "d'", "e'")
  estimate_param <- tryCatch({
    drm(inhibition ~ logconc, data = mat_tbl, 
        fct = LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
        logDose = 10, 
        control = drmc(errorm = FALSE))
    }, 
    warning = function(w) {
      drm(inhibition ~ logconc, data = mat_tbl, 
          fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
          logDose=10)
      },
    error = function(e) {
      drm(inhibition ~ logconc, data = mat_tbl, 
          fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
          logDose=10)
      })
  
  ## Extract and name the coefficients
  coef_estim <- coef(estimate_param)
  ## b is slope; c is min; d is max, and e is IC50
  names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
  
  ## DSS inverts the sign of the slope [see above and compare Yadav et al. supp Eq(1) to Ritz et al. (2)] 
  coef_estim["SLOPE"] <- coef_estim["SLOPE"] * -1 
  
  ## Handle a bunch of corner cases.  This is taken verbatim from Swapnil's curve_Fitting.R.  I don't necessarily follow the logic.
  ## BEGIN CORNER CASES.
  # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  
  # NB: switching IC50 to log10 scale here
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  
  ## Brian: check this--I believe that dss will not calculate dss if IC50 is set to max log concentration, which is why
  ## this is done here.  But, this is not pass directly to DSS, but to nls.
  # similar to previous step but now compare log10(IC50) with log(min. conc.).
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
  # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
  ## Brian: again, may like above?
  coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
  #(Trying to fix curves that need outlier kickout)
  coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
  #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.  
  min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
  min_lower <- ifelse(min_lower >= 100,99,min_lower)
  #similar to previous step but for MAX
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
  #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
  max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
  max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
  max_lower <- ifelse(max_lower < 0,0,max_lower)
  max_lower <- ifelse(max_lower > 100,100,max_lower)
  #(Fix upper maximum for negative slopes)
  run_avg <- runmean(mat_tbl$inhibition, 10)
  max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
  max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
  max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
  max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
  max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
  # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen. 
  mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
  if(mean_inh_last < 60) {
    if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
    else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
  if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
  #add a bit of positive noise to MAX if it is the same as MIN. 
  if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + rnorm(n=1, mean = 0, sd = 0.001)
  
  ## END CORNER CASES.
  
  #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
  nls_result_ic50 <- tryCatch({
    nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", 
        start=list(SLOPE=1,MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
        lower=list(SLOPE=0,MIN=0,MAX=max_lower,IC50=min(mat_tbl$logconc)),
        upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
        control=list(warnOnly=T,minFactor = 1/2048))
    }, 
    error = function(e) {
      minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                        start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                        lower=c(SLOPE=0.5, MIN=min_lower,MAX=100,  IC50=min(mat_tbl$logconc)),
                        upper=c(SLOPE=2.5, MIN=coef_estim["MIN"],MAX=100, IC50=max(mat_tbl$logconc)))
      }
  )
  
  #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
  if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
  {
    if(mean_inh_last > 60)
      coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
    nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",
                           start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                           lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),
                           upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
                           control=list(warnOnly=T,minFactor = 1/2048))
  }
  max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
  
  #############################  
  #############   Final modification & STD error
  
  #prepare final data and convert IC50 back from log scale (inverse)
  coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  #(Fix ic50 for curves in wrong direction)
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
  #(Fix based on MAX)
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
  coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
  #(Fix over sensitive drugs)
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
  
  #Calculate the standard error scores
  sumIC50 = summary(nls_result_ic50); 
  ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]
  ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1)  
  

  # plot IC50
  png(file)
  x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
  y <- predict(nls_result_ic50, data.frame(logconc=x))
  icpl <- ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
    geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
    geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) 
  ## icpl <- icpl + ggtitle(paste0(drug_name," (dss:",round(IC50_dataframe$DSS,1),")\n"))
  icpl <- icpl + theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
    geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background =element_rect(fill = "transparent",colour = NA))

  print(icpl)
  d <- dev.off()
  return(coef_ic50["IC50"])
  
}

library(foreach)

suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

use.synapse <- TRUE

seed <- 1234
set.seed(seed)

if(use.synapse) {
  library(synapseClient)
  synapseLogin()
}

# Pull in the (raw) inhibitor results
raw.inhibitor.data.file <- "inhibitor_data_points_2016_10_10.txt"
if(!file.exists(raw.inhibitor.data.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7440449", downloadFile = TRUE, downloadLocation = ".")
    raw.inhibitor.data.file <- getFileLocation(obj)
  }
}

raw.inhibitor.data <- read.table(raw.inhibitor.data.file, sep="\t", header=TRUE, as.is=TRUE)

## Sort the table so that all entries for a given patient_id and inhibitor are together.  This does
## not affect correctness.  However, if we only do a subset of the table, it will allow us to compare
## what should be replicates.
o <- order(raw.inhibitor.data$patient_id, raw.inhibitor.data$inhibitor, raw.inhibitor.data$replicant, raw.inhibitor.data$lab_id, raw.inhibitor.data$time_of_read)
raw.inhibitor.data <- raw.inhibitor.data[o, ]
flag <- raw.inhibitor.data$inhibitor %in% c("Dasatinib", "JAK Inhibitor I", "AG490", "AST-487", "Erlotinib", "Gefitinib")
raw.inhibitor.data <- raw.inhibitor.data[flag, ]

patient.drug.tbl <- unique(raw.inhibitor.data[,c("patient_id","inhibitor", "replicant", "lab_id", "time_of_read")])

n.iters <- nrow(patient.drug.tbl)
cat(paste0("Total iterations: ", n.iters, "\n"))
## n.iters <- 1
## n.iters <- 20
tbl <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
  # for(i in 1:nrow(patient.drug.tbl)) {
  if( ( i %% 10 ) == 0 ) { cat(paste0("Iter ", i, "\n")) }
  patient.id <- patient.drug.tbl$patient_id[i]
  inhibitor <- patient.drug.tbl$inhibitor[i]
  replicant <- patient.drug.tbl$replicant[i]
  lab_id <- patient.drug.tbl$lab_id[i]
  time_of_read <- patient.drug.tbl$time_of_read[i]
  
  file <- paste0(patient.id, inhibitor, replicant, lab_id, time_of_read)
  file <- paste0(file, ".png")
  
  flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor) & (raw.inhibitor.data$replicant == replicant) & (raw.inhibitor.data$lab_id == lab_id) & (raw.inhibitor.data$time_of_read == time_of_read)
  df <- raw.inhibitor.data[flag,,drop=F]
  ic50 <- calc.ic50(df$well_concentration, df$normalized_viability, file)
  row <- c(patient.id, inhibitor, replicant, lab_id, time_of_read, ic50)
  row
}

rownames(tbl) <- NULL
colnames(tbl) <- c("patient.id", "inhibitor", "replicant", "lab_id", "time_of_read", "ic50")

save.image(file=".RData.tbl")

patient.drug.tbl <- unique(raw.inhibitor.data[,c("patient_id","inhibitor")])

n.iters <- nrow(patient.drug.tbl)
cat(paste0("Total iterations: ", n.iters, "\n"))
## n.iters <- 1
## n.iters <- 20
tbl2 <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
  # for(i in 1:nrow(patient.drug.tbl)) {
  if( ( i %% 10 ) == 0 ) { cat(paste0("Iter ", i, "\n")) }
  patient.id <- patient.drug.tbl$patient_id[i]
  inhibitor <- patient.drug.tbl$inhibitor[i]

  file <- paste(patient.id, inhibitor, sep="-")
  file <- paste(file, ".png", sep="")
  
  flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor)
  df <- raw.inhibitor.data[flag,,drop=F]
  ic50 <- calc.ic50(df$well_concentration, df$normalized_viability, file)
  row <- c(patient.id, inhibitor, ic50)
  row
}


stop("stop")

tbl <- as.data.frame(tbl)
tbl$ic50 <- as.numeric(as.character(tbl$ic50))
flag <- tbl$ic50 < 1
tbl <- tbl[flag, ]

## Do this on a per-patient, per drug basis
## If any replicants, lab ids, or time_of_reads are significant, exclude that drug / patient combo
## l <- lm("ic50 ~ patient.id + replicant + lab_id + time_of_read", data=tbl)
## So fit a model like
## l <- lm("ic50 ~ replicant + lab_id + time_of_read", data = tbl[patient.drug.subset, ])
## 

## How about this: find mean within patient.id + drug and exclude any outside of mean +/- 2 * sd
## This will need to be CV.

## flag <- (tbl$patient.id == "2344") & (tbl$inhibitor=="Dasatinib")
## tbl[flag,]

cvs <- aggregate(tbl$ic50, by=list(tbl$patient.id, tbl$inhibitor), FUN = function(x) sd(x)/mean(x))
cvs <- cvs[!is.na(cvs$x),]
plot(density(cvs$x))