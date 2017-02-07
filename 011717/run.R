suppressPackageStartupMessages(library("nplr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))

suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("caTools"))

suppressPackageStartupMessages(library("scales")) ## for scientific


## This file performs fitting of the Beat AML inhibition data to DSS curves
## It is based heavily on the file curve_Fitting.R provided by
## Swapnil.

output.path <- "output"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

## References to Ritz et al. below are to
## Dose-response analysis using R
## Ritz, Baty, Streibig, and Gerhard (2015) PLoS ONE
## PMID:26717316

## References to Yadav et al. below are to
## Quantitative scoring of differential drug sensitivity for individually optimized anticancer therapies.
## Yadav, Pemovska, Szwajda, Kulesskiy, Kontro, Karjalainen, Majumder, Malani, Murumägi, Knowles, Porkka, Heckman, Kallioniemi, Wennerberg, and Aittokallio (2014) Scientific Reports
## PMID: 24898935

use.nplr <- TRUE

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## This function seems to be invariant of concn_scale, empirically and by looking at the equation.

## y: activity threshold (t) for percent inhibition. This must be greater than 0. Default is 10.
## DSS.type:  type of DSS: 1, 2 or 3. This value corresponds to DSS1, DSS2, DSS3. Default is 3.
## concn_scale:  concentration scale of drugs tested. Set the scale such a way that IC50, Minimum concentration tested and Maximum concentration tested converts to molar (M) unit. Default is \code{NULL}, but, must set if concentration is not in \code{log10} scale. For example, if the drug concentration is in nano molar (nM), then \code{concn.scale = 1e-9}.
##  The activity threshold for percent inhibition is the lower limit in dose response curve (t), the area above which is used for DSS calculation.
  
## BSW add min response (min) as a parameter
dss <- function(ic50,slope,max,min.conc.tested,max.conc.tested,y.arg=10,DSS.type=2,concn_scale=1e-9, min.response=0){
  #rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
  
  ## BSW added:
  ## If the min response is 0, then make the activity threshold = y.arg (e.g,. 10% by default)
  ## If the min response is not 0, then make the activity threshold 10% (by default) of the range (max - min.response) _above_ the min.response
  y <- as.numeric(unname(y.arg))
  if(max < as.numeric(unname(min.response))) { 
    dss <- 0
    return(dss)
  }
  if(as.numeric(unname(min.response)) != 0) {
    range <- max - as.numeric(unname(min.response))
    y <- min(max, as.numeric(unname(min.response)) + (((as.numeric(unname(y.arg)))/100) * range))
  }
  
  a=as.numeric(unname(max))
  
  b=as.numeric(unname(slope))
  d=as.numeric(unname(min.response)) # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)  
  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  ### BSW added
  else if(isTRUE(ic50<=min.conc.tested)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
        else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      
      # This is a logistic function used in Dotmatics.com
      # y = d+(a-d)/(1+10^(b*(c-x)))
      #inverse function
      # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){
        norm_area<-((int_y/total_area)*100)#DSS1
      }
      if(DSS.type==2){
        #       if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){
        #       if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
      }
      if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0} 
  } 
  return (dss)
}

## The inverse of convertToProp.
## i.e., the inverse of p = ( y - T0 ) / ( Ctrl - T0 )
invertProp <- function(p, T0, Ctrl) {
  ( p * ( Ctrl- T0 ) ) + T0
}

## See /home/ubuntu/aml-dss/111816/run.R
nplr.calc.ic50 <- function(dose, viability, file, patient.id, inhibitor, apply.constraints=FALSE, return.fit=TRUE, plot.fit=FALSE) {
  
  ## print(file)
  ##  begin.file <- paste0(file, ".begin")
  ##  end.file <- paste0(file, ".end")
  ##  save(dose, file="dose.Rd")
  ##  save(viability, file="viability.Rd")
  ## Jitter any duplicated data to avoid singularities during fitting
  ## NB: no--don't do this for nplr.
  ##  dupes <- duplicated(viability)
  freq.table <- as.data.frame(table(dose))
  has.replicates <- all(freq.table$Freq >= 2)
  ##  dose.dupes <- duplicated(dose)
  ##  if(any(dupes)) {
  ##    viability[dupes] <- viability[dupes] + rnorm(n = length(which(dupes)), mean = 0, sd = 0.0001)
  ##  }
  
  ## Define inhibition in terms of viability
  inhibition <- 100 - viability
  
  if(length(unique(viability)) == 1) { return(list(tbl=NULL,res=NULL)) }
  
  
  ## Combine the data and sort by dose
  mat_tbl <- data.frame(logconc = log10(dose), dose = dose, viability = viability, inhibition = inhibition)
  mat_tbl <- mat_tbl[order(mat_tbl$logconc), ]
  
  method <- "res"
  if(has.replicates) {
    method <- "sdw"
  }
  ## Calculate the IC50.
  T0 <- min(mat_tbl$inhibition, na.rm = TRUE)
  Ctrl <- max(mat_tbl$inhibition, na.rm = TRUE)
  
  x <- mat_tbl$logconc
  yProp <- convertToProp(mat_tbl$inhibition, T0=T0, Ctrl=Ctrl)
  
  if(any(is.na(yProp)) | any(is.nan(yProp))) { return(list(tbl=NULL,res=NULL)) }
  
  
  ## Constraint vector is over c(bottom, top, xmid=IC50, scal, s)
  lower.limits <- c(-Inf, -Inf, -Inf, -Inf, -Inf)
  upper.limits <- c(Inf, Inf, Inf, Inf, Inf)
  if(apply.constraints == TRUE) {
    ## Constrain xmid (i.e., IC50) to be between the min and max concentrations
    lower.limits[3] <- min(x)
    upper.limits[3] <- max(x)
    
    ## Constrain bottom to be between [min response - 10% of range, max response]
    if(max(yProp) != min(yProp)) {
      lower.limits[1] <- min(yProp) - 0.1 * abs(max(yProp)-min(yProp))
      upper.limits[1] <- max(yProp)
    }
    
    ## Constrain top to be bewteen [min response, max response + 10% of range]
    if(max(yProp) != min(yProp)) {
      lower.limits[2] <- min(yProp)
      upper.limits[2] <- max(yProp) + 0.1 * abs(max(yProp)-min(yProp))
    }
  }
  
  
  result_ic50 <- tryCatch({
    ## 4-parameter model is y = B + (T - B)/[1 + 10^(b*(xmid - x))]
    ## useLog = FALSE --> x is already log10-transformed
    ## T0=0 prevents scaling the bottom of the range to 0
    ## cat("Trying 4 parameter model\n")
    model <- nplr(x = x, y = yProp, npars=4, useLog = FALSE, method = method, apply.constraints = apply.constraints, lower.limits = lower.limits, upper.limits = upper.limits)
    model
  },
  warning = function(w) {
    tryCatch({
      cat("Trying 5 param model\n")
      model <- nplr(x = x, y = yProp, npars=5, useLog = FALSE, method = method, apply.constraints = apply.constraints, lower.limits = lower.limits, upper.limits = upper.limits)
      model
    },
    warning = function(w) { 
      ##      cat(paste0("Returning null for ", file, "\n"))
      return(NULL) 
    },
    error = function(e) { 
      ##      cat(paste0("Returning null for ", file, "\n"))
      return(NULL) 
    })
  },
  error = function(e) {
    tryCatch({
      cat("Trying 5 param model\n")
      model <- nplr(x = x, y = yProp, npars=5, useLog = FALSE, method = method, apply.constraints = apply.constraints, lower.limits = lower.limits, upper.limits = upper.limits)
      model
    },
    warning = function(w) { 
      ##      cat(paste0("Returning null for ", file, "\n"))
      return(NULL) 
    },
    error = function(e) {
      ##      cat(paste0("Returning null for ", file, "\n"))
      return(NULL) 
    })
  }
  )
  
  if(is.null(result_ic50)) { return(list(tbl=NULL,res=NULL)) }
  
  max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
  
  
  #############################  
  #############   Final modification & STD error
  
  ## prepare final data and convert IC50 back from log scale (inverse)
  coef_ic50 <- c(result_ic50@pars$xmid, result_ic50@pars$scal, result_ic50@pars$top, result_ic50@pars$bottom) 
  names(coef_ic50) <- c("IC50", "SLOPE","MAX","MIN"); 
  ## nplr returns IC50 on the log10 scale; whereas dcr returns on the natural scale.  This conversion is also
  ## done by the DSS code based on nls.
  ## Convert IC50 to natural scale.
  coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  
  ## Project the min and max from [0,1] to [0,100]
  coef_ic50["MAX"] <- invertProp(coef_ic50["MAX"], T0, Ctrl)
  coef_ic50["MIN"] <- invertProp(coef_ic50["MIN"], T0, Ctrl)
  if(!is.infinite(lower.limits[1])) { lower.limits[1] <- invertProp(lower.limits[1], T0, Ctrl) }
  if(!is.infinite(upper.limits[1])) { upper.limits[1] <- invertProp(upper.limits[1], T0, Ctrl) }
  if(!is.infinite(lower.limits[2])) { lower.limits[2] <- invertProp(lower.limits[2], T0, Ctrl) }
  if(!is.infinite(upper.limits[2])) { upper.limits[2] <- invertProp(upper.limits[2], T0, Ctrl) }
  
  ## cat("lower limits:\n")
  ## print(lower.limits)
  ## cat("upper limits:\n")
  ## print(upper.limits)
  
  ## NB: nplr (see the first eqn in the R package tutorial by Commo and Bot and/or the 5-parameter logistic regression
  ## eqn under Details in ??nplr) and DSS [see Yadav et al. supp Eq(1)] assume the same form of the logistic
  ## regression equation.  In particular, the exponent in the demoninator is off the form IC50 - x.  Hence, we do
  ## not need to invert the sign of the slope.  We _do_ need to do this when coupling drc with DSS, since dcr
  ## uses the opposite convention in which the exponent is of the form x - IC50 [see eqn (2) of Ritz et al.].
  ## More explicitly, the code below that precedes the calculation of DSS uses nls to fit the same eqn fit by nplr:
  ##     inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc))))
  
  ## coef_ic50["SLOPE"] <- coef_ic50["SLOPE"]*-1 
  
  #  if(method == "sdw") {
  ##  if(abs(log10(coef_ic50["IC50"])) > 10) {
  if(FALSE) {
    if(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50)) {
      base <- paste0(output.path, "/", file, "-inhibition")
      save(coef_ic50, file=paste0(base, ".coef.ic50.Rd"))
      save(mat_tbl, file=paste0(base, ".mat.tbl.Rd"))
      save(result_ic50, file=paste0(base, ".result.ic50.Rd"))
    }
    
    if(abs(log10(coef_ic50["IC50"])) > 10) {
      base <- paste0(output.path, "/", file, "-large")
      save(coef_ic50, file=paste0(base, ".coef.ic50.Rd"))
      save(mat_tbl, file=paste0(base, ".mat.tbl.Rd"))
      save(result_ic50, file=paste0(base, ".result.ic50.Rd"))
    }
  }
  ## Check various error conditions.  dss will return 0 if max <= max_signal.  Hence,
  ## setting max to max_signal will force dss to be zero.
  
  #(Fix ic50 for curves in wrong direction)
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
  #(Fix based on MAX)
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
  ## BSW commented out the condition below
  ##coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
  ## BSW added the condition below
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<coef_ic50["MIN"],max_signal,coef_ic50["IC50"])
  ## BSW commented out the condition below
  ## coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
  #(Fix over sensitive drugs)
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
  
  #Calculate the standard error scores
  estim <- getEstimates(result_ic50, 0.5, conf.level=erf(1/sqrt(2)))
  ic50std_Error <- abs(estim[,4]-estim[,2])/2
  
  residuals <- invertProp(getFitValues(result_ic50), T0, Ctrl) - invertProp(getY(result_ic50), T0, Ctrl)
  ic50std_resid <- round(sqrt(sum((residuals)^2)/(length(residuals)-1)),1)  
  
  ## Get the goodness of fit and the weighted goodness of fit
  gof <- getGoodness(result_ic50)$gof
  wgof <- getGoodness(result_ic50)$wgof
  
  ## The "inhibitor data points" tab of the file "Functional and Clinical Dataset Readme 2-5-16.xlsx"
  ## says of the "well_concentration" column that it is 
  ## "Concentration of drug in a given well, typically in uM."
  ## Hence, concn_scale should be 10^-6
  
  ############################# 
  #############    DSS
  
  dss_score <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, concn_scale=1e-6, min.response=coef_ic50["MIN"]) 
  tmp = coef_ic50;
  coef_ic50 <- c(tmp,Max.Conc.tested=max_signal,Min.Conc.tested=min_signal,DSS=dss_score,IC50_std_error=ic50std_Error,S_est = ic50std_resid)
  
  #inhib. and viab. values in order of dose growth
  ## inhibition_values <- t(matrix(mat_tbl[,"inhibition"],dimnames=list(paste0(rep("D", length(mat_tbl[,"inhibition"])), 1:length(mat_tbl[,"inhibition"])))))
  
  
  #dataframe for IC50
  ## IC50_dataframe <- data.frame(ID=file,PATIENT_ID=patient.id,DRUG_NAME=inhibitor,t(as.matrix(coef_ic50)),inhibition_values)
  IC50_dataframe <- data.frame(ID=file,PATIENT_ID=patient.id,DRUG_NAME=inhibitor,t(as.matrix(coef_ic50)), GOF=gof, WGOF=wgof, MIN.lb=lower.limits[1], MIN.ub=upper.limits[1], MAX.lb=lower.limits[2], MAX.ub=upper.limits[2], IC50.lb=lower.limits[3], IC50.ub=upper.limits[3], SLOPE.lb=lower.limits[4], SLOPE.ub=upper.limits[4])
  
  #round by 2 dex. all the numbers included in range from 4 to N-1 columns
  ## numeric_cols <- sapply(IC50_dataframe, is.numeric)
  ## IC50_dataframe[,numeric_cols] <- round(IC50_dataframe[,numeric_cols],1)
  
  if(FALSE) {
    # plot IC50
    x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
    y <- predict(nls_result_ic50, data.frame(logconc=x))
    icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
      geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
      geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + ggtitle(paste0(drug_name," (dss:",round(IC50_dataframe$DSS,1),")\n")) +
      theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = ifelse(readoutCTX, "% toxicity", "% inhibition"), x = "conc(nM)")  +  ylim(-25, 125) +
      geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
      theme(plot.background = element_rect(fill = "transparent",colour = NA),
            panel.background =element_rect(fill = "transparent",colour = NA))
    
    graphics.off()
    png(filename = file.path(getwd(), "Results", "Curve_fits", "IC50", paste0(product_id, "_IC50_curve_drug.png")),width=190,height=190, bg = "transparent")
    print(icpl)
    dev.off()  
  }
  
  # plot IC50
  if(plot.fit) {
    png(paste0(output.path, "/", file, ".png"))
    showSD <- has.replicates
    plot(result_ic50, pcol="grey40", lcol="skyblue1", showEstim=.5, showInfl=TRUE, cex.main=1.5, showSDerr = showSD)
    title(main = paste0(file, ": DSS = ", dss_score))
    d <- dev.off()
  }
  ## EC50base64 <- gsub("\r?\n|\r", " ", base64::img(paste0(getwd(),"/Results/Curve_fits/EC50/",paste0(product_id, "_EC50_curve_drug.png"))))
  
  #return list with 3 nodes - 1 row for IC50 table and 1 row for EC50 table and EC50 image in base64
  ##cbind(IC50_dataframe, EC50_dataframe, EC50base64)
  if(return.fit) {
    return(list(tbl=IC50_dataframe, res=result_ic50))
  }
  return(list(tbl=IC50_dataframe))
  #c(IC50_dataframe, EC50_dataframe)
  ## END
  
}

## dose: drug dose concentration (in natural, not log10, units)
## viability: cell viability under that dose (in percent--i.e., 100 = 100%)
## NB: viability = 100 - inhibition
## We will fit a log-logistic function, which is scale invariant through the parameter b,
## as noted in Ritz et al.
## Duplicated response (i.e., viability) will be randomly jittered.  Hence, you may want
## to set.seed before calling this function.
calc.ic50 <- function(dose, viability, file, patient.id, inhibitor, return.fit = TRUE, plot.fit = FALSE) {
  
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
  if(plot.fit) {
    png(paste0(output.path, "/", file, ".png"))
    x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
    y <- predict(nls_result_ic50, data.frame(logconc=x))
    icpl <- ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
      geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
      geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) 
    icpl <- icpl + ggtitle(file)
    ## icpl <- icpl + ggtitle(paste0(drug_name," (dss:",round(IC50_dataframe$DSS,1),")\n"))
    icpl <- icpl + theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(min(-25, 1.1*min(inhibition), 1.1*min(y)),max(125,1.1*max(inhibition),1.1*max(y))) +
      geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background =element_rect(fill = "transparent",colour = NA))

    print(icpl)
    d <- dev.off()
  }

  IC50_dataframe <- data.frame(ID=file,PATIENT_ID=patient.id,DRUG_NAME=inhibitor,t(as.matrix(coef_ic50)))
  
  if(return.fit) {
    return(list(tbl=IC50_dataframe, res=nls_result_ic50))
  }
  return(list(tbl=IC50_dataframe))
  
  
}

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
target.drugs <- c("Dasatinib", "JAK Inhibitor I", "AG490", "AST-487", "Erlotinib", "Gefitinib")
flag <- raw.inhibitor.data$inhibitor %in% target.drugs
## flag <- raw.inhibitor.data$inhibitor %in% c("Dasatinib")
## flag <- flag & ( raw.inhibitor.data$patient_id %in% c(1190, 1243, 1285) )
raw.inhibitor.data <- raw.inhibitor.data[flag, ]

## BEGIN debugging
## Pick out a case in which the IC50 exceeds the max

pid <- 1035
inh <- "Gefitinib"
rep <- 2
lid <- "14-00507"
tor <- 24

id <- paste(pid, inh, rep, lid, tor, sep="-")


save.image(".Rdata")

debug.case.with.ic50.outside.range <- FALSE
if(debug.case.with.ic50.outside.range) {
  flag <- (raw.inhibitor.data$patient_id == pid) & (raw.inhibitor.data$inhibitor == inh) & (raw.inhibitor.data$replicant == rep) & (raw.inhibitor.data$lab_id == lid) & (raw.inhibitor.data$time_of_read == tor)
  df <- raw.inhibitor.data[flag,]
  r.constrained <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, "trouble", pid, inh, apply.constraints=TRUE)

  r.unconstrained <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, "trouble", pid, inh, apply.constraints=FALSE)
  
  plot(r.constrained$res)
  plot(r.unconstrained$res)
  
  
  ## STOP debugging
  stop("stop here")
}


## patient.drug.tbl <- patient.drug.tbl[patient.drug.tbl$inhibitor == "Dasatinib",]

do.merged.fits <- function(n.iters, patient.drug.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=FALSE) {
  ## tbl <- as.data.frame(matrix(data=NA, nrow=n.iters, ncol=22))
  
  tbl <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
##  foreach(i=1:n.iters) %dopar% {
  ## for(i in 1:n.iters) {
    ## if( ( i %% 10 ) == 0 ) { cat(paste0("Iter ", i, "\n")) }
    patient.id <- patient.drug.tbl$patient_id[i]
    inhibitor <- patient.drug.tbl$inhibitor[i]
    
    file <- paste(patient.id, inhibitor, sep="-")
    ## file <- paste(file, ".png", sep="")
    
    flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor)
    df <- raw.inhibitor.data[flag,,drop=F]
    ic50 <- NULL
    row <- c(file, rep(NA, 21))
    if(use.nplr) {
      r <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, file, patient.id, inhibitor,apply.constraints=apply.constraints, return.fit=FALSE)
      if(!is.null(r) && !is.null(r[["tbl"]])) {
        row <- r[["tbl"]]
        ## if(i == 1) { colnames(tbl) <- colnames(row) }
        ## results[[i]] <- list(ID=row$ID, result=NA)
        ## results[[i]] <- list(ID=row$ID, result=r[["res"]])
      } else {
##        cat(paste0("Appending row for ", file, "\n"))
##        print(row)
        ## results[[i]] <- list(ID=file, result=NA)
      }
    } else {
      ic50 <- calc.ic50(df$well_concentration, df$normalized_viability, file)
      row <- c(patient.id, inhibitor, ic50)
    }
    ## row
    ## tbl[i,] <- row
    row
  }
  cat("done\n")
  ## return(list(results=results, tbl=tbl))
  return(list(tbl=tbl))
}

do.individual.fits <- function(n.iters, patient.drug.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=FALSE) {
  ## results <- list()
  ## tbl <- as.data.frame(matrix(data=NA, nrow=n.iters, ncol=22))
  
  ## results[[i]] <- list(ID=row$ID, result=NA)
  tbl <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
  ## foreach(i=1:n.iters) %do% {
  ## for(i in 1:n.iters) {
    ## for(i in 1:nrow(patient.drug.tbl)) {
    ## if( ( i %% 10 ) == 0 ) { cat(paste0("Iter ", i, "\n")) }
    patient.id <- patient.drug.tbl$patient_id[i]
    inhibitor <- patient.drug.tbl$inhibitor[i]
    replicant <- patient.drug.tbl$replicant[i]
    lab_id <- patient.drug.tbl$lab_id[i]
    time_of_read <- patient.drug.tbl$time_of_read[i]
  
    file <- paste(patient.id, inhibitor, replicant, lab_id, time_of_read, sep="-")
    ## file <- paste0(file, ".png")
  
    flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor) & (raw.inhibitor.data$replicant == replicant) & (raw.inhibitor.data$lab_id == lab_id) & (raw.inhibitor.data$time_of_read == time_of_read)
    df <- raw.inhibitor.data[flag,,drop=F]
    ## print(df)
    ic50 <- NULL
    
    row <- c(file, rep(NA, 21))
    if(use.nplr) {
      r <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, file, patient.id, inhibitor,apply.constraints=apply.constraints, return.fit=FALSE)
      if(!is.null(r) && !is.null(r[["tbl"]])) {
        row <- r[["tbl"]]
        ## if(i == 1) { colnames(tbl) <- colnames(row) }
        ## results[[i]] <- list(ID=row$ID, result=NA)
        ## cat(paste0("Appending row for ", file, "\n"))
        ## results[[i]] <- list(ID=row$ID, result=r[["res"]])
      } else {
        ## cat(paste0("Appending row for ", file, "\n"))
        ## print(row)
        ## results[[i]] <- list(ID=file, result=NA)
      }
    } else {
      ic50 <- calc.ic50(df$well_concentration, df$normalized_viability, file)
      row <- c(patient.id, inhibitor, replicant, lab_id, time_of_read, ic50)
    }
    ## tbl[i,] <- row
    row
  }
  ## return(list(results=results, tbl=tbl))
  return(list(tbl=tbl))
}

patient.drug.tbl <- unique(raw.inhibitor.data[,c("patient_id","inhibitor")])

n.iters <- nrow(patient.drug.tbl)
cat(paste0("Total iterations: ", n.iters, "\n"))
## n.iters <- 1
## n.iters <- 20
## do -> dopar

## Do the constrained. merged fits
tmp.tbl <- patient.drug.tbl[5460:nrow(patient.drug.tbl),]
tmp.tbl <- patient.drug.tbl
n.iters <- nrow(tmp.tbl)
## n.iters <- 10
cat("Fitting merged with constraints\n")
res.list <- do.merged.fits(n.iters, tmp.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=TRUE)
tbl2.constrained <- res.list[["tbl"]]
## results2.constrained <- res.list[["results"]]

## Do the merged fits
cat("Fitting merged without constraints\n")
res.list <- do.merged.fits(n.iters, patient.drug.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=FALSE)
tbl2 <- res.list[["tbl"]]
## results2 <- res.list[["results"]]


if(use.nplr == FALSE) {
  rownames(tbl2) <- NULL
  colnames(tbl2) <- c("patient.id", "inhibitor", "ic50")
  rownames(tbl2.constrained) <- NULL
  colnames(tbl2.constrained) <- c("patient.id", "inhibitor", "ic50")
}

cat("Saving tbl2\n")
save.image(file=".RData.tbl2")
cat("Saved tbl2\n")

## Fit the individual curves
patient.drug.tbl <- unique(raw.inhibitor.data[,c("patient_id","inhibitor", "replicant", "lab_id", "time_of_read")])

n.iters <- nrow(patient.drug.tbl)
cat(paste0("Total iterations: ", n.iters, "\n"))
## n.iters <- 1
## n.iters <- 20
cat("Fitting individual without constraints\n")
res.list <- do.individual.fits(n.iters, patient.drug.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=FALSE)
tbl <- res.list[["tbl"]]
## results <- res.list[["results"]]

## Do the constrained, individual fits
cat("Fitting individual with constraints\n")
## n.iters <- nrow(patient.drug.tbl)
## n.iters <- 10
res.list <- do.individual.fits(n.iters, patient.drug.tbl, raw.inhibitor.data, use.nplr=TRUE, apply.constraints=TRUE)
tbl.constrained <- res.list[["tbl"]]
## results.constrained <- res.list[["results"]]


if(use.nplr == FALSE) {
  rownames(tbl) <- NULL
  colnames(tbl) <- c("patient.id", "inhibitor", "replicant", "lab_id", "time_of_read", "ic50")
  rownames(tbl.constrained) <- NULL
  colnames(tbl.constrained) <- c("patient.id", "inhibitor", "replicant", "lab_id", "time_of_read", "ic50")
}

save.image(file=".RData.tbl")

stop("stop -- successful analysis")

## Read in Swapnil's data
## Drug Name and replicate number make unique Drug ID and Patient ID and Lab_ID make unique screen ID.
swap <- read.table("OHSU_DSS.txt", sep="\t", header=TRUE)

swap <- swap[swap$DRUG_NAME %in% target.drugs, ]

swap$rep <- unlist(lapply(swap$ID, function(id) {
    r <- regexpr(pattern="_(\\d+)$", text=id,perl=TRUE)
    substring(id, as.numeric(r) + 1, as.numeric(r) + attr(r, "match.length"))
}))

screen.cols <- colnames(swap)[grepl(pattern="^P_", colnames(swap))]

## Extract the patient id and lab id from Swapnil's file
fimm.ohsu <- c()
for(col in screen.cols) {
  str <- col
  r <- regexpr(pattern="P_(\\d+)_(\\d+)\\.(\\d+)$", text=str,perl=TRUE)
  patient.id <- substring(str, attr(r,"capture.start")[1], attr(r,"capture.start")[1] + attr(r, "capture.length")[1] - 1)
  tmp1 <- substring(str, attr(r,"capture.start")[2], attr(r,"capture.start")[2] + attr(r, "capture.length")[2] - 1)
  tmp2 <- substring(str, attr(r,"capture.start")[3], attr(r,"capture.start")[3] + attr(r, "capture.length")[3] - 1)
  lab_id <- paste0(tmp1, "-", tmp2)
  for(i in 1:nrow(swap)) {
    inhibitor <- as.character(swap$DRUG_NAME[i])
    replicant <- swap$rep[i]
    file <- paste(patient.id, inhibitor, replicant, lab_id, sep="-")
    fimm.ohsu <- rbind(fimm.ohsu, c(file, swap[i,col]))
  }
}
colnames(fimm.ohsu) <- c("id", "fimm.dss")
fimm.ohsu <- as.data.frame(fimm.ohsu)
fimm.ohsu$fimm.dss <- as.numeric(as.character(fimm.ohsu$fimm.dss))

## Swapnil does not include the time_of_read in his id, which is the last thing
## I include in mine--strip this off.
tbl$short.id <- unlist(lapply(tbl$ID, function(str) {
  r <- regexpr(pattern="^(.+)-\\d+$", text=as.character(str),perl=TRUE)
  id <- substring(str, attr(r,"capture.start")[1], attr(r,"capture.start")[1] + attr(r, "capture.length")[1] - 1)
  id
}))

## At this point the id in the FIMM DSS values for the OHSU data (fimm.ohsu$id) and 
## in my DSS results for the OHSU data (tbl$short.id) is
## patient_id-inhibitor-replicant-lab_id

## Merge the FIMM and Sage calculations of the DSS values.
## fimm.sage.ohsu <- merge(fimm.ohsu, tbl, by.x="id", by.y="short.id")

save.image(".RData.fimm.sage.merged")

## Compare constrained vs unconstrained individual fits
individ.fits <- merge(tbl, tbl.constrained, by="ID", suffixes=c(".unconstrained", ".constrained"))



## Annotate those that hit a constraint boundary
individ.fits$hit.constraint.bndry <- rep(NA, nrow(individ.fits))
not.na <- !is.na(individ.fits$DSS.constrained)
individ.fits$hit.constraint.bndry[not.na] <- rep(FALSE, length(which(!not.na)))
for(param in c("MIN", "MAX", "IC50", "SLOPE")) {
  col <- paste0(param, ".constrained")
  lb.col <- paste0(param, ".lb.constrained")
  ub.col <- paste0(param, ".ub.constrained")
##  not.na <- unlist(apply(individ.fits[,c(col, lb.col, ub.col)], 1, function(x) !any(is.na(x))))
  individ.fits$hit.constraint.bndry[not.na] <- 
    individ.fits$hit.constraint.bndry[not.na] | unlist(apply(individ.fits[not.na,c(col, lb.col, ub.col)], 1,
                                                        function(row) {
                                                          val <- as.numeric(row[1])
                                                          lb <- as.numeric(row[2])
                                                          ub <- as.numeric(row[2])
                                                          ret <- FALSE
                                                          if(!is.infinite(lb)) {
                                                            if((abs(val-lb)/(1+max(abs(val),abs(lb)))) < 0.05) {ret <- TRUE}
                                                          }
                                                          if(!is.infinite(ub)) {
                                                            if((abs(val-ub)/(1+max(abs(val),abs(ub)))) < 0.05) {ret <- TRUE}
                                                          }
                                                          ret
                                                        }))
}

flag <- individ.fits$hit.constraint.bndry == TRUE
cols <- c()
for(param in c("MIN", "MAX", "IC50", "SLOPE")) {
  cols <- c(cols, paste0(param, ".constrained"), paste0(param, ".lb.constrained"), paste0(param, ".ub.constrained"))
}
head(individ.fits[flag, c(cols, "DSS.constrained")])

## Calculate difference between unconstrained and constrained
individ.fits$diff <- rep(NA, nrow(individ.fits))
## flag <- !is.na(individ.fits$DSS.unconstrained) & !is.na(individ.fits$DSS.constrained) & (individ.fits$DSS.unconstrained != 0) & (individ.fits$DSS.constrained != 0)
flag <- !is.na(individ.fits$DSS.unconstrained) & !is.na(individ.fits$DSS.constrained)
individ.fits$diff[flag] <- as.numeric(individ.fits$DSS.unconstrained[flag]) - as.numeric(individ.fits$DSS.constrained[flag])

## Exclude those along axis and determine sd
flag <- !is.na(individ.fits$DSS.unconstrained) & !is.na(individ.fits$DSS.constrained) & (individ.fits$DSS.unconstrained != 0) & (individ.fits$DSS.constrained != 0)
std.dev <- sd(individ.fits$diff[flag], na.rm=TRUE)

## Annotate those samples for which unconstrained is significantly different than constrained (i.e., > 2 std devs)
flag <- !is.na(individ.fits$DSS.unconstrained) & !is.na(individ.fits$DSS.constrained)
individ.fits$sig.diff <- rep(NA, nrow(individ.fits))
individ.fits$sig.diff[flag] <- unlist(lapply(individ.fits$diff[flag], function(x) ifelse(abs(x) > 2 * std.dev, TRUE, FALSE)))

library(vcd)
tab <- table(!individ.fits$sig.diff, individ.fits$hit.constraint.bndry, dnn = c("Constrained DSS = Unconstrained DSS", "Fit Hit Constraint Boundary"))

## NB: these are for those that are non-zero
## Show me
png("ohsu-nplr-fit-constrained-vs-boundary.png")
mosaic(tab, pop=FALSE)
labeling_cells(text = tab, gp_text = gpar(fontsize = 10), margin=0)(tab)
d <- dev.off()


flag <- individ.fits$hit.constraint.bndry
## x <- individ.fits$DSS.unconstrained
## y <- individ.fits$DSS.constrained
## plot(x[flag], y[flag], xlab="Unconstrained NPLR OHSU DSS2", ylab="Constrained NPLR OHSU DSS2", col = "red", pch = 20)
## points(x[!flag], y[!flag], col = "black", pch = 20)

x <- seq(from=0, to=50, by=0.1)
sd.df <- data.frame(x = x, lb = x - 2 * std.dev, ub = x + 2 * std.dev)

library(gridExtra)
away.from.bndry <- !is.na(individ.fits$hit.constraint.bndry) & (individ.fits$hit.constraint.bndry == FALSE)
hit.bndry <- !is.na(individ.fits$hit.constraint.bndry) & (individ.fits$hit.constraint.bndry == TRUE)
individ.fits$DSS.unconstrained <- as.numeric(individ.fits$DSS.unconstrained)
individ.fits$DSS.constrained <- as.numeric(individ.fits$DSS.constrained)
p1 <- ggplot(data = individ.fits[away.from.bndry,], aes(x = DSS.unconstrained, y = DSS.constrained))
p1 <- p1 + ggtitle(paste0("Fits away from constraint boundary (n = ", length(which(away.from.bndry)), ")"))
p1 <- p1 + geom_point(color = "black")
p1 <- p1 + xlab("Unconstrained NPLR OHSU DSSS2")
p1 <- p1 + ylab("Constrained NPLR OHSU DSSS2")
p1 <- p1 + geom_line(data = sd.df, aes(x = x, y = ub), color = "blue")
p1 <- p1 + geom_line(data = sd.df, aes(x = x, y = lb), color = "blue")

p2 <- ggplot(data = individ.fits[hit.bndry,], aes(x = DSS.unconstrained, y = DSS.constrained))
p2 <- p2 + ggtitle(paste0("Fits hit constraint boundary (n = ", length(which(hit.bndry)), ")"))
p2 <- p2 + geom_point(color = "red")
p2 <- p2 + xlab("Unconstrained NPLR OHSU DSSS2")
p2 <- p2 + ylab("Constrained NPLR OHSU DSSS2")
p2 <- p2 + geom_line(data = sd.df, aes(x = x, y = ub), color = "blue")
p2 <- p2 + geom_line(data = sd.df, aes(x = x, y = lb), color = "blue")

## Show me
png("ohsu-nplr-constrained-vs-unconstrained-dss2.png")
grid.arrange(p1, p2)
d <- dev.off()

## plot(x[!flag], y[!flag], xlab="Unconstrained NPLR OHSU DSS2", ylab="Constrained NPLR OHSU DSS2", col = "black", pch = 20)
## plot(x[flag], y[flag], xlab="Unconstrained NPLR OHSU DSS2", ylab="Constrained NPLR OHSU DSS2", col = "red", pch = 20)
## points(x[flag], y[flag], col = "red", pch = 20)
## x <- seq(from=0, to=50, by=0.1)
## lines(x = x, y = x + 3 * std.dev, lty = 1, col = "blue")
## lines(x = x, y = x - 3 * std.dev, lty = 1, col = "blue")

## Calculate the std dev across replicates based on which of the 4 categories the fits fall in to:
## (1) constrained != unconstrained & hit boundary
## (2) constrained != unconstrained & does not hit boundary
## (3) constrained == unconstrained & hit boundary
## (4) constrained == unconstrained & does not hit boundary

unique.patients.and.drugs <- unique(individ.fits[,c("PATIENT_ID.unconstrained", "DRUG_NAME.unconstrained")])
colnames(unique.patients.and.drugs) <- c("PATIENT_ID", "DRUG_NAME")
unique.patients.and.drugs$std.diff.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$std.diff.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$std.same.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$std.same.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$cv.diff.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$cv.diff.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$cv.same.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$cv.same.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$mean.diff.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$mean.diff.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$mean.same.and.bndry <- rep(NA, nrow(unique.patients.and.drugs))
unique.patients.and.drugs$mean.same.and.no.bndry <- rep(NA, nrow(unique.patients.and.drugs))
individ.fits$DSS.unconstrained <- as.numeric(individ.fits$DSS.unconstrained)
for(i in 1:nrow(unique.patients.and.drugs)) {
  patient <- unique.patients.and.drugs$PATIENT_ID[i]
  drug <- unique.patients.and.drugs$DRUG_NAME[i]
  tbl <- individ.fits[(individ.fits$PATIENT_ID.unconstrained == patient) & (individ.fits$DRUG_NAME.unconstrained == drug) & (!is.na(individ.fits$DSS.unconstrained)) & (!is.na(individ.fits$sig.diff)) & (!is.na(individ.fits$hit.constraint.bndry)),,drop=F]
  if(nrow(tbl) <= 1) { next }
  
  flag <- (tbl$sig.diff == TRUE) & (tbl$hit.constraint.bndry == TRUE)
  if(length(which(flag)) > 1) {
    unique.patients.and.drugs$std.diff.and.bndry[i] <- sd(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$mean.diff.and.bndry[i] <- mean(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$cv.diff.and.bndry[i] <- sd(tbl$DSS.unconstrained[flag])/mean(tbl$DSS.unconstrained[flag])
  }

  flag <- (tbl$sig.diff == TRUE) & (tbl$hit.constraint.bndry == FALSE)
  if(length(which(flag)) > 1) {
    unique.patients.and.drugs$std.diff.and.no.bndry[i] <- sd(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$mean.diff.and.no.bndry[i] <- mean(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$cv.diff.and.no.bndry[i] <- sd(tbl$DSS.unconstrained[flag])/mean(tbl$DSS.unconstrained[flag])
  }

  flag <- (tbl$sig.diff == FALSE) & (tbl$hit.constraint.bndry == TRUE)
  if(length(which(flag)) > 1) {
    unique.patients.and.drugs$std.same.and.bndry[i] <- sd(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$mean.same.and.bndry[i] <- mean(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$cv.same.and.bndry[i] <- sd(tbl$DSS.unconstrained[flag])/mean(tbl$DSS.unconstrained[flag])   
  }
  
  flag <- (tbl$sig.diff == FALSE) & (tbl$hit.constraint.bndry == FALSE)
  if(length(which(flag)) > 1) {
    unique.patients.and.drugs$std.same.and.no.bndry[i] <- sd(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$mean.same.and.no.bndry[i] <- mean(tbl$DSS.unconstrained[flag])
    unique.patients.and.drugs$cv.same.and.no.bndry[i] <- sd(tbl$DSS.unconstrained[flag])/mean(tbl$DSS.unconstrained[flag])    
  }
}

dss.std <- data.frame(condition=rep("DIFF.BNDRY", nrow(unique.patients.and.drugs)), std=unname(unique.patients.and.drugs$std.diff.and.bndry))
dss.std <- rbind(dss.std, data.frame(condition=rep("DIFF.NO.BNDRY", nrow(unique.patients.and.drugs)), std=unname(unique.patients.and.drugs$std.diff.and.no.bndry)))
dss.std <- rbind(dss.std, data.frame(condition=rep("SAME.BNDRY", nrow(unique.patients.and.drugs)), std=unname(unique.patients.and.drugs$std.same.and.bndry)))
dss.std <- rbind(dss.std, data.frame(condition=rep("SAME.NO.BNDRY", nrow(unique.patients.and.drugs)), std=unname(unique.patients.and.drugs$std.same.and.no.bndry)))

library(ggbeeswarm)
## Show me
png("std-beeswarm.png")
g <- ggplot(data = dss.std, aes(x = condition, y = std))
g <- g + geom_boxplot()
g <- g + geom_beeswarm()
g <- g + ylab("DSS Standard Deviation")
print(g)
d <- dev.off()

dss.cv <- data.frame(condition=rep("DIFF.BNDRY", nrow(unique.patients.and.drugs)), cv=unname(unique.patients.and.drugs$cv.diff.and.bndry))
dss.cv <- rbind(dss.cv, data.frame(condition=rep("DIFF.NO.BNDRY", nrow(unique.patients.and.drugs)), cv=unname(unique.patients.and.drugs$cv.diff.and.no.bndry)))
dss.cv <- rbind(dss.cv, data.frame(condition=rep("SAME.BNDRY", nrow(unique.patients.and.drugs)), cv=unname(unique.patients.and.drugs$cv.same.and.bndry)))
dss.cv <- rbind(dss.cv, data.frame(condition=rep("SAME.NO.BNDRY", nrow(unique.patients.and.drugs)), cv=unname(unique.patients.and.drugs$cv.same.and.no.bndry)))

## Show me
png("cv-beeswarm.png")
g <- ggplot(data = dss.cv, aes(x = condition, y = cv))
g <- g + geom_boxplot()
g <- g + geom_beeswarm()
g <- g + ylab("DSS Coefficient of Variation")
print(g)
d <- dev.off()

## TODO
## do above, but plot GOF

## Show me
png("mean-vs-std-same-no-bndry.png")
g <- ggplot(data = unique.patients.and.drugs, aes(x = mean.same.and.no.bndry, y = std.same.and.no.bndry))
g <- g + geom_point()
g <- g + ggtitle("Constrained DSS = Unconstrained DSS and Fit did not hit boundary")
g <- g + xlab("DSS Mean")
g <- g + ylab("DSS Standard Deviation")
print(g)            
d <- dev.off()

## Show me
png("mean-vs-std-same-bndry.png")
g <- ggplot(data = unique.patients.and.drugs, aes(x = mean.same.and.bndry, y = std.same.and.bndry))
g <- g + geom_point()
g <- g + ggtitle("Constrained DSS = Unconstrained DSS and Fit hit boundary")
g <- g + xlab("DSS Mean")
g <- g + ylab("DSS Standard Deviation")
print(g)            
d <- dev.off()


## Merge the FIMM and Sage calculations of the DSS values.
fimm.sage.ohsu <- merge(merge(fimm.ohsu, individ.fits, by.x="id", by.y="short.id"), unique.patients.and.drugs, by.x = c("PATIENT_ID.unconstrained", "DRUG_NAME.unconstrained"), by.y = c("PATIENT_ID", "DRUG_NAME"))

flag <- (fimm.sage.ohsu$sig.diff == TRUE) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g1 <- g1 + ggtitle(paste0("Constrained != Unconstrained\nFit hit constraints (n = ", length(which(flag)), ")"))
g1 <- g1 + geom_point()
g1 <- g1 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- (fimm.sage.ohsu$sig.diff == TRUE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
g2 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g2 <- g2 + ggtitle(paste0("Constrained != Unconstrained\nFit did not hit constraints (n = ", length(which(flag)), ")"))
g2 <- g2 + geom_point()
g2 <- g2 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
g3 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g3 <- g3 + ggtitle(paste0("Constrained = Unconstrained\nFit hit constraints (n = ", length(which(flag)), ")"))
g3 <- g3 + geom_point()
g3 <- g3 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
g4 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g4 <- g4 + ggtitle(paste0("Constrained = Unconstrained\nFit did not hit constraints (n = ", length(which(flag)), ")"))
g4 <- g4 + geom_point()
g4 <- g4 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

## Show me
png("fimm-vs-ohsu-constrained-hit.png")
grid.arrange(grobs = list(g1, g2, g3, g4), layout_matrix = rbind(c(1,2),c(3,4)))
d <- dev.off()

## Sanity check that constrained = unconstrained is true
## Show me
png("ohsu-constrained-vs-unconstrained.png")
flag <- fimm.sage.ohsu$sig.diff == FALSE
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = DSS.unconstrained)) + xlab("Sage Constrained DSS") + ylab("Sage Unconstrained DSS")
g1 <- g1 + ggtitle(paste0("Constrained = Unconstrained (n = ", length(which(flag)), ")"))
g1 <- g1 + geom_point()
g1 <- g1 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- fimm.sage.ohsu$sig.diff == TRUE
g2 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = DSS.unconstrained)) + xlab("Sage Constrained DSS") + ylab("Sage Unconstrained DSS")
g2 <- g2 + ggtitle(paste0("Constrained != Unconstrained (n = ", length(which(flag)), ")"))
g2 <- g2 + geom_point()
g2 <- g2 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

grid.arrange(g1, g2)
d <- dev.off()

flag <- (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
tbl <- fimm.sage.ohsu[flag,]

## Annotate (only the Sage constrained = unconstrained and not hit boundary cases) as:
## (1) FIMM and Sage DSS the same
## Base this on std.dev calculated above for constrained vs unconstrained (std.dev)
tbl$diff <- tbl$fimm.dss - tbl$DSS.unconstrained
tbl$SAME <- rep(FALSE, nrow(tbl))
flag <- !is.na(tbl$diff)
tbl$SAME[flag] <- unlist(lapply(tbl$diff[flag], function(x) ifelse(abs(x) > 2 * std.dev, FALSE, TRUE)))

## (2) !same AND FIMM dss = 0
tbl$FIMM0 <- rep(FALSE, nrow(tbl))
flag <- !is.na(tbl$diff) & !is.na(tbl$SAME) & (tbl$SAME == FALSE)
tbl$FIMM0[flag] <- unlist(lapply(tbl$fimm.dss[flag], function(x) ifelse(x == 0, TRUE, FALSE))) 

## (3) !same AND Sage dss = 0
tbl$SAGE0 <- rep(FALSE, nrow(tbl))
flag <- !is.na(tbl$diff) & !is.na(tbl$SAME) & (tbl$SAME == FALSE)
tbl$SAGE0[flag] <- unlist(lapply(tbl$DSS.unconstrained[flag], function(x) ifelse(x == 0, TRUE, FALSE))) 

## (4) !same AND FIMM dss != 0 AND Sage dss != 0
tbl$DIFF <- rep(FALSE, nrow(tbl))
flag <- !is.na(tbl$diff) & !is.na(tbl$SAME) & (tbl$SAME == FALSE) & (tbl$FIMM0 == FALSE) & (tbl$SAGE0 == FALSE)
tbl$DIFF[flag] <- TRUE

tbl$condition <- unlist(apply(tbl[,c("SAME", "FIMM0", "SAGE0", "DIFF")], 1,
                              function(row) {
                                c("SAME", "FIMM0", "SAGE0", "DIFF")[(row)]
                              }))

save.image(".Rdata.current")

## For the above 4 groups plot: IC50, MIN, MAX, SLOPE, and GOF
cols <- c("IC50.constrained", "MIN.constrained", "MAX.constrained", "SLOPE.constrained", "GOF.constrained")
grobs <- list()
for(i in 1:length(cols)) {
  col <- cols[i]
  prefix <- gsub(x=col, pattern=".constrained", replacement="")
  print(prefix)
  tbl[,col] <- as.numeric(tbl[,col])
##  png(paste0(prefix, "-fimm-vs-sage-beeswarm.png"))
  ## Turn off outliers
  g <- ggplot(data = tbl, aes_string(x = "condition", y = col)) + ylab(prefix)
  g <- g + theme(axis.text.x = element_text(angle = 90))
  g <- g + geom_boxplot(outlier.shape = NA)
##  g <- g + geom_beeswarm()
  g <- g + geom_jitter()
  grobs[[i]] <- g
##  print(g)
##  d <- dev.off()
}

grid.arrange(grobs = grobs, layout_matrix=rbind(c(1,2,3),c(4,5,NA)))


## colors
cols <- c("IC50.constrained", "MIN.constrained", "MAX.constrained", "SLOPE.constrained", "GOF.constrained")
for(i in 1:length(cols)) {
  col <- cols[i]
  prefix <- gsub(x=col, pattern=".constrained", replacement="")
  print(prefix)
  tmp <- tbl
  tmp[,col] <- as.numeric(tmp[,col])
  if(prefix == "IC50") {
    tmp <- tmp[tmp[,col] < 0.25,]
  } else if(prefix == "MAX") {
    tmp <- tmp[tmp[,col] > 0,]
  } else if(prefix == "MIN") {
    tmp <- tmp[tmp[,col] > -100,]
  } else if(prefix == "SLOPE") {
    tmp <- tmp[tmp[,col] < 50,]
  }

  ## show me
  png(paste0(prefix, "-fimm-vs-sage-density.png"))
  ## Turn off outliers
  if(FALSE) {
    g <- ggplot() + ylab(paste0(prefix, " density"))
    for(cond in unique(tmp$condition)) {
      tbl.cond <- tmp[tmp$condition == cond,]
      tbl.cond[,col] <- as.numeric(tbl.cond[,col])
    
      g <- g + geom_density(data = tbl.cond, aes_string(x = col))    
    }
  } else {
    g <- ggplot(data = tmp, aes_string(x = col))
    g <- g + xlab(prefix)
    g <- g + ylab("density")
    g <- g + geom_density(aes(group = condition, colour = condition, fill = condition), alpha=0.3)    
  }
  
  print(g)
  d <- dev.off()
}

## Determine GOF threshold and apply that to 
## (1) FIMM vs Sage
## (2) unconstrained vs constrained

dss.flag <- !is.na(fimm.sage.ohsu$DSS.constrained) & (fimm.sage.ohsu$DSS.constrained > 0) & (fimm.sage.ohsu$sig.diff == FALSE)
## plot(fimm.sage.ohsu$DSS.constrained[dss.flag], fimm.sage.ohsu$DSS.unconstrained[dss.flag])

x.vals <- fimm.sage.ohsu$GOF.constrained[dss.flag]

# Define the density 95% confidence interval around GOF=0.6
delta <- 0.1
dens <- density(x.vals)
## y.subset <- d$y[(d$x < (0.5+delta)) & (d$x > (0.5-delta))]
range.flag <- (dens$x < (0.6+delta)) & (dens$x > (0.6-delta))
x.subset <- dens$x[range.flag]
y.subset <- dens$y[range.flag]
plot(density(y.subset))
## q <- quantile(y.subset, probs=c(0.05,0.5,0.95))
## plot(density(y.subset))

## Define the gof cutoff
mu <- mean(y.subset)
sigma <- sd(y.subset)
indx <- which(dens$y==min(dens$y[dens$y>(mu+3*sigma)]))[1]
gof.threshold <- dens$x[indx]

fimm.sage.ohsu$GOF.constrained <- as.numeric(fimm.sage.ohsu$GOF.constrained)
## Show me
png("nplr-gof-threshold.png")
g <- ggplot(data = fimm.sage.ohsu[dss.flag,], aes(x = GOF.constrained))
g <- g + ggtitle(paste0("GOF threshold = ", round(gof.threshold, digits=2), " (Non-zero DSS fits; constrained = unconstrained)"))
g <- g + xlab("nplr Constrained Goodness of Fit")
g <- g + ylab("density")
g <- g + geom_density()
g <- g + geom_vline(xintercept = gof.threshold)
g <- g + geom_hline(yintercept = (mu - 3 * sigma), linetype="dashed")
g <- g + geom_hline(yintercept = (mu - 0 * sigma))
g <- g + geom_hline(yintercept = (mu + 3 * sigma), linetype="dashed")
print(g)
d <- dev.off()

png("nplr-gof-hist-threshold.png")
g <- ggplot(data = fimm.sage.ohsu[dss.flag,], aes(x = GOF.constrained))
g <- g + ggtitle(paste0("GOF threshold = ", round(gof.threshold, digits=2), " (Non-zero DSS fits; constrained = unconstrained)"))
g <- g + xlab("nplr Constrained Goodness of Fit")
g <- g + ylab("density")
g <- g + geom_histogram(bins=100)
g <- g + geom_vline(xintercept = gof.threshold)
g <- g + geom_hline(yintercept = (mu - 3 * sigma), linetype="dashed")
g <- g + geom_hline(yintercept = (mu - 0 * sigma))
g <- g + geom_hline(yintercept = (mu + 3 * sigma), linetype="dashed")
print(g)
d <- dev.off()


good.dss.flag <- !is.na(fimm.sage.ohsu$DSS.constrained) & (fimm.sage.ohsu$GOF.constrained > gof.threshold)

## Repeat the plots from above
flag <- good.dss.flag & (fimm.sage.ohsu$sig.diff == TRUE) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g1 <- g1 + ggtitle(paste0("Constrained != Unconstrained\nFit hit constraints\nGOF > ", round(gof.threshold, digits=2), " (n = ", length(which(flag))))
g1 <- g1 + geom_point()
g1 <- g1 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- good.dss.flag & (fimm.sage.ohsu$sig.diff == TRUE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
g2 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g2 <- g2 + ggtitle(paste0("Constrained != Unconstrained\nFit did not hit constraints\nGOF > ", round(gof.threshold, digits=2), " (n = ", length(which(flag))))
g2 <- g2 + geom_point()
g2 <- g2 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- good.dss.flag & (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
g3 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g3 <- g3 + ggtitle(paste0("Constrained = Unconstrained\nFit hit constraints\nGOF > ", round(gof.threshold, digits=2), " (n = ", length(which(flag))))
g3 <- g3 + geom_point()
g3 <- g3 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- good.dss.flag & (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
g4 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g4 <- g4 + ggtitle(paste0("Constrained = Unconstrained\nFit did not hit constraints\nGOF > ", round(gof.threshold, digits=2), " (n = ", length(which(flag))))
g4 <- g4 + geom_point()
g4 <- g4 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

## Show me
png("fimm-vs-ohsu-constrained-hit-gof-threshold.png")
grid.arrange(grobs = list(g1, g2, g3, g4), layout_matrix = rbind(c(1,2),c(3,4)))
d <- dev.off()

away.from.bndry <- !is.na(fimm.sage.ohsu$hit.constraint.bndry) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
hit.bndry <- !is.na(fimm.sage.ohsu$hit.constraint.bndry) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
fimm.sage.ohsu$DSS.unconstrained <- as.numeric(fimm.sage.ohsu$DSS.unconstrained)
fimm.sage.ohsu$DSS.constrained <- as.numeric(fimm.sage.ohsu$DSS.constrained)
p1 <- ggplot(data = fimm.sage.ohsu[away.from.bndry,], aes(x = DSS.unconstrained, y = DSS.constrained))
p1 <- p1 + ggtitle(paste0("Fits away from constraint boundary (n = ", length(which(away.from.bndry)), ")"))
p1 <- p1 + geom_point(color = "black")
p1 <- p1 + xlab(paste0("Unconstrained NPLR OHSU DSSS2 (GOF > ", round(gof.threshold, digits=2), ")"))
p1 <- p1 + ylab(paste0("Constrained NPLR OHSU DSSS2 (GOF > ", round(gof.threshold, digits=2), ")"))
p1 <- p1 + geom_line(data = sd.df, aes(x = x, y = ub), color = "blue")
p1 <- p1 + geom_line(data = sd.df, aes(x = x, y = lb), color = "blue")

p2 <- ggplot(data = fimm.sage.ohsu[hit.bndry,], aes(x = DSS.unconstrained, y = DSS.constrained))
p2 <- p2 + ggtitle(paste0("Fits hit constraint boundary (n = ", length(which(hit.bndry)), ")"))
p2 <- p2 + geom_point(color = "red")
p2 <- p2 + xlab(paste0("Unconstrained NPLR OHSU DSSS2 (GOF > ", round(gof.threshold, digits=2), ")"))
p2 <- p2 + ylab(paste0("Constrained NPLR OHSU DSSS2 (GOF > ", round(gof.threshold, digits=2), ")"))
p2 <- p2 + geom_line(data = sd.df, aes(x = x, y = ub), color = "blue")
p2 <- p2 + geom_line(data = sd.df, aes(x = x, y = lb), color = "blue")

## Show me
png("ohsu-nplr-constrained-vs-unconstrained-dss2-gof-threshold.png")
grid.arrange(p1, p2)
d <- dev.off()

## Plot a few of the poorest fits (between FIMM and Sage)
fimm.sage.ohsu$fimm.sage.diff <- as.numeric(fimm.sage.ohsu$fimm.dss) - as.numeric(fimm.sage.ohsu$DSS.constrained)
fimm.sage.ohsu <- fimm.sage.ohsu[order(abs(fimm.sage.ohsu$fimm.sage.diff), decreasing=TRUE),]
head(fimm.sage.ohsu[,c("fimm.dss", "DSS.constrained")])

discrepant.fits <- fimm.sage.ohsu[(fimm.sage.ohsu$fimm.dss > 1) & (fimm.sage.ohsu$DSS.constrained > 1),]
## raw.inhibitor.data$id <- apply(raw.inhibitor.data[,c("patient_id", "inhibitor", "replicant", "lab_id", "time_of_read")], 1, function(row) paste(row, collapse="-"))
raw.inhibitor.data$id <- apply(raw.inhibitor.data[,c("patient_id", "inhibitor", "replicant", "lab_id")], 1, function(row) paste(row, collapse="-"))
raw.inhibitor.data$id <- unlist(lapply(raw.inhibitor.data$id, function(x) gsub(x=x, pattern=" ", replacement="")))
for(i in 1:10) {
  bad.id <- discrepant.fits$id[i]
  print(as.character(bad.id))
  flag <- (as.character(raw.inhibitor.data$id) == as.character(bad.id))  
  df <- raw.inhibitor.data[flag,]
  pid <- unique(df$patient_id)
  inh <- unique(df$inhibitor)
  file <- paste0(bad.id, "-nplr-inconsistent")
  r.constrained <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, file=file, pid, inh, apply.constraints=TRUE, plot.fit=TRUE)
  
  file <- paste0(bad.id, "-drc-inconsistent")
  r.drc <- calc.ic50(df$well_concentration, df$normalized_viability, file=file, pid, inh, plot.fit=TRUE)

  file <- paste0(output.path, "/", bad.id, "-inconsistent.tsv")
  write.table(file=file, df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

## TODO
## 1. can drc fit outside positive rnage
## 2. what happens if we exclude those with values over 100
## Define the max response for each patient/drug/replicant/lab_id
## max.agg <- aggregate(tmp$well_concentration, by=list(raw.inhibitor.data$patient_id, raw.inhibitor.data$inhibitor, raw.inhibitor.data$replicant, raw.inhibitor.data$lab_id), max)
max.agg <- aggregate(raw.inhibitor.data$normalized_viability, by=list(raw.inhibitor.data$id), max)
colnames(max.agg) <- c("id", "max.viability")
ids.outside.range <- max.agg$id[max.agg$max.viability > 100]

flag <- (fimm.sage.ohsu$id %in% ids.outside.range)
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("nplr Constrained DSS") + ylab("drc DSS")
g1 <- g1 + ggtitle(paste0("Response > 100\n(n = ", length(which(flag)), ")"))
g1 <- g1 + geom_point()
g1 <- g1 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- !(fimm.sage.ohsu$id %in% ids.outside.range)
g2 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("nplr Constrained DSS") + ylab("drc DSS")
g2 <- g2 + ggtitle(paste0("Response <= 100\n(n = ", length(which(flag)), ")"))
g2 <- g2 + geom_point()
g2 <- g2 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

png("ohsu-nplr-vs-drc-response.png")
grid.arrange(g1, g2)
d <- dev.off()

base.flag <- !is.na(fimm.sage.ohsu$DSS.constrained) & (fimm.sage.ohsu$GOF.constrained > gof.threshold) & (fimm.sage.ohsu$sig.diff == FALSE) & (fimm.sage.ohsu$hit.constraint.bndry == FALSE)
flag <- base.flag & (fimm.sage.ohsu$id %in% ids.outside.range)
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("nplr Constrained DSS") + ylab("drc DSS")
g1 <- g1 + ggtitle(paste0("Constrained = Unconstrained\nFit did not hit constraints\nGOF > ", round(gof.threshold, digits=2), "\nResponse > 100 (n = ", length(which(flag)), ")"))
g1 <- g1 + geom_point()
g1 <- g1 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

flag <- base.flag & !(fimm.sage.ohsu$id %in% ids.outside.range)
g2 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("nplr Constrained DSS") + ylab("drc DSS")
g2 <- g2 + ggtitle(paste0("Constrained = Unconstrained\nFit did not hit constraints\nGOF > ", round(gof.threshold, digits=2), "\nResponse <= 100 (n = ", length(which(flag)), ")"))
g2 <- g2 + geom_point()
g2 <- g2 + geom_line(data = data.frame(x = c(0, 50), y = c(0, 50)), aes(x = x, y = y), colour="blue", linetype = "dashed")

png("ohsu-nplr-vs-drc-response-no-const-bndry.png")
grid.arrange(g1, g2)
d <- dev.off()

## STOPPED HERE
save.image(".Rdata.current")


## Repeat the plots from above
flag <- good.dss.flag & (fimm.sage.ohsu$sig.diff == TRUE) & (fimm.sage.ohsu$hit.constraint.bndry == TRUE)
g1 <- ggplot(data = fimm.sage.ohsu[flag,], aes(x = DSS.constrained, y = fimm.dss)) + xlab("Sage Constrained DSS") + ylab("FIMM DSS")
g1 <- g1 + ggtitle(paste0("Constrained != Unconstrained\nFit hit constraints\nGOF > ", round(gof.threshold, digits=2), " (n = ", length(which(flag))))


p1 <- p1 + ggtitle(paste0("Fits away from constraint boundary (n = ", length(which(away.from.bndry)), ")"))
p1 <- p1 + geom_point(color = "black")



stop("end")

## TODO: plot GOF vs constrained/unconstrained

## TODO:
## 2.  fimm vs sage dss
##     2a.  of those on axes, how many differed in (1)

## plot(fimm.sage.ohsu$DSS.constrained, fimm.sage.ohsu$fimm.dss, xlab="Sage DSS Calculation", ylab="FIMM DSS Calculation")
## plot(fimm.sage.ohsu$DSS.unconstrained, fimm.sage.ohsu$fimm.dss, xlab="Sage DSS Calculation", ylab="FIMM DSS Calculation")

## str <- "P_4074_16.00950"

## Parse out the replicate number from the ID

## Parse out the patient id and lab id from the column names, like:
## P_118_09.00015

## WORKING HERE

## Plot individual constrained vs unconstrained
t1 <- as.data.frame(tbl)
t2 <- as.data.frame(tbl.constrained)

ind.merged <- merge(t1, t2, by=c("ID"), suffixes=c(".unconstrained", ".constrained"))
for(col in c("DSS.constrained", "DSS.unconstrained")) {
  ind.merged[,col] <- as.numeric(as.character(ind.merged[,col]))
}

plot(ind.merged$DSS.unconstrained, ind.merged$DSS.constrained, xlab="Unconstrained DSS", ylab="Constrained DSS")

## Merge the individual plots and those analyzed following merged
t1 <- as.data.frame(tbl)
t2 <- as.data.frame(tbl2)
m <- merge(t1, t2, by=c("PATIENT_ID", "DRUG_NAME"), suffixes=c(".independent", ".merged"))
m$IC50.independent <- as.numeric(as.character(m$IC50.independent))
m$IC50.merged <- as.numeric(as.character(m$IC50.merged))
## plot(m$ic50.independent, m$ic50.merged)

## plot(log10(m$IC50.independent), log10(m$IC50.merged), xlab="log10 IC50 (Independent Fits of Replicates)", ylab="log10 IC50 (Merged Fits of Replicates)")

## How many are replicated?
## Get the min dose across all drugs/patients.  
tmp <- raw.inhibitor.data

## Define the min concentration for each patient/drug
agg <- aggregate(tmp$well_concentration, by=list(tmp$patient_id, tmp$inhibitor), min)
names(agg) <- c("patient_id", "inhibitor", "min_well_concentration")

tmp.m <- merge(tmp, agg)

## Only keep those with the min concentration
flag <- tmp.m$well_concentration == tmp.m$min_well_concentration
tmp.m <- tmp.m[flag, ]

tmp.m$patient_inhibitor <- apply(tmp.m[,c("patient_id", "inhibitor")], 1, function(row) paste(row[1], row[2], sep="-"))
t <- table(tmp.m$patient_inhibitor)

png("drug-patient-replicate-hist.png")
hist(t[t < 10], xlab="Number of Replicates", main="Histogram of # of Patient/Drug Replicates")
dev.off()

png("independent-vs-merged-ic50-replicates.png")
plot(log10(m$IC50.independent), log10(m$IC50.merged), xlab="log10 IC50 (Independent Fits of Replicates)", ylab="log10 IC50 (Merged Fits of Replicates)")
d <- dev.off()

flag <- (log10(m$IC50.independent) > 10) | (log10(m$IC50.merged) > 10)

## Examine those that have an IC50 > the max concentration tested
ic50.too.large.flag <- (m$IC50.independent > m$Max.Conc.tested.independent)
ic50.too.small.flag <- (m$IC50.independent < m$Min.Conc.tested.independent)
good.ic50.flag <- !ic50.too.large.flag & !ic50.too.small.flag

cat(paste0("Num fits with IC50 > max concentration: ", length(which(ic50.too.large.flag))))
cat(paste0("Num fits with IC50 < min concentration: ", length(which(ic50.too.small.flag))))
cat(paste0("Num fits with min concentration < IC50 < max concentration: ", length(which(good.ic50.flag))))

## Num fits with IC50 > max concentration: 2267
## Num fits with IC50 < min concentration: 1472
## Num fits with min concentration < IC50 < max concentration: 7253

for(indx in which(good.ic50.flag)) {
  id <- m[indx,"ID.independent"]
  cat(paste0(id, " ", m$IC50.independent[indx], "\n"))
  flg <- unlist(lapply(1:length(results), function(i) !is.null(results[[i]]$ID) && (as.character(results[[i]]$ID) == id)))
  if(length(which(flg)) != 1) {
    cat(paste0("Num with matching ids == ", length(which(flg)), "\n"))
    stop("bad id")
  }
  result <- results[[which(flg)]]$result
  png(paste0(output.path, "/", id, "-good-ic50.png"))
  dss <- m[indx,"DSS.independent"]
  ic50 <- m[indx, "IC50.independent"]
  gof <- m[indx, "GOF.independent"]
  plot(result, pcol="grey40", lcol="skyblue1", showInfl=TRUE, cex.main=1.5)
  title(main = paste0(id, "; DSS = ", dss, "; IC50 = ", scientific(ic50, digits=3), "; GOF = ", scientific(gof, digits=3)))
  d <- dev.off()
}

for(indx in which(ic50.too.large.flag)) {
  id <- m[indx,"ID.independent"]
  cat(paste0(id, " ", m$IC50.independent[indx], "\n"))
  flg <- unlist(lapply(1:length(results), function(i) !is.null(results[[i]]$ID) && (as.character(results[[i]]$ID) == id)))
  if(length(which(flg)) != 1) {
    cat(paste0("Num with matching ids == ", length(which(flg)), "\n"))
    stop("bad id")
  }
  result <- results[[which(flg)]]$result
  png(paste0(output.path, "/", id, "-ic50-grt-than-max.png"))
  dss <- m[indx,"DSS.independent"]
  ic50 <- m[indx, "IC50.independent"]
  gof <- m[indx, "GOF.independent"]
  plot(result, pcol="grey40", lcol="skyblue1", showInfl=TRUE, cex.main=1.5)
  title(main = paste0(id, "; DSS = ", dss, "; IC50 = ", scientific(ic50, digits=3), "; GOF = ", scientific(gof, digits=3)))
  d <- dev.off()
}

for(indx in which(ic50.too.small.flag)) {
  id <- m[indx,"ID.independent"]
  cat(paste0(id, " ", m$IC50.independent[indx], "\n"))
  flg <- unlist(lapply(1:length(results), function(i) !is.null(results[[i]]$ID) && (as.character(results[[i]]$ID) == id)))
  if(length(which(flg)) != 1) {
    cat(paste0("Num with matching ids == ", length(which(flg)), "\n"))
    stop("bad id")
  }
  result <- results[[which(flg)]]$result
  png(paste0(output.path, "/", id, "-ic50-less-than-min.png"))
  dss <- m[indx,"DSS.independent"]
  ic50 <- m[indx, "IC50.independent"]
  gof <- m[indx, "GOF.independent"]
  plot(result, pcol="grey40", lcol="skyblue1", showInfl=TRUE, cex.main=1.5)
  title(main = paste0(id, "; DSS = ", dss, "; IC50 = ", scientific(ic50, digits=3), "; GOF = ", scientific(gof, digits=3)))
  d <- dev.off()
}


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

## TODO
## describe GOF on slide
## describe DSS(1-3) on slide
## plot density of GOF (of good guys) -- independent or merged?
## plot density of DSS (of good guys)
f <- good.ic50.flag & (m$DSS.independent > 0)
png(paste0(output.path, "/", "dss-density", ".png"))
plot(density(m$DSS.independent[f]), main="DSS (DSS > 0; min conc < IC50 < max conc)")
d <- dev.off()

## Fit gaussian to peak around 1; exclude those less than 2 stds
png(paste0(output.path, "/", "dof-density", ".png"))
plot(density(m$GOF.independent[f]), main="GOF (DSS > 0; min conc < IC50 < max conc)")
d <- dev.off()

find.density.max <- function(x) {
  d <- density(x)
  y.max <- max(d$y)
  x.max <- d$x[which(d$y == y.max)[1]]
  return(c(x.max, y.max))
}

reflect.data.greater.than <- function(x, reflection.pt) {
  x.grt.or.eq <- x[x >= reflection.pt]
  x.grt <- x.grt.or.eq[x.grt.or.eq > reflection.pt]
  x.less <- reflection.pt - (x.grt - reflection.pt)
  return(c(x.less,x.grt.or.eq))
}

x.vals <- m$GOF.independent[f]

density.max <- find.density.max(x.vals)
x.max <- density.max[1]
y.max <- density.max[2]


if(FALSE) {
  data.reflected <- reflect.data.greater.than(x.vals, x.max)
  mu <- mean(data.reflected)
sigma <- sd(data.reflected)
plot(density(x.vals, n=10000, bw=0.02))
x.reflected <- unique(data.reflected)
uniq.x.vals <- sort(unique(x.vals))
y <- dnorm(x = uniq.x.vals, mean = mu, sd = sigma)
y <- y * y.max / (max(y))
lines(x=uniq.x.vals, y=y)

beta.fit <- fitdistr(x.vals[x.vals > (mu - 3 * sigma)],"beta",list(shape1=1,shape2=1)) 
plot(density(x.vals, n=10000, bw=0.02))
x.reflected <- unique(data.reflected)
uniq.x.vals <- sort(unique(x.vals))
y <- dnorm(x = uniq.x.vals, mean = mu, sd = sigma)
y <- y * y.max / (max(y))
lines(x=uniq.x.vals, y=y)
}

# Define the density 95% confidence interval around GOF=0.5
delta <- 0.1
dens <- density(x.vals)
## y.subset <- d$y[(d$x < (0.5+delta)) & (d$x > (0.5-delta))]
range.flag <- (d$x < (0.6+delta)) & (d$x > (0.6-delta))
x.subset <- dens$x[range.flag]
y.subset <- dens$y[range.flag]
plot(density(y.subset))
## q <- quantile(y.subset, probs=c(0.05,0.5,0.95))
## plot(density(y.subset))

png(paste0(output.path, "/", "dof-density-limits", ".png"))
mu <- mean(y.subset)
sigma <- sd(y.subset)
indx <- which(dens$y==min(dens$y[dens$y>(mu+3*sigma)]))[1]
min.gof <- dens$x[indx]
plot(density(x.vals, n=10000, bw=0.02), main=paste0("Goodness of Fit (threshold = ", round(min.gof, digits=3), ")"))
abline(h=(mu - 3*sigma), lty=3)
abline(h=(mu + 0*sigma), lty=1)
abline(h=(mu + 3*sigma), lty=3)

abline(v=min.gof)
d <- dev.off()

good.gof.flag <- (m$DSS.independent > 0) & (m$IC50.independent < m$Max.Conc.tested.independent) & (m$IC50.independent > m$Min.Conc.tested.independent) & (m$GOF.independent > min.gof)

png(paste0(output.path, "/", "independent-vs-merged-ic50-replicates-filtered.png"))
plot(log10(m$IC50.independent[good.gof.flag]), log10(m$IC50.merged[good.gof.flag]), main=paste0("DSS > 0, min conc < IC50 < max conc; GOF > ", round(min.gof, digits=3)), xlab="log10 IC50 (Independent Fits of Replicates)", ylab="log10 IC50 (Merged Fits of Replicates)")
d <- dev.off()

## Plot some good fits ...

## ... and some bad fits ...

## Now exclude the bad fits and refit the merged
raw.inhibitor.data.good <- raw.inhibitor.data

patient.drug.tbl.good <- unique(raw.inhibitor.data[,c("patient_id","inhibitor", "replicant", "lab_id", "time_of_read")])
patient.drug.tbl.good$ID <- rep(NA, nrow(patient.drug.tbl.good))
n.iters <- nrow(patient.drug.tbl.good)
for(i in 1:n.iters) {
  patient.id <- patient.drug.tbl$patient_id[i]
  inhibitor <- patient.drug.tbl$inhibitor[i]
  replicant <- patient.drug.tbl$replicant[i]
  lab_id <- patient.drug.tbl$lab_id[i]
  time_of_read <- patient.drug.tbl$time_of_read[i]
  
  patient.drug.tbl.good$ID[i] <- paste(patient.id, inhibitor, replicant, lab_id, time_of_read, sep="-")
}
  

patient.drug.tbl.good <- patient.drug.tbl.good[as.character(patient.drug.tbl.good$ID) %in% as.character(m$ID.independent[good.gof.flag]),]
  
n.iters <- nrow(patient.drug.tbl.good)
cat(paste0("Total iterations: ", n.iters, "\n"))
## n.iters <- 1
## n.iters <- 20
## do -> dopar
results2.refit <- list()
tbl2.refit <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
  # for(i in 1:nrow(patient.drug.tbl)) {
  if( ( i %% 10 ) == 0 ) { cat(paste0("Iter ", i, "\n")) }
  patient.id <- patient.drug.tbl.good$patient_id[i]
  inhibitor <- patient.drug.tbl.good$inhibitor[i]
  
  file <- paste(patient.id, inhibitor, sep="-")
  ## file <- paste(file, ".png", sep="")
  
  flag <- (raw.inhibitor.data.good$patient_id == patient.id) & (raw.inhibitor.data.good$inhibitor == inhibitor)
  df <- raw.inhibitor.data.good[flag,,drop=F]
  ic50 <- NULL
  if(use.nplr) {
    r <- nplr.calc.ic50(df$well_concentration, df$normalized_viability, file, patient.id, inhibitor)
    row <- r[["tbl"]]
    results2.refit[[i]] <- list(ID=row$ID, result=r[["res"]])
  } else {
    ic50 <- calc.ic50(df$well_concentration, df$normalized_viability, file)
    row <- c(patient.id, inhibitor, ic50)
  }
  row
}

if(use.nplr == FALSE) {
  rownames(tbl2.refit) <- NULL
  colnames(tbl2.refit) <- c("patient.id", "inhibitor", "ic50")
}

cat("Saving tbl2.good\n")
save.image(file=".RData.tbl2.good")
cat("Saved tbl2.good\n")



## Recalculate merged that exclude the good guys
## (Re-plot) independent vs merged
## Merge the individual plots and those analyzed following merged

## Only keep good fits
tbl.flag <- ( as.character(tbl$ID) %in% as.character(m$ID.independent[good.gof.flag]) ) 
t1 <- as.data.frame(tbl[tbl.flag,])

## Now require there be duplicates
patient.drug <- unlist(apply(t1[,c("PATIENT_ID", "DRUG_NAME")], 1, function(row) paste0(row[1], "-", row[2])))
t1 <- t1[duplicated(patient.drug, fromLast=TRUE) | duplicated(patient.drug, fromLast=FALSE),]

t2 <- as.data.frame(tbl2.refit)
m.refit <- merge(t1, t2, by=c("PATIENT_ID", "DRUG_NAME"), suffixes=c(".independent", ".merged"))
m.refit$IC50.independent <- as.numeric(as.character(m.refit$IC50.independent))
m.refit$IC50.merged <- as.numeric(as.character(m.refit$IC50.merged))
## plot(m$ic50.independent, m$ic50.merged)

png(paste0(output.path, "/", "independent-vs-merged-ic50-replicates-refit.png"))
plot(log10(m.refit$IC50.independent), log10(m.refit$IC50.merged), main=paste0("DSS > 0, min conc < IC50 < max conc; GOF > ", round(min.gof, digits=3)), xlab="log10 IC50 (Independent Fits of Replicates)", ylab="log10 IC50 (Merged Fits of Replicates)")
d <- dev.off()

png(paste0(output.path, "/", "independent-vs-merged-ic50-replicates-refit-zoom.png"))
plot(log10(m.refit$IC50.independent), log10(m.refit$IC50.merged), main=paste0("DSS > 0, min conc < IC50 < max conc; GOF > ", round(min.gof, digits=3)), xlab="log10 IC50 (Independent Fits of Replicates)", ylab="log10 IC50 (Merged Fits of Replicates)", xlim=c(-3,2), ylim=c(-3,2))
d <- dev.off()


## Cases in which independent is poor, but merged is good?

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
k <- kde2d(m.refit$DSS.independent, m.refit$GOF.independent, n=200)
png(paste0(output.path, "/", "dss-gof-density.png"))
image(k, col=r, xlab="DSS", ylab="GOF")
d <- dev.off()

k <- kde2d(m.refit$DSS.merged, m.refit$GOF.merged, n=200)
png(paste0(output.path, "/", "dss-gof-density-merged.png"))
image(k, col=r, xlab="DSS", ylab="GOF")
d <- dev.off()

good.flag <- (m.refit$GOF.independent > 0.8) & (m.refit$GOF.merged > 0.8) & (m.refit$GOF.independent != m.refit$GOF.merged)
indices <- which(good.flag)
indices <- indices[1:min(length(indices),10)]
for(indx in indices) {
  id <- m.refit[indx,"ID.merged"]
  cat(paste0(id, " ", m.refit$IC50.merged[indx], "\n"))
  flg <- unlist(lapply(1:length(results2.refit), function(i) !is.null(results2.refit[[i]]$ID) && (as.character(results2.refit[[i]]$ID) == id)))
  if(length(which(flg)) != 1) {
    cat(paste0("Num with matching ids == ", length(which(flg)), "\n"))
    stop("bad id")
  }
  result <- results2.refit[[which(flg)]]$result
  png(paste0(output.path, "/", id, "-good-refit-ic50.png"))
  dss <- m.refit[indx,"DSS.merged"]
  ic50 <- m.refit[indx, "IC50.merged"]
  gof <- m.refit[indx, "GOF.merged"]
  plot(result, pcol="grey40", lcol="skyblue1", showInfl=TRUE, cex.main=1.5)
  title(main = paste0(id, "; DSS = ", dss, "; IC50 = ", scientific(ic50, digits=3), "; GOF = ", scientific(gof, digits=3)))
  d <- dev.off()
}




library(MASS)
f <- good.ic50.flag & m$DSS.independent > 20 & m$GOF.independent > 0
k <- kde2d(m$DSS.independent[f], m$GOF.independent[f], n=200)
image(k, col=r)

## plot IC50 vs GOF
## plot DSS vs GOF

lab_id <- patient.drug.tbl$lab_id[i]
time_of_read <- patient.drug.tbl$time_of_read[i]

file <- paste(patient.id, inhibitor, replicant, lab_id, time_of_read, sep="-")

## vacuum
## read
## shop
## correlation for gary