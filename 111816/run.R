# This file performs fitting of the Beat AML inhibition data to DSS curves

use.nplr <- FALSE

if(use.nplr) {
  # Do fitting with nplr
  library(nplr)
} else {
  # Do fiiting with drc (for initial values) and minpack.lm
  library(drc)
  library(minpack.lm)
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

patient.drug.tbl <- unique(raw.inhibitor.data[,c("patient_id","inhibitor")])

# The Jacobian
j <- function(x, a, b, c, d){
  expr = expression( d + ( a - d ) / ( 1 + ( 10^(b * ( c - x ) ) ) ) )
  retVal = c(eval(D(expr,'a')),
             eval(D(expr,'b')),
             eval(D(expr,'c')),
             eval(D(expr,'d'))
  )
  return(retVal)
}

f <- function(x, a, b, c, d){
  expr = expression( d + ( a - d ) / ( 1 + ( 10^(b * ( c - x ) ) ) ) )
  eval(expr)
}

eval.f <- function(par, x){
  ( par$d + ( par$a - par$d ) / ( 1 + ( 10^(par$b * ( par$c - x ) ) ) ) )
}

fcn.jac <- function(p, y, x, fcall, jcall){	
  -do.call("jcall",c(list(x=x),as.list(p)));
}

# Residual
fcn <- function(p, y, x, fcall, jcall){	
#  y - do.call("fcall", c(list(x = x, as.list(p))))
  y - eval.f(p, x)
}

# i = 10: patient.id = 118; inhibitor = CYT387 is a good fit
# for(i in 1:6) {  
n.iters <- nrow(patient.drug.tbl)
n.iters <- 6
cat("Fitting logistic curve via nplr\n")
tbl <- foreach(i=1:n.iters, .combine=rbind) %do% {
  # for(i in 1:nrow(patient.drug.tbl)) {
  patient.id <- patient.drug.tbl$patient_id[i]
  inhibitor <- patient.drug.tbl$inhibitor[i]
  flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor)
  df <- raw.inhibitor.data[flag,,drop=F]
  model <- tryCatch({
    if(use.nplr) {
      nplr(x = df$well_concentration, y=convertToProp(df$normalized_viability), npars=4)
    } else {
      # Use drc to get the initial parameters.  
      # FIMM function:
      #   y = d + ( a - d ) / ( 1 + ( 10^( b * ( c - x ) ) ) )
      # Asymptotics: x -> infty --> y -> a
      #              x -> -infty --> y -> d
      # In an inhibition setting, FIMM sets b = 1, d (min asymptote) = 0, and determines a and c via the linear regression self starter
      # For both drug response, we should have x -> infty --> y -> a = 0
      # drc function:
      #   y = c' + ( d' - c' ) / ( 1 + ( exp( b' * ( log(x) - log(e') ) ) )
      # Asymptotics: x -> infty --> y -> c'
      #              x -> -infty --> y -> d'
      # Here, 
      # So, we should set b' = b log 10 =  log 10 and c' = 0
      df$x <- (df$well_concentration)
      df$y <- convertToProp(df$normalized_viability)
      drc <- drm(y ~ x, data = df, fct = LL.4(fixed = c(log(10), 0, NA, NA), names = c("b", "c", "d", "e")))
      b.init <- 1
      a.init <- 0
      d.init <- getInitial(drc)[["d"]]
      c.init <- log(getInitial(drc)[["e"]])
      
      start = list(a = a.init, b = b.init, c = c.init, d = d.init)  
      lower <- c(-Inf, -0.1, -0.1, -Inf)
      upper <- c(+Inf, +1.1, +11, +Inf)
      # nlsLM(y ~ d + ( a - d ) / ( 1 + ( 10^(b * ( c - x ) ) ) ), data = df, start = list(a = a.init, b = b.init, c = c.init, d = d.init), fn = fcn, jac = fcn.jac, fcall = f, jcall = j, x = well_concentration, y = normalized_viability)
      df$x <- log10(df$well_concentration)
      nls.lm(par=start,fn = fcn, jac =fcn.jac,fcall=f, jcall=j, x=df$x, y=df$y, lower=lower, upper=upper)
    }
  }, error = function(e) {
    NA
  }, finally = {
    
  })
  row <- c(patient.id, inhibitor, NA, NA, NA, NA, NA, NA, NA)
  if(class(model) == "nplr") {
    #    dss.tbl$bottom[i] <- model@pars$bottom
    #    dss.tbl$top[i] <- model@pars$top
    #    dss.tbl$xmid[i] <- model@pars$xmid
    #    dss.tbl$scal[i] <- model@pars$scal
    #    dss.tbl$max.viability[i] <- max(df$normalized_viability)
    row <- c(patient.id, inhibitor, model@pars$bottom, model@pars$top, model@pars$xmid, model@pars$scal, max(df$normalized_viability), min(df$well_concentration), max(df$well_concentration))
  } else if(class(model) == "nls.lm") {
    row <- c(patient.id, inhibitor, model$par$a, model$par$b, model$par$c, model$par$d, max(df$normalized_viability), min(df$well_concentration), max(df$well_concentration))
  }
  row
}

dss.tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
rownames(dss.tbl) <- NULL
if(use.nplr) {
  colnames(dss.tbl) <- c("patient.id", "inhibitor", "bottom", "top", "xmid", "scal", "max.viability", "min.concentration", "max.concentration")
  # colnames(dss.tbl) <- c("patient.id", "inhibitor", "bottom", "top", "xmid", "scal", "max.viability")
} else {
  colnames(dss.tbl) <- c("patient.id", "inhibitor", "a", "b", "c", "d", "max.viability", "min.concentration", "max.concentration")
  # colnames(dss.tbl) <- c("patient.id", "inhibitor", "bottom", "top", "xmid", "scal", "max.viability")
  
}
save.image(".Rdata.nplr")

#plot(log10(df$x), df$y)
#x <- seq(from=min(log10(df$x)), to=max(log10(df$x)), by=(max(log10(df$x))-min(log10(df$x)))/1000)
coefs <- coef(model)
x <- seq(from=min(min((df$x)), coefs[["c"]]-(max(df$x)-coefs[["c"]])), to=max((df$x)), by=(max((df$x))-min((df$x)))/1000)
#lines(x, eval.f(as.list(coef(model)), x))
#plot((df$x), df$y, xlab="Log10 Concentration", ylab="Proportional response")
plot(x, eval.f(as.list(coef(model)), x), type="l", xlab="Log10 Concentration", ylab="Proportional response")
points((df$x), df$y)

  

if(use.nplr) {
  cols <- c("bottom", "top", "xmid", "scal", "max.viability", "min.concentration", "max.concentration")
} else {
  cols <- c("a", "b", "c", "d", "max.viability", "min.concentration", "max.concentration")
}

# cols <- c("bottom", "top", "xmid", "scal", "max.viability")
for(col in cols) {
  dss.tbl[,col] <- as.numeric(dss.tbl[,col])
}


# system("svn checkout http://svn.code.sf.net/p/dss-calculation/code/trunk dss-calculation-code")
# system("cd dss-calculation-code; tar xvfz DSS_1.2.tar.gz")
# system("cd dss-calculation-code; R CMD build DSS")
# system("cd dss-calculation-code; R CMD INSTALL DSS_1.2.tar.gz")

library("DSS")

# n.iters <- nrow(patient.drug.tbl)
# conc.sum <- foreach(i=1:n.iters, .combine=rbind) %dopar% {
#  # for(i in 1:nrow(patient.drug.tbl)) {
#  patient.id <- patient.drug.tbl$patient_id[i]
#  inhibitor <- patient.drug.tbl$inhibitor[i]
#  flag <- (raw.inhibitor.data$patient_id == patient.id) & (raw.inhibitor.data$inhibitor == inhibitor)
#  df <- raw.inhibitor.data[flag,,drop=F]
#  row <- c(patient.id, inhibitor, min(df$well_concentration), max(df$well_concentration))
#  row
# }

# conc.tbl <- as.data.frame(conc.sum, stringsAsFactors = FALSE)
# rownames(conc.tbl) <- NULL
# colnames(conc.tbl) <- c("patient.id", "inhibitor", "min.concentration", "max.concentration")
# for(col in c("min.concentration", "max.concentration")) {
#  conc.tbl[,col] <- as.numeric(conc.tbl[,col])
# }


# m <- merge(dss.tbl, conc.tbl)

m <- dss.tbl

#rdata should be in the format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
if(use.nplr) {
  rdata <- data.frame(patient.id = m$patient.id, inhibitor = m$inhibitor, IC50 = m$xmid, SLOPE = m$scal * m$max.viability, MIN = max(0, m$bottom * m$max.viability), MAX = m$top * m$max.viability, MIN.Concentration = m$min.concentration, MAX.Concentration = m$max.concentration)
  rdata <- data.frame(patient.id = m$patient.id, inhibitor = m$inhibitor, IC50 = m$xmid, SLOPE = m$scal * m$max.viability, MIN = unlist(lapply(m$bottom * m$max.viability, function(x) max(x, 0))), MAX = m$top * m$max.viability, MIN.Concentration = m$min.concentration, MAX.Concentration = m$max.concentration)
# rdata <- data.frame(patient.id = m$patient.id, inhibitor = m$inhibitor, IC50 = m$xmid, SLOPE = m$scal * m$max.viability, MIN = unlist(lapply(m$bottom, function(x) max(x, 0))), MAX = m$top, MIN.Concentration = m$min.concentration, MAX.Concentration = m$max.concentration)
} else {
  # NB: this is flipped with respect to Tero's slides--look at the form of the equation: as x -> infinity, our curves
  # go to zero, and the equation goes to a.  Hence, a is the min.
  rdata <- data.frame(patient.id = m$patient.id, inhibitor = m$inhibitor, IC50 = c, SLOPE = b, MIN = a, MAX = d, MIN.Concentration = m$min.concentration, MAX.Concentration = m$max.concentration)  
}
rownames(rdata) <- NULL

concn.scale <- 10^-6
log.transform <- TRUE
# activity.threshold <- 10
activity.threshold <- 0.1
activity.threshold <- 10
concn.scale <- 1e-6

cat("Calculating DSS1\n")
dss.input.cols <- c("IC50", "SLOPE", "MIN", "MAX", "MIN.Concentration", "MAX.Concentration")
# rdata$dss.1 <- DSS(rdata[,dss.input.cols], y = activity.threshold, DSS.type = 1, concn.scale = concn.scale, log = log.transform)
rdata$dss.1 <- DSS(rdata[,dss.input.cols], y = activity.threshold, DSS.type = 1, log = FALSE)$DSS

cat("Calculating DSS2\n")
rdata$dss.2 <- DSS(rdata[,dss.input.cols], y = activity.threshold, DSS.type = 2, concn.scale = concn.scale, log = log.transform)

cat("Calculating DSS3\n")
rdata$dss.3 <- DSS(rdata[,dss.input.cols], y = activity.threshold, DSS.type = 3, concn.scale = concn.scale, log = log.transform)

save.image(".Rdata.dss")

dss.cols <- c("dss.1", "dss.2", "dss.3")
dss.cols <- c("dss.1")
ddss.cols <- paste0("d", dss.cols)
for(ddss.col in ddss.cols) {
  rdata[,ddss.col] <- rep(NA, nrow(rdata))
}

drugs <- unique(patient.drug.tbl$inhibitor)
for(drug in drugs) {
  flag <- rdata$inhibitor == drug
  for(dss.col in dss.cols) {
    ddss.col <- paste0("d", dss.col)
    rdata[flag,ddss.col] <- rdata[flag,dss.col] - median(rdata[flag,dss.col])
  }
}

tmp <- rdata[,c("inhibitor", "patient.id", "ddss.1")]
tmp <- tmp[!is.na(tmp$ddss.1),]
tmp <- tmp[!(tmp$ddss.1 == 0),]

flag <- !grepl(pattern="-", x=tmp$inhibitor)
tmp2 <- tmp[flag,]

drugs <- c("EGFR Inhibitor", "Cytarabine", "Sunitinib")
tmp3 <- tmp2[tmp2$inhibitor %in% drugs,]

ac <- acast(tmp3, inhibitor ~ patient.id ~ ddss.1)
ac <- acast(tmp2, inhibitor ~ patient.id )
for(col in colnames(ac))
# ac[is.na(ac)] <- 0
heatmap(ac)