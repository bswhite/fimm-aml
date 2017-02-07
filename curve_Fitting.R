#.libPaths("");setwd('C:\\projects\\curveHTML'); source("curvefitfunc.R"); library(drc); library(caTools); library(ggplot2); library(plotly); library(gsubfn); library(openxlsx); library(xtable); library(parallel); library(doSNOW); drug_annot_tbl_full = read.xlsx('FA2A annotations.xlsx'); DSS_typ = 2; HTMLreport = T; shortReport = T

require(drc); require(caTools); require(ggplot2); require(gsubfn); library(data.table)
options(stringsAsFactors = F)
dss = compiler::cmpfun(dss)
headerCurvePath <- "/projects/breeze/code/DSRT3/headercurve.txt"


################################################################################################### 
# outlier removal function.

outlier_remove <-  compiler::cmpfun(function(x){
  qq <- unname(quantile(x, probs=c(.25, .75), na.rm = T))
  outlier_detector <- 1.5 * IQR(x, na.rm = T)
  x[x < (qq[1] - outlier_detector) | x > (qq[2] + outlier_detector)] <- NA
  x
})  

################################################################################################### 
# Calculate IC50/EC50/DSS. (...) .LIST VERSION OF FUNC. 

CALC_IC50_EC50_DSS <- compiler::cmpfun(function(i, drug_wells_, xpr_tbl, DSS_typ, readoutCTX = F)
{
  TEC50 = ifelse(readoutCTX, "TC50", "EC50"); drug_wells = drug_wells_[i];
  
  #find indices of wells with drugs 
  idx_filt <- as.character(xpr_tbl$ProductId) %in% drug_wells  
  #extract inhib. and viab. for wells with drugs in current plate
  inhibition <- xpr_tbl$inhibition_percent[idx_filt]
  
  # if there are identical values in inhibition, add a bit noise
  if(any(duplicated(inhibition))) inhibition <- seq(from = 1, length.out = length(inhibition), by = 0.0001) * inhibition; 
  viability = 100-inhibition;
  
  # extract concentrations, unique drug names and product ids for wells with drugs in current plate
  dose <- as.numeric(xpr_tbl$Concentration[idx_filt])
  drug_name <- unique(as.character(xpr_tbl$ProductName)[idx_filt])
  product_id <- unique(as.character(xpr_tbl$ProductId)[idx_filt])
  
  #combine the data and sort by dose.
  mat_tbl <- data.frame(inhibition,dose,logconc = log10(dose),viability)
  mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]  
  
  ############################# 
  #############    IC50
  
  estimate_param <- tryCatch({drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drc::drmc(errorm = F))}, 
                             warning=function(w){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                             error=function(e){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
  # (extract and name coefficients)
  coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
  coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
  
  # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
  # similar to previous step but now compare log10(IC50) with log(min. conc.).
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
  # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
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
  run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
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
  if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
  
  #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
  nls_result_ic50 <- tryCatch({
    nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", 
        start=list(SLOPE=1,MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
        lower=list(SLOPE=0,MIN=0,MAX=max_lower,IC50=min(mat_tbl$logconc)),
        upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
        control=list(warnOnly=T,minFactor = 1/2048))
  }, error = function(e) {
    
    minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                      start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                      lower=c(SLOPE=0.5, MIN=min_lower,MAX=100,  IC50=min(mat_tbl$logconc)),
                      upper=c(SLOPE=2.5, MIN=coef_estim["MIN"],MAX=100, IC50=max(mat_tbl$logconc)))
  })
  
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
  
  
  ############################# 
  #############    DSS
  
  dss_score <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_typ)); tmp = coef_ic50;
  coef_ic50 <- c(tmp,Max.Conc.tested=max_signal,Min.Conc.tested=min_signal,DSS=dss_score,IC50_std_error=ic50std_Error,S_est = ic50std_resid)
  names(tmp) <- c(TEC50, "SLOPE", "MAX", "MIN") # change IC50 to EC/TC
  coef_ec50 <- c(tmp,Max.Conc.tested=max_signal,Min.Conc.tested=min_signal,DSS=dss_score,EC50_std_error=ic50std_Error,S_est = ic50std_resid)
  
  #inhib. and viab. values in order of dose growth
  inhibition_values <- t(matrix(mat_tbl[,"inhibition"],dimnames=
                                  list(paste0(rep("D", length(mat_tbl[,"inhibition"])), 1:length(mat_tbl[,"inhibition"])))))
  
  
  #dataframe for IC50
  IC50_dataframe <- data.frame(
    ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME="IC50",t(as.matrix(coef_ic50)),inhibition_values,GRAPH=NA)
  
  #round by 2 dex. all the numbers included in range from 4 to N-1 columns
  numeric_cols <- sapply(IC50_dataframe, is.numeric)
  IC50_dataframe[,numeric_cols] <- round(IC50_dataframe[,numeric_cols],1)
  
  #dataframe for EC50
  EC50_dataframe <- data.frame(EC50Table=T,
                               ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME=TEC50,t(as.matrix(coef_ec50)),inhibition_values,GRAPH=NA)
  
  #round by 2 dex. all the numeric colums
  numeric_cols <- sapply(EC50_dataframe, is.numeric)
  EC50_dataframe[,numeric_cols] <- round(EC50_dataframe[,numeric_cols],1)
  
  
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
  
  
  # plot EC50
  y = 100 - y; if(readoutCTX) labels = c(100,50,0) else labels = c(0,50,100)
  ecpl <- ggplot2::ggplot(mat_tbl, aes(logconc, viability)) + scale_y_continuous(breaks=c(0,50,100), labels = labels, limits = c(-25, 125)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
    geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
    geom_vline(xintercept = log10(coef_ec50[TEC50]), colour="grey", size = 0.8) + ggtitle(paste0(drug_name," (dss:",round(EC50_dataframe$DSS,1),")\n")) + 
    theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = ifelse(readoutCTX, "% toxicity", "% viability"), x = "conc(nM)")+ 
    geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ec50[TEC50])*0.95, y2=115, text2=ifelse(readoutCTX, "TC50", TEC50)), color="grey", parse=T) +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background =element_rect(fill = "transparent",colour = NA)) 
  
  graphics.off()
  png(filename = file.path(getwd(), "Results", "Curve_fits", "EC50", paste0(product_id, "_EC50_curve_drug.png")),width=190,height=190, bg = "transparent")
  print(ecpl)
  dev.off()  
  
  EC50base64 <- gsub("\r?\n|\r", " ", base64::img(paste0(getwd(),"/Results/Curve_fits/EC50/",paste0(product_id, "_EC50_curve_drug.png"))))

  #return list with 3 nodes - 1 row for IC50 table and 1 row for EC50 table and EC50 image in base64
  cbind(IC50_dataframe, EC50_dataframe, EC50base64)
  #c(IC50_dataframe, EC50_dataframe)
})


################################################################################################### 
# MAIN (BODY)

dirCur = getwd(); load(file.path(dirCur, "RDA", "screen_table.rda")); 
screen_table$ProductId = gsub("/",".",screen_table$ProductId)

# create directories and dictionary
dir.create("./Results/Curve_fits/", recursive=T); dir.create("./Results/DSS/", recursive=T)
dictionary_ = c(); 

# preparation of final tables
final_tbl_conc_drugs <- drug_annot_tbl_full[,c("Batch_ID_Drug", "DRUG_NAME","Mechanism/Targets","Class.explained")]
colnames(final_tbl_conc_drugs)[colnames(final_tbl_conc_drugs) == 'Batch_ID_Drug'] <- 'ID'

for(screen in unique(screen_table$screen_id))
{
  dir.create("./Results/Curve_fits/IC50/", recursive = !0); dir.create("./Results/Curve_fits/EC50/", recursive = !0);
  
  #Load annotations files for this particular cell line.
  xpr_table <- screen_table[as.character(screen_table$screen_id) %in% screen,]
  readoutCTX = toupper(xpr_table$readout[[1]]) %in% c("CTX", "CTXG");   TEC50 = ifelse(readoutCTX, "TC50", "EC50");
  
  #Loop through all plates that make up this screen.
   ic50_ec50_table = data.table:::rbindlist(lapply(unique(xpr_table$Plate), function(plate){ 
    
    #screen table for particular plate
    xpr_tbl <- xpr_table[as.character(xpr_table$Plate) %in% plate,] 
    
    # xpr_tbl$perInh is not NULL in case of PItoDSS pipeline
    if(is.null(xpr_tbl$perInh)){
    # positive and negative data without outliers
    pos_data <- outlier_remove(xpr_tbl$rawIntensity[as.character(xpr_tbl$Content) %in% "pos"])
    neg_data <- outlier_remove(xpr_tbl$rawIntensity[as.character(xpr_tbl$Content) %in% "neg"])
    
    #Calculate percent inhibition and activation
    avg_low <- mean(pos_data,na.rm=T); avg_high <- mean(neg_data,na.rm=T)
    xpr_tbl$inhibition_percent <- ((avg_high-xpr_tbl$rawIntensity)/(avg_high-avg_low))*100
    } else { colnames(xpr_tbl)[colnames(xpr_tbl) == 'perInh'] <- 'inhibition_percent'}

    #average replicates
    xpr_tbl <- xpr_tbl[, c("screen_id", "readout", "Plate", "ProductId", "ProductName", "Concentration", "inhibition_percent")]
    cols_ <- colnames(xpr_tbl)[!grepl("inhibition_percent", colnames(xpr_tbl))] # columns which should be equal to average PI
    X <- as.data.table(xpr_tbl)
    xpr_tbl = as.data.frame(X[,list(inhibition_percent = mean(inhibition_percent)),cols_], stringAsFactors = !1)
    
    #Use only wells with drugs 
    drug_wells <- unique(as.character(xpr_tbl$ProductId)[!(as.character(xpr_tbl$ProductId) %in% c("","BzCl","empty","dmso", "DMSO","cells", NA, "NA"))])
   
    ###########################
    # call parallel calculation of IC50/EC50/DSS
        data.table:::rbindlist(mclapply(seq_along(drug_wells),CALC_IC50_EC50_DSS, drug_wells_ = drug_wells, xpr_tbl = xpr_tbl, 
                                        DSS_typ = DSS_typ, readoutCTX = readoutCTX, mc.cores = parallel::detectCores()), fill = T)
    }), fill = T) 
    
  ic50_ec50_table <- as.data.frame(ic50_ec50_table, stringAsFactors = F)
  
  # extract IC/EC50 dataframes and base64 images for DSS
  IC50_dataframe <- ic50_ec50_table[,c(1:grep("EC50Table",colnames(ic50_ec50_table))-1)]
  EC50_dataframe <- ic50_ec50_table[,c(grep("EC50Table",colnames(ic50_ec50_table))+1):(grep("EC50base64",colnames(ic50_ec50_table))-1)]
  ICDSStable <- cbind(ic50_ec50_table[1], rep(screen,nrow(ic50_ec50_table)),
                      ic50_ec50_table[,c(4, 11, 5, 6, 9, 8 , 10, grep("EC50base64",colnames(ic50_ec50_table)))])
  colnames(ICDSStable)[2] <- "Experiment_id"; colnames(ICDSStable)[10] <- "GRAPH"; 
  final_tbl_conc_drugs = merge(final_tbl_conc_drugs, ICDSStable, by = "ID", all = T, suffixes=c("")); 
  colnames(final_tbl_conc_drugs) <- gsub("NA$", "", colnames(final_tbl_conc_drugs))

  
  # add sDSS columns. and reorder
  EC50_dataframe[,"sDSS"] <- ""; IC50_dataframe[,"sDSS"] <- "";
  
  
  

  #extract data from drug annot. table by ID
  drug_annot_tbl <- drug_annot_tbl_full[drug_annot_tbl_full$Batch_ID_Drug %in% IC50_dataframe$ID,]
  drug_annot_tblIC <- drug_annot_tbl[match(IC50_dataframe$ID, drug_annot_tbl$Batch_ID_Drug),c("Mechanism/Targets","Class.explained","High.phase/Approval.status",
                                                                                            "Label.name","Res..code","Alias","Trade.names","activity.modifier","Active/inactive.in.clinic","Supplier","Supplier.Ref","Solvent","High.conc.(nM)","Plate","DWell")]
  drug_annot_tblEC <- drug_annot_tbl[match(EC50_dataframe$ID, drug_annot_tbl$Batch_ID_Drug),c("Mechanism/Targets","Class.explained","High.phase/Approval.status",
                                                                                              "Label.name","Res..code","Alias","Trade.names","activity.modifier","Active/inactive.in.clinic","Supplier","Supplier.Ref","Solvent","High.conc.(nM)","Plate","DWell")]
  
  ############################# 
  #############    WRITE IC50/EC50 to .xlsx    (_DSRT_analysis_table_Rpipeline.xlsx)
  
  wb = openxlsx::createWorkbook();
  openxlsx::addWorksheet(wb = wb, sheetName = "IC50"); openxlsx::addWorksheet(wb = wb, sheetName = TEC50);
  
  #write IC50 table and add images  
  IC50_dataframe_full <- cbind(IC50_dataframe, drug_annot_tblIC)
  openxlsx::writeDataTable(wb, sheet = "IC50", IC50_dataframe_full, colNames=T)
  G_ = which(colnames(IC50_dataframe_full)=="GRAPH")
  invisible(lapply(1:nrow(IC50_dataframe_full), function (i) {
    img_ <- paste0(dirCur,"/Results/Curve_fits/IC50/",IC50_dataframe_full$ID[i], "_IC50_curve_drug.png");  
    openxlsx::insertImage(wb, sheet = "IC50", file = file.path(img_),
                          width = 700, height = 700, startRow = i+1, startCol = G_, units = "px", dpi = 360)}))
  #change style 
  openxlsx::setColWidths(wb, sheet = "IC50", cols = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, G_, G_+1, G_+2, G_+3, G_+4, G_+5, G_+6, G_+7, G_+8, G_+9, G_+10, G_+11, G_+12, G_+13, G_+14, G_+15), 
                         widths=c(17, 17, 17.5, 7.5, 16, 7, 7, 8, 17.5, 17.5, 6, 6, 6, 6, 6, 26.5, 9, 9, 21, 17, 27, 13, 12, 7, 14, 17, 23, 11, 14, 10, 16), ignoreMergedCells = F)
  openxlsx::setRowHeights(wb, sheet = "IC50", rows = 1:nrow(IC50_dataframe_full)+1, heights = 142.5)
  style_ = openxlsx::createStyle(halign = "center", valign = "center", wrapText = T)
  openxlsx::addStyle(wb, sheet = "IC50", style_, rows = 2:(nrow(IC50_dataframe_full)+1), cols = 1:35, gridExpand = T)
  #add drug annot. table to the right side after images 
  
  if(readoutCTX) colnames(EC50_dataframe)[colnames(EC50_dataframe) == 'EC50_std_error'] <- 'TC50_std_error'
  #write EC50 table and add images
  EC50_dataframe_full <- cbind(EC50_dataframe, drug_annot_tblEC)
  openxlsx::writeDataTable(wb, sheet = TEC50, EC50_dataframe_full, colNames =T)
  G_ = which(colnames(EC50_dataframe_full)=="GRAPH")
  invisible(lapply(1:nrow(EC50_dataframe_full), function (i) {
    img_ <- paste0(dirCur,"/Results/Curve_fits/EC50/",EC50_dataframe_full$ID[i], "_EC50_curve_drug.png");  
    openxlsx::insertImage(wb, sheet = TEC50, file = file.path(img_),
                          width = 700, height = 700, startRow = i+1, startCol = G_, units = "px", dpi = 360)}))
  #change style 
  openxlsx::setColWidths(wb, sheet = TEC50, cols = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, G_, G_+1, G_+2, G_+3, G_+4, G_+5, G_+6, G_+7, G_+8, G_+9, G_+10, G_+11, G_+12, G_+13, G_+14, G_+15), 
                         widths=c(17, 17, 17.5, 7.5, 16, 7, 7, 8, 17.5, 17.5, 6, 6, 6, 6, 6, 26.5, 9, 9, 21, 17, 27, 13, 12, 7, 14, 17, 23, 11, 14, 10, 16), ignoreMergedCells = F)
  openxlsx::setRowHeights(wb, sheet = TEC50, rows = 1:nrow(EC50_dataframe_full)+1, heights = 142.5)
  style_ = openxlsx::createStyle(halign = "center", valign = "center",  wrapText = T)
  openxlsx::addStyle(wb, sheet = TEC50, style_, rows = 2:(nrow(IC50_dataframe_full)+1), cols = 1:35, gridExpand = T)
  openxlsx::saveWorkbook(wb, file = paste0(dirCur,"/Results/Curve_fits/", screen,"_DSRT_analysis_table_Rpipeline.xlsx"), overwrite = T)
  saveRDS(IC50_dataframe_full,file = paste0(dirCur,"/Results/Curve_fits/", screen,"_DSRT_analysis_table_Rpipeline.rds"))
 
   ############################# 
  #############    short report 
  
  if(shortReport) 
  {
    wb = openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb = wb, sheetName = TEC50); openxlsx::addWorksheet(wb = wb, sheetName = "IC50");
    
    #write IC50 table and add images
    IC50_dataframe <- IC50_dataframe[,c("ID","DRUG_NAME","IC50", "GRAPH", "DSS", "sDSS")]
    openxlsx::writeDataTable(wb, sheet = "IC50", cbind(IC50_dataframe, drug_annot_tblIC), colNames =T)
    G_ = which(colnames(IC50_dataframe)=="GRAPH")
    invisible(lapply(1:nrow(IC50_dataframe), function (i) {
      img_ <- paste0(dirCur,"/Results/Curve_fits/IC50/",IC50_dataframe$ID[i], "_IC50_curve_drug.png");  
      openxlsx::insertImage(wb, sheet = "IC50", file = file.path(img_),
                            width = 700, height = 700, startRow = i+1, startCol = G_, units = "px", dpi = 360)}))
    openxlsx::setColWidths(wb, sheet = "IC50", cols = 1:20, 
                           widths=c(17, 17, 10.5, 26.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5), ignoreMergedCells = F)
    openxlsx::setRowHeights(wb, sheet = "IC50", rows = 1, heights = 43)
    openxlsx::setRowHeights(wb, sheet = "IC50", rows = 1:nrow(IC50_dataframe)+1, heights = 142.5)
    
    style_ = openxlsx::createStyle(halign = "center", valign = "center", wrapText = T)
    openxlsx::addStyle(wb, sheet = "IC50", style_, rows = 1:(nrow(IC50_dataframe)+1), cols = 1:21, gridExpand = T, stack = T)
    openxlsx::addStyle(wb, sheet = "IC50", openxlsx::createStyle(fgFill = "#C0C0C0"), rows = 1, cols = 7:21, gridExpand = T, stack = T)
    openxlsx::addStyle(wb, sheet = "IC50", openxlsx::createStyle(valign = "center"), rows = 1, cols = 1:21, gridExpand = T, stack = T)
    
    
    #write EC50 table and add images  and add drug annot. table to the right side after images 
    EC50_dataframe <- EC50_dataframe[,c("ID","DRUG_NAME",TEC50, "GRAPH", "DSS", "sDSS")]
    openxlsx::writeDataTable(wb, sheet = TEC50, cbind(EC50_dataframe, drug_annot_tblEC), colNames =T)
    G_ = which(colnames(EC50_dataframe)=="GRAPH")
    invisible(lapply(1:nrow(EC50_dataframe), function (i) {
      img_ <- paste0(dirCur,"/Results/Curve_fits/EC50/",EC50_dataframe$ID[i], "_EC50_curve_drug.png");  
      openxlsx::insertImage(wb, sheet = TEC50, file = file.path(img_),
                            width = 700, height = 700, startRow = i+1, startCol = G_, units = "px", dpi = 360)}))
    #change style 
    openxlsx::setColWidths(wb, sheet = TEC50, cols = 1:20, 
                           widths=c(17, 17, 10.5, 26.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5), ignoreMergedCells = F)
    openxlsx::setRowHeights(wb, sheet = TEC50, rows = 1, heights = 43)
    openxlsx::setRowHeights(wb, sheet = TEC50, rows = 1:nrow(EC50_dataframe)+1, heights = 142.5)
    style_ = openxlsx::createStyle(halign = "center", valign = "center", wrapText = T)
    openxlsx::addStyle(wb, sheet = TEC50, style_, rows = 1:(nrow(EC50_dataframe)+1), cols = 1:21, gridExpand = T, stack = T)
    openxlsx::addStyle(wb, sheet = TEC50, openxlsx::createStyle(fgFill = "#C0C0C0"), rows = 1, cols = 7:21, gridExpand = T, stack = T)
    openxlsx::addStyle(wb, sheet = TEC50, openxlsx::createStyle(valign = "center"), rows = 1, cols = 1:21, gridExpand = T, stack = T)
    openxlsx::saveWorkbook(wb, file = paste0(dirCur,"/Results/Curve_fits/", screen,"_DSRT_analysis_table_Rpipeline_short.xlsx"), overwrite = T)
  }

  unlink("./Results/Curve_fits/IC50", recursive=T, force = T)
  unlink("./Results/Curve_fits/EC50", recursive=T, force = T)
  
  dictionary_ = rbind(dictionary_, IC50_dataframe[,c("ID","DRUG_NAME")])
}

############################# 
#############    SAVE _dictionary.xlsx  (Merging_IC50_data_tables_dictionary.R in JP's code) and Final table (used for creation of DSS heatmaps and waterfall)

dictionary_ <- dictionary_[!duplicated(dictionary_[,1]),]
openxlsx::write.xlsx(dictionary_, file = paste0(dirCur,"/Results/Curve_fits/", make.names(paste0("IC50_table_DSS_format_",Sys.Date(),
                                                                                                 "_dictionary.xlsx"))), sheetName = "IC50", colNames = T)
save(final_tbl_conc_drugs, file = "final_tbl_conc_drugs.rds")
# write Original_DSS2
i_ = 1:floor(ncol(final_tbl_conc_drugs)/9) # number of screens
final_tbl_conc_drugs_DSS <- final_tbl_conc_drugs[,c(1,2,3+i_*9)]
colnames(final_tbl_conc_drugs_DSS) <- c("ID", "DRUG_NAME", as.character(unique(na.omit(final_tbl_conc_drugs[,c((i_-1) * 9 + 5)]))))
final_tbl_conc_drugs_DSS <- final_tbl_conc_drugs_DSS[ !apply(final_tbl_conc_drugs_DSS[,c(2, 2+i_)], 1, function(x) all(is.na(x))), ]
openxlsx::write.xlsx(final_tbl_conc_drugs_DSS, file = paste0(dirCur,"/Results/DSS/","Original_DSS",as.character(DSS_typ),
                                                             "_",Sys.Date(),".xlsx"), sheetName = "DSS", startRow = 1)

############################# 
#############    HTML report

if(HTMLreport)
{
  mclapply(i_, function(i){
    final_tbl_conc_drugs[,4+i*9] = gsub("\"","\'", final_tbl_conc_drugs[,4+i*9])
    # extract figure table for certain screen comprising figures in base64 format
    figure_table = final_tbl_conc_drugs[,c(1, 4+9*i)];

    EC50_table <- readRDS(paste0(getwd(),"/Results/Curve_fits/", na.omit(final_tbl_conc_drugs[,-4+9*i])[[1]],
                                             "_DSRT_analysis_table_Rpipeline.rds"))

    #exclude GRAPHs and get them from figure_table matching by id.
    EC50_table = EC50_table[ ,!(colnames(EC50_table) %in% c("GRAPH"))]; EC50_table = merge(EC50_table, figure_table, "ID")
    readoutCTX = toupper(as.character(EC50_table$ANALYSIS_NAME[[1]])) %in% c("TC50");
    
    #turning on/off columns...
    columnsEC = "<div class=\"popka\"><div id=\"colsec50\" class=\"col-sm-8\">Columns: <select class=\"selectpicker\" data-style=\"btn-primary\" multiple show-tick data-actions-box=\"T\">" 
    #...creating multiple selects for it
    namescolsEC <- colnames(EC50_table);
    for(j in 2:length(namescolsEC))
      columnsEC = paste0(columnsEC, " <option selected value = \"",j-1,"\">", namescolsEC[j], "</option>")
    columnsEC = paste0(columnsEC, "</select></div>")
    
    # create table and convert it to HTML (not using html.attributes here because removing spaces later and then classes will be merged.)
    htableEC <- print(xtable::xtable(EC50_table), include.rownames = F, caption.placement = "top", type = "html", print.results = F)
        
    # change table alignment, add classes, remove spaces.
    htableEC <- gsub("align=\"right\"", "", htableEC); 
    htableEC <- sub("<table", "<table id=\"ec50table\" class=\"table table-striped table-bordered nowrap hover\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\" ", htableEC)
    
    htableEC <- gsub('<tr>  <th>', '<thead><th>', htableEC); htableEC <- gsub('</th>  </tr>', '</th></thead>', htableEC);
    htableEC <- gsub('<tr> <th>', '<thead><th>', htableEC); htableEC <- gsub('</th> </tr>', '</th></thead>', htableEC);
    htableEC <- gsub('.tested', '', htableEC)
    htableEC <- gsub("&lt;", "<", htableEC); htableEC <- gsub("&gt;", ">", htableEC)
    
    headercurve = base::readChar(headerCurvePath, file.info(headerCurvePath)$size);
    if(readoutCTX) {headercurve <- gsub('href="#EC50">EC50','href="#TC50">TC50',headercurve); headercurve <- gsub('id="EC50"','id="TC50"',headercurve);}
    
    htable <- paste0(headercurve, paste0("<div id=\"EC50\" class=\"tab-pane fade in active\"><br>", columnsEC, htableEC, "</div></div></div>"), 
                     '</body></html>')
    writeChar(htable, paste0(getwd(),"/Results/Curve_fits/",  na.omit(final_tbl_conc_drugs[,-4+9*i])[[1]],"_DSRT_analysis_table_Rpipeline.html"),
              nchar(htable, type = "chars"))
    
  }, mc.cores = parallel::detectCores());
}

####
# pushing to DB

file_name_Db <- list.files(paste0(getwd(),"/Results"),pattern="Screen_Database.xlsx",recursive=!0,full.names = !0)
if(length(file_name_Db)!=0){
  tryCatch(source('/projects/breeze/code/jmpindi/Curve_fitting/Breeze_Db_curvefit_update.R'), 
           error = function(e) print(e))  	
}

####
# JSON for final report, index.html
all_files = list.files('./Results/Curve_fits'); 
# for multiple cond. all_files[Reduce('&', lapply(c(".xlsx", "pipeline"), grepl, all_files))]
listToJSON <- c(doneCF = T, PipelineFiles = list(gsub(".xlsx", "", all_files[grepl("pipeline.xlsx",all_files)])))
jsonL <- rjson::toJSON(listToJSON)
write(jsonL, "./Results/HTMLreport/CF.json")


#remove all objects
unloadNamespace('openxlsx')
rm(list= ls()[!(ls() %in% c('report_name','writeReport','saving_libpath','DSS_typ','DSRT3'))])
print("After removed...")
