# Do this for read.xlsx2
options(java.parameters = "- Xmx1024m")

suppressPackageStartupMessages(library("synapseClient"))

synapseLogin()

cat("Downloading and reading DR1 data\n")
## # Download the AML cpm-transformed expression data
## dr1.obj <- synGet("syn7433879", downloadFile=TRUE)

# Download the AML count expression data
dr1.obj <- synGet("syn7434014", downloadFile=TRUE)
# Load the data
dr1.tbl <- read.table(getFileLocation(dr1.obj), sep=",", header=TRUE, as.is=TRUE)

cat("Downloading and reading DR2 Q1 data\n")
## # Download the AML cpm-transformed expression data
## dr2q1.obj <- synGet("syn7434019", downloadFile=TRUE)
# Download the AML count expression data
dr2q1.obj <- synGet("syn7434036", downloadFile=TRUE)
dr2q1.tbl <- read.table(getFileLocation(dr2q1.obj), sep=",", header=TRUE, as.is=TRUE)

# Merge the data
cat("Merging DR1 and DR2 Q1 data\n")
all.tbl <- merge(dr1.tbl, dr2q1.tbl)

common.cols <- intersect(colnames(dr1.tbl), colnames(dr2q1.tbl))

if(ncol(all.tbl) != (ncol(dr1.tbl) + ncol(dr2q1.tbl) - length(common.cols))) {
  cat("Problem merging columns\n")
}

if(nrow(all.tbl) != nrow(dr1.tbl)) {
  cat("Did not merge some rows from dr1\n")
}

if(nrow(all.tbl) != nrow(dr2q1.tbl)) {
  cat("Did not merge some rows from dr2q1\n")
}

anno.cols <- c("Gene", "Symbol", "Chr", "Exon_Start", "Exon_End", "Strand", "Length", "GeneStart", "GeneEnd")
rownames(all.tbl) <- all.tbl$Gene
expr.cols <- colnames(all.tbl)[!(colnames(all.tbl) %in% anno.cols)]

dr1.cols <- colnames(dr1.tbl)
dr2q1.cols <- colnames(dr2q1.tbl)
dr1.expr.cols <- dr1.cols[!(dr1.cols %in% anno.cols)]
dr2q1.expr.cols <- dr2q1.cols[!(dr2q1.cols %in% anno.cols)]

suppressPackageStartupMessages(library("stats"))

# Free up some memory
rm(dr1.tbl)
rm(dr2q1.tbl)
dr1.tbl <- NULL
dr2q1.tbl <- NULL

# Transform the count data with limma
## suppressPackageStartupMessages(library("limma"))
## vm <- voom(all.tbl[,expr.cols])
## expr.tbl <- vm$E

# Transform the count data with edgeR (and remove batch effects)
# Follows: https://support.bioconductor.org/p/60581/
# NB: the major difference between this approach and voom above is the
# use of normalizeBetweenArrays in zoom and calcNormFactors here.
# Also, voom does a lot of work calculating weights that will not be
# incorporated into our analysis.

suppressPackageStartupMessages(library("edgeR"))
y <- DGEList(counts=all.tbl[,expr.cols])

# Filter non-expressed genes:
A <- aveLogCPM(y)
y2 <- y[A>1,]

# Then normalize and compute log2 counts-per-million with an offset:
# Gordyn uses prior.count = 5 in the above link--I use the default of 0.5
y2 <- calcNormFactors(y2)
expr.tbl <- cpm(y2, log=TRUE, prior.count=0.5)

## biocLite("devtools")
## suppressPackageStartupMessages(library("devtools"))
## install_github("ggbiplot", "vqv")

suppressPackageStartupMessages(library("ggbiplot"))

plot.pca <- function(data.pca, groups, main="") {
  groups <- as.character(groups)
  if(any(is.na(groups))) {
    groups[is.na(groups)] <- "NA" 
  }
  len <- length(unique(groups))
  if(len > 5) { 
    colors <- rainbow(len)
  } else {
    colors <- c("red","black","blue","green","yellow")[1:len]    
  }
  g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, 
                groups = groups, ellipse = FALSE, 
                circle = FALSE, var.axes = FALSE)
  g <- g + ggtitle(main)
  g <- g + scale_colour_manual(values=colors, name = '')
  # g <- g + scale_color_discrete(name = '')
  g <- g + theme_bw()
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  return(g)
}

# Plot the PCA with labels based on the data release.  There is no
# apparent batch effect.
samples <- rownames(t(expr.tbl))
groups <- rep("NULL", length(samples))
groups[samples %in% dr1.expr.cols] <- "DR1"
groups[samples %in% dr2q1.expr.cols] <- "DR2"

cat("Performing PCA\n")
data.pca.all <- prcomp(t(expr.tbl), center=TRUE, scale=TRUE)
out.file <- "aml-dss-pca-all-no-batch.png"
png(out.file)
g <- plot.pca(data.pca.all, groups, main="Avg(CPM) > 1; No Batch Removal")
print(g)
d <- dev.off()

# Then remove batch correct:
expr.tbl.batch <- removeBatchEffect(expr.tbl, groups)

cat("Performing PCA\n")
data.pca.batch <- prcomp(t(expr.tbl.batch), center=TRUE, scale=TRUE)
out.file <- "aml-dss-pca-all-batch.png"
png(out.file)
g <- plot.pca(data.pca.batch, groups, main="Avg(CPM) > 1; Batch Removal")
print(g)
d <- dev.off()

# Pull in the annotations
## biocLite("xlsx")
suppressPackageStartupMessages(library("xlsx"))

# Get the file Diagnosis_Labs_Treatments_Outcomes_2016_07_19.xlsx
diag1.obj <- synGet("syn7124171", downloadFile=TRUE)

# Get the file 	Diagnosis_Labs_Treatments_Outcomes_2016_10_10.xlsx
diag2.obj <- synGet("syn7436986", downloadFile=TRUE)

# Labs is the 3rd panel in each of these files
cat(paste0("Reading diagnostic file ", getFileLocation(diag1.obj), "\n"))
diag1.labs.tbl <- read.xlsx2(getFileLocation(diag1.obj), sheetIndex=3)

cat(paste0("Reading diagnostic file ", getFileLocation(diag2.obj), "\n"))
diag2.labs.tbl <- read.xlsx2(getFileLocation(diag2.obj), sheetIndex=3)

diag1.labs.tbl <- diag1.labs.tbl[!duplicated(diag1.labs.tbl),]
diag2.labs.tbl <- diag2.labs.tbl[!duplicated(diag2.labs.tbl),]

diag.labs.tbl <- diag2.labs.tbl

flag <- duplicated(diag2.labs.tbl[,c("lab_id", "lab_type", "lab_result")])

# Get the rna-seq "dashboard" files from the beat AML project.

# BeatAML_rnaseq_2016_10_03_public_dashboard.xlsx
db1.obj <- synGet("syn7389565", downloadFile=TRUE)

# BeatAML_DR2_rnaseq_2016_10_03_public_dashboard.xlsx
db2.obj <- synGet("syn7389568", downloadFile=TRUE)

# The samples summary is the 4th sheet
cat(paste0("Reading sample summary file ", getFileLocation(db1.obj), "\n"))
db1.sample.summary.tbl <- read.xlsx2(getFileLocation(db1.obj), sheetIndex=4)

cat(paste0("Reading sample summary file ", getFileLocation(db2.obj), "\n"))
db2.sample.summary.tbl <- read.xlsx2(getFileLocation(db2.obj), sheetIndex=4)

# It appears that these tables only differ in that SpecimenHasSeqCap and PatientHasSeqCap have
# been updated in db2.  Let's drop these columns and merge.
db1.sample.summary.tbl <- db1.sample.summary.tbl[,!(colnames(db1.sample.summary.tbl) %in% c("SpecimenHasSeqCap", "PatientHasSeqCap"))]
db2.sample.summary.tbl <- db2.sample.summary.tbl[,!(colnames(db2.sample.summary.tbl) %in% c("SpecimenHasSeqCap", "PatientHasSeqCap"))]
sample.summary.tbl <- merge(db1.sample.summary.tbl, db2.sample.summary.tbl, all=TRUE)

# Make sure the number of rows are unchanged from the first file and that there are no duplicates.
if(any(duplicated(sample.summary.tbl$SeqID))) {
  cat("Duplicated rows in merged sample summary\n")
}

if(nrow(sample.summary.tbl) != nrow(db1.sample.summary.tbl)) {
  cat("Merged sample summary has more rows than original file\n")
}

# This sample.summary table provides a map from patient id to sequence id.
# Expect that since the sample sequence ids where column headers in our expression
# data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
# So change that.
sample.summary.tbl$transformed.seqid <- unlist(lapply(sample.summary.tbl$SeqID, function(id) gsub(x=paste0("X",id), pattern="-", replacement=".")))

# There are a handful of expression samples for which we don't have clinical annotations.
# Drop samples for which we don't have annotations
flag <- (colnames(expr.tbl) %in% sample.summary.tbl$transformed.seqid)
expr.tbl <- expr.tbl[,flag]

# Read in the drug sensitivity data: interpreted_inhibitor_results_2016_07_19.txt
# inh1.obj <- synGet("syn7124173", downloadFile=TRUE)

# Read in the drug sensitivty data: interpreted_inhibitor_results_2016_10_10.txt
# inh2.obj <- synGet("syn7436989", downloadFile=TRUE)

# cat(paste0("Reading inhibitor file ", getFileLocation(inh1.obj), "\n"))
# inh1.tbl <- read.table(getFileLocation(inh1.obj), sep="\t", header=TRUE, as.is=TRUE)

# cat(paste0("Reading inhibitor file ", getFileLocation(inh2.obj), "\n"))
# inh2.tbl <- read.table(getFileLocation(inh2.obj), sep="\t", header=TRUE, as.is=TRUE)

# Read in the drug sensitivity data interpreted_inhibitor_results_2016_10_10.txt
inh.obj <- synGet("syn7440451", downloadFile=TRUE)
cat(paste0("Reading inhibitor file ", getFileLocation(inh.obj), "\n"))
inh.tbl <- read.table(getFileLocation(inh.obj), sep="\t", header=TRUE, as.is=TRUE)
inh.tbl <- inh.tbl[!duplicated(inh.tbl),]

# Drop any expression columns that are not used in at least one inhibitor study
patients.in.inhibitor.studies <- unique(inh.tbl$patient_id)
flag <- !(patients.in.inhibitor.studies %in% sample.summary.tbl$PatientID)
if(any(flag)) {
  cat(paste0("Could not find annotations for ", length(which(flag)), " patients used in inhibitor studies\n"))
}

flag <- !(sample.summary.tbl$PatientID %in% patients.in.inhibitor.studies)
if(any(flag)) {
  cat(paste0("Could not find inhibitor results for ", length(which(flag)), " patients with annotations\n"))
}

tmp <- unique(sample.summary.tbl[,c("PatientID","transformed.seqid")])
tmp <- tmp[!is.na(tmp$PatientID) & !is.na(tmp$transformed.seqid),]
inh.tbl.expr <- merge(inh.tbl, tmp, by.x="patient_id", by.y="PatientID", all=FALSE)
inh.tbl.expr <- inh.tbl.expr[inh.tbl.expr$transformed.seqid %in% colnames(expr.tbl),]

# Drop expression samples without any drug response data
flag <- (colnames(expr.tbl) %in% inh.tbl.expr$transformed.seqid)
expr.tbl <- expr.tbl[,flag]

# Plot the PCA with labels based on the data release.  There is no
# apparent batch effect.
samples <- rownames(t(expr.tbl))
groups <- rep("NULL", length(samples))
groups[samples %in% dr1.expr.cols] <- "DR1"
groups[samples %in% dr2q1.expr.cols] <- "DR2"

# Drop any genes that have the same expression across all samples
flag <- unlist(apply(expr.tbl, 1, function(row) all(row==row[1])))
expr.tbl <- expr.tbl[!flag,]

# NB: transpose the data so that columns are the features/genes
cat("Performing PCA\n")
data.pca <- prcomp(t(expr.tbl), center=TRUE, scale=TRUE)
out.file <- "aml-dss-pca-harmonized.png"
png(out.file)
g <- plot.pca(data.pca, groups, main="Avg(CPM) > 1; Harmonized with Clinical")
print(g)
d <- dev.off()

# To identify outliers
# biocLite("plotly")
# ggplotly()
# But, really, this doesn't tell us anything more than their positions.

# outliers <- names(which(data.pca$x[,"PC2"] > 220))

# cat("Performing PCA\n")
# expr.tbl.no.outliers <- expr.tbl[,!(colnames(expr.tbl) %in% outliers)]
# if(ncol(expr.tbl) != (ncol(expr.tbl.no.outliers) + length(outliers))) {
#   cat("Did not properly drop outliers")
# }

# # Drop any genes that have the same expression across all samples
# flag <- unlist(apply(expr.tbl.no.outliers, 1, function(row) all(row==row[1])))
# expr.tbl.no.outliers <- expr.tbl.no.outliers[!flag,]
# data.pca.no.outliers <- prcomp(t(expr.tbl.no.outliers), center=TRUE, scale=TRUE)
# samples.no.outliers <- rownames(t(expr.tbl.no.outliers))
# groups.no.outliers <- rep("NULL", length(samples.no.outliers))
# groups.no.outliers[samples.no.outliers %in% dr1.expr.cols] <- "DR1"
# groups.no.outliers[samples.no.outliers %in% dr2q1.expr.cols] <- "DR2"

# out.file <- "aml-dss-pca-rm-out1.png"
# png(out.file)
# g <- plot.pca(data.pca.no.outliers, groups.no.outliers)
# print(g)
# d <- dev.off()

out.file <- "aml-dss-pca-harmonized-var.png"
png(out.file)
sum <- summary(data.pca)
plot(sum$importance[2,1:10], type="l", ylab="Proportion of variance", xlab="Principal Component Index", main="Avg(CPM) > 1; Harmonized with Clinical")
d <- dev.off()

# Do a PCA labeled by specimen type
sample.names <- rownames(data.pca$x)
tbl <- data.frame(sample.name=sample.names)
tbl <- merge(tbl, sample.summary.tbl, by.x="sample.name", by.y="transformed.seqid", all.x=TRUE)
rownames(tbl) <- tbl$sample.name
tbl <- tbl[sample.names,]
if(any(rownames(tbl) != sample.names)) {
  cat("Rows are misaligned with PCA data")
}
groups.specimen.type <- tbl$SpecimenType
print(table(groups.specimen.type))

out.file <- "aml-dss-pca-specimen-type.png"
png(out.file)
g <- plot.pca(data.pca, groups.specimen.type, main="Avg(CPM) > 1; Harmonized with Clinical")
print(g)
d <- dev.off()

# Do a PCA labeled by diagnosis
groups.diagnosis <- tbl$Diagnosis
groups.diagnosis <- as.character(groups.diagnosis)
flag <- grepl(x=groups.diagnosis, pattern="AMBIGUOUS")
groups.diagnosis[flag] <- "ambiguous"
flag <- grepl(x=groups.diagnosis, pattern="AML")
groups.diagnosis[flag] <- "AML"
flag <- grepl(x=groups.diagnosis, pattern="SYNDROMES")
groups.diagnosis[flag] <- "MDS"

print(table(groups.diagnosis))

out.file <- "aml-dss-pca-diagnosis.png"
png(out.file)
g <- plot.pca(data.pca, groups.diagnosis, main="Avg(CPM) > 1; Harmonized with Clinical")
print(g)
d <- dev.off()

# Plot PCA vs inferred sex

# Plot PCA vs FAB
tmp <- unique(sample.summary.tbl[,c("PatientID","transformed.seqid")])
tmp <- tmp[!is.na(tmp$PatientID) & !is.na(tmp$transformed.seqid),]
diag.labs.expr <- merge(diag.labs.tbl, tmp, by.x="patient_id", by.y="PatientID", all=FALSE)

fab.df <- unique(diag.labs.expr[grepl(x=diag.labs.expr$lab_type, pattern="FAB"),c("transformed.seqid","lab_result")])
# Remove any samples with conflicting FAB status
fab.df <- fab.df[!duplicated(fab.df$transformed.seqid, fromLast=TRUE) & !duplicated(fab.df$transformed.seqid, fromLast=FALSE),]

fab.df <- fab.df[(fab.df$lab_result %in% c("NOS","M1","M2","M3","M4","M5","M6","M7")),]
fab.df <- fab.df[(fab.df$transformed.seqid %in% colnames(expr.tbl)),]

cat("Performing PCA\n")
data.pca.fab <- prcomp(t(expr.tbl[,fab.df$transformed.seqid]), center=TRUE, scale=TRUE)
groups = fab.df$lab_result
out.file <- "aml-dss-pca-fab.png"
png(out.file)
g <- plot.pca(data.pca.fab, groups, main="Avg(CPM) > 1; Harmonized with Clinical")
print(g)
d <- dev.off()


cols.to.numeric <- function(x,y) {
  res <- as.character(y)
  flag <- is.na(res) | (res == "n/a")
  x <- x[!flag]
  res <- as.character(res[!flag])
  res <- gsub(x=res, pattern="%", replacement="")
  res <- gsub(x=res, pattern=" ", replacement="")
  
  flag <- suppressWarnings(!is.na(as.numeric(res)))
  res <- as.numeric(res[flag])
  x <- x[flag]
  return(list(x=x, y=res))
}

test.correlation <- function(clin.var, drugs, inh.tbl, diag.labs.tbl) {
  for(drug in drugs) {

    inh.tbl.drug <- inh.tbl[inh.tbl$drug == drug,]
    
    diag.sub <- diag.labs.tbl[grepl(pattern=clin.var, diag.labs.tbl$lab_type),]
    
    m <- merge(inh.tbl.drug, diag.sub, by="patient_id")
    
    aucs <- m$Area_under_the_curve
    l <- cols.to.numeric(m$Area_under_the_curve, m$lab_result)
    
    # plot(l$x, l$y, xlab=clin.var, ylab="AUC")
    lm.fit <- summary(lm(l$y ~ l$x))
    # Is this one coefficient significant?
    r2 <- lm.fit$r.squared
    f <- lm.fit$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    if((p < 0.05) && (r2 > 0.4)) { 
      cat(paste0(drug, " AUC correlated with ", clin.var, ": R^2 = ", r2, " p: ", p, "\n")) 
    }
  }
  
}

test.continuous.correlation <- function(clin.var, drugs, inh.tbl, diag.labs.tbl) {
  for(drug in drugs) {
    
    inh.tbl.drug <- inh.tbl[inh.tbl$drug == drug,]
    
    diag.sub <- diag.labs.tbl[diag.labs.tbl$lab_type == clin.var,]
    
    m <- merge(inh.tbl.drug, diag.sub, by="patient_id")
    
    aucs <- m$Area_under_the_curve
    l <- cols.to.numeric(m$Area_under_the_curve, m$lab_result)
    if(length(l$x) > 5) {
      # plot(l$x, l$y, xlab=clin.var, ylab="AUC")
      lm.obj <- lm(l$y ~ l$x)
      lm.fit <- summary(lm.obj)
      # Is this one coefficient significant?
      r2 <- lm.fit$r.squared
      f <- lm.fit$fstatistic
      if(!is.null(f) && !is.nan(f[1])) {
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        if((p < 0.05) && (r2 > 0.4)) { 
          cat(paste0(drug, " AUC correlated with ", clin.var, ": R^2 = ", r2, " p: ", p, "\n")) 
          main <- paste0(drug, " vs ", clin.var, ": r2 = ", format(r2, digits=2))
          plot(l$x, l$y, xlab=clin.var, ylab="AUC", main=main)
          abline(lm.obj)
        }
      }
    }
  }
}

clin.vars <- unique(diag.labs.tbl$lab_type)
drugs <- unique(inh.tbl$drug)

for(clin.var in clin.vars) {
  skip <- c("Karyotype", "Other Cytogenetics", "Surface Antigens (Immunohistochemical Stains)")
  if(clin.var %in% skip) { next }
  res <- diag.labs.tbl$lab_result[diag.labs.tbl$lab_type==clin.var]
  flag <- suppressWarnings(!is.na(as.numeric(as.character(res))))
  if(any(flag)) {
    cat(paste0("Testing continuous variable ", clin.var, "\n"))
    test.continuous.correlation(clin.var, drugs, inh.tbl, diag.labs.tbl)
  }
}

# Correlate % blast in PB with drug AUC
clin.var <- "Blasts in PB"
test.correlation(clin.var, drugs, inh.tbl, diag.labs.tbl)

# Correlate % blast in BM with drug AUC
clin.var <- "Blasts in BM"
test.correlation(clin.var, drugs, inh.tbl, diag.labs.tbl)

# Correlate WBC with drug AUC
clin.var <- "WBC Count"
test.correlation(clin.var, drugs, inh.tbl, diag.labs.tbl)

# Correlate specimen type with AUC
for(drug in drugs[1:5]) {
  
  inh.tbl.drug <- inh.tbl[inh.tbl$drug == drug,]
  
  lm.obj <- lm(inh.tbl.drug$Area_under_the_curve ~ inh.tbl.drug$specimen_type)
  
  # Is this one coefficient significant?
  r2 <- lm.fit$r.squared
  f <- lm.fit$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  if((p < 0.05) && (r2 > 0.4)) { 
    cat(paste0(drug, " AUC correlated with ", clin.var, ": R^2 = ", r2, " p: ", p, "\n")) 
  }
}


