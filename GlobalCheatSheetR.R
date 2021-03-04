#######!/usr/bin/Rscript
options(echo=FALSE)
message("############################################################################################################")
message("##                           R-script to record R commands, syntax and notes                              ##")
message("##                                                                                                        ##")
message("##                                                                                                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                                      ##")
message("##    Date: 2021-03-01                                                                                    ##")
message("##                                                                                                        ##")
message("## Version: 1.0 No additional comments (Test OK on Windows 10 platform)                                   ##")
message("## Require: R v4.0.3                                                                                      ##")
message("##      OS: CentOS7, HMS O2 cluster                                                                       ##")
message("##                                                                                                        ##")
message("## To install & update Bioconductor:                                                                      ##")
message("##    if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")               ##")
message("##    BiocManager::install()                                                                              ##")
message("## To install a package, for example the \"GenomicFeatures\". Run the following commands:                   ##")
message("##    BiocManager::install(c(\"GenomicFeatures\", \"AnnotationDbi\"))                                         ##")
message("## To search available packages programmatically, use the following:                                      ##")                                                                    ##")
message("##    BiocManager::available()                                                                            ##")
message("##                                                                                                        ##")
message("##                                                                                                        ##")
message("## There are several ways to re-direct all outputs into to a file during running R script:                ##")
message("## Methods One: (Using R CMD BATCH)                                                                       ##")
message("## 1.Run your R scripts (test.R as an example here) under Unix or Windows command mode (or under the      ##")
message("##   Terminal environment of RStudio) by giving the following command line                                ##")
message("## 2.$ R CMD BATCH test.R &                                                                               ##")
message("## 3.It will redirect the whole output to a file named \"YourScriptFullNameout\"                            ##")
message("##   (\"test.Rout\" as a result here) and the output file contains all input commands and running results   ##")
message("## Methods Two: (Using sink function)                                                                     ##")
message("## 1.Run the following commands inside the R console                                                      ##")
message("##       con <- file(\"test.log\")                                                                          ##")
message("##       sink(con, append=TRUE)                                                                           ##")
message("##       sink(con, append=TRUE, type=\"message\")                                                           ##")
message("##       source(\"./test.R\", echo=FALSE, max.deparse.length=10000) # if echo=T will add all input commands ##")
message("##       sink()                                                                                           ##")
message("##       sink(type=\"message\")                                                                             ##")
message("##       cat(readLines(\"test.log\"), sep=\"\\n\")                                                             ##")
message("## 2.The output file will or will not contain all input commands depending on your ECHO setting           ##")
message("## Methods Three: (Using Rscript command)                                                                 ##")
message("##  (Rscript ./TMA_TCell_TumorCellSubset.AsFunction.R \"D:/BWHMPE/TMA_TCellPanel/\"                         ##")
message("##   \"D:/BWHMPE/2019-11-11/\") >& ./TMA_TCell_TumorCellSubset.Function.R.log                               ##")
message("##                                                                                                        ##")
message("##                                                                                                        ##")
message("##                                                                                                        ##")
message("## Diagnostic Functions in R:                                                                             ##")
message("## length()       Retrieve or set the dimension of an object                                              ##")
message("## names()        Names of elements within an object                                                      ##")
message("## class()        Retrieves the internal class of an object                                               ##")
message("## mode()         Get or set the type or storage mode of an object                                        ##")
message("## str()          Compactly display the internal structure of an R object                                 ##")
message("## dim()          Retrieve or set the dimension of an object                                              ##")
message("## typeof()       Returns a character string that corresponds to the internal type of the object          ##")
message("## attributes()   Names & Dimensions of Matrices and arrays                                               ##")
message("## sessionInfo()  Print version information about R and attached or loaded packages                       ##")
message("## options()      To set and examine a variety of global options which affect the way of R's computations ##")
message("## storage.mode() *The same with typeof()                                                                 ##")
message("## try() tryCatch() debug()                                                                               ##")
message("##                                                                                                        ##")
message("## rm(list = ls()) Remove all objects existed in the Environment at one time                              ##")
message("## require(parallel, quiet=TRUE)                                                                          ##")
message("## mcores <- detectCores()                                                                                ##")
message("##                                                                                                        ##")
message("############################################################################################################")
message("")
message("# Using print() to print out the outputs while running script in background mode otherwise run the command  ")
message("# directly without putting it inside the print()                                                            ")
options(echo=TRUE)  # see commands in output file
print(Sys.time())
print(R.version)
print(R.home())
print(sessionInfo())
#Using the following 4 command lines to overwrite the system environmental variable
#path <- Sys.getenv("PATH")
#path <- c("C:\\rtools40\\usr\\bin", "C:\\rtools40\\usr\\bin\make.exe", path)
#path <- paste(path,collapse=";")
#Sys.setenv(PATH=path)
#system(sprintf("mkdir %s", snps))
d  <- Sys.Date()
set.seed(12345)
st <- proc.time()
timing <- ( proc.time() - st )
d <- Sys.time()
d <- word(d, 1:2) # function word() needs "stringr"

############################################################################################################
## Packages installation and managements                                                        
############################################################################################################
print(.libPaths())
# .libPaths()
# .libPaths("C:/Program Files/R/R-4.0.3/library")
.libPaths()
# installed.packages()
# install.packages("name-of-your-package", lib="~/R/library")

# The followings show you a special way to install package through "devtools" but need rtools compiler
# library(devtools)
# install_github("vqv/ggbiplot")
# devtools::install_git('http://gitwise.pharma.aventis.com/R-PACKAGER-SANOFI/SRtools.git', ref = '1.3.2', upgrade="never")

# update.packages("spatstat")
# remove.packages("packagename", lib=NULL)
# search()
# nrow(installed.packages())         #Tell you how many packages installed currently
# (.packages())                      #Will list all loaded packages
# packageVersion("spatstat")         #Will list the version of package "spatstat"
# .libPaths("/udd/retwc/R/library/3.1")
# .libPaths( c( .libPaths(), "/udd/retwc/R/library/3.1") )
# detach("package:GGtools" , unload=TRUE)

############################################################################################################
## Working directory and paths                                                         
############################################################################################################
# setwd("D:/BWHMPE/2020-10-19/IF_TMA_CellCompac932/RDA/")
# RDAList <- list.files( pattern="\\.rda$" , full.names = FALSE, recursive = F)
# setwd("C:/Users/E0475408/Documents/Projects/MGH_OLINK_COVID_Oct_12_2020")
# setwd("D:/BWHMPE/2020-10-19")
# dir("../../../data/imputed", full=TRUE) #paste0(dir("../fromPauline", full=T), "TEST")
dir()         # To list the files under the working directory
ls()          # To list the R objects and variables in the current working session
rm(list=ls()) # Clear the workspace
              # CTRL+L to clear the screen of Console

############################################################################################################
## Reading & Writing files; Loading & Saving objects/libraries                                                         
############################################################################################################
# CellSegDataRDA.Path <- "D:/BWHMPE/2019-11-11/LocalLaptop/TotalCellSegData_2019-11-11.rda"
# LinkingPath         <- "D:/BWHMPE/2019-12-23/NHSHPFS TMA linking file final 20181220 MC.csv"
# require("cluster", lib.loc="~/R/library")
# library("stringr")

# source("umap.R") # load and execute a script of R commands
# get(load(object))
# load("D:/BWHMPE/2020-10-19/IF_TMA_CellCompac932Clinical_4178614x89_2020-10-18.rda")
# load("D:/BWHMPE/2019-11-11/LocalLaptop/TotalCellSegData_2019-11-11.rda")
# load("D:/BWHMPE/2020-01-06/TumorID_added_on_TissueTumor_2020-01-06.rda")
# long_model <- readRDS("~/Projects/fromPauline/ResultsR/long_model_INF_pval.Rds")
# ClinicData <- readRDS("./ClinicalData_966x40_2020-01-29.Rds") 
# loading RDS object will return both the name and content but RDA only return name
# readRDS() will directly return content to the assigned name variable

# save(x, file=paste("./Cell932/RDA/Cell932_TumorID_", i, "_", d, ".rda", sep=""))
# save(TissueOther, file=paste("./TissueOther_206_", d, ".rda", sep=""))
# save(TotalTumor,  file=paste("./TotalTumorC_2981441x57_", d, ".rda", sep=""))
# save(ClinicData, file=paste("./ClinicalData_966x40_", d, ".rda", sep="")) 
# save() will save R-objects and their environment information
# saveRDs() only save the R-object without environment information, but it will give the variable name directly
# saveRDS(ClinicData, file=paste("./ClinicalData_966x40_", d, ".Rds", sep="")) 
# savehistory(file="makeSE.hist.txt")

exprPre  <- read.csv("./InputFiles/LV135preX.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
# It means the read.csv function won't convert those special characters such as space or others of column_name into dot character
# (for example "Col 11" -> "Col.11") when you assigned "check.names=FALSE" 
TID <- read.csv("./TumorID_List.csv", header=TRUE, check.names=FALSE, as.is=T)                #--1
TID <- read.csv("./TumorID_List.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE) #--2
   # The 1st & 2nd statements are identical. when "as.is=T" means the entire columns are not treated as factor type
TID <- read.csv("./TumorID_List.csv", header=TRUE, check.names=FALSE, as.is=2)                #--3
   # The 3th statement will treat only the 2nd column as not factor but all other columns as factor type

options(stringsAsFactors = F)                                         #--4
TID <- read.csv("./TumorID_List.csv", header=TRUE, check.names=FALSE) #--5
   # The 4th and 5th statements together will treat all columns as not factor type

write.csv(TotalWeights,    file=paste("./TotalWeights(928Tumors)_", d, ".csv", sep=""), row.names=F)
write.csv(TotalTumorCells, file=paste("./TotalTumorCe_2981441x57_", d, ".csv", sep=""), row.names=FALSE)
write.csv(ClinicalDMeta, file=paste("./ClinicalDataMeta_35x28_", d, ".csv", sep=""), row.names=FALSE)

NPX.original <- read.table("./MGH_COVID_OLINK_NPX.txt", header=T, sep=";")
TempFile <- read.table(filepath, sep="\t", header = TRUE, stringsAsFactors = FALSE, na.strings = "#N/A", comment.char = "")
# The space existed in the column name would be replaced by period symbol. 
# for example, "Tissue Category" would become "Tissue.Category"

write.table(x, file=paste("./Cell932/TAB/Cell932_TumorID_", i, "_", d, ".txt", sep=""), 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(TotalTumorCells, file=paste("./TotalTumorCells_2981441x57_", d, ".txt", sep=""), 
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x, file=paste("./TAB/TTuOther_PatholDiff01_TumorID_", i, "_", d, ".txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

VariablesDes <- read.xlsx("./variable_descriptions.xlsx", 1) # need "rJava,xlsx,xlsxjars" to read the first sheet 


############################################################################################################
## Checking objects, data and types .....                                                        
############################################################################################################
summary(NPX.original)
class(NPX.original)
str(NPX.original)
dim(NPX.original)
table(Clinfo$COVID)
table(Clinfo$COVID, Clinfo$HTN)
print(table(dat$genotype, useNA="ifany"))

length(levels(factor(NPX.original$UniProt))) # 1420 unique proteins
length(levels(factor(NPX.original$SampleID))) # 786
length(levels(factor(NPX.original$subject_id))) # 383
length(levels(factor(NPX.original$Index)))
length(levels(factor(NPX.original$OlinkID)))

head(table(NPX.original$subject_id, NPX.original$Timepoint))
tail(table(NPX.original$subject_id, NPX.original$Timepoint))
SID_Time <- table(NPX.original$subject_id, NPX.original$Timepoint)
SID_Panel <- table(NPX.original$subject_id, NPX.original$Panel)

nrow(NPX.original) / length(levels(factor(NPX.original$SampleID))) #1461.725
nrow(NPX.original) / length(levels(factor(NPX.original$subject_id))) #2999.781

levels(factor(NPX.original$Timepoint))

tail(table(NPX.original$subject_id, NPX.original$Panel))
tail(table(NPX.original$SampleID, NPX.original$Panel))

with(Clinfo, table(COVID, HTN, dnn = c("CO", "pre"))) # Give names "co" to Y-axis and "pre" to X-axis of the table
with(Clinfo, table(COVID, DIABETES, dnn = c("COVID-19", "Pre_existed_Diabetes")))
with(Clinfo, table(COVID, HEART, dnn = c("COVID-19", "Pre_existed_HeartProblem")))

TID_109 <- scale(TID_109, center = TRUE, scale = TRUE) #Perform the Z-score normalization, value-colMeans/colSD
TID_109 <- as.data.frame(TID_109)                      #Convert the dtype from matrix back to data frame




############################################################################################################
## Loops & Iterations                                                         
############################################################################################################
# x <- 10
# while (x > 0) {print(x); x <- x-1;}            #R used ";" to separate the codes
# repeat { print(x); x <- x+1; if (x>10) break}  #R used "{}" to wrap the code block


############################################################################################################
## Statistics                                                    
############################################################################################################
x <- rnorm(10, mean=0, sd=0.5)                 #Randomly pick 10 points from Normal Distribution N(0, 0.5)
x <- seq(from=1, to=10, by=0.2)
fivenum(a)
fivenum(tmp1$columnX)

t.test(tmp1$columnX,tmp2$columnX)

summary(lm(breaks ~ wool + tension, data = warpbreaks))

message("# Setting the seed = 12345 now ")
message("# Setting the MAF threshold as MAF >= 0.15 ")
message("# Setting cis-radius as 100kb ")
message("# Setting the permutation times as 6 to calculate plug-in FDR values ")
message("# Setting the lower bound genotype frequency (lbgtf) as 0.023 to filter out the one-patient minor allele homozygous issue ")
workerFunc <- function(n) { 
  print(table(seqnames(expr.list[[n]][[1]][1:10,])))
  # print(cisCount(summex=expr.list[[n]][[1]], vcf.tf=expr.list[[n]][[2]], lbmaf=0.05, lbgtf=0.01, cisradius=100000))
  # run = cisAssoc(summex=expr.list[[n]][[1]], vcf.tf=expr.list[[n]][[2]], nperm=6, lbmaf=0.05, lbgtf=0.01, cisradius=100000)
  run = cisEsts(summex=expr.list[[n]][[1]], vcf.tf=expr.list[[n]][[2]], nperm=6, lbmaf=0.15, lbgtf=0.023, cisradius=100000)
  run$piFDR = gQTLstats:::pifdr(run$chisq, c(run$permScore_1, run$permScore_2, run$permScore_3, run$permScore_4, run$permScore_5, run$permScore_6))
  message("The amount of tested SNP-GENE association pairs belong to this chromosome is ", length(run))
  message("The amount of successful-FDR-calculations on association pairs belong to this chromosome is ", length(run$piFDR))
  message("Finishing the gQTLstats cis-eQTL on this Chromosome ")
  return(df_fdr.1 <- run[which(run$piFDR<= 1)])
}
res <- mclapply(seq_along(expr.list), workerFunc, mc.cores = mcores)
exprSeQTL_fdr.1  <- do.call(rbind, lapply(res, as.data.frame))
exprSeQTL_fdr.1  <- exprSeQTL_fdr.1[order(exprSeQTL_fdr.1$piFDR), ]
exprSeQTL_fdr.1  <- exprSeQTL_fdr.1[order(exprSeQTL_fdr.1$chisq, decreasing = TRUE), ]
exprSeQTL_fdr.05 <- exprSeQTL_fdr.1[which(exprSeQTL_fdr.1$piFDR<=0.05), ]


############################################################################################################
## Data subset and  selection                                                   
############################################################################################################
dataFrame1$class <- ifelse(dataFrame$class == "A", "A", "B") #Change those values are not "A" to "B"

chosen <- sample(unique(TotalTumorCells$Entire.Cell.DAPI), 500000)
TTumor739844_167  <- subset(TotalTumorCells, Entire.Cell.DAPI %in% chosen)
TTumor2241597_167 <- subset(TotalTumorCells, Entire.Cell.DAPI %nin% chosen)

TissueOther <- subset(TotalCellSegData, Tissue.Category=="Other")
TissueTumor <- subset(TissueTumor, select=-c(col1, col2, col3)) # col1~col3 will be removed
TissueTumorOther <- subset(TissueTumor, Phenotype=="Other")
exprPre  <- subset(exprPre, chr=="chr20" | chr=="chr21" | chr=="chr22" | chr=="chr23")  
phenotMale <- subset(phenotypes, exclude==0 & Sex=="F", select=-c(HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
# The value of column "exclude" equal to zero and value of column "Sex" equal to "F" will be selected 
y1 <- subset(x, Tumor.ID==4 & Phenotype=="CD3- CD45RO+" & Tissue.Category=="Other")
y2 <- subset(IF_TMA_CellCompac932Clinical, Tumor.ID==4 & Phenotype=="CD3- CD45RO+" & Tissue.Category=="Other")
IF_TMA_932_CellProportionNoBigNucleusSize   <- subset(IF_TMA_932_TCellNucBigger748Pixels, BigTProportion == 0.0)
IF_TMA_932_CellProportionWithBigNucleusSize <- subset(IF_TMA_932_TCellNucBigger748Pixels, BigTProportion > 0.0)

IF_TMA_OnlyImmunoCells_932Clinical <- IF_TMA_CellCompac932Clinical[!(IF_TMA_CellCompac932Clinical$Phenotype=="Other"), ]
exprPre  <- exprPre[,  intersect(rownames(phenotypesFeMale), colnames(exprPre))]

TID <- levels(factor(ClinicData$TumorID))
nlevels(factor(TTumorIntersection$Tumor.ID))     #nlevel() return number rather than strings
sum(is.na(TTumorIntersection))
CRCSL_T1N1 <- as.character(na.omit(CRCSurvLTID$TID_T1N1))

tmp  <- a[which(a[,1] == "xxx"), ]   # method 1, 選取第一個column裏頭為xxx的所有rows
tmp1 <-        a[(a[,1]=="xxx"), ]   # method 2, 不加which括號也行
tmp2 <-        a[(a[,1]!="xxx"), ]
tmp$column1 <- NULL # will remove column1 from tmp data frame
TotalCellSegData <- TotalCellSegData[-c(1),] # To remove the first row which is filled with "NA"
identical(x, IF_TMA_CellCompac932Clinical)
TissueTumor <- TissueTumor[,  c(1,41,46:48,2:3,42:44,4:5,45,28:34,6:27,35:40)] #To re-arrange the column order  
ClinicData <- ClinicData[order(ClinicData$TumorID), ] #Sorting the TumorID from small to big (4 -> 3730)

############################################################################################################
## MISC                                                    
############################################################################################################
 grep("xyz", a[,1])    # 列出第一個column中含有"xyz"字元的 row number
grepl("xyz", a[,1])    # 傳回的是 TRUE FALSE, 不再是 row number
TissueTumor$Confidence <- as.numeric(sub("%", "",TissueTumor$Confidence, fixed=TRUE))/100
TID_Dim_List$Dimension <- gsub(', 42)', '', TID_Dim_List$Dimension) #To remove the pattern ", 42)"
TID_Dim_List$Dimension <- gsub('\\(', '', TID_Dim_List$Dimension)   #You need \\ because \ itself is a special char also
colnames(TotalWei.cor) <- gsub('V', 'W', colnames(TotalWeights.cor))
colnames(TotalWeights) <- gsub('TumorID_', '', colnames(TotalWeights))
                          gsub(" ", "", gsub("[^0-9]*:", "", summary(C1)[4,2]), fixed=TRUE)
   GENE$ugene_id  <- paste(GENE$ugene_id, "X", sep="")
     df$ugene_id  <- sub('X$', '', df$ugene_id) # remove the postfixed "X" protection
   rownames(Covs) <- gsub("(B.*V)", "\\1_\\1", rownames(Covs)) # change BxxxxV back to BxxxxV_BxxxxV
colnames(exprPre) <- gsub("(B.*V)_(B.*V)", "\\1_Pre",  colnames(exprPre)) # change BxxxxV_BxxxxV to BxxxxV_Pre , 135 subjects
colnames(exprPos) <- gsub("(B.*V)_(B.*V)", "\\1_Post", colnames(exprPost)) # change BxxxxV_BxxxxV to BxxxxV_Post, 133 subjects

exprZero <- exprPre[(rowSums(exprPre==0.0)  == ncol(exprPre)),]
exprPre  <- exprPre[!(rowSums(exprPre==0.0)  == ncol(exprPre)),]
exprPre  <- exprPre[setdiff(rownames(exprPre),  rownames(exprPreZero)),] 
exprPre  <- exprPre[, intersect(rownames(phenotypesPre),colnames(exprPre))]
exprPre  <- exprPre[rowSums(exprPre>0.1)  > 10,] # 18212 genes left
exprPre  <- log2(exprPre+1.1)


colnames(TotalWeights) <- str_pad(colnames(TotalWeights), 4, pad = "0") #To add leading zero on each tumor ID number in 4 digits length
colnames(TotalWeights)[2:37] <- str_pad(colnames(TotalWeights)[2:37], 2, pad = "0")
colnames(TotalWeights)[2:37] <- paste("W", colnames(TotalWeights)[2:37], sep="")
TotalWeights <- TotalWeights[, order(names(TotalWeights))] #To sort the order of column name
colnames(TotalWeights) <- paste("TumorID_", colnames(TotalWeights), sep="") #To add prefix TumorID back to column name
TotalWeights <- t(TotalWeights) #Transpose the TotalWeights matrix
fudge <- max(TotalWeights)-min(TotalWeights)
sum(TotalWeights < 0)
sum(is.na(TID_Dim_List$Tstage)) #Has 67 NA
phenotypesPre  <- merge(GenoCov120, phenotypesPre,  by="row.names", all.x=TRUE)



message("# The amount of NA (Not Available; missing values) in PRE dataset after log2 transformation is ", sum(is.na(exprPre)))
message("# The amount of NaN (Not a Number; 0/0) in PRE dataset after log2 transformation is ", sum(is.nan(as.matrix(exprPre))))
message("# The amount of Inf (Infinity; caused by N/0) in PRE dataset after log2 transformation is ", sum(is.infinite(as.matrix(exprPre))))
message("# The amount of zero in PRE dataset after log2 transformation is ", sum(exprPre==0, na.rm=TRUE))
message("# The amount of negative values in PRE dataset after log2 transformation is ", sum(exprPre < 0, na.rm=TRUE))
message("# The amount of positive values in PRE dataset after log2 transformation is ", sum(exprPre > 0, na.rm=TRUE))
message("# The Normalized PRE dataset has ", sum(exprPre==max(exprPre), na.rm=TRUE), " Max values and the Max value is ", max(exprPre, na.rm=TRUE))
message("# The Normalized PRE dataset has ", sum(exprPre==min(exprPre), na.rm=TRUE), " min values and the min value is ", min(exprPre, na.rm=TRUE))


categorical_varaibles = c("Sex", "SeqCenter", "Reads")
message("# The identified categorical variables include SEX, Sequencing-Center and Reads-length ")
for(x in categorical_varaibles) {CovsPre  = cbind(CovsPre,  value=1);  CovsPre[,x]=paste0(x,  CovsPre[,x]); CovsPre  = dcast(CovsPre,  as.formula(paste0("... ~ ", x)), fill=0);}

# if(file.exists(".RData")) load(".RData") else{
#     message("# Didn't find the RData file so please load RData ...")}

getMode <- function(x) { #Computing mode of a given vector
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

args <- commandArgs(TRUE)
message("#The first command-line argument will be the path of TMA_TCell source files ")
message("#The second command-line argument will be the path of working directory ")
if (length(args)==0) {
  SourcePath <- "D:/BWHMPE/TMA_TCellPanel/"
  WorkingDIR <- "D:/BWHMPE/2019-11-11/"
  message("#No command-line arguments found, use the default paths #")
  # stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  SourcePath  <- args[1]
  WorkingDIR  <- args[2]
}

SDIR <- "~/Projects/RIPK1i(SAR443122_DNL758)/fromPauline/ResultsR/"
WDIR <- "C:/Users/E0475408/Documents/Projects/RIPK1i(SAR443122_DNL758)/MergingPDY16879_MGH/"
setwd(WDIR)
getwd()
RDsList <- list.files(path=SDIR, pattern="\\.Rds$", include.dirs = TRUE, # pattern is case sensitive, ".rds$" won't get data
                      full.names = TRUE, recursive = F)
for (i in RDsList) {
  assign((tools::file_path_sans_ext(basename(i))), readRDS(i)) # Using assign() to give data to a dynamic variable
  rm(i); rm(f)
}

EXP <- list(df1=exprPre, df2=exprPost)
EXP <- lapply(EXP, function(df) {    
  df$ugene_id <- sub('X$', '', df$ugene_id) # remove the postfixed "X" protection
  df <- df[,-c(2:5)] # remove columns of gene coordinates information
  rownames(df) <- df[,1]
  df <- df[, -1]
  df
})
exprPre  <- EXP$df1
exprPost <- EXP$df2
rm(EXP)

addingTumorID_fromLinking <- function(CellSegData, Linking, FileNameForSaved) {
  #CellSegData <- get(load(CellSegDataRDA.Path))
  #Linking     <- get(read.csv(LinkingPath, stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE))
  #Linking     <- Linking[, -1] #Remove the first column
  #message("#The path of Cell Segment Dataset (in R_Object format): ", CellSegDataRDA.Path)
  #message("#The path of Linking File: ", LinkingPath)
  #message("#The current working directory is ", getwd())
  #CellSegData$TMA.ID <- CellSegData$HTMA.ID <- CellSegData$Tumor.ID <- NULL
  inPutName <- deparse(substitute(CellSegData))
  message("Start running function with input file ", inPutName)
  print(Sys.time())
  
  ND <- CellSegData                                                               #Copy to NewData(ND) file
  for (i in 1:nrow(ND)) {
    ND$TMA.ID[i]  <- i
    ND$HTMA.ID[i] <- i
    firstWord <- word(ND$Sample.Name[i], 1:3)                                     #To get the first 3 words ("HTMA" "251" "Variant")
    firstWord <- na.omit(firstWord)                                               #To remove the hidden missing value (NA) in the second or the third word position
    firstWord <- str_extract(str_flatten(firstWord), "^[[:alpha:]]+[[:digit:]]+") #To merge these words into one word without space in between
    fs <- substr(firstWord,1,1)                                                   #then extract the first word ending with number such as "HTMA251"
    if ( fs != "T") {
      ND$HTMA.ID[i] <- firstWord                                                  #If the value is not started with "T" then save the value to HTMA.ID column
    } else {
      ND$TMA.ID[i]  <- firstWord                                                  #If the value is started with "T" then save the value to TMA.ID column
    }
  }
  ND$TMA.ID  <- gsub("^[[:digit:]]+", "NA", ND$TMA.ID)                            #To convert those serial number in TMA.ID column into NA missing vlaue
  ND$HTMA.ID <- gsub("^[[:digit:]]+", "NA", ND$HTMA.ID)                           #To convert those serial number in HTMA.ID column into NA missing vlaue
  message("Finished creating HTMA.ID & TMA.ID columns from Sample.Name column to dataset ", inPutName)
  print(Sys.time())
  
  for (i in 1:nrow(ND)) {
    if ( ND$HTMA.ID[i] != "NA" ) {                                                #If the HTMA.ID is available then use it to lookup the linking file ....
      coord <- substring(as.character(ND$Sample.Name[i]), regexpr("\\[",as.character(ND$Sample.Name[i]))+1, regexpr("\\]",as.character(ND$Sample.Name[i]))-1)
      coord <- str_c(word(coord, 2:3, sep=fixed(",")), collapse=",")             #To extract the coordinates (such as "1,10,18") and remove the first digit (became "10,18")
      ND$Tumor.ID[i] <- Linking[which((Linking$coordinates==coord) & (Linking$HTMA.ID==ND$HTMA.ID[i])), ][1] #If both coordinates and HTMA.ID matched then 
      ND$TMA.ID[i]   <- Linking[which((Linking$coordinates==coord) & (Linking$HTMA.ID==ND$HTMA.ID[i])), ][6] #add both Tumor.ID and TMA.ID to CellSegDataset
    } else {                                                                   #If the HTMA.ID is not available then use TMA.ID to lookup the linking file ....
      coord <- substring(as.character(ND$Sample.Name[i]), regexpr("\\[",as.character(ND$Sample.Name[i]))+1, regexpr("\\]",as.character(ND$Sample.Name[i]))-1)
      coord <- str_c(word(coord, 2:3, sep=fixed(",")), collapse=",")             #To extract the coordinates (such as "1,10,18") and remove the first digit (became "10,18")
      ND$Tumor.ID[i] <- Linking[which((Linking$coordinates==coord) & (Linking$TMA==ND$TMA.ID[i])), ][1]      #If both coordinates and TMA.ID matched then
      ND$HTMA.ID[i]  <- Linking[which((Linking$coordinates==coord) & (Linking$TMA==ND$TMA.ID[i])), ][4]      #add both Tumor.ID and HTMA.ID to CellSegDataset
    }
  }
  message("Finished adding Tumor.ID columns onto input Cell-Seg-Dataset ", inPutName)
  print(Sys.time())
  ND <- as.data.frame(lapply(ND, unlist))
  ND <- ND[order(ND$Tumor.ID), ]
  tumorsN <- nlevels(factor(ND$Tumor.ID))
  totalcells <- nrow(ND)
  message("There are totally ", totalcells, " cells belong to ", tumorsN, " tumor specimens")
  message("There are totally ", sum(is.na(ND$Tumor.ID)), " missing values in Tumor.ID column")
  #return(ND)
  print(ND[1:20, c(2:3, 201:203, 207:209)])
  save(ND, file=paste("./TumorID_added_on_", FileNameForSaved, "_", d, ".rda", sep=""))
  message("Finished adding Tumor.ID to each cell record and saving results to file postfixed ", FileNameForSaved)
  message("Finished running the R function with input cell-seg-dataset ", inPutName)
  print(Sys.time())
}

message("#Here is the current working directory ", getwd())
message("#Loading input files now ")
print(Sys.time())
args <- commandArgs(TRUE)
message("#The 1st command-line argument will be the linking file ")
message("#The 2nd command-line argument will be the 1st input ")
message("#The 3rd command-line argument will be the 2nd input ")
message("#The 4th command-line argument will be the 3rd input ")
if (length(args)==0) {
  LinkingPath <- "/home/tac15/Temp/2019_12_27/NHSHPFS TMA linking file final 20181220 MC.csv"
  input_1     <- "/home/tac15/Temp/2019_12_27/TissueOther_206_2019-12-26.rda"
  input_2     <- "/home/tac15/Temp/2019_12_27/TissueStroma_206_2019-12-26.rda"
  input_3     <- "/home/tac15/Temp/2019_12_27/TissueTumor_206_2019-12-26.rda"
  message("#No command-line arguments found, use the default paths #")
  # stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==4) {
  LinkingPath  <- args[1]
  input_1      <- args[2]
  input_2      <- args[3]
  input_3      <- args[4]
}

message("#Finished loading input files now ")
print(Sys.time())
LinkingFile  <- read.csv(LinkingPath, stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
LinkingFile  <- LinkingFile[, -1]
load(input_1)
load(input_2)
load(input_3)

message("#Start running the adding_Tumor.ID_function with 3 input files respectively now")
print(Sys.time())
addingTumorID_fromLinking(TissueOther,  LinkingFile, FileNameForSaved="TissueOther")
addingTumorID_fromLinking(TissueStroma, LinkingFile, FileNameForSaved="TissueStroma")
addingTumorID_fromLinking(TissueTumor,  LinkingFile, FileNameForSaved="TissueTumor")

message("#Finished running the R function with all input files now ")
print(Sys.time())
timing <- ( proc.time() - st )
message("#Finishing The entire R Script with running time (in unit of seconds) listed in the following lines ")
print(timing)
timestamp()





############################################################################################################
## Plot & Graph                                                       
############################################################################################################
?par ?pch                                      #Set or Query Graphical Parameters
colors()                                       #To list all color names adopted by R
par()                                          #To list all parameters used in par() function
par("mar")                                     #To display the default setting of margin size of the graph
grep("blue", colors(), value=T)                #Find those values contain "blue" string and return these strings. It returns number only when set value=F
x <- rnorm(10, mean=0, sd=0.5)                 #Randomly pick 10 points from Normal Distribution N(0, 0.5)
y <- rep(1,10)                                 #pch (point character) cex (point size)

dev.new(width=5, height=4, unit="in")          #Will open a new window (inches) to show R plot()
plot(x=1:20, y=-1:-20, pch=19, cex=3, col=rgb(0.75,0.10,0.80))
dev.off()                                      #Will close this new window

dev.new(width = 550, height = 330, unit = "px")#Will open a new window (pixels) to show R plot()
plot(x=1:15, y=5:20, pch=19, cex=3, col=rgb(1.0,0.2,0.3))
graphics.off()                                 #Will close all opened windows even those inside plots window of RStudio

plot(x=1:10, y=rep(1,10), pch=19, cex=2, col="dark red")                             #The following codes show you different ways to use rgb function
points(x=1:10, y=rep(2,10), pch=19, cex=2, col="557799")                             #Color in HEX code
points(x=1:10, y=rep(3,10), pch=19, cex=2, col=rgb(.25, .5, .3))                     #RGB(0~1)
points(x=1:10, y=rep(4,10), pch=19, cex=2, col=rgb(10, 100, 100, maxColorValue=255)) #RGB(0~255)
plot(x=1:5, y=rep(5,5), pch=19, cex=8, col=rgb(.25, .5, .3, alpha=.5), xlim=c(0,6))  #alpha(0~1) means opacity/transparency
                                                                                     #xlim() means to set the limits of the X axis.
par(bg="gray40")                                                 #To set the background color using Graphical Parameters Function called "par()"
col.tr <- grDevices::adjustcolor("557799", alpha=0.7)            #To call a function "adjustcolor" from a package named "grDevices"
plot(x=1:5, y=rep(5,5), pch=19, cex=10, col=col.tr, xlim=c(0,6)) #Plot the color changed dynamically with "col.tr"
palette1 <- heat.colors(5, alpha=1)                              # 5 colors from the heat palette, opaque
palette2 <- rainbow(5, alpha=.5)                                 # 5 colors from the heat palette, transparent
plot(x=1:10, y=1:10, pch=19, cex=5, col=palette1)

PaletteOurOwnFunction <- colorRampPalette(c("gray80", "dark red"))      #To call a function named "colorRampPalette" to create our own color gradient (returns a function)
plot(x=10:1, y=1:10, pch=19, cex=5, col=PaletteOurOwnFunction(5))       #The color gradient recycle for every  5 circles 
plot(x=10:1, y=1:10, pch=19, cex=5, col=PaletteOurOwnFunction(10))      #The color gradient recycle for every 10 circles
palf <- colorRampPalette(c(rgb(1,1,1, .2),rgb(.8,0,0, .7)), alpha=TRUE) 
plot(x=10:1, y=1:10, pch=19, cex=5, col=palf(10))

plot(x=c(1:7), y=c( gsub(" ", "", gsub("[^0-9]*:", "", summary(C1)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C2)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C3)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C4)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C5)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C6)[4,2]), fixed=TRUE),
                    gsub(" ", "", gsub("[^0-9]*:", "", summary(C7)[4,2]), fixed=TRUE) ), 
     type="l", lwd=2, cex=3.0, col="blue", cex.main=3.5, cex.lab=1.5, cex.axis=1.5,
     xlab = "cluster ID", ylab = "CRC_Survival Average years", pch = 19, frame = TRUE)

pdf("Stdistribution.pdf",paper="usr", width = 0, height = 0) #3x3 (Each page contains 3 rows, each row contains 3 plots)
par(mfrow=c(3,3), mar=c(2,2,2,1))
sapply(colnames(TotalTumorCells),function(x) hist(TotalTumorCells[,x],main=x,breaks=10000, cex.main=1.0))
replicate(ceiling(ncol(TotalTumorCells)/3)*3-ncol(TotalTumorCells),plot.new())
# sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"),function(x) hist(apply(TTumorCells5000[,grep(x,colnames(TTumorCells5000))],1,median),main=x,breaks=100, cex.main=0.9))
dev.off()

pdf("BoxPlot167.pdf", paper="usr", width=0, height=0, useDingbats=FALSE)
boxplot(TotalTumorCells, cex=0.4, axes=FALSE)
axis(1, cex.axis=0.6)
axis(2, cex.axis=0.6)
mtext('167 columns (features)', side=1, line=2)
mtext('Quantity', side=2, line=2)
mtext('Data Distribution of 167 Columns', side=3, line=0)
dev.off()

png(filetag, height=3000, width=3000)
corrplot(TotalWeights.cor, tl.pos = NULL,
         upper.col=col4, lower.col=col4, tl.cex=4.0, cl.cex=4.0, number.cex=2.0, tl.col="red")
dev.off()

png(filename, height=800, width=800)
plot(x, y, main = maintitle, col="dark red",
     xlab = "TumorCell.X.Position", ylab = "TumorCell.Y.Position",
     pch = 19, frame = FALSE)
dev.off()

#barplot(TID_Dim_List$Dimension, main="Tumor size distribution (# cells)", xlab="932 Tumor Specimens")
hist(TID_Dim_List$Dimension, main="Tumor size distribution (# cells)", breaks=100, cex.main=1.0, xlab="932 Tumor Specimens", ylab="Frequency")
#The above hist() divided the 932 points into 100 groups (9.32/group, 100 groups so totally 932 points)

pheatmap(x, scale="row", cluster_rows=TRUE, clustering_distance_rows="euclidean", cluster_cols=TRUE, 
         clustering_distance_cols = "correlation", main="T cell nucleus size clustering with clinincal features", 
         color=colorRampPalette(c("blue", "red"))(81920), fontsize=9, fontsize_row=6)


#==================================================================================================================================================================================#
# Response variable is dependent variable (denote as Y: Outcome)
# Explanatory variable is independent variable (denote as X: Cause), also known as "predictors, covariates, features, columns …"
#==================================================================================================================================================================================#
### Using lm() to perform "multiple linear regression"
#multi.fit = lm(mpg~., data=mtcars)     # "." dot notation means all columns except "mpg" column here
#summary(multi.fit)                     # Y (dependent variable DV: to be predicted) is "mpg" ; X (independent variable IV: the predictor) is "cyl", "disp" ....."carb"
#> multi.fit = lm(mpg~., data=mtcars)   # To test if those columns affected the mpg performance
#> summary(multi.fit)                   # DV = f(IV+CV), CV:covariates
#
#Call:
#  lm(formula = mpg ~ ., data = mtcars)
#
#Residuals:                                 #The difference between the real value and the predicted value and it follows the next 4 rules
#  Min       1Q   Median      3Q     Max    #The mean of the errors is zero (and the sum of the errors is zero)
#-3.4506  -1.6044  -0.1196  1.2193  4.6271  #The distribution of the errors are normal
##All of the errors are independent
##Variance of errors is constant (Homoscedastic)
##Pr(>|t|) means statistical significance
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 12.30337   18.71788   0.657   0.5181  
#cyl         -0.11144    1.04502  -0.107   0.9161  
#disp         0.01334    0.01786   0.747   0.4635  
#hp          -0.02148    0.02177  -0.987   0.3350  
#drat         0.78711    1.63537   0.481   0.6353  
#wt          -3.71530    1.89441  -1.961   0.0633 . #the only p-value which is close to 0.05 and it means the mpg performance decrease -3.715 unit when wt increase one unit 
#qsec         0.82104    0.73084   1.123   0.2739  
#vs           0.31776    2.10451   0.151   0.8814  
#am           2.52023    2.05665   1.225   0.2340  
#gear         0.65541    1.49326   0.439   0.6652  
#carb        -0.19942    0.82875  -0.241   0.8122  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.65 on 21 degrees of freedom
#Multiple R-squared:  0.869,	Adjusted R-squared:  0.8066 
#F-statistic: 13.93 on 10 and 21 DF,  p-value: 3.793e-07
#
#Performance Measures: Three sets of measurements are provided.
#Residual Standard Error: This is the standard deviation of the residuals.  Smaller is better.
#Multiple / Adjusted R-Square: For one variable, the distinction doesn’t really matter.  
#
#R-squared shows the amount of variance explained by the model
#
#Adjusted R-Square takes into account the number of variables and is most useful for multiple-regression.
#F-Statistic: The F-test checks if at least one variable’s weight is significantly different than zero.  
#This is a global test to help asses a model.  If the p-value is not significant 
#(e.g. greater than 0.05) than your model is essentially not doing anything.
#==================================================================================================================================================================================#
#==================================================================================================================================================================================#

#Cox Proportional-Hazards Model	
#==================================================================================================================================================================================#
#==================================================================================================================================================================================#
#=======================================================================================
#   Description
# Survival in patients with advanced lung cancer from the North Central Cancer Treatment 
# Group. Performance scores rate how well the patient can perform usual daily activities.
#
#
# Usage
# lung
# cancer
# Format
# inst:	Institution code
# time:	Survival time in days
# status:	censoring status 1=censored, 2=dead
# age:	Age in years
# sex:	Male=1 Female=2
# ph.ecog:	ECOG performance score as rated by the physician. 0=asymptomatic, 
# 1= symptomatic but completely ambulatory, 2= in bed <50% of the day, 
# 3= in bed > 50% of the day but not bedbound, 4 = bedbound
# ph.karno:	Karnofsky performance score (bad=0-good=100) rated by physician
# pat.karno:	Karnofsky performance score as rated by patient
# meal.cal:	Calories consumed at meals
# wt.loss:	Weight loss in last six months
#======================================================================================
  
library("survival")
library("survminer")
library(lubridate)
library(MASS)
#status: censoring status 1=censored, 2=dead
#sex:	Male=1 Female=2
#data("lung")
#res.cox <- coxph(Surv(time, status) ~ sex, data = lung)

#> summary(IF_TMA_932_CellProportionWithBigNucleusSize)
#summary(IF_TMA_932_CellProportionWithBigNucleusSize$BigTProportion)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0001097 0.0010839 0.0021438 0.0050848 0.0052636 0.0540540 
#> d <- Sys.Date()
#
#x <- IF_TMA_932_CellProportionWithBigNucleusSize # 118 x 8 

#x$BigTProportion[which(x$BigTProportion >= 0.0052636)] <- 3
#x$BigTProportion[which((x$BigTProportion >= 0.0010839) & (x$BigTProportion < 3.0))] <- 2
#x$BigTProportion[which(x$BigTProportion < 2.0)] <- 1

#The status column has to be neumeric format
#x$BigTProportion[which(x$BigTProportion == 3.0)] <- "High"
#x$BigTProportion[which(x$BigTProportion == 2.0)] <- "Moderate"
#x$BigTProportion[which(x$BigTProportion == 1.0)] <- "Low"

#x$Gender[which(x$Gender == 1)] <- 2 #Female
#x$Gender[which(x$Gender == 0)] <- 1

#x$CRCensored[which(x$CRCensored == 1)] <- 2 #dead
#x$CRCensored[which(x$CRCensored == 0)] <- 1 #alive(censored)

#res.cox <- coxph(Surv(CRCSurv, CRCensored) ~ BigTProportion + Age + Gender, data =  x)
#summary(res.cox)
#> summary(res.cox)
#Call:
#  coxph(formula = Surv(CRCSurv, CRCensored) ~ BigTProportion + 
#          Age + Gender, data = x)

#n= 932, number of events= 284 

#coef exp(coef)  se(coef)      z   Pr(>|z|)
#BigTProportion -0.082725  0.920604  0.073742 -1.122    0.262  #beta=-0.08 means protection, HR = e^(-0.08)=0.9 means BigerT nucleus resuced 92% risk, but p=0.262 so the conclusion is not significant
#Age             0.011127  1.011189  0.007153  1.556    0.120  #HR=1 no affaect, HR <1 means protection, HR > 1 means hazardous
#Gender         -0.059827  0.941927  0.123335 -0.485    0.628  #beta = slope, when slope is negative means decreasing so it means protection because of decreasing risk

#exp(coef) exp(-coef) lower .95 upper .95
#BigTProportion    0.9206     1.0862    0.7967     1.064
#Age               1.0112     0.9889    0.9971     1.025
#Gender            0.9419     1.0617    0.7397     1.199

#Concordance= 0.53  (se = 0.018 )
#Likelihood ratio test= 3.9  on 3 df,   p=0.3
#Wald test            = 3.9  on 3 df,   p=0.3
#Score (logrank) test = 3.9  on 3 df,   p=0.3
#======================================================================
#  Call:
#  coxph(formula = Surv(CRCSurv, CRCensored) ~ BigTProportion + 
#          Age + Gender, data = x)

#n= 118, number of events= 36 

#coef exp(coef) se(coef)      z Pr(>|z|)   
#BigTProportion  0.82414   2.27993  0.25658  3.212  0.00132 **
#  Age            -0.01870   0.98148  0.02018 -0.926  0.35426   
#Gender          0.21163   1.23569  0.35366  0.598  0.54957   
#---
#  Signif. codes:  0 �***� 0.001 �**� 0.01 �*� 0.05 �.� 0.1 � � 1
#
#exp(coef) exp(-coef) lower .95 upper .95
#BigTProportion    2.2799     0.4386    1.3789     3.770
#Age               0.9815     1.0189    0.9434     1.021
#Gender            1.2357     0.8093    0.6178     2.471

#Concordance= 0.664  (se = 0.046 )
#Likelihood ratio test= 11.98  on 3 df,   p=0.007
#Wald test            = 11.46  on 3 df,   p=0.009
#Score (logrank) test = 12.02  on 3 df,   p=0.007
#=======================================================================
  # broom::tidy(
  #  coxph(Surv(CRCSurv, CRCensored) ~ BigTProportion + Age + Gender, data = x), 
  #  exp = TRUE
  #  ) %>% 
  #  kable()
  
#  coxph(Surv(CRCSurv, CRCensored) ~ BigTProportion + Age + Gender, data = x) %>% 
#  gtsummary::tbl_regression(exp = TRUE) 


# Create the new data  
#BiggerTcell <- with(x, data.frame(BigTProportion = c(1, 2, 3), Age = rep(mean(Age, na.rm = TRUE), 3), Gender = c(2, 2, 2))) #Female
#BiggerTcell <- with(x, data.frame(BigTProportion = c(1, 2, 3), Age = rep(mean(Age, na.rm = TRUE), 3), Gender = c(1, 1, 1))) #Male

#fit <- survfit(res.cox, newdata = BiggerTcell)
#ggsurvplot(fit, data = BiggerTcell, conf.int = TRUE, 
#           legend.labs=c("BigT=1", "BigT=2", "BigT=3"),
#           ggtheme = theme_minimal())


#ggsurvplot(survfit(res.cox), data = x, palette = "#2E9FDF", 
#           ggtheme = theme_minimal(), legend = "top")


#plot(survfit(Surv(CRCSurv, CRCensored) ~ 1, data = x), 
#     xlab = "Years", 
#     ylab = "CRC survival probability")

#ggsurvplot(
#  fit = survfit(Surv(CRCSurv, CRCensored) ~ 1, data = x), 
#  xlab = "Years",
#  ylab = "CRC survival probability")


#Drawing heat map for clustering analyses
#library(pheatmap)
#x <- IF_TMA_932_CellProportionWithBigNucleusSize # 118 x 8 
#rownames(x) <- x[, 1]
#x <- x[, -1]
#x <- x[, -4] #removed column CRCensored

#x <- data.matrix(x)
#pheatmap(x, scale="row", cluster_rows=TRUE, clustering_distance_rows="euclidean", cluster_cols=TRUE, 
#         clustering_distance_cols = "correlation", main="T cell nucleus size clustering with clinincal features", 
#         color=colorRampPalette(c("blue", "red"))(81920), fontsize=9, fontsize_row=6)

#==================================================================================================================================================================================#
#==================================================================================================================================================================================#




