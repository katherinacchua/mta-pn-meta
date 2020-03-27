### This script is to pull post-imputation data in 100,000 snp chunks and run through 
### the CoxPH analysis. 

### Inputs for this pipeline are:
### PARAM1 = this script and PARAM2 = vcf.gz 
### .info file
### C40502phen2data-imp-adjdose.RData

### Outputs from this pipeline are:
### "NonSNPsites_chrX_C40502.txt" = writes out snps that were removed from post-imputation data and the reason for so
### "C40502_chrX_infofilt.txt" = writes out filter info file, excluding REF/ALT != ATCG, >1; REF==ALT, Rsq > 0.3
### "chrX_dsdat[0-9].txt" = writes out the dosage data that is pulled from vcf.gz files
### "chrX_dsdat.filt[0-9].txt" = writes out the dosage data that is filtered by infofilt file
### "CALGB40502_imp.dsdatX-PN2-g_euro-adjdose.RData" = saved phenotype and dosage data merged
### "CALGB40502_Cox_PN2_results_adjdose_chrX_imp.dsdatX.txt" = saved results from each chunk at chrX

basedir <- "~/katherinachua/C40502/"
libdir <- "~/R/x86_64-redhat-linux-gnu-library/3.3/"
setwd(basedir)

args <- commandArgs(T)

## Files to read in 
info_Name <- args[1]


library(VariantAnnotation, lib.loc = libdir)
library(survival)
library(gdata)
library(data.table)
library(parallel)

sessionInfo()

## Load params for dsdata
myfile <- info_Name
myvcf <- paste0(basedir, myfile)
mygenome <- "hg19"
tab <- TabixFile(myvcf, yieldSize = 100000)
mychr <- strsplit(myfile, ".", fixed = T)[[1]][1]
chr.no <- strsplit(mychr, "[^[:digit:]]")[[1]][4]
print(mychr)

## Load info data
imp.info <- read.table(paste0(basedir, "chr", chr.no, ".info", sep = ""), header = T, stringsAsFactors = F)
names(imp.info)[1:4] <- c("SNP", "REF", "ALT", "ALT_Frq")

#################################################################################################################
## 1. Filter 'Non-snps' (as checkVCF.py does) meaning remove those that have length of ref or alt allele > 1,
## 		not composed of ATCG, and ref == alt. Also filter by Rsq > 0.3.
#################################################################################################################

ATCG = list('A', 'T', 'C', 'G')

remove.indexAll <- remove.index <- remove.index2 <- NULL ##Filter SNPs that are not in ATCG format

for (i in 1:length(imp.info$SNP)){
    if ( nchar(imp.info$REF[i]) != 1 | nchar(imp.info$ALT[i]) != 1 | !(imp.info$REF[i] %in% ATCG) | !(imp.info$ALT[i] %in% ATCG) ){
        remove.index <- c(i, imp.info$SNP[i], "REF/ALT > 1 or REF/ALT != ATCG")
        remove.indexAll <- rbind(remove.indexAll, remove.index)
    }
    if ( imp.info$REF[i] == imp.info$ALT[i] ){
        remove.index2 <- c(i, imp.info$SNP[i], "REF == ALT")
        remove.indexAll <- rbind(remove.indexAll, remove.index2)
    }else{
        next
    }
}

remove.indexAll <- data.frame(remove.indexAll, stringsAsFactors = F, row.names = NULL)
names(remove.indexAll) <- c("index", "SNP", "removereason")
remove.indexAll$filterindex <- paste0(remove.indexAll$index, "_", remove.indexAll$SNP, sep = "")
write.table(remove.indexAll, file = paste0("NonSNPsites_chr", chr.no, "_C40502.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

remove.indexAll.nodups <- remove.indexAll[!duplicated(remove.indexAll$filterindex),]
duplicates <- remove.indexAll[duplicated(remove.indexAll$filterindex),] ##Are there duplicates? There shouldn't be

if (length(duplicates$index) > 0) { ##If there are duplicates, filter by remove.indexAll.nodups but if not filter by remove.indexAll
    imp.info.filt <- subset(imp.info, !(row.names(imp.info) %in% remove.indexAll.nodups$index))
} else{
    imp.info.filt <- subset(imp.info, !(row.names(imp.info) %in% remove.indexAll$index))
}

imp.info.filt <- subset(imp.info.filt, Rsq > 0.3 | Rsq == "-")
imp.info.filt$snp.index <- paste0(imp.info.filt$SNP, "_", imp.info.filt$REF, "_", imp.info.filt$ALT)
write.table(imp.info.filt, file = paste0("C40502_chr", chr.no, "_infofilt.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names =T)

#################################################################################################################
## 2.  Use this imputed info file to filter out the dosage information. Since analysis is only on dosage data.
#################################################################################################################

open(tab)
i <- 0
param <- ScanVcfParam(fixed = "ALT", geno=c("GT", "DS"), info = c("AF", "MAF", "R2"))
while(nrow(vcf_test <- readVcf(tab, mygenome, param = param))){
    i <- i + 1
    #    gdat <- paste0(mychr, "_", "gdat", i, ".txt", sep = "")
    #    write.table(geno(vcf_test)$GT, file = gdat, col.names = T, quote = F, sep = "\t", row.names = T)
    
    snpids<-vcf_test@rowRanges@ranges@NAMES
    index <- c(1:vcf_test@rowRanges@seqnames@lengths)
    refalt <- data.frame(fixed(vcf_test), snpids, index, stringsAsFactors = F)
    refalt$snp.index <- paste0(snpids, "_", refalt$REF, "_", refalt$ALT)
    maf <- data.frame(info(vcf_test), snpids, index, stringsAsFactors = F)
    snpinfo <- merge(refalt, maf, by = c("index", "snpids"))
    #    write.table(snpinfo, file = snpdat, col.names = T, quote = F, sep = "\t", row.names = F)
    
    dsdat <- paste0(mychr, "_", "dsdat", i, ".txt", sep = "")
    dinfo <- data.frame(index, row.names(geno(vcf_test)$DS), geno(vcf_test)$DS, stringsAsFactors = F)
    names(dinfo)[1:2] <- c("index", "snpids")
    final.dinfo <- merge(snpinfo, dinfo, by = c("index", "snpids")) 
    write.table(final.dinfo, file = dsdat, col.names = T, quote = F, sep = "\t", row.names = T)
    
    cat("vcf gen dim:", dim(vcf_test), "\n")
}
close(tab)
warnings()

rm(dinfo, final.dinfo, imp.info, maf, refalt, snpinfo, remove.indexAll, remove.indexAll.nodups, ATCG, dsdat, i, index,
   param, remove.index, remove.index2, tab, snpids, vcf_test)


files <- list.files(pattern = paste0(mychr, "_dsdat[0-9]"))
for (i in 1:length(files)){
    name <- paste(mychr, "_dsdat.filt", i, sep = "")
    data <- read.table(files[i], stringsAsFactors = F, header = T, sep = "\t")
    data <- subset(data, data$snp.index %in% imp.info.filt$snp.index)
    #    data <- conv.gt(data) #for gtdata only; do not need this line for dsdata
    
    write.table(data, file = paste0(name, ".txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)
}

rm(data, files, name)

#################################################################################################################
## 3.  Load the phenotype data into R. Read in all the filtered dosage data. Merge phenotype and dosage data.
##			Run the Cox model and write out the results.		
#################################################################################################################
load(paste0(basedir, "C40502phen2data-imp-adjdose.RData")) ##variable = PN2data_imp

files <- list.files(pattern = paste0(mychr, "_dsdat.filt[0-9]"))

for (i in 1:length(files)){
    imp.dsdat <- read.table(files[i], stringsAsFactors = F, header = T, sep = "\t", row.names = NULL)
    imp.dsdat$snp.index <- gsub(paste0(chr.no,":"), paste0(mychr,"."), imp.dsdat$snp.index)
    imp.dsdat <- subset(imp.dsdat, select = -(c(index, snpids, REF, ALT, AF, MAF, R2)))
    imp.dsdat <- t(data.frame(imp.dsdat, stringsAsFactors = F, row.names = NULL))
    colnames(imp.dsdat) <- imp.dsdat["snp.index",]
    imp.dsdat <- imp.dsdat[c(2:length(row.names(imp.dsdat))),]
    all.imp.dsdat <- merge(PN2data_imp, imp.dsdat, by.x = "SAMPID", by.y = "row.names")
    save(all.imp.dsdat, file = paste0("CALGB40502_", mychr, "_imp.dsdat", i, "-PN2-g_euro-adjdose.RData"))
}
rm(files, imp.dsdat, all.imp.dsdat)


## Baseline formula for cox model
pnfirstgr2.mod0 <- as.formula('Surv(all.imp.dsdat$adj.PN2cum, all.imp.dsdat$pn2plus_event) ~ 
                            strata(all.imp.dsdat$Arm) + log10(all.imp.dsdat$Age)')
## Baseline formula for cox model without covariates
pnfirstgr2.mod <- as.formula('Surv(all.imp.dsdat$adj.PN2cum, all.imp.dsdat$pn2plus_event) ~ 
                            strata(all.imp.dsdat$Arm)')

## Function for fitting Cox model with snp
## Function for fitting Cox model with snp
func.pnfirstgr2 <- function(x){
    rsid <- names(pndata)[x]
    input <- paste0("as.numeric(as.character(pndata$",rsid,"))")
    pnfirstgr2.mod1 <- update(pnfirstgr2.mod0, paste0('~.+',input))
    pnfirstgr2.mod2 <- update(pnfirstgr2.mod, paste0('~.+',input))
    if(length(unique(pndata[[rsid]])) > 1) { # check if rsid has more than one value, otherwise, skip to else and set pval = NA
        mod <- try(summary(coxph(pnfirstgr2.mod1)))
        if (class(mod) != "try-error"){ # model will return try-error if there is not big enough difference in fit
            n <- mod$n
            beta <- coef(mod)[input, "coef"]
            wpval <- coef(mod)[input,"Pr(>|z|)"] # wpval: p-value
            hr <- coef(mod)[input,"exp(coef)"] # hr: hazards ratio
            se <- coef(mod)[input, "se(coef)"] # se: standard error
            sc_pval <- mod$sc[["pvalue"]] # score test : p-value
            sc_df <- mod$sctest["df"] 
        } else{
            n <- beta <- wpval <- hr <- se <- sc_pval <- sc_df <- NA 
        } # set NA to those that model will not fit appropriately (those with only one value or those that have little difference in dosage)
        
        mod2 <- try(summary(coxph(pnfirstgr2.mod2)))
        if (class(mod2) != "try-error"){	
            n2 <- mod2$n
            beta2 <- coef(mod2)[input, "coef"]
            wpval2 <- coef(mod2)[input,"Pr(>|z|)"] # wpval: p-value
            hr2 <- coef(mod2)[input,"exp(coef)"] # hr: hazards ratio
            se2 <- coef(mod2)[input, "se(coef)"] # se: standard error
            sc_pval2 <- mod2$sctest[["pvalue"]] # score test: p-value
            sc_df2 <- mod2$sctest["df"] # score test: degrees of freedom
        } else{ n2 <- beta2 <- wpval2 <- hr2 <- se2 <- sc_pval2 <- sc_df2 <- NA }
    } else{  
        n <- beta <- wpval <- hr <- se <- sc_pval <- sc_df <- n2 <- beta2 <- wpval2 <- hr2 <- se2 <- sc_pval2 <- sc_df2 <- NA
    }
    
    coll <- data.frame(rsid, n, beta, hr, se, wpval, sc_pval, sc_df, n2, beta2, hr2, se2, wpval2, sc_pval2, sc_df2, stringsAsFactors = F)
    return(coll)
}

files <- list.files(pattern = paste0("CALGB40502_", mychr, "_imp.dsdat[0-9]"))

for (i in 1:length(files)){
    load(paste0(basedir, "CALGB40502_", mychr, "_imp.dsdat", i, "-PN2-g_euro-adjdose.RData"))
    print(paste("Running cox model on imputation data set", i, "out of", length(files), sep = " "))
    
    ## Separate snp data into multiple lists so that you can run in parallel
    pndata.snps <- all.imp.dsdat[, c((length(PN2data_imp)+1):length(all.imp.dsdat))]
    numsets <- floor(length(pndata.snps)/1000)
    
    for (j in 1:numsets){
        name <- paste("pndata.snps",j, sep="")
        assign(name, pndata.snps[,c(1:1000)])
        pndata.snps <- pndata.snps[,-c(1:1000)]
    }
    pndata.snps.list <- lapply(ls(pattern = "pndata.snps"), get)
    
    ## Use function to run gwas in for loop using snp list
    for (k in 1:length(pndata.snps.list)){
        pndata <- pndata.snps.list[[k]]
        out <- mclapply(1:length(names(pndata)), func.pnfirstgr2, mc.cores = 3)
        results <- paste("results.pndata", k, sep = "")
        print(results)
        assign(results, do.call(rbind, out))
    }
    
    ## Save results
    results.pnfirstgr2.list <- lapply(ls(pattern = "results.pndata[0-9]"), get)
    results.pnfirstgr2.all <- rbindlist(lapply(results.pnfirstgr2.list, as.data.frame), fill = T)
    write.table(results.pnfirstgr2.all, file = paste0("CALGB40502_Cox_PN2_results_adjdose_ageonly", mychr, "_imp.dsdat",i, ".txt"), sep ="\t", row.names = F, col.names = T, quote = F) 
    
    rm(list = ls(pattern = "pndata"))
    rm(all.imp.dsdat)
    rm(list = ls(pattern = "results.pndata[0-9]"))
    rm(numsets, name, out, results, results.pnfirstgr2.list, results.pnfirstgr2.all, i, j, k)
}

rm(files)

results.files <- list.files(pattern = paste0("CALGB40502_Cox_PN2_results_adjdose_ageonly", mychr, "_imp.dsdat[0-9]"))
for (i in 1:length(results.files)){
    name <- paste("GWAS.results", i, sep = "")
    data <- read.table(results.files[i], stringsAsFactors = F, header = T, sep = "\t", fill = T)
    assign(name, data)
}
results.list <- lapply(ls(pattern = "GWAS.results[0-9]"), get)
results.all <- rbindlist(lapply(results.list, as.data.frame), fill = T)
write.table(results.all, file = paste0("CALGB40502_Cox_PN2_results_adjdose_ageonly", mychr, "_imp.dsdat.ALL.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

rm(results.files, name, data, results.list, results.all)
rm(list = ls(pattern = "GWAS.results[0-9]"))
warnings()