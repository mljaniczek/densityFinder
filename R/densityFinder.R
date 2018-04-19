#############################
# densityFindR
# Author: Margaret Hannum
# 12/1/17
# To find position before which contains highest density of nonsense mutations
# Then test if that density is significantly different than the background mutation rate (calculated using synonymous mutations)
#
library(tidyverse)
setwd("/Users/TinyDragon/Desktop/Data")
load("NSM_syn_onco.RData") #Amr's subset of nonsense and synonymous mutations, excluding indels, and categorized via Chang et al paper
#Split into nonsense and synonymous
NSM_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Nonsense_Mutation")
syn_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Silent")

#We need to make sure amino acid position is numeric and complete
nsmProt <- subset(NSM_syn_onco, complete.cases(NSM_syn_onco$Amino_Acid_Position))
nsmProt$Amino_Acid_Position <- as.numeric(nsmProt$Amino_Acid_Position)

#Now make subset of columns so loop runs faster
features = c("Hugo_Symbol", "Amino_Acid_Length", "Amino_Acid_Position", "Variant_Classification")
nsmselect <- select(nsmProt, features)
#Subsetting just genes that have both nonsense and silent mutations
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% syn_onco$Hugo_Symbol)
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% NSM_onco$Hugo_Symbol)
#637232 observations

genename <-c(unique(nsmselect$Hugo_Symbol)) #16411 common genes between Nonsense and silent
#Fisher exact test on nonsense vs silent mutations for each gene after selecting position using our densityFindR algorithm
start.time = Sys.time()
cutoff = c()
p = c()
or = c()
length = c()
hugo = c()
for (i in 100:200){   #test with a few
#for (i in 1:length(genename)){  #for full set : took 49 minutes
  tryCatch({
  gene = genename[i]
  sub <- subset(NSM_syn_onco, NSM_syn_onco$Hugo_Symbol == gene)
  subnon <- subset(sub, sub$Variant_Classification == "Nonsense_Mutation")
  if (length(subnon$Hugo_Symbol) >=15){
  ###Initial length and pos if using protein length and Amino acid position
  #Get initial lengths and positions - using PROTEIN instead of gene length
  protlength <- sub$Amino_Acid_Length[1]
  begin <- 1
  center <- round((protlength - begin)/2, digits = 0)
  up <- round((center - begin)/2, digits = 0)
  down <- round((center + (protlength - center)/2), digits = 0)
  while ((down - center)>15){
    #Calculate first densities
    udense <- sum(c(subnon$Amino_Acid_Position < up), na.rm = TRUE)/up
    cdense <- sum(c(subnon$Amino_Acid_Position < center), na.rm = TRUE)/center
    ddense <- sum(c(subnon$Amino_Acid_Position < down), na.rm = TRUE)/down
    #Choose bound with largest density upstream of point
    newc <- c(up, center, down)[which.max(c(udense, cdense, ddense))] #Which max provides index of which one is max
    #Need to make new upper and lower bounds
    newu <- round(newc - (down - center)/2, digits = 0)
    newd <- round(newc + (down - center)/2, digits = 0)
    #redefine center, up, down
    center = newc
    up = newu
    down = newd
  } #second while end curly
  sub$before.cutoff <- as.numeric(sub$Amino_Acid_Position <= center)
  genesubsetTest <- fisher.test(table(sub$Variant_Classification, sub$before.cutoff)[,c(2,1)])
  p <- c(p, genesubsetTest$p.value)
  or <- c(or, genesubsetTest$estimate)
  length <- c(length, protlength)
  cutoff <- c(cutoff, center)
  hugo <- c(hugo, gene)
  print(i)
  } #end first while curly
  else {
    print(i)}
  }, #Trycatch end curly
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
   } #end loop
end.time = Sys.time()
end.time - start.time #Took 49 minutes

#1278 genes had more than 15 nonsense mutations
densityResults <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
densityResults$gene <- hugo
densityResults$pval <- p
densityResults$odds <- or
densityResults$Amino_Acid_Pos <- cutoff 
densityResults$Amino_Acid_Length <- length

save(densityResults, file = "densityfinder15mutAll.RData")
#Now make subset of p-val under FDR cutoff of 0.1: 64 observations in our sample
densityResultsFDR <- subset(densityResults,densityResults$p.adjust.p..method....BH.. <= 0.1)

#Can now compare against COSMIC census to see cancer genes
cosmicGenes <- read.csv("Census_allMon Nov 20 16_36_34 2017.csv")
cosmicSelect <- subset(cosmicGenes, cosmicGenes$Gene.Symbol %in% densityResultsFDR$gene)





###################
# OLD CODE BELOW  #
###################

###################
# Try Z-score test using density of non/syn #
#############################################
### Normal approx to binom
start.time = Sys.time()
results <- c()
pval = c()
AA_Length  = c()
Cutoff = c()
RegionLength = c()
Syn_ctRegion = c()
BG_Density = c()
sigma = c()
pvalue = c()
bgVar = c()
#or (i in 10:20){
for (i in 1:length(densityResults$gene)){
  tryCatch({
  sub <- subset(NSM_syn_onco, NSM_syn_onco$Hugo_Symbol == densityResults$gene[i])
  subnon <- subset(sub, sub$Variant_Classification == "Nonsense_Mutation")
  subsyn <- subset(sub, sub$Variant_Classification == "Silent")
  n = sub$Amino_Acid_Length[1] - densityResults$cut[i]
  mu = length(subset(subsyn, subsyn$Amino_Acid_Position < densityResults$cut[i])[,1])
  p = mu/n
  q = 1-p
  sigma = sqrt(n*p*q)
  #check if approx is valid?
  #Do normal approximation p-val
  pval = 1 - pnorm(length(subnon$Hugo_Symbol), mean = mu, sd = sigma)
  gene = c(gene, densityResults$gene[i])
  AA_Length  = c(AA_Length, sub$Amino_Acid_Length[1])
  Cutoff = c(Cutoff, densityResults$cut[i])
  RegionLength = c(RegionLength, n)
  Syn_ctRegion = c(Syn_ctRegion, mu)
  BG_Density = c(BG_Density, p)
  bgVar = c(bgVar, sigma)
  pvalue = c(pvalue, pval)
print(i) 
 }, #Trycatch end curly
error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
end.time = Sys.time()
end.time - start.time #Took 49 minutes
colnames(results) <- c("Gene", "Amino_Acid_Length", "Cutoff", "RegionLength", "Syn_ctRegion", "BG_Density", "sigma", "pval")

result.df$padj <- c(p.adjust(result.df$pval, method = "BH"))

results_fdr1 <- subset(results, results$padj < 0.1)


results <- as.data.frame(results)
rownames(results) <- results[,1]
result.df <- data.frame(results[,-1])
i=6
sub <- subset(NSM_syn_onco, NSM_syn_onco$Hugo_Symbol == densityResults$gene[i])
subnon <- subset(sub, sub$Variant_Classification == "Nonsense_Mutation")
subsyn <- subset(sub, sub$Variant_Classification == "Silent")
mu = length(subset(subsyn, subsyn$Amino_Acid_Position < densityResults$cut[i])[,1])

##################
#Individual gene #
#Used below prior to making loop above
#
#Subset for gene
sub <- subset(NSM_syn_onco, NSM_syn_onco$Hugo_Symbol == "DST")
sub$Amino_Acid_Position <- as.numeric(sub$Amino_Acid_Position)
subnon <- subset(sub, sub$Variant_Classification == "Nonsense_Mutation")
###Initial length and pos if using protein length and Amino acid position
#Get initial lengths and positions - using PROTEIN instead of gene length
protlength <- sub$Amino_Acid_Length[1]
begin <- 1
center <- round((protlength - begin)/2, digits = 0)
up <- round((center - begin)/2, digits = 0)
down <- round((center + (protlength - center)/2), digits = 0)

#Calculate first densities
udense <- sum(c(subnon$Amino_Acid_Position < up), na.rm = TRUE)/up
cdense <- sum(c(subnon$Amino_Acid_Position < center), na.rm = TRUE)/center
ddense <- sum(c(subnon$Amino_Acid_Position < down), na.rm = TRUE)/down

#Choose bound with largest density upstream of point
newc <- c(up, center, down)[which.max(c(udense, cdense, ddense))] #Which max provides index of which one is max
#Need to make new upper and lower bounds
newu <- round(newc - (down - center)/2, digits = 0)
newd <- round(newc + (down - center)/2, digits = 0)

#redefine center, up, down
center = newc
up = newu
down = newd

#Maybe stop loop here if center, up, and down are within 10 bp of each other? and choose center?
#Otherwise go back and calculate densities, find new center, repeat. 

