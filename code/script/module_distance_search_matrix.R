library(openxlsx)
library(igraph)

options(stringsAsFactors = F)

setwd("C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/SWMmer vs NM/code")
######################################
source("script/computeDegreeDistribution.R")
source("script/computeNetworkDistance.R")
source("script/computeNetworkSeparation.R")
source("script/computeMinimum.R")
source("script/selectRandomNodes.R")
source("script/performIterationSeparation.R")
source("script/computeStatistics.R")
######################################
perc_thr <- 5
iter <- 100
######################################
path <- "C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/code/files/"

interactome <- read.table(paste0(path,"Supplementary_data1.txt"),header = T, sep = '\t', check.names = F, quote = "")

filename1 <- "C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/code/files/disease_gene.xlsx"

filename2 <- filename1 #both the files are the same and we take xlxs file which has 7 sheets one for each disease to measure the separation measure

check <- identical(filename1,filename2) #TRUE 

# esempio: 
# filename1 - disease genes per 6 malattie tra cui breast e lung cancer
# filename2 - sono gli switch genes per le stesse malattie
# voglio calcolare se switch e disease genes si overlappano
# ovvero la matrice delle separation tra switch genes e disease genes
# sulle righe avremo gli switch genes e sulle colonne i disease genes
# la matrice non e' simmetrica infatti
# la separation tra gli switch genes del breast cancer e i disease genes del lung cacer e' diversa
# dalla separation dei disease genes del breast cancer e gli switch genes del lung cancer
######################################
disease1 <- getSheetNames(filename1) #the sheet nemaes are crossponds to the diseases to be tested
disease2 <- getSheetNames(filename2)

n1 <- length(disease1) #it is 7 total number of sheets or diseases
n2 <- length(disease2) #it is 7 
######################################
graph <- graph_from_data_frame(interactome, directed = F) #create graph from interactome network

graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                  edge.attr.comb = igraph_opt("edge.attr.comb"))

node <- names(V(graph))

degree <- degree(graph, v = V(graph))

rm(interactome)
######################################

separation <- matrix(NA, nrow = n1, ncol = n2) #separation value matrix  , if the value is negative the diseases overlapp
pval <- matrix(NA, nrow = n1, ncol = n2) # p value matrix , statistically significant interaction if less than 0.05

rownames(separation) <- disease1
colnames(separation) <- disease2

rownames(pval) <- disease1
colnames(pval) <- disease2

# compute separation values and p values for each pair of diseases and build the matrices
for(i in 1:n1){
  
  list1 <- read.xlsx(filename1, sheet = i)
  
  start <- ifelse(check,i,1)
  for(j in start:n2){
    
    list2 <- read.xlsx(filename2, sheet = j)
    
    separation[i,j] <- computeNetworkSeparation(list1$ID,list2$ID)
  
    degree_distribution_1 <- computeDegreeDistribution(list1$ID) 
    degree_distribution_2 <- computeDegreeDistribution(list2$ID)
    random_distribution <- performIterationSeparation(degree_distribution_1,degree_distribution_2,iter)
    par <- computeStatistics(random_distribution,separation[i,j],"Network separation",plot=F)
    
    pval[i,j] <- par$pval
    
    if(check){
      separation[j,i] <- separation[i,j] 
      pval[j,i] <- pval[i,j]
    }
    
  }}

write.table(separation,"matrix_separation.txt",row.names = T, quote=F, col.names = NA, sep ="\t")
write.table(pval,"matrix_separation_pval.txt",row.names = T, quote=F, col.names = NA, sep ="\t")
