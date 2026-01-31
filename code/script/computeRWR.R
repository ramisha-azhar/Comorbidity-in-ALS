# https://medium.com/@chaitanya_bhatia/random-walk-with-restart-and-its-applications-f53d7c98cb9

library(RANKS)

options(stringsAsFactors = F)
##############################
setwd("C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/code")

source("script/createAdjMatrix.R")
##############################
disease_sel <- "Amyotrophic Lateral Sclerosis" #the selected disease
##############################
path <- "C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/code/files/"
disease_gene <- read.table(paste0(path,"Phenopedia.txt"), header = T, sep = '\t', check.names = F, quote = "")
disease_gene$disease <- gsub(pattern = "\"", x = disease_gene$disease, replacement = "")

nodes_sel <- disease_gene[which(disease_gene$disease %in% disease_sel),"GeneID"] #initial set of the starting nodes for the random walk
##############################
#path <- "C:/Users/ramis/OneDrive/Desktop/BIOINFORMETICS 2/Comorbidity/code/files/"
#interactome <- read.table(paste0(path,"Supplementary_data1.txt"),header = T, sep = '\t', check.names = F, quote = "")
#interactome <- interactome[,-3]
#colnames(interactome) <- c("source","target")

#W <- createAdjMatrix(interactome)
#save(W, file = "Interactome_adj.RData")
 
#rm(interactome)
##############################
load("Interactome_adj.RData") # we can run the above script to create the matrix 
#once the matrix is created we can comment the above and just load the matrix 
#there is no need to create the matrix again and again

ind.positives <- which(rownames(W) %in% nodes_sel) #initail vector of the disease genes
names(ind.positives) <- rownames(W)[ind.positives]

res <- RWR(W, ind.positives, gamma = 0.8, tmax = 1000, eps = 1e-10, norm = TRUE) #random walk is computes

p <- res$p #probabilities for each node in the network
node_prob <- data.frame(GeneID = names(p), probability = p)

df <- merge(node_prob, disease_gene, by = "GeneID", all.x = T) 
#build a data frame with the disease associated to each node and compute the number of reached genes for each disease

summary <- data.frame(table(df$disease)) #the number of times the related disease genes is reached
View(summary)
##############################
list <- split(df$probability,df$disease)
mean <- data.frame(mean_prob=sapply(list,mean))

results <- merge(summary,mean,by.x="Var1",by.y = 0)
colnames(results) <- c("disease","node_num","prob_mean")

ind <- which(results$node_num > 30) #threshold for the number of reached disease genes
results <- results[ind,]

results$zscore <- (results$prob_mean - mean(results$prob_mean)) / sd(results$prob_mean)# it is to normalize the probability mean
results$zscore_mod <- (results$prob_mean - median(results$prob_mean)) / mad(results$prob_mean)
results <- results[order(results$zscore, decreasing = T),]
##############################
observation <- results[results$disease == "Frontotemporal Dementia","zscore_mod"]
observation1 <- results[results$disease == "Polyneuropathies","zscore_mod"]
observation2 <- results[results$disease == "Dementia, Vascular","zscore_mod"]
observation3 <- results[results$disease == "Muscle Weakness","zscore_mod"]
observation4 <- results[results$disease == "Neuromuscular Diseases","zscore_mod"]
observation5 <- results[results$disease == "Diabetic Neuropathies","zscore_mod"]

xmin <- min(-3,observation)
xmax = max(3,observation)
fun <- function(x) dnorm(x,m=0,sd=1)
plot <- curve(fun, xlim = c(xmin,xmax), xlab = "ALS closeness", ylab = "Probability density", n = 1000)
polygon(plot, col = "grey")

abline(v = observation, lty = 2, lwd = 2, col = "red")
abline(v = observation1, lty = 2, lwd = 2, col = "green")
abline(v = observation2, lty = 2, lwd = 2, col = "blue")
abline(v = observation3, lty = 2, lwd = 2, col = "violet")
abline(v = observation4, lty = 2, lwd = 2, col = "orange")
abline(v = observation5, lty = 2, lwd = 2, col = "turquoise")

write.table(results, "results_RWR.txt", sep = "\t", row.names = F, col.names = T, quote = F )

