if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#Other packages
install.packages("ape")
install.packages("phangorn")
install.packages("VennDiagram")

library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)

setwd("~/School_Work_PU/ANSC595/Project")

#OTU table (shared file)
OTU = read.table("stability.opti_mcc.shared", header=TRUE, sep="\t")

#Taxonomy of each OTU
tax = read.table("stability.taxonomy", header=TRUE, sep="\t")

#Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta = read.table("Metadata.txt", header=TRUE, row.names=1, sep="\t")


#OTU Table
row.names(OTU) = OTU$Group
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]

#Taxonomy Table
row.names(tax) = tax$OTU
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean = separate(tax, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "Strain", "OTU"))]

#Order the data
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
meta = meta[order(row.names(meta)),]


#Set seed
set.seed(8765)

##ALPHA DIVERSITY
#Create 2x2 plot environment so that we can see all 4 metrics at once. 
par(mfrow = c(2, 2))

#Then plot each metric.
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$simpson, main="Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)

#Change to inverse Simpson
#Create 2x2 plot environment 
par(mfrow = c(2, 2))

#Plots
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(1/meta$simpson, main="Inverse Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)

#Test for normalcy
shapiro.test(meta$shannon)
shapiro.test(1/meta$simpson)
shapiro.test(meta$chao)
shapiro.test(meta$ace)

#Run the ANOVA and save it as an object
aov.shannon.age = aov(shannon ~ Feed.id, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.age)
TukeyHSD(aov.shannon.age)

#Re-order the groups because the default is 1yr-2w-8w
meta$Feed.id = factor(meta$Feed.id, c("CTM","MR1","MR2"))
#Return the plot area to 1x1
par(mfrow = c(1, 1))
#Plot
boxplot(shannon ~ Feed.id, data=meta, ylab="Shannon's diversity")


##Non-normally distributed metrics - Categorical variables
kruskal.test(chao ~ Feed.id, data=meta)
pairwise.wilcox.test(meta$chao, meta$Feed.id, p.adjust.method="fdr")

#Create 1x1 plot environment
par(mfrow = c(1, 1))
#Plot
boxplot(chao ~ Feed.id, data=meta, ylab="Chao richness")





##BETA DIVERSITY
#Dot plots
BC.nmds = metaMDS(OTU.clean, distance="bray", k=2, trymax=1000)
levels(meta$Feed.id)
#Plot nMDS with different color
par(mfrow = c(1, 1))
#Create a blank plot for the nmds
plot(BC.nmds, type="n", main="Bray-Curtis")
#Add the points colored by age
points(BC.nmds, display="sites", pch=20, col=c("blue", "green", "red")[meta$Feed.id])
#Add a legend
legend(-5.5, 2.5, legend=c("MR1","MR2","CTM"), col=c("green","red","blue"), pch=20)

BC.nmds$stress

# a lower stress means the spatial distances in your NMDS more accurately represent your calculated similarites. <0.2 is basically required, <0.1 is better but uncommon. If its bad, maybe you need to recalculate the mds line above and use transformed data 
#autotransform = TRUE

# I like to merge my NMDS coordinates in together with my metadata to make one big dataframe, I think this makes plotting easier later on

nmds <-as.data.frame(BC.nmds$points)
metanmds <- merge(meta, nmds, by.x = 0, by.y = 0)
row.names(metanmds) <- metanmds[,1]
metanmds <- metanmds[,-1]
str(metanmds)
metanmds$Feed.id <- factor(metanmds$Feed.id)

NMDS.mean <- aggregate(metanmds[,7:8], list(group=metanmds$Feed.id), mean)
colnames(NMDS.mean) <- c('CTM','MR1','MR2')

metanmds <- merge(metanmds, NMDS.mean, by.x = "Feed.id", by.y="CTM")
str(metanmds)

ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Feed.id)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(BC.nmds$stress, digits = 2))) +
  stat_ellipse(aes(color=Feed.id), level = 0.95) +
  theme(legend.title = element_blank()) 


#Jaccard Index
J.nmds = metaMDS(OTU.clean, distance="jaccard", k=2, trymax=1000)

plot(J.nmds, type="n", main="Jaccard")
points(J.nmds, display="sites", pch=20, col=c("blue", "green", "red")[meta$Feed.id])
legend(-3, 1.5, legend=c("MR1","MR2","CTM"), col=c("green","red","blue"), pch=20)

##plot standard error ellipses 
plot(BC.nmds, type="n", main="Bray-Curtis")
legend(-5.5, 2.5, legend=c("MR1","MR2","CTM"), col=c("green","red","blue"), pch=20)
#Add an ellipse for MR1
ordiellipse(BC.nmds, groups=meta$Feed.id, display="sites", kind="se", conf=0.99, label=FALSE, col="green", draw="polygon", alpha=200, show.groups = c("MR1"), border=FALSE)
#Add an ellipse for MR2
ordiellipse(BC.nmds, groups=meta$Feed.id, display="sites", kind="se", conf=0.99, label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("MR2"), border=FALSE)
#Add an ellipse for CTM
ordiellipse(BC.nmds, groups=meta$Feed.id, display="sites", kind="se", conf=0.99, label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("CTM"), border=FALSE)

J.nmds$stress

nmds2 <-as.data.frame(J.nmds$points)
metanmds2 <- merge(meta, nmds2, by.x = 0, by.y = 0)
row.names(metanmds2) <- metanmds2[,1]
metanmds2 <- metanmds2[,-1]
str(metanmds2)
metanmds2$Feed.id <- factor(metanmds2$Feed.id)

NMDS2.mean <- aggregate(metanmds2[,7:8], list(group=metanmds2$Feed.id), mean)
colnames(NMDS2.mean) <- c('CTM','MR1','MR2')

metanmds2 <- merge(metanmds2, NMDS2.mean, by.x = "Feed.id", by.y="CTM")
str(metanmds2)

ggplot(metanmds2, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Feed.id)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(J.nmds$stress, digits = 2))) +
  stat_ellipse(aes(color=Feed.id), level = 0.95) +
  theme(legend.title = element_blank()) 

