## Code for streamplots: Fig. 2B, S2
## Anna M. Seekatz
## 4.8.16


# this is the complete code to recreate the PCoA comparisons done for Figures 2B, S2
# Overall description:
	# use relative abundance file created in Streamplot code
	# For exp 1/2 comparisons (Fig. S2), combine sample relative abundances and recalculate distance matrix
			
## files used:
	# experiment 1 (Fig. 2B):
		# t.m3_phylofrac_filtered.txt: relative abundance of phylotypes (fecal samples only)
	# experiment 1 and 2 (Fig. S2):
		# m3_metadata.txt
		# citro_metadata.txt
		# t.m3_phylofrac_filtered.txt
		# t.citro2_phylofrac.txt
		# exp1.2.comp.data.txt

###
#-----------
###

# Part I:  Create Figure 2B: PCoA of samples, by day of diet

library(labdsv)
library(plyr)
library(vegan)

# color palette:
diets.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "FF" ) {
colorvec[i] = "red"
}
if (n[i] == "FR" ) {
colorvec[i] = "chartreuse4"
}
if( n[i] == "pre") {
colorvec[i] = "grey67"
}
if( n[i] == "1_FR_FF") {
colorvec[i] = "blue3"
}
if( n[i] == "baseline") {
colorvec[i] = "chartreuse4"
}
}
c(colorvec)
}

# read files and add new matadata line:
phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)

# select only groups FF, FR, pre, 1day:
phylofrac<-phylofrac[phylofrac$group_ID %in% c("FF", "FR", "pre", "1_FR_FF"), ]
phylofrac<-droplevels(phylofrac)

# make a matrix out of the data to get the
phylo<-phylofrac[, c(10:23)]						#chooses only your relative abundance columns
rownames(phylo)<-phylofrac$sampleID
phylo<-as.matrix(phylo)

phylo.dist.bray<-vegdist(phylo, method="bray")

# pcoa (Fig2B):

bray.pco<-pco(phylo.dist.bray, k=10)
percent.variance<-round(100*bray.pco$eig/sum(bray.pco$eig))
percent.variance

#diet group (used in manuscript):
plot(bray.pco$points[,1], bray.pco$points[,2], bg=diets.col(phylofrac$group_ID), pch=21, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), xlab="PCOA 1 (75%)", ylab="PCOA 2 (9%)")
#text(bray.pco$points[,1], bray.pco$points[,2], phylofrac$mouse, cex=0.75)
legend("bottomright",legend=c("Fiber-free (FF)", "Fiber-rich (FR)", "prebiotic", "1-day FR/FF"), pt.bg=c("red", "chartreuse4","grey67", "blue3"), pch=21, col="black")

#day of diet (if wanted by diet actually fed that day):
plot(bray.pco$points[,1], bray.pco$points[,2], bg=diets.col(phylofrac$day_of_diet), pch=21, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), xlab="PCOA 1 (75%)", ylab="PCOA 2 (9%)")
#text(bray.pco$points[,1], bray.pco$points[,2], phylofrac$mouse, cex=0.75)
legend("bottomright",legend=c("Fiber-free (FF)", "Fiber-rich (FR)", "prebiotic"), pt.bg=c("red", "chartreuse4","grey67"), pch=21, col="black")
	

###
#-----------
###

# Part II:  Comparison of experiments and other experimental variables (sex, cages) - Fig. S2

# Since the PCOA was done according to phylotype, this should be pretty straight forward
	# (if OTUs had been used instead, would have had to regenerate a .shared file with ALL of the sequences)
	
# caveats:
	#1. Cage information might be different than noted in the meta data files before switching the mice from Fiber rich (FR) to Fiber free (FF) diet. 
		#However, for simplicity, we could only consider the cage information post diet switch (that is the info given in the meta data files).
	#2. We do not need Citrobacter infection samples in this PCOA plot. Thus, please exclude all samples after the infection (this refers only to Exp#2).
	#3. The first three batches (that is, first three dates) of fecal samples in both Exp#1 and Exp#2 for the FF group are actually FR, as the mice were switched to FF group just after collecting the third batch of fecal samples for both experiments (Exp#1 and Exp#2). 
		#Therefore, in the PCOA plot, please consider the first three FF batches as FR (as done in the old PCOA plot).
	# Mahesh also sent over newly designated information, added to the older metadata files:
		# citro_metadata_MD.xlsx
		# m3_metadata_MD.xlsx
		
library(data.table)
library(vegan)
library(labdsv)


		
# read in files, add new data, and merge:

exp1.meta<-read.table(file="m3_metadata.txt", header=TRUE)
exp2.meta<-read.table(file="citro_metadata.txt", header=TRUE)
exp1.phylo<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)
exp2.phylo<-read.table(file="t.citro2_phylofrac.txt", header=TRUE)

# based on information given:
exp1.meta$experiment<-"exp1"
exp2.meta$experiment<-"exp2"
exp1.meta$cage[exp1.meta$mouse %in% c(1,2)] <- 'c1'
exp1.meta$cage[exp1.meta$mouse %in% c(3,4)] <- 'c2'
exp1.meta$cage[exp1.meta$mouse %in% c(5,6)] <- 'c3'
exp1.meta$cage[exp1.meta$mouse %in% c(7,8)] <- 'c4'
exp1.meta$sex[exp1.meta$mouse %in% c(1,2,5,6)] <- 'male'
exp1.meta$sex[exp1.meta$mouse %in% c(3,4,7,8)] <- 'female'
exp2.meta$cage[exp2.meta$mouse %in% c(1,2,3)] <- 'c5'
exp2.meta$cage[exp2.meta$mouse %in% c(4,5,6,7)] <- 'c6'
exp2.meta$cage[exp2.meta$mouse %in% c(8)] <- 'c7'
exp2.meta$cage[exp2.meta$mouse %in% c(9)] <- 'c8'
exp2.meta$cage[exp2.meta$mouse %in% c(10)] <- 'c9'
exp2.meta$cage[exp2.meta$mouse %in% c(11,12,13)] <- 'c10'
exp2.meta$cage[exp2.meta$mouse %in% c(14)] <- 'c11'
exp2.meta$sex<-"male"

# merge metadata, but only for FR/FF groups in both experiments
exp1<-exp1.meta[exp1.meta$group_ID %in% c("FF", "FR"), c("sampleID", "mouse", "day_of_diet", "day", "experiment", "cage", "sex")]
exp2<-exp2.meta[exp2.meta$day < 53, c("sampleID", "mouse", "day_of_diet", "day", "experiment", "cage", "sex")]
exp2$mouse<-as.factor(exp2$mouse)
meta<-rbind(exp1, exp2)
#write.table(meta, file="mahesh_exp.comp.meta.txt", sep="\t", quote=FALSE, col.names=NA)
	# this file should have all the samples you want to compare

#organize/merge phylo, then merge with meta:
exp1<-exp1.phylo[, c("sampleID", "B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis", "E_rectale", "C_symbiosum","B_formatexigens", "R_intestinalis", "F_prausnitzii", "E_coli", "D_piger", "C_aerofaciens", "A_muciniphila")]
exp2<-exp2.phylo[, c("sampleID", "B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis", "E_rectale", "C_symbiosum","B_formatexigens", "R_intestinalis", "F_prausnitzii", "E_coli", "D_piger", "C_aerofaciens", "A_muciniphila")]
#setnames(exp1, "E_rectal", "E_rectale")
phylo<-rbind(exp1, exp2)
combined<-merge(meta, phylo, by.x=c("sampleID"), by.y=c("sampleID"))
combined$cagegroup<-paste(combined$day_of_diet, combined$cage, sep="_")
#write.table(combined, file="exp1.2.comp.data.txt", sep="\t", quote=FALSE, col.names=NA)
	# this file has the metadata plus phylotype information

--
## PCOA plots, using phylotype data instead of OTUs:
# this will give a better representation of the actual phylotypes present
# differentiation by regular OTUs is not as clean, since some Firmicutes do not split at 97% cutoff

combined<-read.table(file="exp1.2.comp.data.txt", header=TRUE)

# make a matrix out of the data to get the
phylo<-combined[, c(8:21)]						#chooses only your relative abundance columns
rownames(phylo)<-combined$sampleID
phylo<-as.matrix(phylo)

phylo.dist.bray<-vegdist(phylo, method="bray")

bray.pco<-pco(phylo.dist.bray, k=10)
percent.variance<-round(100*bray.pco$eig/sum(bray.pco$eig))
percent.variance
plot(bray.pco$points)
	# this looks good too

#let's add some color
group.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "FF" ) {
colorvec[i] = "red"
}
if (n[i] == "FR" ) {
colorvec[i] = "chartreuse4"
}
}
c(colorvec)
}

#display.brewer.pal(n = 5, name = 'Reds')
#display.brewer.pal(n = 5, name = 'Greens')
greens<-brewer.pal(n = 8, name = "Greens")
reds<-brewer.pal(n = 4, name = "Reds")
col.dist <- function(inp, comp) sum( abs(inp - col2rgb(comp) ) )
colors()[ apply(col2rgb(greens), 2, 
             function(z) which.min( sapply(colors(), 
                           function(x) col.dist(inp=z, comp=x) ) ) ) ]

cage.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "FR_c1" ) {
colorvec[i] = "gray97"
}
if (n[i] == "FR_c10" ) {
colorvec[i] = "honeydew2"
}
if ( n[i] == "FR_c11" ) {
colorvec[i] = "darkseagreen1"
}
if (n[i] == "FR_c2" ) {
colorvec[i] = "darkseagreen3"
}
if ( n[i] == "FR_c3" ) {
colorvec[i] = "palegreen3"
}
if (n[i] == "FR_c4" ) {
colorvec[i] = "mediumseagreen"
}
if ( n[i] == "FR_c7" ) {
colorvec[i] = "seagreen"
}
if (n[i] == "FR_c8" ) {
colorvec[i] = "darkgreen"
}
if (n[i] == "FR_c9" ) {
colorvec[i] = "darkslategrey"
}
if (n[i] == "FF_c1" ) {
colorvec[i] = "mistyrose"
}
if ( n[i] == "FF_c2" ) {
colorvec[i] = "burlywood2"
}
if (n[i] == "FF_c5" ) {
colorvec[i] = "tomato"
}
if (n[i] == "FF_c6" ) {
colorvec[i] = "firebrick3"
}
}
c(colorvec)
}
#c("FR_c1", "FR_c10", "FR_c11", "FR_c2", "FR_c3", "FR_c4", "FR_c7", "FR_c8", "FR_c9", "FF_c1", "FF_c2", "FF_c5", "FF_c6")
#c("gray97", "honeydew2", "darkseagreen1", "darkseagreen3", "palegreen3", "mediumseagreen", "seagreen", "darkgreen", "darkslategrey", "mistyrose", "burlywood2", "tomato", "firebrick3")

type.pch <- function(n) {
pchvec <- vector(mode="numeric", length=length(n))
for (i in 1:length(n)) {
pchvec[i] = 1
if ( n[i] == "exp1" ) {
pchvec[i] = 21
}
if (n[i] == "exp2" ) {
pchvec[i] = 24
}
}
c(pchvec)
}

sex.pch <- function(n) {
pchvec <- vector(mode="numeric", length=length(n))
for (i in 1:length(n)) {
pchvec[i] = 1
if ( n[i] == "male" ) {
pchvec[i] = 21
}
if (n[i] == "female" ) {
pchvec[i] = 24
}
}
c(pchvec)
}

# by experiment:
par(mfrow=c(1,2))
plot(bray.pco$points, col="black", bg=cage.col(combined$cagegroup), pch=type.pch(combined$experiment), xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), xlab="PCO 1 (80%)", ylab="PCO 2 (7%)", main="By Experiment")
legend("bottomright",legend=c("FR: cage 1", "FR: cage 10", "FR: cage 11", "FR: cage 2", "FR: cage 3", "FR: cage 4", "FR: cage 7", "FR: cage 8", "FR: cage 9", "FF: cage 1", "FF: cage 2", "FF: cage 5", "FF: cage 6"), col=c("gray97", "honeydew2", "darkseagreen1", "darkseagreen3", "palegreen3", "mediumseagreen", "seagreen", "darkgreen", "darkslategrey", "mistyrose", "burlywood2", "tomato", "firebrick3"), pch=19, cex=0.8)
legend("bottomleft",legend=c("Exp. 1", "Exp. 2"), pt.bg=c("white"), col=c("black"), pch=c(21, 24), cex=0.8)

# by sex:
plot(bray.pco$points, col="black", bg=cage.col(combined$cagegroup), pch=sex.pch(combined$sex), xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), xlab="PCO 1 (80%)", ylab="PCO 2 (7%)", main="By Sex")
legend("bottomright",legend=c("FR: cage 1", "FR: cage 10", "FR: cage 11", "FR: cage 2", "FR: cage 3", "FR: cage 4", "FR: cage 7", "FR: cage 8", "FR: cage 9", "FF: cage 1", "FF: cage 2", "FF: cage 5", "FF: cage 6"), col=c("gray97", "honeydew2", "darkseagreen1", "darkseagreen3", "palegreen3", "mediumseagreen", "seagreen", "darkgreen", "darkslategrey", "mistyrose", "burlywood2", "tomato", "firebrick3"), pch=19, cex=0.8)
legend("bottomleft",legend=c("Male", "Female"), pt.bg=c("white"), col=c("black"), pch=c(21, 24), cex=0.8)

	# a couple of the FF samples are hanging out in the FR group; let's check on their dates
text(bray.pco$points, rownames(bray.pco$points), cex=0.5)
	# looks like these are 18_May, 23_May samples...


### if wanted to do AMOVA:
library(pegas)
library(ape)

# follow the example below:
# http://www.inside-r.org/packages/cran/pegas/docs/amova
require(ape)
data(woodmouse)
d <- dist.dna(woodmouse)
g <- factor(c(rep("A", 7), rep("B", 8)))
p <- factor(c(rep(1, 3), rep(2, 4), rep(3, 4), rep(4, 4)))
amova(d ~ g/p, nperm = 100) # 2 levels
amova(d ~ p, nperm = 100) # 1 level
amova(d ~ g, nperm = 100)

a<-phylo.dist.bray
f1<-combined$cage
f2<-combined$day_of_diet
amova(a ~ f1, nperm = 100)
amova(a ~ f1/f2, nperm = 100)
