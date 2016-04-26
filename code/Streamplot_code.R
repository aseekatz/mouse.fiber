## Code for streamplots: Fig. 2A, S6
## Anna M. Seekatz
## 4.7.16


# this is the complete code to recreate community abundance streamplots (excuse the beastliness)
# Overall description:
	# generate a relative abundance file from the mothur-generated phylotype file (classified using a unique classifier for the community)
		# clean up samples < 5,000 seqs
		# condense community to 14 species (or 15, for C. rodentium data)
			# all classified sequences (including unknowns) were confirmed (BLAST)
			# note: some 'unknown' sequences pop up here and there
			# in general, these were either totally unknown (small portion), too close to one of the Bacteroides species (unresolved), or found to belong to Lactococcus
			# these sequences were not included in the analysis
			# Lactococcus sequences were only observed in days where __ mice were fed __ diet
			# sequences belonging to this classification did not persist past __ diet
			# it was determined that these sequences represented something found in the (autoclaved) diet, and were not viable
			# for more info, contact A. seekatz or E. Martens
			
## files used:
	# experiment 1 (Fig. 4):
		# martens3.trim.contigs.good.unique.good.filter.unique.precluster.pick.Mahesh_14_species_taxonomy.wang.tax.summary: created in mothur
		# m3_metadata.txt: contains all sample-related information
		# t.m3_phylofrac_filtered.txt (created): relative abundance of phylotypes (fecal samples only)
		# t.m3.cecal_phylofrac_filtered.txt (created): relative abundance of phylotypes (cecal samples only)
	# experiment 2 (Fig. S6):
		# citro2.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.Mahesh_15_species_taxonomy.wang.tax.summary: created in mothur
		# citro_metadata.txt: metadata of samples
		# t.citro2_phylofrac.txt (created)

###
#-----------
###

# part I: generate an abundance file for Experiment 1 data

library(plyr)

# read in file and condense to nonredundant phylotypes
# note: naming from classifier is a little off, so you might need to be creative when pulling the final list (can't just use taxlevel=7)

phylotype<-read.table(file="martens3.trim.contigs.good.unique.good.filter.unique.precluster.pick.Mahesh_14_species_taxonomy.wang.tax.summary.txt", header=TRUE)
phylo.final<-phylotype[phylotype$taxon %in% c("aerofaciens", "caccae", "ovatus", "thetaiotaomicron", "uniformis", "intestinihomis", "symbiosum", 
												"rectale", "formatexigens", "intestinalis", "prausnitzii", "Desulfovibrio", "Escherichia", "muciniphila"), ]
rownames(phylo.final)<-phylo.final$taxon
phylos<-phylo.final[ , -c(1:4)]

# add phylotype info, colors, and organize by phyla and abundance:
phylos$phylotype<-mapvalues(rownames(phylos), from = c("aerofaciens", "caccae", "ovatus", "thetaiotaomicron", "uniformis", "intestinihomis", "symbiosum", "rectale", "formatexigens", "intestinalis", "prausnitzii", "Desulfovibrio", "Escherichia", "muciniphila"), 
												to = c("C_aerofaciens","B_caccae", "B_ovatus", "B_theta", "B_uniformis", "B_intestihominis", "C_symbiosum", "E_rectale", "B_formatexigens", "R_intestinalis", "F_prausnitzii", "D_piger", "E_coli", "A_muciniphila"))

phylos$phylum[phylos$phylotype %in% c("B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis")]<-"1_Bacteroides"
phylos$phylum[phylos$phylotype %in% c("E_rectale", "C_symbiosum", "B_formatexigens", "R_intestinalis", "F_prausnitzii")]<-"2_Firmicutes"
phylos$phylum[phylos$phylotype %in% c("E_coli", "D_piger")]<-"3_Proteobacteria"
phylos$phylum[phylos$phylotype %in% c("C_aerofaciens")]<-"4_Actinobacteria"
phylos$phylum[phylos$phylotype %in% c("A_muciniphila")]<-"5_Verrucomicrobia"
phylos<-phylos[order(phylos$phylum, -phylos$total), ]
phylos<-phylos[order(factor(phylos$phylotype,levels=c(c("B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis", "E_rectale", "C_symbiosum", "B_formatexigens", "R_intestinalis", "F_prausnitzii", "E_coli", "D_piger", "C_aerofaciens", "A_muciniphila")))),]
phylos$color<-c("green4","chartreuse3","lawngreen","darkolivegreen2","palegreen3","midnightblue","dodgerblue4","blue","deepskyblue","slateblue1","gold","orange1","purple1","red")
rownames(phylos)<-phylos$phylotype

# transpose, clean up, create matrix, and transform into relative abundance
t.phylos<-t(phylos[,(2:680)])
f.phylos<-t.phylos[which(rowSums(t.phylos)>=1000),]	
phylo<-f.phylos/rowSums(f.phylos)*100	

#write.table(phylo, file="t.m3_phylofrac.txt", quote=FALSE, sep="\t", col.names=NA)

# merge with metadata:
phylo<-as.data.frame(phylo)
phylo$sampleID<-rownames(phylo)											#adds a column with sampleID
var<-read.table(file="m3_metadata.txt", header=TRUE)
phylo$sampleID<-gsub("X", "", phylo$sampleID)							# since seq. names started with a number, R added an X
m<-merge(var, phylo, by.x=c("sampleID"), by.y=c("sampleID")) 			#this merges the data based on the sampleID/group1 match

# check on number of samples per day since some samples were repeated:
num.samples <- ddply(m, c("group_ID", "day"), summarise,
               N    = length(B_theta))
sorted.n<-num.samples[order(num.samples$N),]
	# FF and FR and pre groups should have at most 4; others should have at most 3
	# unfortunately, must do this manually, as each case will be different
	
# groups that should only have 3:
oversampled3<-num.samples[num.samples$N>3 & num.samples$group_ID %in% c("1_FR_FF", "1_pre_FF", "4_FR_FF", "4_pre_FF", "pre"),]
	# some have definitely been oversampled
	# let's take a look at get rid of these first
FR1_samples<-m[m$group_ID==c("1_FR_FF") & m$day %in% c(18, 54), ]
	# d54 has both fecal and cecal, but Cecal12_pl2 is repeated
	# d18: 27OCT_12_pl2, 27OCT_13_pl2, 27OCT_14_pl2 are repeated
pre1_samples<-m[m$group_ID==c("1_pre_FF") & m$day %in% c(21, 54), ]
	# d54 ok
	# d21: 30OCT_18_pl2 repeated
FR4_samples<-m[m$group_ID==c("4_FR_FF") & m$day %in% c(18, 20, 54), ]
	# d18: 27OCT_15_pl2, , 27OCT_17_pl2, 27OCT_16_pl2 repeated
	# d20: 29OCT_17_pl2 repeated
pre4_samples<-m[m$group_ID==c("4_pre_FF") & m$day %in% c(22, 54), ]
	# d22: 31OCT_23_pl2 repeated
pre_samples<-m[m$group_ID==c("pre") & m$day %in% c(18, 54), ]
	# d54, cecal: Cecal10_pl2, Cecal11_pl2 repeated
	# d18: 27OCT_10_pl2, 27OCT_11_pl2 repeated
	
# groups that should have 4:
oversampled4<-num.samples[num.samples$N>4 & num.samples$group_ID %in% c("FR", "FF"),]
FR_samples<-m[m$group_ID==c("FR") & m$day %in% c(21, 54), ]
	# d21: 30OCT_6_pl2 (not sure what sample this is)
	# d54: ok
FF_samples<-m[m$group_ID==c("FF") & m$day %in% c(54), ]  #ok

# remove repeated samples so that only one representative sample for that animal is included (not pretty code)
# note: some samples had already been removed during mothur
# kept plate 1 sequences
# afterwards, can repeat above code to make sure all samples were removed
filtered.m<-subset(m, !(sampleID %in% c("Cecal12_pl2", "27OCT_12_pl2", "27OCT_13_pl2", "27OCT_14_pl2", "30OCT_18_pl2", "31OCT_23_pl2", "27OCT_15_pl2", "27OCT_16_pl2", "27OCT_17_pl2", "27OCT_17_plr", "Cecal10_pl2", "Cecal11_pl2", "27OCT_10_pl2", "27OCT_11_pl2", "30OCT_6_pl2", "29OCT_17_pl2")))

# remove other control samples:
levels(filtered.m$group_ID)											#this lists groups left--we also want to remove the control samples
f2.m <- filtered.m[!grepl("10",filtered.m$group_ID),]
f3.m <- f2.m[!grepl("other",f2.m$group_ID),]
phylofrac<-f3.m[order(f3.m$group_ID, as.numeric_version(f3.m$day)),]	# sorts the table by group, then by day
phylofrac[, c("sampleID", "mouse", "day")]
	# scanning through, it looks like there are no doubly-sampled mice per day

# separate fecal samples:
phylofrac.cecal <- phylofrac[grepl("Cecal",phylofrac$sampleID),]
phylofrac.fecal <- phylofrac[!grepl("Cecal",phylofrac$sampleID),]

# create tables for future use:
write.table(phylofrac.fecal, file="t.m3_phylofrac_filtered.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(phylofrac.cecal, file="t.m3.cecal_phylofrac_filtered.txt", quote=FALSE, sep="\t", col.names=NA)

#### You now have a data.frame (phylo.frac) that has been filtered of all other samples



###
#-----------
###

# Part II:  Creating Fig. 2A: Streamplots
# note: This code is pretty beastly--any help on making it easier/more efficient is appreciated!

# Stream plots of FF, FR and pre:
		
# in R:

library(plyr)

# to get initial definitions:
phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)

mean.phylofrac<-ddply(phylofrac, c("group_ID", "day"), colwise(mean, is.numeric))

m.split<-split(mean.phylofrac, mean.phylofrac$group_ID)
mean.FF<-m.split$'FF'
mean.FR<-m.split$'FR'
mean.pre<-m.split$'pre'
mean.1_FR_FF<-m.split$'1_FR_FF'
mean.4_FR_FF<-m.split$'4_FR_FF'
mean.4_pre_FF<-m.split$'4_pre_FF'
mean.1_pre_FF<-m.split$'1_pre_FF'

ylim <- c(0, 100)
par(xpd=TRUE)
layout(matrix(1:4, ncol=4), widths=1, respect=FALSE)
par(mfrow=c(2,3), mai = c(0.5, 0.2, 0.1, 0.1), oma=c(1,4,0,4))

## for FR:

mean.FR<-mean.FR[order(as.numeric_version(mean.FR$day)),]
mean.FR$day <- as.numeric(as.character(mean.FR$day))

x1<-mean.FR$day[c(1:15)]			#the first 15 rows represent days 6-26
x2<-mean.FR$day[c(16:28)]			#the last 13 rows represent the end days
xx.p1<-c(x1, rev(x1))
xx.p2<-c(x2, rev(x2))
yyB_theta.p1 <- c(rep(0, nrow(mean.FR[1:15,])), rev(mean.FR$B_theta[1:15]))
yyB_theta.p2 <- c(rep(0, nrow(mean.FR[16:28,])), rev(mean.FR$B_theta[16:28]))
plot(x=mean.FR$day, y=mean.FR$B_theta, ylim=ylim, col='darkgreen', type='l', xaxt='n',
ylab='Relative abundance', xlab='Day', main='FR group', cex.lab=0.8, cex.axis=0.8)
polygon(xx.p1, yyB_theta.p1, col='darkgreen')
polygon(xx.p2, yyB_theta.p2, col='darkgreen')
yyB_ovatus.p1 <- c(mean.FR$B_theta[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]))
yyB_ovatus.p2 <- c(mean.FR$B_theta[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]))
polygon(xx.p1, yyB_ovatus.p1, col='chartreuse2')
polygon(xx.p2, yyB_ovatus.p2, col='chartreuse2')
yyB_uniformis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]))
yyB_uniformis.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]))
polygon(xx.p1, yyB_uniformis.p1, col='green')
polygon(xx.p2, yyB_uniformis.p2, col='green')
yyB_caccae.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]))
yyB_caccae.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]))
polygon(xx.p1, yyB_caccae.p1, col='limegreen')
polygon(xx.p2, yyB_caccae.p2, col='limegreen')
yyB_intestihominis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]))
yyB_intestihominis.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]))
polygon(xx.p1, yyB_intestihominis.p1, col='seagreen')
polygon(xx.p2, yyB_intestihominis.p2, col='seagreen')
yyE_rectale.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]))
yyE_rectale.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]))
polygon(xx.p1, yyE_rectale.p1, col='midnightblue')
polygon(xx.p2, yyE_rectale.p2, col='midnightblue')
yyC_symbiosum.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]))
yyC_symbiosum.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]))
polygon(xx.p1, yyC_symbiosum.p1, col='dodgerblue4')
polygon(xx.p2, yyC_symbiosum.p2, col='dodgerblue4')
yyB_formatexigens.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]))
yyB_formatexigens.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]))
polygon(xx.p1, yyB_formatexigens.p1, col='blue')
polygon(xx.p2, yyB_formatexigens.p2, col='blue')
yyR_intestinalis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]))
yyR_intestinalis.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]))
polygon(xx.p1, yyR_intestinalis.p1, col='deepskyblue')
polygon(xx.p2, yyR_intestinalis.p2, col='deepskyblue')
yyF_prausnitzii.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]))
yyF_prausnitzii.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28] + mean.FR$R_intestinalis[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]) + rev(mean.FR$F_prausnitzii[16:28]))
polygon(xx.p1, yyF_prausnitzii.p1, col='slateblue1')
polygon(xx.p2, yyF_prausnitzii.p2, col='slateblue1')
yyE_coli.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]))
yyE_coli.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28] + mean.FR$R_intestinalis[16:28] + mean.FR$F_prausnitzii[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]) + rev(mean.FR$F_prausnitzii[16:28]) + rev( mean.FR$E_coli[16:28]))
polygon(xx.p1, yyE_coli.p1, col='gold')
polygon(xx.p2, yyE_coli.p2, col='gold')
yyD_piger.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]))
yyD_piger.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28] + mean.FR$R_intestinalis[16:28] + mean.FR$F_prausnitzii[16:28] + mean.FR$E_coli[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]) + rev(mean.FR$F_prausnitzii[16:28]) + rev( mean.FR$E_coli[16:28]) + rev(mean.FR$D_piger[16:28]))
polygon(xx.p1, yyD_piger.p1, col='orange1')
polygon(xx.p2, yyD_piger.p2, col='orange1')
yyC_aerofaciens.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15] + mean.FR$D_piger[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]) + rev(mean.FR$C_aerofaciens[1:15]))
yyC_aerofaciens.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28] + mean.FR$R_intestinalis[16:28] + mean.FR$F_prausnitzii[16:28] + mean.FR$E_coli[16:28] + mean.FR$D_piger[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]) + rev(mean.FR$F_prausnitzii[16:28]) + rev( mean.FR$E_coli[16:28]) + rev(mean.FR$D_piger[16:28]) + rev(mean.FR$C_aerofaciens[16:28]))
polygon(xx.p1, yyC_aerofaciens.p1, col='purple1')
polygon(xx.p2, yyC_aerofaciens.p2, col='purple1')
yyA_muciniphila.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15] + mean.FR$D_piger[1:15] + mean.FR$C_aerofaciens[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]) + rev(mean.FR$C_aerofaciens[1:15]) + rev(mean.FR$A_muciniphila[1:15]))
yyA_muciniphila.p2 <- c(mean.FR$B_theta[16:28] + mean.FR$B_ovatus[16:28] + mean.FR$B_uniformis[16:28] + mean.FR$B_caccae[16:28] +mean.FR$B_intestihominis[16:28] + mean.FR$E_rectale[16:28] + mean.FR$C_symbiosum[16:28] +mean.FR$B_formatexigens[16:28] + mean.FR$R_intestinalis[16:28] + mean.FR$F_prausnitzii[16:28] + mean.FR$E_coli[16:28] + mean.FR$D_piger[16:28] + mean.FR$C_aerofaciens[16:28], rev(mean.FR$B_theta[16:28]) + rev(mean.FR$B_ovatus[16:28]) + rev(mean.FR$B_uniformis[16:28]) + rev(mean.FR$B_caccae[16:28]) + rev(mean.FR$B_intestihominis[16:28]) + rev(mean.FR$E_rectale[16:28]) + rev(mean.FR$C_symbiosum[16:28]) + rev(mean.FR$B_formatexigens[16:28]) + rev(mean.FR$R_intestinalis[16:28]) + rev(mean.FR$F_prausnitzii[16:28]) + rev( mean.FR$E_coli[16:28]) + rev(mean.FR$D_piger[16:28]) + rev(mean.FR$C_aerofaciens[16:28]) + rev(mean.FR$A_muciniphila[16:28]))
polygon(xx.p1, yyA_muciniphila.p1, col='red3')
polygon(xx.p2, yyA_muciniphila.p2, col='red3')
#legend(x=max(mean.FR$day), y=100, c(colnames(mean.FR[,5:18])), fill=c('darkgreen', 'chartreuse2', 'green', 'limegreen', 'seagreen', 'midnightblue', 'dodgerblue4', 'blue', 'deepskyblue', 'slateblue1', 'gold', 'orange1', 'purple1', 'red3'), cex=0.6)
axis(1, at=mean.FR$day, labels=mean.FR$day, cex.lab=0.8, cex.axis=0.8)


## for FF:

mean.FF<-mean.FF[order(as.numeric_version(mean.FF$day)),]
mean.FF$day <- as.numeric(as.character(mean.FF$day))

x1<-mean.FF$day[c(1:15)]			#the first 15 rows represent days 6-26
x2<-mean.FF$day[c(16:28)]			#the last 13 rows represent the end days
xx.p1<-c(x1, rev(x1))
xx.p2<-c(x2, rev(x2))
yyB_theta.p1 <- c(rep(0, nrow(mean.FF[1:15,])), rev(mean.FF$B_theta[1:15]))
yyB_theta.p2 <- c(rep(0, nrow(mean.FF[16:28,])), rev(mean.FF$B_theta[16:28]))
plot(x=mean.FF$day, y=mean.FF$B_theta, ylim=ylim, col='darkgreen', type='l', xaxt='n',
ylab='Relative abundance', xlab='Day', main='FF group', cex.lab=0.8, cex.axis=0.8)
polygon(xx.p1, yyB_theta.p1, col='darkgreen')
polygon(xx.p2, yyB_theta.p2, col='darkgreen')
yyB_ovatus.p1 <- c(mean.FF$B_theta[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]))
yyB_ovatus.p2 <- c(mean.FF$B_theta[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]))
polygon(xx.p1, yyB_ovatus.p1, col='chartreuse2')
polygon(xx.p2, yyB_ovatus.p2, col='chartreuse2')
yyB_uniformis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]))
yyB_uniformis.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]))
polygon(xx.p1, yyB_uniformis.p1, col='green')
polygon(xx.p2, yyB_uniformis.p2, col='green')
yyB_caccae.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]))
yyB_caccae.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]))
polygon(xx.p1, yyB_caccae.p1, col='limegreen')
polygon(xx.p2, yyB_caccae.p2, col='limegreen')
yyB_intestihominis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]))
yyB_intestihominis.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]))
polygon(xx.p1, yyB_intestihominis.p1, col='seagreen')
polygon(xx.p2, yyB_intestihominis.p2, col='seagreen')
yyE_rectale.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]))
yyE_rectale.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]))
polygon(xx.p1, yyE_rectale.p1, col='midnightblue')
polygon(xx.p2, yyE_rectale.p2, col='midnightblue')
yyC_symbiosum.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]))
yyC_symbiosum.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]))
polygon(xx.p1, yyC_symbiosum.p1, col='dodgerblue4')
polygon(xx.p2, yyC_symbiosum.p2, col='dodgerblue4')
yyB_formatexigens.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]))
yyB_formatexigens.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]))
polygon(xx.p1, yyB_formatexigens.p1, col='blue')
polygon(xx.p2, yyB_formatexigens.p2, col='blue')
yyR_intestinalis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]))
yyR_intestinalis.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]))
polygon(xx.p1, yyR_intestinalis.p1, col='deepskyblue')
polygon(xx.p2, yyR_intestinalis.p2, col='deepskyblue')
yyF_prausnitzii.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]))
yyF_prausnitzii.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28] + mean.FF$R_intestinalis[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]) + rev(mean.FF$F_prausnitzii[16:28]))
polygon(xx.p1, yyF_prausnitzii.p1, col='slateblue1')
polygon(xx.p2, yyF_prausnitzii.p2, col='slateblue1')
yyE_coli.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]))
yyE_coli.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28] + mean.FF$R_intestinalis[16:28] + mean.FF$F_prausnitzii[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]) + rev(mean.FF$F_prausnitzii[16:28]) + rev( mean.FF$E_coli[16:28]))
polygon(xx.p1, yyE_coli.p1, col='gold')
polygon(xx.p2, yyE_coli.p2, col='gold')
yyD_piger.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]))
yyD_piger.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28] + mean.FF$R_intestinalis[16:28] + mean.FF$F_prausnitzii[16:28] + mean.FF$E_coli[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]) + rev(mean.FF$F_prausnitzii[16:28]) + rev( mean.FF$E_coli[16:28]) + rev(mean.FF$D_piger[16:28]))
polygon(xx.p1, yyD_piger.p1, col='orange1')
polygon(xx.p2, yyD_piger.p2, col='orange1')
yyC_aerofaciens.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15] + mean.FF$D_piger[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]) + rev(mean.FF$C_aerofaciens[1:15]))
yyC_aerofaciens.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28] + mean.FF$R_intestinalis[16:28] + mean.FF$F_prausnitzii[16:28] + mean.FF$E_coli[16:28] + mean.FF$D_piger[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]) + rev(mean.FF$F_prausnitzii[16:28]) + rev( mean.FF$E_coli[16:28]) + rev(mean.FF$D_piger[16:28]) + rev(mean.FF$C_aerofaciens[16:28]))
polygon(xx.p1, yyC_aerofaciens.p1, col='purple1')
polygon(xx.p2, yyC_aerofaciens.p2, col='purple1')
yyA_muciniphila.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15] + mean.FF$D_piger[1:15] + mean.FF$C_aerofaciens[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]) + rev(mean.FF$C_aerofaciens[1:15]) + rev(mean.FF$A_muciniphila[1:15]))
yyA_muciniphila.p2 <- c(mean.FF$B_theta[16:28] + mean.FF$B_ovatus[16:28] + mean.FF$B_uniformis[16:28] + mean.FF$B_caccae[16:28] +mean.FF$B_intestihominis[16:28] + mean.FF$E_rectale[16:28] + mean.FF$C_symbiosum[16:28] +mean.FF$B_formatexigens[16:28] + mean.FF$R_intestinalis[16:28] + mean.FF$F_prausnitzii[16:28] + mean.FF$E_coli[16:28] + mean.FF$D_piger[16:28] + mean.FF$C_aerofaciens[16:28], rev(mean.FF$B_theta[16:28]) + rev(mean.FF$B_ovatus[16:28]) + rev(mean.FF$B_uniformis[16:28]) + rev(mean.FF$B_caccae[16:28]) + rev(mean.FF$B_intestihominis[16:28]) + rev(mean.FF$E_rectale[16:28]) + rev(mean.FF$C_symbiosum[16:28]) + rev(mean.FF$B_formatexigens[16:28]) + rev(mean.FF$R_intestinalis[16:28]) + rev(mean.FF$F_prausnitzii[16:28]) + rev( mean.FF$E_coli[16:28]) + rev(mean.FF$D_piger[16:28]) + rev(mean.FF$C_aerofaciens[16:28]) + rev(mean.FF$A_muciniphila[16:28]))
polygon(xx.p1, yyA_muciniphila.p1, col='red3')
polygon(xx.p2, yyA_muciniphila.p2, col='red3')
#legend(x=max(mean.FF$day), y=100, c(colnames(mean.FF[,5:18])), fill=c('darkgreen', 'chartreuse2', 'green', 'limegreen', 'seagreen', 'midnightblue', 'dodgerblue4', 'blue', 'deepskyblue', 'slateblue1', 'gold', 'orange1', 'purple1', 'red3'), cex=0.6)
axis(1, at=mean.FF$day, labels=mean.FF$day, cex.lab=0.8, cex.axis=0.8)

## for pre:

mean.pre<-mean.pre[order(as.numeric_version(mean.pre$day)),]
mean.pre$day <- as.numeric(as.character(mean.pre$day))

x1<-mean.pre$day[c(1:15)]			#the first 15 rows represent days 6-26
x2<-mean.pre$day[c(16:28)]			#the last 13 rows represent the end days
xx.p1<-c(x1, rev(x1))
xx.p2<-c(x2, rev(x2))
yyB_theta.p1 <- c(rep(0, nrow(mean.pre[1:15,])), rev(mean.pre$B_theta[1:15]))
yyB_theta.p2 <- c(rep(0, nrow(mean.pre[16:28,])), rev(mean.pre$B_theta[16:28]))
plot(x=mean.pre$day, y=mean.pre$B_theta, ylim=ylim, col='darkgreen', type='l', xaxt='n',
ylab='Relative abundance', xlab='Day', main='pre group', cex.lab=0.8, cex.axis=0.8)
polygon(xx.p1, yyB_theta.p1, col='darkgreen')
polygon(xx.p2, yyB_theta.p2, col='darkgreen')
yyB_ovatus.p1 <- c(mean.pre$B_theta[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]))
yyB_ovatus.p2 <- c(mean.pre$B_theta[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]))
polygon(xx.p1, yyB_ovatus.p1, col='chartreuse2')
polygon(xx.p2, yyB_ovatus.p2, col='chartreuse2')
yyB_uniformis.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]))
yyB_uniformis.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]))
polygon(xx.p1, yyB_uniformis.p1, col='green')
polygon(xx.p2, yyB_uniformis.p2, col='green')
yyB_caccae.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]))
yyB_caccae.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]))
polygon(xx.p1, yyB_caccae.p1, col='limegreen')
polygon(xx.p2, yyB_caccae.p2, col='limegreen')
yyB_intestihominis.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]))
yyB_intestihominis.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]))
polygon(xx.p1, yyB_intestihominis.p1, col='seagreen')
polygon(xx.p2, yyB_intestihominis.p2, col='seagreen')
yyE_rectale.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]))
yyE_rectale.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]))
polygon(xx.p1, yyE_rectale.p1, col='midnightblue')
polygon(xx.p2, yyE_rectale.p2, col='midnightblue')
yyC_symbiosum.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]))
yyC_symbiosum.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]))
polygon(xx.p1, yyC_symbiosum.p1, col='dodgerblue4')
polygon(xx.p2, yyC_symbiosum.p2, col='dodgerblue4')
yyB_formatexigens.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]))
yyB_formatexigens.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]))
polygon(xx.p1, yyB_formatexigens.p1, col='blue')
polygon(xx.p2, yyB_formatexigens.p2, col='blue')
yyR_intestinalis.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]))
yyR_intestinalis.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]))
polygon(xx.p1, yyR_intestinalis.p1, col='deepskyblue')
polygon(xx.p2, yyR_intestinalis.p2, col='deepskyblue')
yyF_prausnitzii.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15] + mean.pre$R_intestinalis[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]) + rev(mean.pre$F_prausnitzii[1:15]))
yyF_prausnitzii.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28] + mean.pre$R_intestinalis[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]) + rev(mean.pre$F_prausnitzii[16:28]))
polygon(xx.p1, yyF_prausnitzii.p1, col='slateblue1')
polygon(xx.p2, yyF_prausnitzii.p2, col='slateblue1')
yyE_coli.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15] + mean.pre$R_intestinalis[1:15] + mean.pre$F_prausnitzii[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]) + rev(mean.pre$F_prausnitzii[1:15]) + rev( mean.pre$E_coli[1:15]))
yyE_coli.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28] + mean.pre$R_intestinalis[16:28] + mean.pre$F_prausnitzii[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]) + rev(mean.pre$F_prausnitzii[16:28]) + rev( mean.pre$E_coli[16:28]))
polygon(xx.p1, yyE_coli.p1, col='gold')
polygon(xx.p2, yyE_coli.p2, col='gold')
yyD_piger.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15] + mean.pre$R_intestinalis[1:15] + mean.pre$F_prausnitzii[1:15] + mean.pre$E_coli[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]) + rev(mean.pre$F_prausnitzii[1:15]) + rev( mean.pre$E_coli[1:15]) + rev(mean.pre$D_piger[1:15]))
yyD_piger.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28] + mean.pre$R_intestinalis[16:28] + mean.pre$F_prausnitzii[16:28] + mean.pre$E_coli[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]) + rev(mean.pre$F_prausnitzii[16:28]) + rev( mean.pre$E_coli[16:28]) + rev(mean.pre$D_piger[16:28]))
polygon(xx.p1, yyD_piger.p1, col='orange1')
polygon(xx.p2, yyD_piger.p2, col='orange1')
yyC_aerofaciens.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15] + mean.pre$R_intestinalis[1:15] + mean.pre$F_prausnitzii[1:15] + mean.pre$E_coli[1:15] + mean.pre$D_piger[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]) + rev(mean.pre$F_prausnitzii[1:15]) + rev( mean.pre$E_coli[1:15]) + rev(mean.pre$D_piger[1:15]) + rev(mean.pre$C_aerofaciens[1:15]))
yyC_aerofaciens.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28] + mean.pre$R_intestinalis[16:28] + mean.pre$F_prausnitzii[16:28] + mean.pre$E_coli[16:28] + mean.pre$D_piger[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]) + rev(mean.pre$F_prausnitzii[16:28]) + rev( mean.pre$E_coli[16:28]) + rev(mean.pre$D_piger[16:28]) + rev(mean.pre$C_aerofaciens[16:28]))
polygon(xx.p1, yyC_aerofaciens.p1, col='purple1')
polygon(xx.p2, yyC_aerofaciens.p2, col='purple1')
yyA_muciniphila.p1 <- c(mean.pre$B_theta[1:15] + mean.pre$B_ovatus[1:15] + mean.pre$B_uniformis[1:15] + mean.pre$B_caccae[1:15] +mean.pre$B_intestihominis[1:15] + mean.pre$E_rectale[1:15] + mean.pre$C_symbiosum[1:15] +mean.pre$B_formatexigens[1:15] + mean.pre$R_intestinalis[1:15] + mean.pre$F_prausnitzii[1:15] + mean.pre$E_coli[1:15] + mean.pre$D_piger[1:15] + mean.pre$C_aerofaciens[1:15], rev(mean.pre$B_theta[1:15]) + rev(mean.pre$B_ovatus[1:15]) + rev(mean.pre$B_uniformis[1:15]) + rev(mean.pre$B_caccae[1:15]) + rev(mean.pre$B_intestihominis[1:15]) + rev(mean.pre$E_rectale[1:15]) + rev(mean.pre$C_symbiosum[1:15]) + rev(mean.pre$B_formatexigens[1:15]) + rev(mean.pre$R_intestinalis[1:15]) + rev(mean.pre$F_prausnitzii[1:15]) + rev( mean.pre$E_coli[1:15]) + rev(mean.pre$D_piger[1:15]) + rev(mean.pre$C_aerofaciens[1:15]) + rev(mean.pre$A_muciniphila[1:15]))
yyA_muciniphila.p2 <- c(mean.pre$B_theta[16:28] + mean.pre$B_ovatus[16:28] + mean.pre$B_uniformis[16:28] + mean.pre$B_caccae[16:28] +mean.pre$B_intestihominis[16:28] + mean.pre$E_rectale[16:28] + mean.pre$C_symbiosum[16:28] +mean.pre$B_formatexigens[16:28] + mean.pre$R_intestinalis[16:28] + mean.pre$F_prausnitzii[16:28] + mean.pre$E_coli[16:28] + mean.pre$D_piger[16:28] + mean.pre$C_aerofaciens[16:28], rev(mean.pre$B_theta[16:28]) + rev(mean.pre$B_ovatus[16:28]) + rev(mean.pre$B_uniformis[16:28]) + rev(mean.pre$B_caccae[16:28]) + rev(mean.pre$B_intestihominis[16:28]) + rev(mean.pre$E_rectale[16:28]) + rev(mean.pre$C_symbiosum[16:28]) + rev(mean.pre$B_formatexigens[16:28]) + rev(mean.pre$R_intestinalis[16:28]) + rev(mean.pre$F_prausnitzii[16:28]) + rev( mean.pre$E_coli[16:28]) + rev(mean.pre$D_piger[16:28]) + rev(mean.pre$C_aerofaciens[16:28]) + rev(mean.pre$A_muciniphila[16:28]))
polygon(xx.p1, yyA_muciniphila.p1, col='red3')
polygon(xx.p2, yyA_muciniphila.p2, col='red3')
#legend(x=max(mean.pre$day), y=100, c(colnames(mean.pre[,5:18])), fill=c('darkgreen', 'chartreuse2', 'green', 'limegreen', 'seagreen', 'midnightblue', 'dodgerblue4', 'blue', 'deepskyblue', 'slateblue1', 'gold', 'orange1', 'purple1', 'red3'), cex=0.6)
axis(1, at=mean.pre$day, labels=mean.pre$day, cex.lab=0.8, cex.axis=0.8)

#axis(1, at=mean.pre$day, labels=mean.pre$day)

# if you want line segments of the diet timeline:
#segments(6, -2, 15, -2, lwd=4, col="black")
#segments(15, -2, 54, -2, lwd=4, col="purple4")

plot(1, type="n", axes=FALSE, xlab="", ylab="")			#blank spot for legend...
legend("topleft", c(colnames(mean.FF[,5:18])), fill=c('darkgreen', 'chartreuse2', 'green', 'limegreen', 'seagreen', 'midnightblue', 'dodgerblue4', 'blue', 'deepskyblue', 'slateblue1', 'gold', 'orange1', 'purple1', 'red3'), cex=0.8)
#legend("bottomleft", c("Fiber-rich", "Fiber-free", "Prebiotic"), col=c("black", "darkred", "purple4"), lty=1, lwd=5)
mtext(side=2, c("Relative abundance"), outer=TRUE, line=2, cex=0.8, adj=0.68)
mtext(side=1, c("day"), outer=TRUE, line=-1, cex=0.8, adj=0.4)


###
#-----------
###

# Part III:  Create an abundance file for Experiment 2 data

# read in file and condense to nonredundant phylotypes
# note: naming from classifier is a little off, so you might need to be creative when pulling the final list (can't just use taxlevel=7)

phylotype<-read.table(file="citro2.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.Mahesh_15_species_taxonomy.wang.tax.summary", header=TRUE)
phylo.final<-phylotype[phylotype$taxon %in% c("aerofaciens", "caccae", "ovatus", "thetaiotaomicron", "uniformis", "intestinihomis", "symbiosum", 
												"rectale", "formatexigens", "intestinalis", "prausnitzii", "Desulfovibrio", "Escherichia", "rodentium", "muciniphila"), ]
rownames(phylo.final)<-phylo.final$taxon
phylos<-phylo.final[ , -c(1:4)]

# add phylotype info, colors, and organize by phyla and abundance:
phylos$phylotype<-mapvalues(rownames(phylos), from = c("aerofaciens", "caccae", "ovatus", "thetaiotaomicron", "uniformis", "intestinihomis", "symbiosum", "rectale", "formatexigens", "intestinalis", "prausnitzii", "Desulfovibrio", "rodentium", "Escherichia", "muciniphila"), 
												to = c("C_aerofaciens","B_caccae", "B_ovatus", "B_theta", "B_uniformis", "B_intestihominis", "C_symbiosum", "E_rectale", "B_formatexigens", "R_intestinalis", "F_prausnitzii", "D_piger", "C_rodentium", "E_coli", "A_muciniphila"))

phylos$phylum[phylos$phylotype %in% c("B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis")]<-"1_Bacteroides"
phylos$phylum[phylos$phylotype %in% c("E_rectale", "C_symbiosum", "B_formatexigens", "R_intestinalis", "F_prausnitzii")]<-"2_Firmicutes"
phylos$phylum[phylos$phylotype %in% c("E_coli", "D_piger", "C_rodentium")]<-"3_Proteobacteria"
phylos$phylum[phylos$phylotype %in% c("C_aerofaciens")]<-"4_Actinobacteria"
phylos$phylum[phylos$phylotype %in% c("A_muciniphila")]<-"5_Verrucomicrobia"
#phylos<-phylos[order(phylos$phylum, -phylos$total), ]
phylos<-phylos[order(factor(phylos$phylotype,levels=c("B_theta", "B_ovatus", "B_uniformis", "B_caccae", "B_intestihominis", "E_rectale", "C_symbiosum", "B_formatexigens", "R_intestinalis", "F_prausnitzii", "E_coli", "D_piger", "C_aerofaciens", "A_muciniphila", "C_rodentium"))),]
phylos$color<-c("green4","chartreuse3","lawngreen","darkolivegreen2","palegreen3","midnightblue","dodgerblue4","blue","deepskyblue","slateblue1","gold","orange1","purple1","red", "pink")
rownames(phylos)<-phylos$phylotype

# convert to relative abundance:
t.phylos<-t(phylos[,(2:466)])
f.phylo<-t.phylos[which(rowSums(t.phylos)>=1000),]			
phylo<-f.phylo/rowSums(f.phylo)*100

# merge with metadata:
phylo<-as.data.frame(phylo)
phylo$sampleID<-rownames(phylo)											#adds a column with sampleID
var<-read.table(file="citro_metadata.txt", header=TRUE)
phylo$sampleID<-gsub("X", "", phylo$sampleID)	
m<-merge(var, phylo, by.x=c("sampleID"), by.y=c("sampleID")) 			#this merges the data based on the sampleID/group1 match
	# this should eliminate all control samples
	# none of the samples in this experiment were re-sequenced, so there should be no duplicates

#write.table(m, file="t.citro2_phylofrac.txt", quote=FALSE, sep="\t", col.names=NA)

#### You now have a data.frame (phylo.frac) that has been filtered of all other samples

###
#-----------
###

# Part IV:  Creating Fig. 2A: Streamplots

library(plyr)

phylofrac<-read.table(file="t.citro2_phylofrac.txt", header=TRUE)
mean.phylofrac<-ddply(phylofrac, c("day_of_diet", "day"), colwise(mean, is.numeric))

# then, you can split these into their respective groups:
m.split<-split(mean.phylofrac, mean.phylofrac$day_of_diet)
mean.FF<-m.split$'FF'
mean.FR<-m.split$'FR'

# order as numeric
mean.FF<-mean.FF[order(as.numeric_version(mean.FF$day)),]
mean.FF$day <- as.numeric(as.character(mean.FF$day))

# graph it:
ylim <- c(0, 100)
par(xpd=TRUE)
layout(matrix(1:4, ncol=4), widths=1, respect=FALSE)
par(mfrow=c(1,3), mai = c(0.5, 0.1, 0.1, 0.1), oma=c(1,4,0,4))

x1<-mean.FF$day[c(1:15)]			#the first 15 rows represent days 6-26
x2<-mean.FF$day[c(16:38)]			#the last 13 rows represent the end days
xx.p1<-c(x1, rev(x1))
xx.p2<-c(x2, rev(x2))

# FR:
mean.FR<-mean.FR[order(as.numeric_version(mean.FR$day)),]
mean.FR$day <- as.numeric(as.character(mean.FR$day))
yyB_theta.p1 <- c(rep(0, nrow(mean.FR[1:15,])), rev(mean.FR$B_theta[1:15]))
yyB_theta.p2 <- c(rep(0, nrow(mean.FR[16:38,])), rev(mean.FR$B_theta[16:38]))
plot(x=mean.FR$day, y=mean.FR$B_theta, ylim=ylim, col='darkgreen', type='l', xaxt='n',
ylab=NA, xlab=NA, main='FR group')
polygon(xx.p1, yyB_theta.p1, col='darkgreen')
polygon(xx.p2, yyB_theta.p2, col='darkgreen')
	# then, repeat the layering of the phylotypes individually:
yyB_ovatus.p1 <- c(mean.FR$B_theta[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]))
yyB_ovatus.p2 <- c(mean.FR$B_theta[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]))
polygon(xx.p1, yyB_ovatus.p1, col='chartreuse4')
polygon(xx.p2, yyB_ovatus.p2, col='chartreuse4')
yyB_uniformis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]))
yyB_uniformis.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]))
polygon(xx.p1, yyB_uniformis.p1, col='green')
polygon(xx.p2, yyB_uniformis.p2, col='green')
yyB_caccae.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]))
yyB_caccae.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]))
polygon(xx.p1, yyB_caccae.p1, col='limegreen')
polygon(xx.p2, yyB_caccae.p2, col='limegreen')
yyB_intestihominis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]))
yyB_intestihominis.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]))
polygon(xx.p1, yyB_intestihominis.p1, col='seagreen')
polygon(xx.p2, yyB_intestihominis.p2, col='seagreen')
yyE_rectale.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]))
yyE_rectale.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]))
polygon(xx.p1, yyE_rectale.p1, col='midnightblue')
polygon(xx.p2, yyE_rectale.p2, col='midnightblue')
yyC_symbiosum.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]))
yyC_symbiosum.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]))
polygon(xx.p1, yyC_symbiosum.p1, col='dodgerblue4')
polygon(xx.p2, yyC_symbiosum.p2, col='dodgerblue4')
yyB_formatexigens.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]))
yyB_formatexigens.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]))
polygon(xx.p1, yyB_formatexigens.p1, col='blue')
polygon(xx.p2, yyB_formatexigens.p2, col='blue')
yyR_intestinalis.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]))
yyR_intestinalis.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]))
polygon(xx.p1, yyR_intestinalis.p1, col='deepskyblue')
polygon(xx.p2, yyR_intestinalis.p2, col='deepskyblue')
yyF_prausnitzii.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]))
yyF_prausnitzii.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]))
polygon(xx.p1, yyF_prausnitzii.p1, col='slateblue1')
polygon(xx.p2, yyF_prausnitzii.p2, col='slateblue1')
yyE_coli.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]))
yyE_coli.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38] + mean.FR$F_prausnitzii[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]) + rev( mean.FR$E_coli[16:38]))
polygon(xx.p1, yyE_coli.p1, col='gold')
polygon(xx.p2, yyE_coli.p2, col='gold')
yyD_piger.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]))
yyD_piger.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38] + mean.FR$F_prausnitzii[16:38] + mean.FR$E_coli[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]) + rev( mean.FR$E_coli[16:38]) + rev(mean.FR$D_piger[16:38]))
polygon(xx.p1, yyD_piger.p1, col='orange1')
polygon(xx.p2, yyD_piger.p2, col='orange1')
yyC_aerofaciens.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15] + mean.FR$D_piger[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]) + rev(mean.FR$C_aerofaciens[1:15]))
yyC_aerofaciens.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38] + mean.FR$F_prausnitzii[16:38] + mean.FR$E_coli[16:38] + mean.FR$D_piger[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]) + rev( mean.FR$E_coli[16:38]) + rev(mean.FR$D_piger[16:38]) + rev(mean.FR$C_aerofaciens[16:38]))
polygon(xx.p1, yyC_aerofaciens.p1, col='purple1')
polygon(xx.p2, yyC_aerofaciens.p2, col='purple1')
yyA_muciniphila.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15] + mean.FR$D_piger[1:15] + mean.FR$C_aerofaciens[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]) + rev(mean.FR$C_aerofaciens[1:15]) + rev(mean.FR$A_muciniphila[1:15]))
yyA_muciniphila.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38] + mean.FR$F_prausnitzii[16:38] + mean.FR$E_coli[16:38] + mean.FR$D_piger[16:38] + mean.FR$C_aerofaciens[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]) + rev( mean.FR$E_coli[16:38]) + rev(mean.FR$D_piger[16:38]) + rev(mean.FR$C_aerofaciens[16:38]) + rev(mean.FR$A_muciniphila[16:38]))
polygon(xx.p1, yyA_muciniphila.p1, col='red')
polygon(xx.p2, yyA_muciniphila.p2, col='red')
yyC_rodentium.p1 <- c(mean.FR$B_theta[1:15] + mean.FR$B_ovatus[1:15] + mean.FR$B_uniformis[1:15] + mean.FR$B_caccae[1:15] +mean.FR$B_intestihominis[1:15] + mean.FR$E_rectale[1:15] + mean.FR$C_symbiosum[1:15] +mean.FR$B_formatexigens[1:15] + mean.FR$R_intestinalis[1:15] + mean.FR$F_prausnitzii[1:15] + mean.FR$E_coli[1:15] + mean.FR$D_piger[1:15] + mean.FR$C_aerofaciens[1:15] + mean.FR$A_muciniphila[1:15], rev(mean.FR$B_theta[1:15]) + rev(mean.FR$B_ovatus[1:15]) + rev(mean.FR$B_uniformis[1:15]) + rev(mean.FR$B_caccae[1:15]) + rev(mean.FR$B_intestihominis[1:15]) + rev(mean.FR$E_rectale[1:15]) + rev(mean.FR$C_symbiosum[1:15]) + rev(mean.FR$B_formatexigens[1:15]) + rev(mean.FR$R_intestinalis[1:15]) + rev(mean.FR$F_prausnitzii[1:15]) + rev( mean.FR$E_coli[1:15]) + rev(mean.FR$D_piger[1:15]) + rev(mean.FR$C_aerofaciens[1:15]) + rev(mean.FR$A_muciniphila[1:15]) + rev(mean.FR$C_rodentium[1:15]))
yyC_rodentium.p2 <- c(mean.FR$B_theta[16:38] + mean.FR$B_ovatus[16:38] + mean.FR$B_uniformis[16:38] + mean.FR$B_caccae[16:38] +mean.FR$B_intestihominis[16:38] + mean.FR$E_rectale[16:38] + mean.FR$C_symbiosum[16:38] +mean.FR$B_formatexigens[16:38] + mean.FR$R_intestinalis[16:38] + mean.FR$F_prausnitzii[16:38] + mean.FR$E_coli[16:38] + mean.FR$D_piger[16:38] + mean.FR$C_aerofaciens[16:38] + mean.FR$A_muciniphila[16:38], rev(mean.FR$B_theta[16:38]) + rev(mean.FR$B_ovatus[16:38]) + rev(mean.FR$B_uniformis[16:38]) + rev(mean.FR$B_caccae[16:38]) + rev(mean.FR$B_intestihominis[16:38]) + rev(mean.FR$E_rectale[16:38]) + rev(mean.FR$C_symbiosum[16:38]) + rev(mean.FR$B_formatexigens[16:38]) + rev(mean.FR$R_intestinalis[16:38]) + rev(mean.FR$F_prausnitzii[16:38]) + rev( mean.FR$E_coli[16:38]) + rev(mean.FR$D_piger[16:38]) + rev(mean.FR$C_aerofaciens[16:38]) + rev(mean.FR$A_muciniphila[16:38])+ rev(mean.FR$C_rodentium[16:38]))
polygon(xx.p1, yyC_rodentium.p1, col='pink')
polygon(xx.p2, yyC_rodentium.p2, col='pink')

axis(1, at=mean.FR$day, labels=mean.FR$day)
#segments(6, -2, 54, -2, lwd=4, col="black")

# FF:
mean.FF<-mean.FF[order(as.numeric_version(mean.FF$day)),]
mean.FF$day <- as.numeric(as.character(mean.FF$day))
yyB_theta.p1 <- c(rep(0, nrow(mean.FF[1:15,])), rev(mean.FF$B_theta[1:15]))
yyB_theta.p2 <- c(rep(0, nrow(mean.FF[16:38,])), rev(mean.FF$B_theta[16:38]))
plot(x=mean.FF$day, y=mean.FF$B_theta, ylim=ylim, col='darkgreen', type='l', xaxt='n', yaxt='n',
ylab=NA, xlab=NA, main='FF group')
polygon(xx.p1, yyB_theta.p1, col='darkgreen')
polygon(xx.p2, yyB_theta.p2, col='darkgreen')
	# then, repeat the layering of the phylotypes individually:
yyB_ovatus.p1 <- c(mean.FF$B_theta[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]))
yyB_ovatus.p2 <- c(mean.FF$B_theta[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]))
polygon(xx.p1, yyB_ovatus.p1, col='chartreuse4')
polygon(xx.p2, yyB_ovatus.p2, col='chartreuse4')
yyB_uniformis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]))
yyB_uniformis.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]))
polygon(xx.p1, yyB_uniformis.p1, col='green')
polygon(xx.p2, yyB_uniformis.p2, col='green')
yyB_caccae.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]))
yyB_caccae.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]))
polygon(xx.p1, yyB_caccae.p1, col='limegreen')
polygon(xx.p2, yyB_caccae.p2, col='limegreen')
yyB_intestihominis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]))
yyB_intestihominis.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]))
polygon(xx.p1, yyB_intestihominis.p1, col='seagreen')
polygon(xx.p2, yyB_intestihominis.p2, col='seagreen')
yyE_rectale.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]))
yyE_rectale.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]))
polygon(xx.p1, yyE_rectale.p1, col='midnightblue')
polygon(xx.p2, yyE_rectale.p2, col='midnightblue')
yyC_symbiosum.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]))
yyC_symbiosum.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]))
polygon(xx.p1, yyC_symbiosum.p1, col='dodgerblue4')
polygon(xx.p2, yyC_symbiosum.p2, col='dodgerblue4')
yyB_formatexigens.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]))
yyB_formatexigens.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]))
polygon(xx.p1, yyB_formatexigens.p1, col='blue')
polygon(xx.p2, yyB_formatexigens.p2, col='blue')
yyR_intestinalis.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]))
yyR_intestinalis.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]))
polygon(xx.p1, yyR_intestinalis.p1, col='deepskyblue')
polygon(xx.p2, yyR_intestinalis.p2, col='deepskyblue')
yyF_prausnitzii.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]))
yyF_prausnitzii.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]))
polygon(xx.p1, yyF_prausnitzii.p1, col='slateblue1')
polygon(xx.p2, yyF_prausnitzii.p2, col='slateblue1')
yyE_coli.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]))
yyE_coli.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38] + mean.FF$F_prausnitzii[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]) + rev( mean.FF$E_coli[16:38]))
polygon(xx.p1, yyE_coli.p1, col='gold')
polygon(xx.p2, yyE_coli.p2, col='gold')
yyD_piger.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]))
yyD_piger.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38] + mean.FF$F_prausnitzii[16:38] + mean.FF$E_coli[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]) + rev( mean.FF$E_coli[16:38]) + rev(mean.FF$D_piger[16:38]))
polygon(xx.p1, yyD_piger.p1, col='orange1')
polygon(xx.p2, yyD_piger.p2, col='orange1')
yyC_aerofaciens.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15] + mean.FF$D_piger[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]) + rev(mean.FF$C_aerofaciens[1:15]))
yyC_aerofaciens.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38] + mean.FF$F_prausnitzii[16:38] + mean.FF$E_coli[16:38] + mean.FF$D_piger[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]) + rev( mean.FF$E_coli[16:38]) + rev(mean.FF$D_piger[16:38]) + rev(mean.FF$C_aerofaciens[16:38]))
polygon(xx.p1, yyC_aerofaciens.p1, col='purple1')
polygon(xx.p2, yyC_aerofaciens.p2, col='purple1')
yyA_muciniphila.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15] + mean.FF$D_piger[1:15] + mean.FF$C_aerofaciens[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]) + rev(mean.FF$C_aerofaciens[1:15]) + rev(mean.FF$A_muciniphila[1:15]))
yyA_muciniphila.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38] + mean.FF$F_prausnitzii[16:38] + mean.FF$E_coli[16:38] + mean.FF$D_piger[16:38] + mean.FF$C_aerofaciens[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]) + rev( mean.FF$E_coli[16:38]) + rev(mean.FF$D_piger[16:38]) + rev(mean.FF$C_aerofaciens[16:38]) + rev(mean.FF$A_muciniphila[16:38]))
polygon(xx.p1, yyA_muciniphila.p1, col='red')
polygon(xx.p2, yyA_muciniphila.p2, col='red')
yyC_rodentium.p1 <- c(mean.FF$B_theta[1:15] + mean.FF$B_ovatus[1:15] + mean.FF$B_uniformis[1:15] + mean.FF$B_caccae[1:15] +mean.FF$B_intestihominis[1:15] + mean.FF$E_rectale[1:15] + mean.FF$C_symbiosum[1:15] +mean.FF$B_formatexigens[1:15] + mean.FF$R_intestinalis[1:15] + mean.FF$F_prausnitzii[1:15] + mean.FF$E_coli[1:15] + mean.FF$D_piger[1:15] + mean.FF$C_aerofaciens[1:15] + mean.FF$A_muciniphila[1:15], rev(mean.FF$B_theta[1:15]) + rev(mean.FF$B_ovatus[1:15]) + rev(mean.FF$B_uniformis[1:15]) + rev(mean.FF$B_caccae[1:15]) + rev(mean.FF$B_intestihominis[1:15]) + rev(mean.FF$E_rectale[1:15]) + rev(mean.FF$C_symbiosum[1:15]) + rev(mean.FF$B_formatexigens[1:15]) + rev(mean.FF$R_intestinalis[1:15]) + rev(mean.FF$F_prausnitzii[1:15]) + rev( mean.FF$E_coli[1:15]) + rev(mean.FF$D_piger[1:15]) + rev(mean.FF$C_aerofaciens[1:15]) + rev(mean.FF$A_muciniphila[1:15]) + rev(mean.FF$C_rodentium[1:15]))
yyC_rodentium.p2 <- c(mean.FF$B_theta[16:38] + mean.FF$B_ovatus[16:38] + mean.FF$B_uniformis[16:38] + mean.FF$B_caccae[16:38] +mean.FF$B_intestihominis[16:38] + mean.FF$E_rectale[16:38] + mean.FF$C_symbiosum[16:38] +mean.FF$B_formatexigens[16:38] + mean.FF$R_intestinalis[16:38] + mean.FF$F_prausnitzii[16:38] + mean.FF$E_coli[16:38] + mean.FF$D_piger[16:38] + mean.FF$C_aerofaciens[16:38] + mean.FF$A_muciniphila[16:38], rev(mean.FF$B_theta[16:38]) + rev(mean.FF$B_ovatus[16:38]) + rev(mean.FF$B_uniformis[16:38]) + rev(mean.FF$B_caccae[16:38]) + rev(mean.FF$B_intestihominis[16:38]) + rev(mean.FF$E_rectale[16:38]) + rev(mean.FF$C_symbiosum[16:38]) + rev(mean.FF$B_formatexigens[16:38]) + rev(mean.FF$R_intestinalis[16:38]) + rev(mean.FF$F_prausnitzii[16:38]) + rev( mean.FF$E_coli[16:38]) + rev(mean.FF$D_piger[16:38]) + rev(mean.FF$C_aerofaciens[16:38]) + rev(mean.FF$A_muciniphila[16:38])+ rev(mean.FF$C_rodentium[16:38]))
polygon(xx.p1, yyC_rodentium.p1, col='pink')
polygon(xx.p2, yyC_rodentium.p2, col='pink')

axis(1, at=mean.FF$day, labels=mean.FF$day)
#segments(6, -2, 15, -2, lwd=4, col="black")
#segments(15, -2, 54, -2, lwd=4, col="darkred")

plot(1, type="n", axes=FALSE, xlab="", ylab="")			#blank spot for legend...
legend("topleft", c("B. theta", "B. ovatus", "B. uniformis", "B. caccae", "B. intestihominis", "E. rectale", "C. symbiosum", "B. formatexigens", "R. intestinalis", "F. prausnitzii", "E. coli", "D. piger", "C. aerofaciens", "A. muciniphila", "C. rodentium"), fill=c('darkgreen', 'chartreuse4', 'green', 'limegreen', 'seagreen', 'midnightblue', 'dodgerblue4', 'blue', 'deepskyblue', 'slateblue1', 'gold', 'orange1', 'purple1', 'red', 'pink'), cex=1)
#legend("bottomleft", c("Fiber-rich", "Fiber-free", "Prebiotic"), col=c("black", "darkred", "purple4"), lty=1, lwd=5)
mtext(side=2, c("Relative abundance"), outer=TRUE, line=2, cex=0.8, adj=0.68)
mtext(side=1, c("day"), outer=TRUE, line=-1, cex=0.8, adj=0.4)


###
#-----------
###

# Part V:  Bargraphs of avg. relative abundance in cecal samples (part of Fig 2A)

# read in unfiltered data:

phylo<-read.table(file="t.m3.cecal_phylofrac_filtered.txt", header=TRUE)

# get mean for the groups:
mean.cecal<-ddply(phylo, c("group_ID"), colwise(mean, is.numeric))

m.split<-split(mean.cecal, mean.cecal$group_ID)
mean.FF<-m.split$'FF'
mean.FR<-m.split$'FR'
mean.pre<-m.split$'pre'
means.all3<-rbind(mean.FF, mean.FR, mean.pre)

col.cecal<- c("midnightblue","blue","dodgerblue4","deepskyblue","slateblue1","green4","chartreuse3","lawngreen","darkolivegreen2","palegreen3","gold","orange1","purple1","red","lightskyblue","pink","darkgreen","white","magenta","magenta","black","purple","darkblue","darkblue")

par(mar=c(5,4,4,10))
par(xpd=T)
cecal.matrix<-t(means.all3[,5:18])
colnames(cecal.matrix)<-means.all3$group_ID
barplot(cecal.matrix, las=2, ylab="Relative abundance-genera (%)", cex.names=0.5, ylim=c(0,100), col=col.cecal)
legend(9,100,legend=rownames(cecal.matrix),col=col.cecal,fill=col.cecal,cex=0.6)

# note: colors were changed for final manuscript
