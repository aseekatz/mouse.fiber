## Code for relative abundance over time (median, change in relative abundance, and heatmaps of PoMA): Fig. 2C, S3
## Anna M. Seekatz
## 4.22.16

# this is the complete code to recreate community abundance streamplots
# Overall description:
	# generate relative abundance over time for specified phylotypes (linegraphs, Fig. S3)
	# generate CHANGE in relative abundance (from previous day) for specified phylotypes (Fig. 2C)
	# generate Percent of Maximum Abundance per phylotype (heatmaps, Fig. S3)
					
## files used:
	# t.m3_phylofrac_filtered.txt (created in earlier Streamplot code)
	# relabund.summary.stats.txt (summary stats of the relative abundance, created in this code)
	# delta.stats_allgroups.txt (change in relative abundance per phylotype over time, created in this code)
	# delta.sumstats_allgroups.txt (summary stats of the change in relative abundance per group, created in this code)
	# phylofrac_prop.by.max.txt (poma calculations for all groups, created in this code)
	# poma_summary.stats.txt (summary stats of poma calculations, created in this code)
	
	
###
#-----------
###

# part I: generate relative abundance abundance over time -> linegraphs in Fig. S3

## first, let's create a file of the summary stats per phylotype:
phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)

# non-loop format (works, but not really high throughput)
sum.stats <- ddply(phylofrac, c("group_ID", "day"), summarise,
               N    = length(B_theta),
               avg.B_theta = mean(B_theta, na.rm=TRUE),
               sd.B_theta   = sd(B_theta, na.rm=TRUE),
               med.B_theta   = median(B_theta, na.rm=TRUE),
               lq.B_theta	 = quantile(B_theta, 0.25, na.rm=TRUE),
			   hq.B_theta	 = quantile(B_theta, 0.75, na.rm=TRUE),
			   avg.B_intestihominis = mean(B_intestihominis, na.rm=TRUE),
               sd.B_intestihominis   = sd(B_intestihominis, na.rm=TRUE),
               med.B_intestihominis   = median(B_intestihominis, na.rm=TRUE),
               lq.B_intestihominis	 = quantile(B_intestihominis, 0.25, na.rm=TRUE),
			   hq.B_intestihominis	 = quantile(B_intestihominis, 0.75, na.rm=TRUE),
			   avg.R_intestinalis = mean(R_intestinalis, na.rm=TRUE),
               sd.R_intestinalis   = sd(R_intestinalis, na.rm=TRUE),
               med.R_intestinalis   = median(R_intestinalis, na.rm=TRUE),
               lq.R_intestinalis	 = quantile(R_intestinalis, 0.25, na.rm=TRUE),
			   hq.R_intestinalis	 = quantile(R_intestinalis, 0.75, na.rm=TRUE),
			   avg.C_aerofaciens = mean(C_aerofaciens, na.rm=TRUE),
               sd.C_aerofaciens   = sd(C_aerofaciens, na.rm=TRUE),
               med.C_aerofaciens   = median(C_aerofaciens, na.rm=TRUE),
               lq.C_aerofaciens	 = quantile(C_aerofaciens, 0.25, na.rm=TRUE),
			   hq.C_aerofaciens	 = quantile(C_aerofaciens, 0.75, na.rm=TRUE),
			   avg.B_caccae = mean(B_caccae, na.rm=TRUE),
               sd.B_caccae   = sd(B_caccae, na.rm=TRUE),
               med.B_caccae   = median(B_caccae, na.rm=TRUE),
               lq.B_caccae	 = quantile(B_caccae, 0.25, na.rm=TRUE),
			   hq.B_caccae	 = quantile(B_caccae, 0.75, na.rm=TRUE),
			   avg.B_formatexigens = mean(B_formatexigens, na.rm=TRUE),
               sd.B_formatexigens   = sd(B_formatexigens, na.rm=TRUE),
               med.B_formatexigens   = median(B_formatexigens, na.rm=TRUE),
               lq.B_formatexigens	 = quantile(B_formatexigens, 0.25, na.rm=TRUE),
			   hq.B_formatexigens	 = quantile(B_formatexigens, 0.75, na.rm=TRUE),
			   avg.D_piger = mean(D_piger, na.rm=TRUE),
               sd.D_piger   = sd(D_piger, na.rm=TRUE),
               med.D_piger   = median(D_piger, na.rm=TRUE),
               lq.D_piger	 = quantile(D_piger, 0.25, na.rm=TRUE),
			   hq.D_piger	 = quantile(D_piger, 0.75, na.rm=TRUE),
			   avg.B_ovatus = mean(B_ovatus, na.rm=TRUE),
               sd.B_ovatus   = sd(B_ovatus, na.rm=TRUE),
               med.B_ovatus   = median(B_ovatus, na.rm=TRUE),
               lq.B_ovatus	 = quantile(B_ovatus, 0.25, na.rm=TRUE),
			   hq.B_ovatus	 = quantile(B_ovatus, 0.75, na.rm=TRUE),
			   avg.F_prausnitzii = mean(F_prausnitzii, na.rm=TRUE),
               sd.F_prausnitzii   = sd(F_prausnitzii, na.rm=TRUE),
               med.F_prausnitzii   = median(F_prausnitzii, na.rm=TRUE),
               lq.F_prausnitzii	 = quantile(F_prausnitzii, 0.25, na.rm=TRUE),
			   hq.F_prausnitzii	 = quantile(F_prausnitzii, 0.75, na.rm=TRUE),
			   avg.A_muciniphila = mean(A_muciniphila, na.rm=TRUE),
               sd.A_muciniphila   = sd(A_muciniphila, na.rm=TRUE),
               med.A_muciniphila   = median(A_muciniphila, na.rm=TRUE),
               lq.A_muciniphila	 = quantile(A_muciniphila, 0.25, na.rm=TRUE),
			   hq.A_muciniphila	 = quantile(A_muciniphila, 0.75, na.rm=TRUE),
			   avg.B_uniformis = mean(B_uniformis, na.rm=TRUE),
               sd.B_uniformis   = sd(B_uniformis, na.rm=TRUE),
               med.B_uniformis   = median(B_uniformis, na.rm=TRUE),
               lq.B_uniformis	 = quantile(B_uniformis, 0.25, na.rm=TRUE),
			   hq.B_uniformis	 = quantile(B_uniformis, 0.75, na.rm=TRUE),
			   avg.C_symbiosum = mean(C_symbiosum, na.rm=TRUE),
               sd.C_symbiosum   = sd(C_symbiosum, na.rm=TRUE),
               med.C_symbiosum   = median(C_symbiosum, na.rm=TRUE),
               lq.C_symbiosum	 = quantile(C_symbiosum, 0.25, na.rm=TRUE),
			   hq.C_symbiosum	 = quantile(C_symbiosum, 0.75, na.rm=TRUE),
			   avg.E_coli = mean(E_coli, na.rm=TRUE),
               sd.E_coli   = sd(E_coli, na.rm=TRUE),
               med.E_coli   = median(E_coli, na.rm=TRUE),
               lq.E_coli	 = quantile(E_coli, 0.25, na.rm=TRUE),
			   hq.E_coli	 = quantile(E_coli, 0.75, na.rm=TRUE),
			   avg.E_rectalee = mean(E_rectalee, na.rm=TRUE),
               sd.E_rectalee   = sd(E_rectalee, na.rm=TRUE),
               med.E_rectalee   = median(E_rectalee, na.rm=TRUE),
               lq.E_rectalee	 = quantile(E_rectalee, 0.25, na.rm=TRUE),
			   hq.E_rectalee	 = quantile(E_rectalee, 0.75, na.rm=TRUE)
			   )
#write.table(sum.stats, file="relabund.summary.stats.txt", quote=FALSE, sep="\t", col.names=NA)

sum.stats<-read.table(file="relabund.summary.stats.txt", header=TRUE)

# then, define this script's stats to be graphed:
m.sum.split<-split(sum.stats, sum.stats$group_ID)
sum.FF<-m.sum.split$'FF'
sum.FR<-m.sum.split$'FR'
sum.pre<-m.sum.split$'pre'
sum.1_FR_FF<-m.sum.split$'1_FR_FF'
sum.4_FR_FF<-m.sum.split$'4_FR_FF'
sum.4_pre_FF<-m.sum.split$'4_pre_FF'
sum.1_pre_FF<-m.sum.split$'1_pre_FF'

# each of these must be ordered by day, now:
sum.FF<-sum.FF[order(as.numeric_version(sum.FF$day)),]
sum.FF$day <- as.numeric(as.character(sum.FF$day))
sum.FR<-sum.FR[order(as.numeric_version(sum.FR$day)),]
sum.FR$day <- as.numeric(as.character(sum.FR$day))
sum.pre<-sum.pre[order(as.numeric_version(sum.pre$day)),]
sum.pre$day <- as.numeric(as.character(sum.pre$day))
sum.1_FR_FF<-sum.1_FR_FF[order(as.numeric_version(sum.1_FR_FF$day)),]
sum.1_FR_FF$day <- as.numeric(as.character(sum.1_FR_FF$day))
sum.4_FR_FF<-sum.4_FR_FF[order(as.numeric_version(sum.4_FR_FF$day)),]
sum.4_FR_FF$day <- as.numeric(as.character(sum.4_FR_FF$day))
sum.4_pre_FF<-sum.4_pre_FF[order(as.numeric_version(sum.4_pre_FF$day)),]
sum.4_pre_FF$day <- as.numeric(as.character(sum.4_pre_FF$day))
sum.1_pre_FF<-sum.1_pre_FF[order(as.numeric_version(sum.1_pre_FF$day)),]
sum.1_pre_FF$day <- as.numeric(as.character(sum.1_pre_FF$day))

# now graph it (Fig. S3):
par(xpd=TRUE)
layout(matrix(1:4, ncol=4), widths=1, respect=FALSE)
par(mfrow=c(2,4))
par(mai = c(0.3, 0.2, 0.5, 0.1), oma=c(1,4,0,4))

#A_muciniphila:
Am<-plot(sum.FR$day, sum.FR$med.A_muciniphila, 
	pch=".", xlim=c(5,55), ylim=c(10,50), xaxt='n', col="chartreuse4", cex=0.6, 
	main="A. muciniphila", ylab="", xlab="")
axis(1, at=sum.FF$day, labels=sum.FF$day, cex.axis=0.8)
lines(sum.FR$day, sum.FR$med.A_muciniphila, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.A_muciniphila, sum.FR$day, sum.FR$hq.A_muciniphila, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.A_muciniphila, col="red")
arrows(sum.FF$day, sum.FF$lq.A_muciniphila, sum.FF$day, sum.FF$hq.A_muciniphila, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.A_muciniphila, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.A_muciniphila, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.A_muciniphila, sum.pre$day, sum.pre$hq.A_muciniphila, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.A_muciniphila, pch=".", col="grey67", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.A_muciniphila, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.A_muciniphila, sum.1_FR_FF$day, sum.1_FR_FF$hq.A_muciniphila, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.A_muciniphila, pch=".", col="blue3", cex=0.6)

#B_caccae:
Bc<-plot(sum.FR$day, sum.FR$med.B_caccae, 
	pch=".", xlim=c(5,55), ylim=c(0,40), xaxt='n', col="chartreuse4", cex=0.6, 
	main="B. caccae", ylab="", xlab="")
axis(1, at=sum.FR$day, labels=sum.FR$day)
lines(sum.FR$day, sum.FR$med.B_caccae, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.B_caccae, sum.FR$day, sum.FR$hq.B_caccae, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.B_caccae, col="red")
arrows(sum.FF$day, sum.FF$lq.B_caccae, sum.FF$day, sum.FF$hq.B_caccae, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.B_caccae, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.B_caccae, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.B_caccae, sum.pre$day, sum.pre$hq.B_caccae, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.B_caccae, pch=".", col="grey67", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.B_caccae, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.B_caccae, sum.1_FR_FF$day, sum.1_FR_FF$hq.B_caccae, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.B_caccae, pch=".", col="blue3", cex=0.6)

#D_piger:
Dp<-plot(sum.FR$day, sum.FR$med.D_piger, 
	pch=".", xlim=c(5,55), ylim=c(0,10), xaxt='n', col="chartreuse4", cex=0.6, 
	main="D. piger", ylab="", xlab="")
axis(1, at=sum.FR$day, labels=sum.FR$day)
lines(sum.FR$day, sum.FR$med.D_piger, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.D_piger, sum.FR$day, sum.FR$hq.D_piger, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.D_piger, col="red")
arrows(sum.FF$day, sum.FF$lq.D_piger, sum.FF$day, sum.FF$hq.D_piger, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.D_piger, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.D_piger, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.D_piger, sum.pre$day, sum.pre$hq.D_piger, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.D_piger, pch=".", col="grey67", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.D_piger, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.D_piger, sum.1_FR_FF$day, sum.1_FR_FF$hq.D_piger, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.D_piger, pch=".", col="blue3", cex=0.6)

# this '4th plot is a ghost plot to keep the legend in it:
plot(1, type="n", axes=FALSE, xlab="", ylab="")			#blank spot for legend...
legend("topleft", c('Fiber-rich (FR)', 'Fiber-free (FF)', 'prebiotic', '1-day FR/FF'), col=c('chartreuse4', 'red', 'grey67', 'blue4'), cex=0.8, lty=1, lwd=2)

# E. rectale
Er<-plot(sum.FR$day, sum.FR$med.E_rectalee, 
	pch=".", xlim=c(5,55), ylim=c(0,20), xaxt='n', col="chartreuse4", cex=0.6, 
	main="E. rectale", ylab="", xlab="")
axis(1, at=sum.FR$day, labels=sum.FR$day)
lines(sum.FR$day, sum.FR$med.E_rectalee, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.E_rectalee, sum.FR$day, sum.FR$hq.E_rectalee, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.E_rectalee, col="red")
arrows(sum.FF$day, sum.FF$lq.E_rectalee, sum.FF$day, sum.FF$hq.E_rectalee, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.E_rectalee, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.E_rectalee, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.E_rectalee, sum.pre$day, sum.pre$hq.E_rectalee, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.E_rectalee, pch=".", col="grey67", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.E_rectalee, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.E_rectalee, sum.1_FR_FF$day, sum.1_FR_FF$hq.E_rectalee, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.E_rectalee, pch=".", col="blue3", cex=0.6)

#B_ovatus:
Bo<-plot(sum.FR$day, sum.FR$med.B_ovatus, 
	pch=".", xlim=c(5,55), ylim=c(0,30), xaxt='n', col="chartreuse4", cex=0.6, 
	main="B. ovatus", ylab="", xlab="")
axis(1, at=sum.FR$day, labels=sum.FR$day)
lines(sum.FR$day, sum.FR$med.B_ovatus, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.B_ovatus, sum.FR$day, sum.FR$hq.B_ovatus, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.B_ovatus, col="red")
arrows(sum.FF$day, sum.FF$lq.B_ovatus, sum.FF$day, sum.FF$hq.B_ovatus, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.B_ovatus, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.B_ovatus, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.B_ovatus, sum.pre$day, sum.pre$hq.B_ovatus, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.B_ovatus, pch=".", col="grey67", cex=0.6)#E_rectalee:
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.B_ovatus, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.B_ovatus, sum.1_FR_FF$day, sum.1_FR_FF$hq.B_ovatus, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.B_ovatus, pch=".", col="blue3", cex=0.6)

#M_formatexigens:
Mf<-plot(sum.FR$day, sum.FR$med.B_formatexigens, 
	pch=".", xlim=c(5,55), ylim=c(0,15), xaxt='n', col="chartreuse4", cex=0.6, 
	main="M. formatexigens", ylab="", xlab="")
axis(1, at=sum.FR$day, labels=sum.FR$day)
lines(sum.FR$day, sum.FR$med.B_formatexigens, col="chartreuse4")
arrows(sum.FR$day, sum.FR$lq.B_formatexigens, sum.FR$day, sum.FR$hq.B_formatexigens, length=0.025, angle=90, code=3, col="chartreuse4")
lines(sum.FF$day, sum.FF$med.B_formatexigens, col="red")
arrows(sum.FF$day, sum.FF$lq.B_formatexigens, sum.FF$day, sum.FF$hq.B_formatexigens, length=0.025, angle=90, code=3, col="red")
points(sum.FF$day, sum.FF$med.B_formatexigens, pch=".", col="red", cex=0.6)
lines(sum.pre$day, sum.pre$med.B_formatexigens, col="grey67", lty=1)
arrows(sum.pre$day, sum.pre$lq.B_formatexigens, sum.pre$day, sum.pre$hq.B_formatexigens, length=0.025, angle=90, code=3, col="grey67")
points(sum.pre$day, sum.pre$med.B_formatexigens, pch=".", col="grey67", cex=0.6)#E_rectalee:
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.B_formatexigens, col="blue3", lty=1)
arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.B_formatexigens, sum.1_FR_FF$day, sum.1_FR_FF$hq.B_formatexigens, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.B_formatexigens, pch=".", col="blue3", cex=0.6)

mtext(side=2, c("Relative abundance (%)"), outer=TRUE, line=2, cex=0.8, adj=0.68)

###
#-----------
###

# part II: generate CHANGE in relative abundance abundance over time -> Fig. 2C
# To calculate change over time (difference between consecutive samplings...)

phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)

phylofrac<- phylpfrac[order(phylofrac$day, phylofrac$mouse) ,]
	# must order by mouse and day!!

# subset by mouse:
mouse.split<-split(phylofrac, phylofrac$mouse)
# FF mice = 1-4::
m1<-mouse.split$'1'
m2<-mouse.split$'2'
m3<-mouse.split$'3'
m4<-mouse.split$'4'
# FR mice = 5-8:
m5<-mouse.split$'5'
m6<-mouse.split$'6'
m7<-mouse.split$'7'
m8<-mouse.split$'8'
# pre mice = 9-11:
m9<-mouse.split$'9'
m10<-mouse.split$'10'
m11<-mouse.split$'11'
# 1-day FR/FF mice = 12-14:
m12<-mouse.split$'12'
m13<-mouse.split$'13'
m14<-mouse.split$'14'
# 4-day FF/FR mice = 15-17
m15<-mouse.split$'15'
m16<-mouse.split$'16'
m17<-mouse.split$'17'
# 1-day pre/FF = 18-20:
m18<-mouse.split$'18'
m19<-mouse.split$'19'
m20<-mouse.split$'20'
# 4-day pre/FF = 21-23:
m21<-mouse.split$'21'
m22<-mouse.split$'22'
m23<-mouse.split$'23'

# this command takes the difference between consecutive rows:
#diff(m5$B_theta)

# note: time points are not necessarily equal between all the samples/days
	# this might impact naming the rows for some of the ind. mice
	
## For each group, do the following:

## FF group:
# make each of these a matrix to apply the 'diff' command to multiple column variables:
m1.matrix<-m1[, c(10:23)]						#chooses only your relative abundance columns
m2.matrix<-m2[, c(10:23)]
m3.matrix<-m3[, c(10:23)]
m4.matrix<-m4[, c(10:23)]
	#rownames(m4.matrix)<-m4$day
						
# apply diff function and rename rowname as day after change occurred (also, added pertinent metadata)
change.m1<-as.data.frame(lapply(m1.matrix, diff))			#calculates difference between consecutive days--woo hoo!
rownames(change.m1)<-m1$day[2:28]
change.m1$group_ID<-c("FF")
change.m1$mouse<-c(1)
change.m1$day<-rownames(change.m1)

change.m2<-as.data.frame(lapply(m2.matrix, diff))			
rownames(change.m2)<-m2$day[2:28]
change.m2$group_ID<-c("FF")
change.m2$mouse<-c(2)
change.m2$day<-rownames(change.m2)

change.m3<-as.data.frame(lapply(m3.matrix, diff))			
rownames(change.m3)<-m3$day[2:28]
change.m3$group_ID<-c("FF")
change.m3$mouse<-c(3)
change.m3$day<-rownames(change.m3)

change.m4<-as.data.frame(lapply(m4.matrix, diff))
rownames(change.m4)<-m4$day[2:26]
change.m4$group_ID<-c("FF")
change.m4$mouse<-c(1)
change.m4$day<-rownames(change.m4)

## FR group:
m5.matrix<-m5[, c(10:23)]
m6.matrix<-m6[, c(10:23)]
m7.matrix<-m7[, c(10:23)]
m8.matrix<-m8[, c(10:23)]
						
change.m5<-as.data.frame(lapply(m5.matrix, diff))
rownames(change.m5)<-m5$day[2:28]
change.m5$group_ID<-c("FR")
change.m5$mouse<-c(5)
change.m5$day<-rownames(change.m5)

change.m6<-as.data.frame(lapply(m6.matrix, diff))			
rownames(change.m6)<-m6$day[2:28]
change.m6$group_ID<-c("FR")
change.m6$mouse<-c(6)
change.m6$day<-rownames(change.m6)

change.m7<-as.data.frame(lapply(m7.matrix, diff))			
rownames(change.m7)<-m7$day[2:28]
change.m7$group_ID<-c("FR")
change.m7$mouse<-c(7)
change.m7$day<-rownames(change.m7)

change.m8<-as.data.frame(lapply(m8.matrix, diff))		
rownames(change.m8)<-m8$day[2:28]
change.m8$group_ID<-c("FR")
change.m8$mouse<-c(8)
change.m8$day<-rownames(change.m8)

## pre group:
m9.matrix<-m9[, c(10:23)]
m10.matrix<-m10[, c(10:23)]
m11.matrix<-m11[, c(10:23)]
						
change.m9<-as.data.frame(lapply(m9.matrix, diff))
rownames(change.m9)<-m9$day[2:27]
change.m9$group_ID<-c("pre")
change.m9$mouse<-c(9)
change.m9$day<-rownames(change.m9)

change.m10<-as.data.frame(lapply(m10.matrix, diff))			
rownames(change.m10)<-m10$day[2:28]
change.m10$group_ID<-c("pre")
change.m10$mouse<-c(10)
change.m10$day<-rownames(change.m10)

change.m11<-as.data.frame(lapply(m11.matrix, diff))			
rownames(change.m11)<-m11$day[2:28]
change.m11$group_ID<-c("pre")
change.m11$mouse<-c(11)
change.m11$day<-rownames(change.m11)

## 1-day FR/FF group:
m12.matrix<-m12[, c(10:23)]
m13.matrix<-m13[, c(10:23)]
m14.matrix<-m14[, c(10:23)]
						
change.m12<-as.data.frame(lapply(m12.matrix, diff))
rownames(change.m12)<-m12$day[2:28]
change.m12$group_ID<-c("1_FR_FF")
change.m12$mouse<-c(12)
change.m12$day<-rownames(change.m12)

change.m13<-as.data.frame(lapply(m13.matrix, diff))			
rownames(change.m13)<-m13$day[2:28]
change.m13$group_ID<-c("1_FR_FF")
change.m13$mouse<-c(13)
change.m13$day<-rownames(change.m13)

change.m14<-as.data.frame(lapply(m14.matrix, diff))			
rownames(change.m14)<-m14$day[2:28]
change.m14$group_ID<-c("1_FR_FF")
change.m14$mouse<-c(14)
change.m14$day<-rownames(change.m14)

## 4-day FR/FF group:
# make each of these a matrix to apply the 'diff' command to multiple column variables:
m15.matrix<-m15[, c(10:23)]
m16.matrix<-m16[, c(10:23)]
m17.matrix<-m17[, c(10:23)]
						
change.m15<-as.data.frame(lapply(m15.matrix, diff))
rownames(change.m15)<-m15$day[2:27]
change.m15$group_ID<-c("4_FR_FF")
change.m15$mouse<-c(15)
change.m15$day<-rownames(change.m15)

change.m16<-as.data.frame(lapply(m16.matrix, diff))			
rownames(change.m16)<-m16$day[2:28]
change.m16$group_ID<-c("4_FR_FF")
change.m16$mouse<-c(16)
change.m16$day<-rownames(change.m16)

change.m17<-as.data.frame(lapply(m17.matrix, diff))			
rownames(change.m17)<-m17$day[2:27]
change.m17$group_ID<-c("4_FR_FF")
change.m17$mouse<-c(17)
change.m17$day<-rownames(change.m17)

## 1-day pre/FF group:
m18.matrix<-m18[, c(10:23)]
m19.matrix<-m19[, c(10:23)]
m20.matrix<-m20[, c(10:23)]
						
change.m18<-as.data.frame(lapply(m18.matrix, diff))
rownames(change.m18)<-m18$day[2:28]
change.m18$group_ID<-c("1_pre_FF")
change.m18$mouse<-c(18)
change.m18$day<-rownames(change.m18)

change.m19<-as.data.frame(lapply(m19.matrix, diff))			
rownames(change.m19)<-m19$day[2:28]
change.m19$group_ID<-c("1_pre_FF")
change.m19$mouse<-c(19)
change.m19$day<-rownames(change.m19)

change.m20<-as.data.frame(lapply(m20.matrix, diff))			
rownames(change.m20)<-m20$day[2:27]
change.m20$group_ID<-c("1_pre_FF")
change.m20$mouse<-c(20)
change.m20$day<-rownames(change.m20)

## 4-day pre/FF group:
m21.matrix<-m21[, c(10:23)]
m22.matrix<-m22[, c(10:23)]
m23.matrix<-m23[, c(10:23)]
						
change.m21<-as.data.frame(lapply(m21.matrix, diff))
rownames(change.m21)<-m21$day[2:26]
change.m21$group_ID<-c("4_pre_FF")
change.m21$mouse<-c(21)
change.m21$day<-rownames(change.m21)

change.m22<-as.data.frame(lapply(m22.matrix, diff))			
rownames(change.m22)<-m22$day[2:27]
change.m22$group_ID<-c("4_pre_FF")
change.m22$mouse<-c(22)
change.m22$day<-rownames(change.m22)

change.m23<-as.data.frame(lapply(m23.matrix, diff))			
rownames(change.m23)<-m23$day[2:28]
change.m23$group_ID<-c("4_pre_FF")
change.m23$mouse<-c(23)
change.m23$day<-rownames(change.m23)

# then bind them all together:

# then bind them all together:
change.all<-rbind(change.m1, change.m2, change.m3, change.m4, change.m5, change.m6, change.m7, change.m8, change.m9, change.m10, change.m11, change.m12, change.m13, change.m14, change.m15, change.m16, change.m17, change.m18, change.m19, change.m20, change.m21, change.m22, change.m23)
write.table(change.all, file="delta.stats_allgroups.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	# this file can be used to calculate stats of change over time for one treatment group (d13 compared to d14, etc)

# use sum.stats script to get summary statistics on this:
sum.all <- ddply(change.all, c("group_ID", "day"), summarise,
               N    = length(B_theta),
               avg.B_theta = mean(B_theta, na.rm=TRUE),
               sd.B_theta   = sd(B_theta, na.rm=TRUE),
               med.B_theta   = median(B_theta, na.rm=TRUE),
               lq.B_theta	 = quantile(B_theta, 0.25, na.rm=TRUE),
			   hq.B_theta	 = quantile(B_theta, 0.75, na.rm=TRUE),
			   avg.B_intestihominis = mean(B_intestihominis, na.rm=TRUE),
               sd.B_intestihominis   = sd(B_intestihominis, na.rm=TRUE),
               med.B_intestihominis   = median(B_intestihominis, na.rm=TRUE),
               lq.B_intestihominis	 = quantile(B_intestihominis, 0.25, na.rm=TRUE),
			   hq.B_intestihominis	 = quantile(B_intestihominis, 0.75, na.rm=TRUE),
			   avg.R_intestinalis = mean(R_intestinalis, na.rm=TRUE),
               sd.R_intestinalis   = sd(R_intestinalis, na.rm=TRUE),
               med.R_intestinalis   = median(R_intestinalis, na.rm=TRUE),
               lq.R_intestinalis	 = quantile(R_intestinalis, 0.25, na.rm=TRUE),
			   hq.R_intestinalis	 = quantile(R_intestinalis, 0.75, na.rm=TRUE),
			   avg.C_aerofaciens = mean(C_aerofaciens, na.rm=TRUE),
               sd.C_aerofaciens   = sd(C_aerofaciens, na.rm=TRUE),
               med.C_aerofaciens   = median(C_aerofaciens, na.rm=TRUE),
               lq.C_aerofaciens	 = quantile(C_aerofaciens, 0.25, na.rm=TRUE),
			   hq.C_aerofaciens	 = quantile(C_aerofaciens, 0.75, na.rm=TRUE),
			   avg.B_caccae = mean(B_caccae, na.rm=TRUE),
               sd.B_caccae   = sd(B_caccae, na.rm=TRUE),
               med.B_caccae   = median(B_caccae, na.rm=TRUE),
               lq.B_caccae	 = quantile(B_caccae, 0.25, na.rm=TRUE),
			   hq.B_caccae	 = quantile(B_caccae, 0.75, na.rm=TRUE),
			   avg.B_formatexigens = mean(B_formatexigens, na.rm=TRUE),
               sd.B_formatexigens   = sd(B_formatexigens, na.rm=TRUE),
               med.B_formatexigens   = median(B_formatexigens, na.rm=TRUE),
               lq.B_formatexigens	 = quantile(B_formatexigens, 0.25, na.rm=TRUE),
			   hq.B_formatexigens	 = quantile(B_formatexigens, 0.75, na.rm=TRUE),
			   avg.D_piger = mean(D_piger, na.rm=TRUE),
               sd.D_piger   = sd(D_piger, na.rm=TRUE),
               med.D_piger   = median(D_piger, na.rm=TRUE),
               lq.D_piger	 = quantile(D_piger, 0.25, na.rm=TRUE),
			   hq.D_piger	 = quantile(D_piger, 0.75, na.rm=TRUE),
			   avg.B_ovatus = mean(B_ovatus, na.rm=TRUE),
               sd.B_ovatus   = sd(B_ovatus, na.rm=TRUE),
               med.B_ovatus   = median(B_ovatus, na.rm=TRUE),
               lq.B_ovatus	 = quantile(B_ovatus, 0.25, na.rm=TRUE),
			   hq.B_ovatus	 = quantile(B_ovatus, 0.75, na.rm=TRUE),
			   avg.F_prausnitzii = mean(F_prausnitzii, na.rm=TRUE),
               sd.F_prausnitzii   = sd(F_prausnitzii, na.rm=TRUE),
               med.F_prausnitzii   = median(F_prausnitzii, na.rm=TRUE),
               lq.F_prausnitzii	 = quantile(F_prausnitzii, 0.25, na.rm=TRUE),
			   hq.F_prausnitzii	 = quantile(F_prausnitzii, 0.75, na.rm=TRUE),
			   avg.A_muciniphila = mean(A_muciniphila, na.rm=TRUE),
               sd.A_muciniphila   = sd(A_muciniphila, na.rm=TRUE),
               med.A_muciniphila   = median(A_muciniphila, na.rm=TRUE),
               lq.A_muciniphila	 = quantile(A_muciniphila, 0.25, na.rm=TRUE),
			   hq.A_muciniphila	 = quantile(A_muciniphila, 0.75, na.rm=TRUE),
			   avg.B_uniformis = mean(B_uniformis, na.rm=TRUE),
               sd.B_uniformis   = sd(B_uniformis, na.rm=TRUE),
               med.B_uniformis   = median(B_uniformis, na.rm=TRUE),
               lq.B_uniformis	 = quantile(B_uniformis, 0.25, na.rm=TRUE),
			   hq.B_uniformis	 = quantile(B_uniformis, 0.75, na.rm=TRUE),
			   avg.C_symbiosum = mean(C_symbiosum, na.rm=TRUE),
               sd.C_symbiosum   = sd(C_symbiosum, na.rm=TRUE),
               med.C_symbiosum   = median(C_symbiosum, na.rm=TRUE),
               lq.C_symbiosum	 = quantile(C_symbiosum, 0.25, na.rm=TRUE),
			   hq.C_symbiosum	 = quantile(C_symbiosum, 0.75, na.rm=TRUE),
			   avg.E_coli = mean(E_coli, na.rm=TRUE),
               sd.E_coli   = sd(E_coli, na.rm=TRUE),
               med.E_coli   = median(E_coli, na.rm=TRUE),
               lq.E_coli	 = quantile(E_coli, 0.25, na.rm=TRUE),
			   hq.E_coli	 = quantile(E_coli, 0.75, na.rm=TRUE),
			   avg.E_rectalee = mean(E_rectalee, na.rm=TRUE),
               sd.E_rectalee   = sd(E_rectalee, na.rm=TRUE),
               med.E_rectalee   = median(E_rectalee, na.rm=TRUE),
               lq.E_rectalee	 = quantile(E_rectalee, 0.25, na.rm=TRUE),
			   hq.E_rectalee	 = quantile(E_rectalee, 0.75, na.rm=TRUE)
			   )
sum.all<-sum.all[order(as.numeric_version(sum.all$day)),]
sum.all$day <- as.numeric(as.character(sum.all$day))
write.table(sum.all, file="delta.sumstats_allgroups.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# to plot the change in relative abundance of the following phylotypes for specific diet group, do the following (Fig. 2C):

# read in files and subset as different groups:
deltasum<-read.table(file="delta.sumstats_allgroups.txt", header=TRUE)
sum.1_FR_FF<-deltasum[deltasum$group_ID==c("1_FR_FF"), ]
sum.FR<-deltasum[deltasum$group_ID==c("FR"), ]
sum.FF<-deltasum[deltasum$group_ID==c("FF"), ]

# plotting each of these graphs:

# FF, FR, and 1_FR_FF:

par(xpd=TRUE)
layout(matrix(1:4, ncol=4), widths=1, respect=FALSE)
par(mfrow=c(2,2))
par(mai = c(0.3, 0.2, 0.5, 0.1), oma=c(1,4,0,4))

# FR, FF, 1FR/FF:
# note: can add error bars, if want (removed in figure)
Am<-plot(sum.FR$day, sum.FR$med.A_muciniphila, 
	pch=".", xlim=c(5,55), ylim=c(-15,25), xaxt='n', col="green3", cex=0.6, 
	main="A. muciniphila", ylab="", xlab="")
axis(1, at=sum.FF$day, labels=sum.FF$day, cex.axis=0.8)
lines(sum.FR$day, sum.FR$med.A_muciniphila, col="green3")
#arrows(sum.FR$day, sum.FR$lq.A_muciniphila, sum.FR$day, sum.FR$hq.A_muciniphila, length=0.025, angle=90, code=3, col="green3")
lines(sum.FF$day, sum.FF$med.A_muciniphila, col="red3")
#arrows(sum.FF$day, sum.FF$lq.A_muciniphila, sum.FF$day, sum.FF$hq.A_muciniphila, length=0.025, angle=90, code=3, col="red3")
points(sum.FF$day, sum.FF$med.A_muciniphila, pch=".", col="red3", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.A_muciniphila, col="blue3", lty=1)
#arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.A_muciniphila, sum.1_FR_FF$day, sum.1_FR_FF$hq.A_muciniphila, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.A_muciniphila, pch=".", col="blue3", cex=0.6)

Bc<-plot(sum.FR$day, sum.FR$med.B_caccae, 
	pch=".", xlim=c(5,55), ylim=c(-15,10), xaxt='n', col="green3", cex=0.6, 
	main="B. caccae", ylab="", xlab="")
axis(1, at=sum.FF$day, labels=sum.FF$day, cex.axis=0.8)
lines(sum.FR$day, sum.FR$med.B_caccae, col="green3")
#arrows(sum.FR$day, sum.FR$lq.B_caccae, sum.FR$day, sum.FR$hq.B_caccae, length=0.025, angle=90, code=3, col="green3")
lines(sum.FF$day, sum.FF$med.B_caccae, col="red3")
#arrows(sum.FF$day, sum.FF$lq.B_caccae, sum.FF$day, sum.FF$hq.B_caccae, length=0.025, angle=90, code=3, col="red3")
points(sum.FF$day, sum.FF$med.B_caccae, pch=".", col="red3", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.B_caccae, col="blue3", lty=1)
#arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.B_caccae, sum.1_FR_FF$day, sum.1_FR_FF$hq.B_caccae, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.B_caccae, pch=".", col="blue3", cex=0.6)
legend("bottomright", c('Fiber-rich (FR)', 'Fiber-free (FF)', '1-day FR/FF'), col=c('green3', 'red3', 'blue3'), cex=0.8, lty=1, lwd=2)

Er<-plot(sum.FR$day, sum.FR$med.E_rectalee, 
	pch=".", xlim=c(5,55), ylim=c(-10,15), xaxt='n', col="green3", cex=0.6, 
	main="E. rectale", ylab="", xlab="")
axis(1, at=sum.FF$day, labels=sum.FF$day, cex.axis=0.8)
lines(sum.FR$day, sum.FR$med.E_rectalee, col="green3")
#arrows(sum.FR$day, sum.FR$lq.E_rectalee, sum.FR$day, sum.FR$hq.E_rectalee, length=0.025, angle=90, code=3, col="green3")
lines(sum.FF$day, sum.FF$med.E_rectalee, col="red3")
#arrows(sum.FF$day, sum.FF$lq.E_rectalee, sum.FF$day, sum.FF$hq.E_rectalee, length=0.025, angle=90, code=3, col="red3")
points(sum.FF$day, sum.FF$med.E_rectalee, pch=".", col="red3", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.E_rectalee, col="blue3", lty=1)
#arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.E_rectalee, sum.1_FR_FF$day, sum.1_FR_FF$hq.E_rectalee, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.E_rectalee, pch=".", col="blue3", cex=0.6)

Bo<-plot(sum.FR$day, sum.FR$med.B_ovatus, 
	pch=".", xlim=c(5,55), ylim=c(-15,10), xaxt='n', col="green3", cex=0.6, 
	main="B. ovatus", ylab="", xlab="")
axis(1, at=sum.FF$day, labels=sum.FF$day, cex.axis=0.8)
lines(sum.FR$day, sum.FR$med.B_ovatus, col="green3")
#arrows(sum.FR$day, sum.FR$lq.B_ovatus, sum.FR$day, sum.FR$hq.B_ovatus, length=0.025, angle=90, code=3, col="green3")
lines(sum.FF$day, sum.FF$med.B_ovatus, col="red3")
#arrows(sum.FF$day, sum.FF$lq.B_ovatus, sum.FF$day, sum.FF$hq.B_ovatus, length=0.025, angle=90, code=3, col="red3")
points(sum.FF$day, sum.FF$med.B_ovatus, pch=".", col="red3", cex=0.6)
lines(sum.1_FR_FF$day, sum.1_FR_FF$med.B_ovatus, col="blue3", lty=1)
#arrows(sum.1_FR_FF$day, sum.1_FR_FF$lq.B_ovatus, sum.1_FR_FF$day, sum.1_FR_FF$hq.B_ovatus, length=0.025, angle=90, code=3, col="blue3")
points(sum.1_FR_FF$day, sum.1_FR_FF$med.B_ovatus, pch=".", col="blue3", cex=0.6)


###
#-----------
###

# part III: generate PoMA heatmaps (Percent of Maximum Abundance) -> Fig. S3, heatmaps

# proportions: max proportion
# trying to get proportion of each species' relative abundance at each time point of the max abundance per mouse for that species...

phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)

# first, define the max (or mean) of all time points within each mouse:
bymouse <- ddply(phylofrac, c("mouse"), summarise,
               N    = length(B_theta),
			   mean.B_theta  =mean(B_theta, na.rm=TRUE),
			   max.B_theta  =max(B_theta, na.rm=TRUE),
			   mean.B_intestihominis = mean(B_intestihominis, na.rm=TRUE),
               max.B_intestihominis   = max(B_intestihominis, na.rm=TRUE),
               max.R_intestinalis   = max(R_intestinalis, na.rm=TRUE),
               mean.R_intestinalis   = mean(R_intestinalis, na.rm=TRUE),
               max.C_aerofaciens   = max(C_aerofaciens, na.rm=TRUE),
               mean.C_aerofaciens   = mean(C_aerofaciens, na.rm=TRUE),
               max.B_caccae   = max(B_caccae, na.rm=TRUE),
               mean.B_caccae   = mean(B_caccae, na.rm=TRUE),
               max.B_formatexigens   = max(B_formatexigens, na.rm=TRUE),
               mean.B_formatexigens   = mean(B_formatexigens, na.rm=TRUE),
               max.D_piger   = max(D_piger, na.rm=TRUE),
               mean.D_piger   = mean(D_piger, na.rm=TRUE),
               max.B_ovatus   = max(B_ovatus, na.rm=TRUE),
               mean.B_ovatus   = mean(B_ovatus, na.rm=TRUE),
               max.F_prausnitzii   = max(F_prausnitzii, na.rm=TRUE),
               mean.F_prausnitzii   = mean(F_prausnitzii, na.rm=TRUE),
               max.A_muciniphila   = max(A_muciniphila, na.rm=TRUE),
               mean.A_muciniphila   = mean(A_muciniphila, na.rm=TRUE),
               max.B_uniformis   = max(B_uniformis, na.rm=TRUE),
               mean.B_uniformis   = mean(B_uniformis, na.rm=TRUE),
               max.C_symbiosum   = max(C_symbiosum, na.rm=TRUE),
               mean.C_symbiosum   = mean(C_symbiosum, na.rm=TRUE),
               max.E_coli   = max(E_coli, na.rm=TRUE),
               mean.E_coli   = mean(E_coli, na.rm=TRUE),
               max.E_rectale   = max(E_rectale, na.rm=TRUE),
               mean.E_rectale   = mean(E_rectale, na.rm=TRUE)
               )
bymouse<- bymouse[ order(bymouse$mouse), ]

# add mean/max colums to phylofrac:
phylomax<-merge(phylofrac, bymouse, by.x=c("mouse"), by.y=c("mouse"), all.x=TRUE)
#write.table(phylomax, file="phylofrac_w.max.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


# then, divide each species max abundance by its relative abundance on that day:
phylomax$prop.B_theta<-phylomax$B_theta/phylomax$max.B_theta
phylomax$prop.B_intestihominis<-phylomax$B_intestihominis/phylomax$max.B_intestihominis
phylomax$prop.R_intestinalis<-phylomax$R_intestinalis/phylomax$max.R_intestinalis
phylomax$prop.C_aerofaciens<-phylomax$C_aerofaciens/phylomax$max.C_aerofaciens
phylomax$prop.B_caccae<-phylomax$B_caccae/phylomax$max.B_caccae
phylomax$prop.B_formatexigens<-phylomax$B_formatexigens/phylomax$max.B_formatexigens
phylomax$prop.D_piger<-phylomax$D_piger/phylomax$max.D_piger
phylomax$prop.B_ovatus<-phylomax$B_ovatus/phylomax$max.B_ovatus
phylomax$prop.F_prausnitzii<-phylomax$F_prausnitzii/phylomax$max.F_prausnitzii
phylomax$prop.A_muciniphila<-phylomax$A_muciniphila/phylomax$max.A_muciniphila
phylomax$prop.B_uniformis<-phylomax$B_uniformis/phylomax$max.B_uniformis
phylomax$prop.C_symbiosum<-phylomax$C_symbiosum/phylomax$max.C_symbiosum
phylomax$prop.E_coli<-phylomax$E_coli/phylomax$max.E_coli
phylomax$prop.E_rectale<-phylomax$E_rectale/phylomax$max.E_rectale

# some of the values are 0, so you get an 'inf' value
#replace these:
phylomax[phylomax==Inf] <- 0

# then, select the columns that have been adjusted to create a new data fram of the phylofractions for further analysis:
phylo.prop<-phylomax[,c("sampleID", "mouse", "group_ID", "day", "diet", "day_of_diet", "prop.B_theta", "prop.B_intestihominis", "prop.R_intestinalis", "prop.C_aerofaciens", "prop.B_caccae", "prop.B_formatexigens", "prop.D_piger", "prop.B_ovatus", "prop.F_prausnitzii", "prop.A_muciniphila", "prop.B_uniformis", "prop.C_symbiosum", "prop.E_coli", "prop.E_rectale")] 
#write.table(phylo.prop, file="phylofrac_prop.by.max.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

# in future, can just upload this file to get the new groupings:
phylo.prop<-read.table(file="phylofrac_prop.by.max.txt", header=TRUE)

# test out summary stats and graphing on this:

prop.sum.stats <- ddply(phylo.prop, c("group_ID", "day"), summarise,
               N    = length(prop.B_theta),
               avg.B_theta = mean(prop.B_theta, na.rm=TRUE),
               sd.B_theta   = sd(prop.B_theta, na.rm=TRUE),
               med.B_theta   = median(prop.B_theta, na.rm=TRUE),
               lq.B_theta	 = quantile(prop.B_theta, 0.25, na.rm=TRUE),
			   hq.B_theta	 = quantile(prop.B_theta, 0.75, na.rm=TRUE),
			   avg.B_intestihominis = mean(prop.B_intestihominis, na.rm=TRUE),
               sd.B_intestihominis   = sd(prop.B_intestihominis, na.rm=TRUE),
               med.B_intestihominis   = median(prop.B_intestihominis, na.rm=TRUE),
               lq.B_intestihominis	 = quantile(prop.B_intestihominis, 0.25, na.rm=TRUE),
			   hq.B_intestihominis	 = quantile(prop.B_intestihominis, 0.75, na.rm=TRUE),
			   avg.R_intestinalis = mean(prop.R_intestinalis, na.rm=TRUE),
               sd.R_intestinalis   = sd(prop.R_intestinalis, na.rm=TRUE),
               med.R_intestinalis   = median(prop.R_intestinalis, na.rm=TRUE),
               lq.R_intestinalis	 = quantile(prop.R_intestinalis, 0.25, na.rm=TRUE),
			   hq.R_intestinalis	 = quantile(prop.R_intestinalis, 0.75, na.rm=TRUE),
			   avg.C_aerofaciens = mean(prop.C_aerofaciens, na.rm=TRUE),
               sd.C_aerofaciens   = sd(prop.C_aerofaciens, na.rm=TRUE),
               med.C_aerofaciens   = median(prop.C_aerofaciens, na.rm=TRUE),
               lq.C_aerofaciens	 = quantile(prop.C_aerofaciens, 0.25, na.rm=TRUE),
			   hq.C_aerofaciens	 = quantile(prop.C_aerofaciens, 0.75, na.rm=TRUE),
			   avg.B_caccae = mean(prop.B_caccae, na.rm=TRUE),
               sd.B_caccae   = sd(prop.B_caccae, na.rm=TRUE),
               med.B_caccae   = median(prop.B_caccae, na.rm=TRUE),
               lq.B_caccae	 = quantile(prop.B_caccae, 0.25, na.rm=TRUE),
			   hq.B_caccae	 = quantile(prop.B_caccae, 0.75, na.rm=TRUE),
			   avg.B_formatexigens = mean(prop.B_formatexigens, na.rm=TRUE),
               sd.B_formatexigens   = sd(prop.B_formatexigens, na.rm=TRUE),
               med.B_formatexigens   = median(prop.B_formatexigens, na.rm=TRUE),
               lq.B_formatexigens	 = quantile(prop.B_formatexigens, 0.25, na.rm=TRUE),
			   hq.B_formatexigens	 = quantile(prop.B_formatexigens, 0.75, na.rm=TRUE),
			   avg.D_piger = mean(prop.D_piger, na.rm=TRUE),
               sd.D_piger   = sd(prop.D_piger, na.rm=TRUE),
               med.D_piger   = median(prop.D_piger, na.rm=TRUE),
               lq.D_piger	 = quantile(prop.D_piger, 0.25, na.rm=TRUE),
			   hq.D_piger	 = quantile(prop.D_piger, 0.75, na.rm=TRUE),
			   avg.B_ovatus = mean(prop.B_ovatus, na.rm=TRUE),
               sd.B_ovatus   = sd(prop.B_ovatus, na.rm=TRUE),
               med.B_ovatus   = median(prop.B_ovatus, na.rm=TRUE),
               lq.B_ovatus	 = quantile(prop.B_ovatus, 0.25, na.rm=TRUE),
			   hq.B_ovatus	 = quantile(prop.B_ovatus, 0.75, na.rm=TRUE),
			   avg.F_prausnitzii = mean(prop.F_prausnitzii, na.rm=TRUE),
               sd.F_prausnitzii   = sd(prop.F_prausnitzii, na.rm=TRUE),
               med.F_prausnitzii   = median(prop.F_prausnitzii, na.rm=TRUE),
               lq.F_prausnitzii	 = quantile(prop.F_prausnitzii, 0.25, na.rm=TRUE),
			   hq.F_prausnitzii	 = quantile(prop.F_prausnitzii, 0.75, na.rm=TRUE),
			   avg.A_muciniphila = mean(prop.A_muciniphila, na.rm=TRUE),
               sd.A_muciniphila   = sd(prop.A_muciniphila, na.rm=TRUE),
               med.A_muciniphila   = median(prop.A_muciniphila, na.rm=TRUE),
               lq.A_muciniphila	 = quantile(prop.A_muciniphila, 0.25, na.rm=TRUE),
			   hq.A_muciniphila	 = quantile(prop.A_muciniphila, 0.75, na.rm=TRUE),
			   avg.B_uniformis = mean(prop.B_uniformis, na.rm=TRUE),
               sd.B_uniformis   = sd(prop.B_uniformis, na.rm=TRUE),
               med.B_uniformis   = median(prop.B_uniformis, na.rm=TRUE),
               lq.B_uniformis	 = quantile(prop.B_uniformis, 0.25, na.rm=TRUE),
			   hq.B_uniformis	 = quantile(prop.B_uniformis, 0.75, na.rm=TRUE),
			   avg.C_symbiosum = mean(prop.C_symbiosum, na.rm=TRUE),
               sd.C_symbiosum   = sd(prop.C_symbiosum, na.rm=TRUE),
               med.C_symbiosum   = median(prop.C_symbiosum, na.rm=TRUE),
               lq.C_symbiosum	 = quantile(prop.C_symbiosum, 0.25, na.rm=TRUE),
			   hq.C_symbiosum	 = quantile(prop.C_symbiosum, 0.75, na.rm=TRUE),
			   avg.E_coli = mean(prop.E_coli, na.rm=TRUE),
               sd.E_coli   = sd(prop.E_coli, na.rm=TRUE),
               med.E_coli   = median(prop.E_coli, na.rm=TRUE),
               lq.E_coli	 = quantile(prop.E_coli, 0.25, na.rm=TRUE),
			   hq.E_coli	 = quantile(prop.E_coli, 0.75, na.rm=TRUE),
			   avg.E_rectale = mean(prop.E_rectale, na.rm=TRUE),
               sd.E_rectale   = sd(prop.E_rectale, na.rm=TRUE),
               med.E_rectale   = median(prop.E_rectale, na.rm=TRUE),
               lq.E_rectale	 = quantile(prop.E_rectale, 0.25, na.rm=TRUE),
			   hq.E_rectale	 = quantile(prop.E_rectale, 0.75, na.rm=TRUE)
			   )
#write.table(prop.sum.stats, file="poma_summary.stats.txt", quote=FALSE, sep="\t", col.names=NA)

# can read in the produced file, and split according to group before graphing:

prop.sum.stats<-read.table(file="poma_summary.stats.txt", header=TRUE)		   

# then, you can split these into their respective groups:
m.sum.split<-split(prop.sum.stats, prop.sum.stats$group_ID)
sum.FF<-m.sum.split$'FF'
sum.FR<-m.sum.split$'FR'
sum.pre<-m.sum.split$'pre'
sum.1_FR_FF<-m.sum.split$'1_FR_FF'
sum.4_FR_FF<-m.sum.split$'4_FR_FF'
sum.4_pre_FF<-m.sum.split$'4_pre_FF'
sum.1_pre_FF<-m.sum.split$'1_pre_FF'

# each of these must be ordered by day, now:
sum.FF<-sum.FF[order(as.numeric_version(sum.FF$day)),]
sum.FF$day <- as.numeric(as.character(sum.FF$day))
sum.FR<-sum.FR[order(as.numeric_version(sum.FR$day)),]
sum.FR$day <- as.numeric(as.character(sum.FR$day))
sum.pre<-sum.pre[order(as.numeric_version(sum.pre$day)),]
sum.pre$day <- as.numeric(as.character(sum.pre$day))
sum.1_FR_FF<-sum.1_FR_FF[order(as.numeric_version(sum.1_FR_FF$day)),]
sum.1_FR_FF$day <- as.numeric(as.character(sum.1_FR_FF$day))
sum.4_FR_FF<-sum.4_FR_FF[order(as.numeric_version(sum.4_FR_FF$day)),]
sum.4_FR_FF$day <- as.numeric(as.character(sum.4_FR_FF$day))
sum.4_pre_FF<-sum.4_pre_FF[order(as.numeric_version(sum.4_pre_FF$day)),]
sum.4_pre_FF$day <- as.numeric(as.character(sum.4_pre_FF$day))
sum.1_pre_FF<-sum.1_pre_FF[order(as.numeric_version(sum.1_pre_FF$day)),]
sum.1_pre_FF$day <- as.numeric(as.character(sum.1_pre_FF$day))

# graphing median POMA
## make a heatmap of the results (as in Fig. S3):

library(RColotBrewer)
library(ggplots)

phylofrac<-read.table(file="t.m3_phylofrac_filtered.txt", header=TRUE)
	# need some metadata from this file
	
# for 1_FR_FF:
# select columns for matrix:
pomap<-subset(sum.1_FR_FF, select=c(avg.A_muciniphila, avg.B_caccae, avg.D_piger, avg.B_uniformis, avg.E_rectale, avg.B_ovatus, avg.B_theta, avg.B_formatexigens, avg.F_prausnitzii, avg.C_aerofaciens, avg.C_symbiosum, avg.R_intestinalis, avg.B_intestihominis, avg.E_coli) )
###test <- sum.1_FR_FF[, grepl(sum.1_FR_FF$day, "avg")]			#this did not work...
rownames(pomap)<-sum.1_FR_FF$day
pomatrix<-as.matrix(pomap)

# get list of day of diet:
sub <- subset(phylofrac, group_ID=="1_FR_FF", select = c(day, day_of_diet))		#subset your group of interest
diet1<-sub[!duplicated(sub$day), ]								#select the unique value to get the list of diet on a particular day
diet1$color<-mapvalues(diet1$day_of_diet, from = c("FF","FR"), to = c("red4", "chartreuse3"))
cbind(row.names(pomatrix), diet1)												#check that diet lines up

# define colors:
myCol <- c("white", "lightblue1", "skyblue1", "steelblue1", "dodgerblue", "royalblue3", "blue", "navyblue")
myBreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 1)

# better
myCol2 <- colorRampPalette(brewer.pal(9,"Blues"))(10)
myBreaks2 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# not as good
my.col <- colorRampPalette(c("aliceblue","blue4"))(n = 249)
my.breaks = c(seq(0,.1,length=50),  
               seq(.1,.2,length=50),
               seq(0.2,0.3,length=50),
               seq(0.3,0.5,length=50),
               seq(0.5,0.8,length=25),
               seq(0.8,1,length=25))

result <- heatmap.2(pomatrix,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(8,8),    
                    col=myCol2,        
                    breaks=myBreaks2,    
                    dendrogram="column",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=0.25,
                    lwid = c(1.2,5),
                    lhei = c(1.2,5),
                    main ="1-day FR/FF",
                    labRow=rownames(pomatrix),
                    labCol=colnames(pomatrix),
                    RowSideColors = as.character(diet1$color)
)

# for 4-day FR/FF:
pomap.4_FR_FF<-subset(sum.4_FR_FF, select=c(avg.A_muciniphila, avg.B_caccae, avg.D_piger, avg.B_uniformis, avg.E_rectale, avg.B_ovatus, avg.B_theta, avg.B_formatexigens, avg.F_prausnitzii, avg.C_aerofaciens, avg.C_symbiosum, avg.R_intestinalis, avg.B_intestihominis, avg.E_coli) )
rownames(pomap.4_FR_FF)<-sum.4_FR_FF$day
pomatrix.4_FR_FF<-as.matrix(pomap.4_FR_FF)

sub <- subset(phylofrac, group_ID=="4_FR_FF", select = c(day, day_of_diet))		#subset your group of interest
diet.4_FR_FF<-sub[!duplicated(sub$day), ]										#select the unique value to get the list of diet on a particular day
diet.4_FR_FF$color<-mapvalues(diet.4_FR_FF$day_of_diet, from = c("FF", "FR","pre"), to = c("red4", "Chartreuse3", "grey47"))
cbind(row.names(pomatrix.4_FR_FF), diet.4_FR_FF)												#check that diet lines up

result <- heatmap.2(pomatrix.4_FR_FF,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(8,8),    
                    col=myCol2,        
                    breaks=myBreaks2,    
                    dendrogram="column",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=0.25,
                    lwid = c(1.2,5),
                    lhei = c(1.2,5),
                    main ="4-day FR/FF",
                    labRow=rownames(pomatrix.4_FR_FF),
                    labCol=colnames(pomatrix.4_FR_FF),
                    RowSideColors = as.character(diet.4_FR_FF$color)
)
