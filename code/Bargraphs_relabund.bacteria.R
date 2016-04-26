## Bargraphs of combined mucus-degrading bacteria: Fig. 2D, 2E
## Anna M. Seekatz
## 4.22.16

# This code was used to generate bargraphs of the relative abundance of mucus-degrading bacteria (individual or combined)
	# Fig 2D: % relative abundance of mucus-degrading bacteria in cecum (combined )
	# Fig. 2E: comparison of the relative abundance of different bacteria in lumen or mucus samples

# files used:
	# t.m3.cecal_phylofrac_filtered.txt
	
	
###
#-----------
###

# Fig. 2D:

f2.m<-read.table(file="t.m3.cecal_phylofrac_filtered.txt", header=TRUE)

# now add relative abundance of mucin degraders together (A. muciniphila, B. caccae, B. intestihominis, B. theta)
f2.m$mucin<-f2.m$B_theta+f2.m$B_caccae+f2.m$B_intestihominis+f2.m$A_muciniphila
f2.m$mucin2<-f2.m$A_muciniphila+f2.m$B_caccae
f2.m$fiber<-f2.m$B_ovatus+f2.m$E_rectale
f2.m<-droplevels(f2.m)

# get median/interquartiles, etc
#note: this script does not work on a single variable (but you can graph them as single later)
median<- sapply(levels(f2.m$group_ID), function(class) sapply(f2.m[f2.m$group_ID == class, c("mucin", "mucin2", "fiber", "A_muciniphila", "B_theta", "B_caccae", "B_intestihominis")], median))
min<- sapply(levels(f2.m$group_ID), function(class) sapply(f2.m[f2.m$group_ID == class, c("mucin", "mucin2", "fiber", "A_muciniphila", "B_theta", "B_caccae", "B_intestihominis")], min))
max<- sapply(levels(f2.m$group_ID), function(class) sapply(f2.m[f2.m$group_ID == class, c("mucin", "mucin2", "fiber", "A_muciniphila", "B_theta", "B_caccae", "B_intestihominis")], max))
lq<- sapply(levels(f2.m$group_ID), function(class) sapply(f2.m[f2.m$group_ID == class, c("mucin", "mucin2", "fiber", "A_muciniphila", "B_theta", "B_caccae", "B_intestihominis")], function(x) quantile(x)[2]))
hq<- sapply(levels(f2.m$group_ID), function(class) sapply(f2.m[f2.m$group_ID == class, c("mucin", "mucin2", "fiber", "A_muciniphila", "B_theta", "B_caccae", "B_intestihominis")], function(x) quantile(x)[4]))

# reorder to desired x-axis order:
median<-median[, c("FR", "FF", "pre", "1_FR_FF", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
min<-min[, c("FR", "FF", "pre", "1_FR_FF", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
max<-max[, c("FR", "FF", "pre", "1_FR_FF", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
lq<-lq[, c("FR", "FF", "pre", "1_FR_FF", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
hq<-hq[, c("FR", "FF", "pre", "1_FR_FF", "4_FR_FF", "1_pre_FF", "4_pre_FF")]

## Fig. 2D: to graph total mucin degraders for all groups (option1):
gr.col<-c("chartreuse4", "red", "grey67", "blue3", "black", "brown", "purple4")
mp1<-barplot(median[1,], beside=TRUE, col=gr.col, las=2, ylim=c(0,80), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,35))
segments(mp1, max[1,], mp1, min[1,])

# to graph total mucin degraders for only the groups FF, FR, pre, and 1day oscillation (option2):
sub.median<-median[, select=c("FR", "FF", "pre", "1_FR_FF")]
sub.max<-max[, select=c("FR", "FF", "pre", "1_FR_FF")]
sub.min<-min[, select=c("FR", "FF", "pre", "1_FR_FF")]

# graph plots together:
gr2.col<-c("chartreuse4", "red", "grey67", "blue3")
par(mfrow=c(1,2))
mp2<-barplot(sub.median[2,], beside=TRUE, col=gr2.col, las=2, ylim=c(0,80), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,5), names.arg="")
segments(mp2, sub.max[2,], mp2, sub.min[2,])
text(x =  mp2, y = par("usr")[3]-3, srt = 45, adj = 1, labels = c("FR", "FF", "pre", "1-day FR/FF"), xpd = TRUE, cex=0.8)
title('Mucin degraders', cex.main=0.75)

mp3<-barplot(sub.median[3,], beside=TRUE, col=gr2.col, ylim=c(0,40), ylab="Relative Abundance (%)", xlim=c(0,5), names.arg="")
segments(mp3, sub.max[3,], mp3, sub.min[3,])
title('Fiber degraders', cex.main=0.75)
text(x =  mp2, y = par("usr")[3]-3, srt = 45, adj = 1, labels = c("FR", "FF", "pre", "1-day FR/FF"), xpd = TRUE, cex=0.8)

#supplemental figures:
# graph the remaining three groups (with FR as a control):
sub.median.sup<-median[, select=c("FR", "pre", "1_pre_FF", "4_pre_FF")]
sub.max.sup<-max[, select=c("FR", "pre", "1_pre_FF", "4_pre_FF")]
sub.min.sup<-min[, select=c("FR", "pre", "1_pre_FF", "4_pre_FF")]
gr2.sup.col<-c("chartreuse4", "grey67", "purple", "purple4")
par(mfrow=c(1,2))
mp2<-barplot(sub.median.sup[2,], beside=TRUE, col=gr2.sup.col, las=2, ylim=c(0,80), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,5), names.arg="")
segments(mp2, sub.max.sup[2,], mp2, sub.min.sup[2,])
text(x =  mp2, y = par("usr")[3]-3, srt = 45, adj = 1, labels = c("FR", "pre", "1-day pre/FF", "4-day pre/FF"), xpd = TRUE, cex=0.8)
title('Mucin degraders', cex.main=0.75)
	#with others:
sub.median.sup2<-median[, select=c("FR", "FF", "pre", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
sub.max.sup2<-max[, select=c("FR", "FF", "pre", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
sub.min.sup2<-min[, select=c("FR", "FF", "pre", "4_FR_FF", "1_pre_FF", "4_pre_FF")]
gr2.sup2.col<-c("chartreuse4", "red", "grey67", "blue", "purple", "purple4")
mp2<-barplot(sub.median.sup2[2,], beside=TRUE, col=gr2.sup2.col, las=2, ylim=c(0,80), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,5), names.arg="")
segments(mp2, sub.max.sup2[2,], mp2, sub.min.sup2[2,])
text(x =  mp2, y = par("usr")[3]-3, srt = 45, adj = 1, labels = c("FR", "FF", "pre", "4-day FR/FF", "1-day pre/FF", "4-day pre/FF"), xpd = TRUE, cex=0.8)
title('Mucin degraders', cex.main=0.75)

## stats:
# kruskal-wallis test:
kruskal.test(mucin2~group_ID, data=f2.m)
	#data:  mucin2 by group_ID
	#Kruskal-Wallis chi-squared = 16.7708, df = 6, p-value = 0.01016
kruskal.test(fiber~group_ID, data=f2.m)
	#data:  fiber by group_ID
	#Kruskal-Wallis chi-squared = 19.3004, df = 6, p-value = 0.003685
kruskal.test(mucin~group_ID, data=f2.m)
	#data:  mucin by group_ID
	#Kruskal-Wallis chi-squared = 14.083, df = 6, p-value = 0.02872
	
# anova:
fit<-aov(mucin~group_ID, data=f2.m)
summary(fit)
            Df Sum Sq Mean Sq F value   Pr(>F)    
#group_ID     6 1442.7  240.46   15.71 1.06e-05 ***
#Residuals   15  229.6   15.31 
TukeyHSD(fit)
#Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = mucin ~ group_ID, data = f2.m)
#
#$group_ID
#                         diff        lwr        upr     p adj
#1_pre_FF-1_FR_FF    0.7292792 -10.071624  11.530182 0.9999834
#4_FR_FF-1_FR_FF     5.4804554  -5.320448  16.281358 0.6168819
#4_pre_FF-1_FR_FF    7.5691100  -3.231793  18.370013 0.2755756
#FF-1_FR_FF          4.6295144  -6.171389  15.430417 0.7683439
#FR-1_FR_FF        -16.6012317 -26.704551  -6.497912 0.0008564*
#pre-1_FR_FF         2.1228773  -8.678026  12.923780 0.9927975
#4_FR_FF-1_pre_FF    4.7511762  -6.049727  15.552079 0.7478872
#4_pre_FF-1_pre_FF   6.8398308  -3.961072  17.640734 0.3789550
#FF-1_pre_FF         3.9002352  -6.900668  14.701138 0.8754464
#FR-1_pre_FF       -17.3305109 -27.433830  -7.227191 0.0005521*
#pre-1_pre_FF        1.3935981  -9.407305  12.194501 0.9992864
#4_pre_FF-4_FR_FF    2.0886546  -8.712248  12.889558 0.9933895
#FF-4_FR_FF         -0.8509410 -11.651844   9.949962 0.9999589
#FR-4_FR_FF        -22.0816871 -32.185007 -11.978368 0.0000375*
#pre-4_FR_FF        -3.3575781 -14.158481   7.443325 0.9330303
#FF-4_pre_FF        -2.9395955 -13.740498   7.861307 0.9633154
#FR-4_pre_FF       -24.1703417 -34.273661 -14.067022 0.0000127*
#pre-4_pre_FF       -5.4462327 -16.247136   5.354670 0.6231996
#FR-FF             -21.2307461 -31.334066 -11.127427 0.0000594*
#pre-FF             -2.5066371 -13.307540   8.294266 0.9830713
#pre-FR             18.7241090   8.620789  28.827429 0.0002431*

###
#-----------
###

# Fig. 2E:

library(plyr)
	
# received following sheet from Mahesh (of previous data):
lcm<-read.table(file="LCM_data.txt", header=TRUE)

# one way:
med.lcm<-ddply(lcm, c("site", "diet"), colwise(median, is.numeric))
max.lcm<-ddply(lcm, c("site", "diet"), colwise(max, is.numeric))
min.lcm<-ddply(lcm, c("site", "diet"), colwise(min, is.numeric))
	
# better way (this way was used):
lcm.stats <- ddply(lcm, c("site", "diet"), summarise,
               N    = length(E_rectale),
               med.B_caccae = mean(B_caccae, na.rm=TRUE),
               lq.B_caccae	 = quantile(B_caccae, 0.25, na.rm=TRUE),
			   hq.B_caccae	 = quantile(B_caccae, 0.75, na.rm=TRUE),
			   med.A_municiphila = mean(A_municiphila, na.rm=TRUE),
               lq.A_municiphila	 = quantile(A_municiphila, 0.25, na.rm=TRUE),
			   hq.A_municiphila	 = quantile(A_municiphila, 0.75, na.rm=TRUE),
			   med.E_rectale = mean(E_rectale, na.rm=TRUE),
               lq.E_rectale	 = quantile(E_rectale, 0.25, na.rm=TRUE),
			   hq.E_rectale	 = quantile(E_rectale, 0.75, na.rm=TRUE),
			   med.B_ovatus = mean(B_ovatus, na.rm=TRUE),
               lq.B_ovatus	 = quantile(B_ovatus, 0.25, na.rm=TRUE),
			   hq.B_ovatus	 = quantile(B_ovatus, 0.75, na.rm=TRUE)
			   )
lcm.stats<-lcm.stats[order(lcm.stats$diet, lcm.stats$site),]
mp1<-barplot(c(lcm.stats$med.B_caccae, lcm.stats$med.A_municiphila, lcm.stats$med.E_rectale, lcm.stats$med.B_ovatus), beside=TRUE, col=c("red", "red", "chartreuse4", "chartreuse4"), las=2, ylim=c(0,60), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,20), space=c(.6,.2,.2,.2))
mp1<-barplot(c(lcm.stats$med.B_caccae, lcm.stats$med.A_municiphila, lcm.stats$med.E_rectale, lcm.stats$med.B_ovatus), beside=TRUE, col=c("black"), las=2, ylim=c(0,60), ylab="Relative Abundance (%)", cex.names=0.8, xlim=c(0,20), add=TRUE, density=c(0,20,0,20), space=c(.6,.2,.2,.2))
segments(mp1, c(lcm.stats$lq.B_caccae, lcm.stats$lq.A_municiphila, lcm.stats$lq.E_rectale, lcm.stats$lq.B_ovatus), mp1, c(lcm.stats$hq.B_caccae, lcm.stats$hq.A_municiphila, lcm.stats$hq.E_rectale, lcm.stats$hq.B_ovatus))
#text(x =  1.5, y = -4, srt = 45, adj = 1, labels = c("lumen"), xpd = TRUE, cex=0.8)
#text(x =  3.5, y = -4, srt = 45, adj = 1, labels = c("mucus"), xpd = TRUE, cex=0.8)
text(x =  6.5, y = -4, adj = 2, labels = c("B. caccae"), xpd = TRUE, cex=0.8)
text(x =  13.5, y = -4, adj = 2, labels = c("A. municiphila"), xpd = TRUE, cex=0.8)
text(x =  16.5, y = -4, adj = 2, labels = c("E. rectale"), xpd = TRUE, cex=0.8)
text(x =  21.5, y = -4, adj = 2, labels = c("B. ovatus"), xpd = TRUE, cex=0.8)
legend(18, 40, c("FF", "FR"), col=c("red", "chartreuse4"), pch=15, cex=0.8)
legend(18, 60, c("lumen", "mucus"), col=c("black"), density=c(0,20), cex=0.8)

# stats for these:
# using wilcoxon test as before:

# by site (mucus vs. lumen):
# must regroup based on lumen/mucus status:
lcm<-read.table(file="LCM_data.txt", header=TRUE)
diet.split<-split(lcm, lcm$diet)
sum.FF<-diet.split$'FF'
sum.FR<-diet.split$'FR'

# for each comparison (within each diet group):
wilcox.test(B_caccae~site, data=sum.FF)
	#data:  B_caccae by site
	#W = 2, p-value = 0.1143

