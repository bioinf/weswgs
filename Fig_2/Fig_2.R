library(ggplot2)
library(reshape2)
library(cowplot)

# Color scheme
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

# BIE analysis - Figure 2a

curvefr = read.table('./WES_pfrag_curves.tsv', sep='\t')
wgscurvefr = read.table('WGS_pfrag_curve.tsv', sep='\t')
curvefr = rbind(curvefr, wgscurvefr)
colnames(curvefr) = c('COV', 'FRAC', 'BAIT_SET')
curvefr[curvefr$COV == 0.00, 2] = 1.0
curvefr$COV = as.factor(curvefr$COV)

meancurves = aggregate(curvefr$FRAC, list(curvefr$COV, curvefr$BAIT_SET), mean)
sdcurves = aggregate(curvefr$FRAC, list(curvefr$COV, curvefr$BAIT_SET), sd)
colnames(meancurves) = c('COV', 'BAIT_SET', 'FRAC')
colnames(sdcurves) = c('COV', 'BAIT_SET', 'FRAC')

meancurves$LOWER = meancurves$FRAC - (sdcurves$FRAC/sqrt(10))
meancurves$UPPER = meancurves$FRAC + (sdcurves$FRAC/sqrt(10))
meancurves$COV = as.numeric(meancurves$COV) * 0.01

panel_a <- ggplot(meancurves, aes(x=COV, y=FRAC, col=BAIT_SET)) +
  scale_color_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6),
                     name="Capture\ntechnology") +
  geom_line(lwd=1) + theme_bw() + guides(col=F) +
  geom_ribbon(aes(ymin=LOWER, ymax=UPPER, fill=BAIT_SET), 
              alpha=0.3, col=NA) +
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6),
                    name="Capture\ntechnology") +
  guides(fill=F, col=F) + xlab('Normalized coverage') + ylab('Fraction of intervals')

print(panel_a)

# BIE scores boxplot 0 Figure 2b (!
# IMPORTANT
# all samples used (only 10 of those in the main paper text!!!))

toevenness = curvefr
toevenness$COV = as.numeric(toevenness$COV)
toevenness = toevenness[toevenness$COV <= 100, ]
toevenness$BINNED = toevenness$FRAC * 0.01
toevenness$SAMPLENO = rep(1:177, each=100)
evsc = aggregate(toevenness$BINNED, list(toevenness$BAIT_SET, toevenness$SAMPLENO), sum)

panel_b <- ggplot(evsc, aes(x=Group.1, y=x, fill=Group.1)) + 
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('BIE') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0.65, 0.925))

print(panel_b)

# WIE profiles - Fig 2d
pbase_roche = read.table('roche.hist',
                         sep = ' ', dec = '.', stringsAsFactors = F)[, 1:100]
se_roche = apply(pbase_roche, 2, function(elt) sd(elt)/sqrt(length(elt)))


pbase_illumina = read.table('illumina.hist', 
                            sep = ' ', dec = '.', stringsAsFactors = F)[, 1:100]
se_illumina = apply(pbase_illumina, 2, function(elt) sd(elt)/sqrt(length(elt)))


pbase_truseq = read.table('truseq.hist', 
                          sep = ' ', dec = '.', stringsAsFactors = F)[, 1:100]
se_truseq = apply(pbase_truseq, 2, function(elt) sd(elt)/sqrt(length(elt)))

pbase_agilent = read.table('agilent.hist', 
                           sep = ' ', dec = '.', stringsAsFactors = F)[, 1:100]
se_agilent = apply(pbase_agilent, 2, function(elt) sd(elt)/sqrt(length(elt)))


pbase_wgs = read.table('wgs.pbase.hist', sep = ' ', dec = '.', stringsAsFactors = F)[, 2:101]
colnames(pbase_wgs) = colnames(pbase_truseq)
se_wgs = apply(pbase_wgs, 2, function(elt) sd(elt)/sqrt(length(elt)))

myline = as.data.frame(colMeans(pbase_roche))
myline$reldist = seq(0, 100, length.out = 100)
myline$nextera = colMeans(pbase_illumina)
myline$truseq = colMeans(pbase_truseq)
myline$agilent = colMeans(pbase_agilent)
myline$wgs = colMeans(pbase_wgs)
colnames(myline) = c('roche', 'reldist', 'nextera', 'truseq', 'agilent', 'wgs')

topl = melt(myline, id.vars = 'reldist')
by(topl$value, topl$variable, function(elt) sum(elt*0.01))

sds = c(se_roche, se_illumina, se_truseq, se_agilent, se_wgs)
topl$y_min = topl$value - sds
topl$y_max = topl$value + sds
panel_d <- ggplot(topl, aes(x = reldist, y = value, col = variable)) + 
  geom_line(lwd=0.6) + geom_ribbon(aes(x=reldist, ymin=y_min, ymax=y_max, fill=variable), 
                                   alpha=0.6, col=NA) +
  theme_bw() + xlab('Relative distance') + ylab('Relative coverage') +
  scale_fill_manual(values=c(mycol3, mycol2, mycol4, mycol1, mycol6), name='Capture\ntechnology') + 
  scale_color_manual(values=c(mycol3, mycol2, mycol4, mycol1, mycol6), name='Capture\ntechnology') + 
  guides(fill=F, col=F) +
  scale_y_continuous(limits = c(0.5, 1.5))

print(panel_d)

# WIE score calculation and plotting - figure 2e

pbase_all = rbind(pbase_agilent, pbase_illumina, pbase_roche, pbase_truseq, pbase_wgs)
pbase_all[pbase_all > 1.0] = 1.0
meanvals = apply(pbase_all, 1, sum)
toboxpl = data.frame(BAIT_SET = c(rep(c('agilent', 'nextera', 'roche', 'truseq'), each=10), rep('wgs', 30)), 
                     SMOOTH = meanvals)
toboxpl$SMOOTH = toboxpl$SMOOTH * 0.01

panel_e <- ggplot(toboxpl, aes(x=BAIT_SET, y=SMOOTH, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('Smoothness') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0.85, 1))

print(panel_e)

# Statistical evaluation
pairwise.wilcox.test(toboxpl$SMOOTH, toboxpl$BAIT_SET, p.adjust.method = 'bonf')
# All differences well significant