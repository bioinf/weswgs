library(ggplot2)
library(reshape2)
library(cowplot)

# Color scheme
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol1a = rgb(85, 138, 221, alpha=120, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol2a = rgb(255, 192, 78, alpha=120, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol3a = rgb(255, 121, 177, alpha=120, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol4a = rgb(221, 221, 221, alpha=120, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

# Boxplots for numbers of variants
hist2 = read.table('NUMVARS_RESTR_WGS.tsv', header=F, sep='\t')
colnames(hist2) = c('NUMBER', 'BAIT_SET', 'TYPE')
#hist2 = hist2[c(1:6,8:80), ]
hist3 = read.table('NUM_NR_WGS.tsv', header=F, sep='\t')
colnames(hist3) = c('NUMBER', 'BAIT_SET', 'TYPE')
pooled_hist = rbind(hist2, hist3)
pooled_hist$RESTR = rep(c('Targeted CDS', 'All CDS'), each=94)

hist_all = pooled_hist[hist2$TYPE == 'ALL', ]
hist_low = hist2[hist2$TYPE == 'LOWGQ', 1:2]

hist_all$NUMBER = hist_all$NUMBER/1000
panel_a <- ggplot(hist_all, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('Number of variants (k)') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(22, 27)) + facet_wrap(~RESTR, nrow=2)
print(panel_a)

panel_b <- ggplot(hist_low, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('lowGQ variants') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0, 1400))
print(panel_b)

#+========================================
# Indel analysis - Figure 4c
hist4 = read.table('INDELS_NEW_WGS.tsv', header=F, sep='\t')
colnames(hist4) = c('NUMBER', 'BAIT_SET', 'TYPE')

hist_all = hist4[hist4$TYPE == 'ALL', ]
hist_low = hist4[hist4$TYPE == 'LOWGQ', 1:2]

hist_all$NUMBER = hist_all$NUMBER/1000
panel_c1 <- ggplot(hist_all, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('Number of variants (k)') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0.75, 1.25))# + facet_wrap(~RESTR, nrow=2)

panel_c2 <- ggplot(hist_low, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('lowGQ variants') + guides(fill=FALSE)
#  scale_y_continuous(limits=c(0, 2000))

plot_grid(panel_c1, panel_c2, nrow=1, labels=c('C', ''))

# Allele ratio comparison = Figure 4d
hist = read.table('SNP_HETS_WGS.tsv', header=F, sep='\t', dec='.')
colnames(hist) = c('AR', 'BAIT_SET')

panel_d <- ggplot(hist, aes(x=BAIT_SET, y=AR, fill=BAIT_SET)) + geom_violin(lwd=0.5) + 
  geom_boxplot(fill='white', width=0.2, lwd=0.5, outlier.shape=NA) +
  theme_bw() + theme(axis.text.x=element_blank()) + scale_y_continuous(limits=c(0, 1)) +
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) +
  xlab('Technology') + ylab('Allele ratio') + guides(fill=F)

print(panel_d)
by(hist$AR, hist$BAIT_SET, function(elt) 1 - median(elt))

# Mappability bootstrapping 
# Density calculations

dvars = read.table('./ExAC_All_Counts.bed', sep='\t', header=F)
cols = 1:3
rownames(dvars) = do.call(paste, c(dvars[cols], sep="-"))
tregs = read.table('./target_regions.bed', sep='\t', header=F)
tset = do.call(paste, c(tregs[cols], sep="-"))
mean(dvars[tset, 'V5'])
mean(dvars[, 'V5'])

non_mf = dvars[!(rownames(dvars) %in% tset), ]
bs_vals = c()
for (i in 1:10000){
  bs_vals = c(bs_vals, mean(sample(non_mf[, 'V5'], size=length(tset))))
}
exp = (1 - mean(tregs[,4])) * mean(bs_vals)
print(exp)

tdf = data.frame(bs = bs_vals)
ggplot(tdf, aes(x=bs)) + geom_density(fill='red') + scale_x_continuous(limits=c(0, 0.16)) + 
  geom_vline(xintercept = mean(dvars[tset, 'V5'])) + theme_bw()

mfs = read.table('all_MF_sorted.tsv', header=T, sep='\t')
mfs$DV = dvars[, 'V5']
ggplot(mfs, aes(x=MF_agilent, y=DV)) + geom_point() + theme_bw() +
  scale_y_continuous(limits = c(0,1)) + xlab('Multimapper coverage fraction') + 
  ylab('ExAC variant site density per nucl.')

#==============================================
# Variants inside bait regions - SI
hist2 = read.table('bait_vars.hist', header=F, sep='\t')
colnames(hist2) = c('NUMBER', 'BAIT_SET', 'TYPE')
hist2 = hist2[c(1:6,8:80), ]

hist_all = hist2[hist2$TYPE == 'ALL', 1:2]
hist_low = hist2[hist2$TYPE == 'LOWGQ', 1:2]

panel_3a <- ggplot(hist_all, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) + 
  xlab('Technology') + ylab('Number of variants') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0, 65000))

panel_3b <- ggplot(hist_low, aes(x=BAIT_SET, y=NUMBER, fill=BAIT_SET)) +
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) + 
  xlab('Technology') + ylab('lowGQ variants') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0, 3500))

theme_set(theme_bw())
plot_grid(panel_3a, panel_3b, nrow=1, labels=c('A', 'B'))