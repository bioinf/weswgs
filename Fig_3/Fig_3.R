library(ggplot2)
library(cowplot)
library(reshape2)
library(lattice)
library(RColorBrewer)
library(colorRamps)

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

# Creating a dataset for modeling
rv <- read.table('NEW_VALUES_WPLS.tsv', sep='\t', header=F)
coverages = read.table('q10_filtered_pfrag.hist', sep='\t')
#coverages[209538, ] = rep(0, 197)
exones <- read.table('exones_plusgc.bed', sep='\t', header=F)
exlens <- read.table('exon_sizes.tab')
inclusion <- read.table('all.inclusion_075.tab', sep='\t')
colnames(inclusion) = c('INC_agilent', 'INC_illumina', 'INC_roche', 'INC_truseq', 'INC_wgs')
colnames(exones) = c('CHR', 'START', 'END', 'GC')
exones$LEN = exlens$V1
exones <- cbind(exones, inclusion)
norm_coverages <- t(t(coverages)/rv$V2)
dim(norm_coverages)
for (i in c('agilent', 'illumina', 'roche', 'truseq', 'wgs')) {
  exones[, paste0('MCOV_', i)] = rowMeans(norm_coverages[, rv$V4 == i])
}
for (i in c('agilent', 'illumina', 'roche', 'truseq', 'wgs')) {
  exones[, paste0('SDCOV_', i)] = apply(norm_coverages[, rv$V4 == i], 1, sd)
}
mfs = read.table('all_MF_sorted.tsv', header=T, sep='\t')
exones_tt = cbind(exones, mfs)
pls = c('agilent', 'illumina', 'roche', 'truseq', 'wgs')

# Simulating values based on BIE * WIE
simvals = data.frame(coverages = seq(20, 200, 10), agilent = rep(NA, 19), illumina = rep(NA, 19), 
                     roche = rep(NA, 19), truseq = rep(NA, 19), wgs = rep(NA, 19))

for (j in 1:19){
  mcov = simvals$coverages[j]
  ncov = 10/mcov
  for (i in 1:5){
    print(paste0('Processing j ', as.character(j), ' i ', as.character(i)))
    smprof = read.table(paste0(pls[i], '.finalsmprof'), sep='\t')
    exx = as.data.frame(cbind(exones, smprof))
    simvals[j, pls[i]] = as.integer(sum(apply(exx, 1, function(elt) sum((rnorm(1, 
                                                                               mean = as.numeric(elt[10 + i]), 
                                                                               sd = as.numeric(elt[15 + i])) * as.numeric(elt[21:120])) < ncov) * as.numeric(elt[5]) * 0.01)))
  }
}

write.table(simvals, file='simulated_new.tsv', sep='\t', row.names=F)
sv = melt(simvals, id.vars='coverages')
colnames(sv) = c('V2', 'V4', 'MB')
sv$MB = sv$MB/1000000
write.table(sv, file='simvals_new.tsv', sep='\t', quote=F, row.names = F)

sv[c(1, 22, 39, 59 ), 'MB'] = 3.5
sv[96, ] = c(12, 'wgs', 3.5)
sv$V2 = as.numeric(sv$V2)
sv$MB = as.numeric(sv$MB)

# Making figure 3a

rv$MB = rv$V3/1000000
p_a <- ggplot(rv, aes(x=V2, y=MB, col=V4, fill=V4)) + geom_point() + 
  geom_line(data=sv, aes(x=V2, y=MB, col=V4), lwd=1) + 
  facet_wrap(~V4, nrow=1) + theme_bw() +
  scale_color_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) +
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + guides(col=F, fill=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(0, 230)) + scale_y_continuous(limits=c(0, 3.5))

print(p_a)

# Making figure 3b
#Correlation plots
# Distributions above heatmap

topl = melt(exones[, c(11:15)])
panel_top <- ggplot(topl, aes(x=value, fill=variable)) + 
  geom_histogram(stat='bin', binwidth = 0.01) +  
  scale_x_continuous(limits=c(0, 2)) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  geom_vline(xintercept = 0.155) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank()) + xlab('Normalized coverage') + ylab('Interval count') +
  facet_wrap(~variable, nrow=1) + guides(fill=F)

theme_set(theme_bw())
plot_grid(panel_top, panel_top, ncol=1)
print(panel_top)

# Heatmap
cors = matrix(rep(NA, 25), nrow=5, ncol=5)
rownames(cors) = pls
colnames(cors) = pls
for (i in pls){
  for (j in pls){
    cors[i, j] = mean(cor(coverages[, rv$V4 == i], coverages[, rv$V4 == j]))
  }
}

drawCoolHM = function(df){
  e = round(df, digits=2)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(-0.1, 1, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=0)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)
drawCoolHM(cors)

# Making figure 3c - GC-bias evaluation

gcbias = read.table('ALL_NEW_WGS_GCBIAS.hist', sep='\t')
colnames(gcbias) = c('V1', 'V2', 'BAIT_SET')
gcsamps = read.table('samples.list')
rownames(rv) = rv$V1
mcovs = rv[gcsamps$V1, 2]
#mcovs = sapply(1:50, function(elt) mean(gcbias[(1 + 191678 * (elt - 1)):(191678 * elt), 2], na.rm = T))
gcbias$MCOV = rep(mcovs, each=191678)
gcbias$NORMCOV = gcbias$V2 / gcbias$MCOV

panel_3d <- ggplot(gcbias, aes(x=as.factor(V1), y=NORMCOV, fill=BAIT_SET)) + 
  #  stat_boxplot(geom='errorbar', lwd=0.5, width=0.65) + 
  geom_violin() + 
  geom_boxplot(outlier.shape=NA, width=0.35, lwd=0.5, fill='white')  + 
  theme_bw() +
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) +
  scale_y_continuous(limits=c(0,5)) +
  xlab('Decile') + ylab('Coverage') + guides(fill=F) + facet_wrap(~BAIT_SET, nrow=1)
print(panel_3d)


# Figure 3e and related - multimappers and linear analysis


dataset = read.table('./all_WES_filt.tsv', sep='\t', header=F)
rlens = read.table('./read_lengths.tab', sep='\t', row.names=1)
inserts = read.table('./insert_sizes.tab', sep='\t', row.names=1)
qc_data = read.csv('./WES.CDS.UPD.tsv', sep = '\t', header=T, stringsAsFactors = F, row.names=1)
dataset$BAIT_SET = qc_data[as.character(dataset[,1]), 1]
dataset$RLEN = rlens[as.character(dataset[,1]), 1]
dataset$INSSIZE = inserts[as.character(dataset[,1]), 1]
dataset = na.omit(dataset)

genomes = read.table('./genomes.tsv', sep='\t', header=F)
genomes$BAIT_SET = rep('wgs', nrow(genomes))
genomes_lens = read.table('./wgs_lens.tsv', sep='\t', row.names = 1)
genomes_ins = read.table('./wgs_inserts.tsv', sep='\t', row.names = 1)
genomes$RLEN = genomes_lens[as.character(genomes[,1]), 1]
genomes$INSSIZE = genomes_ins[as.character(genomes[,1]), 1]
dat2 <- as.data.frame(rbind(dataset, genomes))

mf = dataset[dataset$V2 == 1, ]

# bases with MF >= Xvs insert 
ggplot(dat2, aes(x=INSSIZE, y=V3, col=BAIT_SET)) + geom_point(size=2) + 
  scale_color_manual(values = c(mycol1, mycol2, mycol3, mycol4, mycol6)) +
  theme_bw() + facet_wrap(~V2, nrow=2) + guides(col=F) + ylab('Bases with MF > N') +
  xlab('Insert size')

# bases with MF = 1 by insert by technology
ggplot(mf, aes(x=INSSIZE, y=V3, col=BAIT_SET)) + geom_point(size=2) + 
  scale_color_manual(values = c(mycol1, mycol2, mycol3, mycol4)) +
  theme_bw() + facet_wrap(~BAIT_SET, nrow=2) + guides(col=F) + ylab('Bases with MF = 1') +
  xlab('Insert size')

# bases with MF = 1 over insert by read length and technology
ggplot(mf, aes(x=INSSIZE, y=V3, col=BAIT_SET)) + geom_point(size=2) + 
  scale_color_manual(values = c(mycol1, mycol2, mycol3, mycol4)) +
  theme_bw() + facet_wrap(~BAIT_SET + RLEN, nrow=3) + guides(col=F) + ylab('Bases with MF = 1') +
  xlab('Insert size')


# Figure 3e
ggplot(dat2, aes(x=BAIT_SET, y=V3, fill=BAIT_SET)) + 
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  facet_wrap(~V2, nrow=2) + guides(fill=FALSE) + xlab('Technology') + 
  ylab('Bases with MF > N')

# Over read length
ggplot(dataset, aes(x=as.factor(RLEN), y=V3, fill=as.factor(RLEN))) + 
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  facet_wrap(~V2, nrow=2)# + guides(fill=FALSE)

fit_1 <- lm(V3 ~ RLEN, mf)
summary(fit_1)
fit_2 <- lm(V3 ~ INSSIZE, mf)
summary(fit_2)
fit_3 <- lm(V3 ~ RLEN + INSSIZE, mf)
summary(fit_3)
fit_4 <- lm(V3 ~ RLEN + INSSIZE + BAIT_SET, mf)
summary(fit_4)
 
# Dumbbells

lccds = data.frame(BAIT_SET = pls, dfs = as.numeric(table(rv$V4)) - 1, at1 = rep(NA, 5), at2 = rep(NA, 5))

for (i in 1:5){
  lccds$at1[i] = sum(apply(exones, 1, function(elt) ifelse(pt((0.1 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95, 
                                                           as.numeric(elt[3]) - as.numeric(elt[2]), 0)))
}

for (i in 1:5){
  lccds$at2[i] = sum(apply(exones, 1, function(elt) ifelse(pt((0.2 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95, 
                                                           as.numeric(elt[3]) - as.numeric(elt[2]), 0)))
}

# To inclusion
for (i in 1:5){
  vecind = apply(exones, 1, function(elt) pt((0.1 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95)
  cool_exones = exones[vecind, 1:3]
  write.table(cool_exones, file=paste0(pls[i], '_low_cov_01_CDS.bed'), quote=F, col.names=F, sep='\t', row.names = F)
}

for (i in 1:5){
  vecind = apply(exones, 1, function(elt) pt((0.2 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95)
  cool_exones = exones[vecind, 1:3]
  write.table(cool_exones, file=paste0(pls[i], '_low_cov_02_CDS.bed'), quote=F, col.names=F, sep='\t', row.names = F)
}

# Not inclusion
lccds$nincat1 = c(-210005, -211561, -347255, -216457, 0)
lccds$nincat2 = c(-227617, -223198, -395114, -223653, 0)

for (i in 1:5){
  lccds$nincat1[i] = -sum(apply(exones, 
                                1, 
                                function(elt) ifelse(pt((0.1 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95 & elt[5 + i] == 0, 
                                                     as.numeric(elt[3]) - as.numeric(elt[2]), 0)))
}

for (i in 1:5){
  lccds$nincat2[i] = -sum(apply(exones, 
                                1, 
                                function(elt) ifelse(pt((0.2 - as.numeric(elt[10 + i]))/as.numeric(elt[15 + i]), lccds$dfs[i]) > 0.95 & elt[5 + i] == 0, 
                                                     as.numeric(elt[3]) - as.numeric(elt[2]), 0)))
}

lccds$at1 = lccds$at1 + lccds$nincat1
lccds$at2 = lccds$at2 + lccds$nincat2

stuff = melt(lccds, id.vars = c('BAIT_SET', 'dfs'))
stuff$COVERAGE = rep(c(rep('10', 5), rep('20', 5)), 2)
stuff$value = stuff$value/1000
ggplot(stuff, aes(x=BAIT_SET, y=value, size=COVERAGE, col=BAIT_SET)) + geom_point(position=position_dodge(width=0.64)) + 
  scale_color_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) +
  scale_size_manual(values=c(5, 3)) + geom_hline(yintercept = 0) +
  geom_linerange(aes(x=BAIT_SET, ymax=value, ymin=0, group=COVERAGE), size=1.5, position=position_dodge(width=0.64)) +
  theme_bw() + coord_flip() + scale_y_reverse(limits=c(800, -800)) + guides(fill=F, col=F, size=F) + 
  ylab('CDS with coverage less that X') + xlab('Technology')



# Writing dataset to make fits
exones_tt = exones_tt[exones$LEN > 0, ]
exones_tt = exones_tt[exones$CHR != 'Y', ]
write.table(x = exones_tt, 'all_exo_q10_mf.tsv', sep='\t', quote = F)



