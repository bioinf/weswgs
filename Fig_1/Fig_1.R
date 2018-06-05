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

# Read data and filter non-ExAC samples and other outliers
qc_data = read.csv('./WES.CDS.UPD.tsv', sep = '\t', header=T, stringsAsFactors = F)
head(qc_data)
qc_data = qc_data[qc_data$TOTAL_READS < 150000000, ]
qc_data$TOTAL_READS = qc_data$TOTAL_READS/1000000
qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), 
                              function(elt) elt[1]) == "wes_4"),]
qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), 
                           function(elt) elt[1]) == "wes_6"),]
qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), 
                           function(elt) elt[1]) == "wes_7"),]
qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), 
                           function(elt) elt[1]) == "wes_8"),]
qc_data$BAIT_SET[qc_data$BAIT_SET == 'illumina'] = 'nextera'

rownames(qc_data) = qc_data$SAMPLE
#qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), function(elt) elt[1]) == "wes_7"),]
#qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), function(elt) elt[1]) == "wes_6"),]
#qc_data = qc_data[!(sapply(strsplit(qc_data$SAMPLE, '[.]'), function(elt) elt[1]) == "wes_8"),]

qc_data_exac = qc_data[qc_data$X20X > 0.8, ]
qc_merged = rbind(qc_data, qc_data_exac)
qc_merged$filter = c(rep('All', nrow(qc_data)), rep('ExAC', nrow(qc_data_exac)))

# Function to plot boxes with jitters
#=================================================
cool_box_plot <- function(varname, xlabel='', ylabel='', 
                          lower=min(qc_merged[, varname]), 
                          upper=max(qc_merged[, varname])) {
  p <- ggplot(qc_data_exac, aes(x=BAIT_SET, y=get(varname), fill=BAIT_SET)) + 
    stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
    geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
    geom_jitter(width=0.25, size=0.65, col=mycol5) + 
    theme_bw() + theme(axis.text.x=element_blank()) + 
    scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) + 
    xlab(xlabel) + ylab(ylabel) + guides(fill=FALSE) +
    scale_y_continuous(limits=c(lower, upper))# + facet_wrap(~filter)
  return(p)
}

# Total reads and selection
panel_a = cool_box_plot('TOTAL_READS', 'Technology', 'Millions of reads', 0, 150)
panel_b = cool_box_plot('PCT_SELECTED_BASES', 'Technology', 'Fraction of CDS bases', 
                        0, 1)
panels_ab <- plot_grid(panel_a, panel_b, nrow=2, ncol=1, labels = c('A', 'B'))

# Overall fit of selection to read depth
dpth_select <- ggplot(qc_data_exac, aes(x=TOTAL_READS, y=PCT_SELECTED_BASES)) + 
  geom_point(aes(fill=BAIT_SET), colour="black", pch = 21, stroke=0.2, size=2) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) +
  geom_smooth(method='lm') + xlab('Millions of reads') + ylab('Fraction of CDS bases') +
  scale_y_continuous(limits=c(0.325, 0.825)) + guides(fill=F) +
  theme_bw()
print(dpth_select)
theme_set(theme_bw())


# The same per platform - panel C
fig_c <- ggplot(qc_data_exac, aes(x=TOTAL_READS, y=PCT_SELECTED_BASES)) + 
  geom_point(aes(fill=BAIT_SET), colour="black", pch = 21, stroke=0.8, size=3) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) +
  geom_smooth(method='lm') + xlab('Millions of reads') + ylab('Fraction of CDS bases') +
  facet_wrap(~BAIT_SET, ncol=2, nrow=2) +
  theme_bw()
print(fig_c)

# Some per-base coverage boxplots
panel_s2a <- cool_box_plot('X10X', 'Technology', 'Fraction of bases 10x', 0, 1)
panel_s2b <- cool_box_plot('X20X', 'Technology', 'Fraction of bases 20x', 0, 1)
panel_s2c <- cool_box_plot('X30X', 'Technology', 'Fraction of bases 30x', 0, 1)
panel_s2d <- cool_box_plot('X100X', 'Technology', 'Fraction of bases 100x',0, 1)
theme_set(theme_bw())
plot_grid(panel_s2a, panel_s2b, panel_s2c, panel_s2d, ncol=2, nrow=2,
          labels=c('A', 'B', 'C', 'D'))

# Mean bait coverage boxes - 1D
panel_d <- cool_box_plot('MEAN_BAIT_COVERAGE', 'Technology', 'Mean bait coverage', 0, 200)
print(panel_d)

# Normalized coverage (curves) 1E

curvefr = read.table('./complete_evenness.tab', sep='\t')
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

# Plot curves
fig_e <- ggplot(meancurves, aes(x=COV, y=FRAC, col=BAIT_SET)) +
  scale_color_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6),
                     name="Capture\ntechnology") +
  geom_line(lwd=1) + theme_bw() + guides(col=F) +
  geom_ribbon(aes(ymin=LOWER, ymax=UPPER, fill=BAIT_SET), 
              alpha=0.3, col=NA) +
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6),
                    name="Capture\ntechnology") +
  guides(fill=F, col=F) + xlab('Normalized coverage') + ylab('Fraction of bases')

# Evenness calculations - for figure 1f
toevenness = curvefr
toevenness$COV = as.numeric(toevenness$COV)
toevenness = toevenness[toevenness$COV <= 100, ]
toevenness$BINNED = toevenness$FRAC * 0.01
toevenness$SAMPLENO = rep(1:70, each=100)
evsc = aggregate(toevenness$BINNED, list(toevenness$BAIT_SET, toevenness$SAMPLENO), sum)

fig_f <- ggplot(evsc, aes(x=Group.1, y=x, fill=Group.1)) + 
  stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
  geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
  geom_jitter(width=0.25, size=0.65, col=mycol5) + 
  theme_bw() + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4, mycol6)) + 
  xlab('Technology') + ylab('Evenness') + guides(fill=FALSE) +
  scale_y_continuous(limits=c(0.65, 0.9))

theme_set(theme_bw())

# Some related SI figures
qc_basedp = qc_data[, c('SAMPLE', 'BAIT_SET', 'TOTAL_READS', 'BAIT_SET', 
                        'X2X', 'X10X', 'X20X', 'X30X', 'X40X', 'X50X', 'X100X')]

qc_basedp$X0X = 1
qc_basedp = qc_basedp[, c('SAMPLE', 'BAIT_SET', 'X0X', 'X2X',
                          'X10X', 'X20X', 'X30X', 'X40X', 'X50X', 'X100X')]

toplot = melt(qc_basedp, id.vars = c('SAMPLE', 'BAIT_SET'))
toplot$coverage = rep(c(0, 2, 10, 20, 30, 40, 50, 100), each=250)
toplot$constant = qc_data[toplot$SAMPLE, 'X20X']
toplot <- transform(toplot, SAMPLE = reorder(SAMPLE, constant))

ggplot(toplot, aes(x=SAMPLE, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='identity', width=1) +
  theme_bw() + theme(axis.text.x = element_blank()) + 
  facet_wrap(~BAIT_SET, scales='free_x')


#_================================================

ggplot(qc_data, aes(x=MEAN_BAIT_COVERAGE, y=X20X))  + theme_bw() +
  geom_point(aes(fill=BAIT_SET), size=2, pch=21, col='black') + 
  scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) +
  guides(fill=F) + geom_hline(yintercept = 0.8, lwd=1, lty=2) +
  geom_vline(xintercept = 65) + xlab('Mean CDS coverage') + 
  ylab('Fraction of CDS bases 20x') +
  scale_x_continuous(breaks=c(50, 65, 100, 150, 200))# + facet_wrap(~BAIT_SET)

#==================

cool_box_plot('FOLD_ENRICHMENT', 'Technology', 'Fold enrichment', lower=0)

# UTR_Coverage
utrqc = read.table('UTR_Coverage.tsv', header=T, row.names=1, sep='\t')
qc_utrs = qc_data_exac[rownames(qc_data_exac) %in% rownames(utrqc), ]
qc_utrs$UTR10X = utrqc[rownames(qc_utrs), 29]
qc_utrs$UTR20X = utrqc[rownames(qc_utrs), 30]
qc_utrs$UTR30X = utrqc[rownames(qc_utrs), 31]
qc_utrs$UTR100X = utrqc[rownames(qc_utrs), 34]

cool_box_plot_utr <- function(varname, xlabel='', ylabel='', 
                          lower=min(qc_merged[, varname]), 
                          upper=max(qc_merged[, varname])) {
  p <- ggplot(qc_utrs, aes(x=BAIT_SET, y=get(varname), fill=BAIT_SET)) + 
    stat_boxplot(geom='errorbar', lwd=0.5, width=0.45) + 
    geom_boxplot(width=0.6, lwd=0.5, outlier.shape = NA) + 
    geom_jitter(width=0.25, size=0.65, col=mycol5) + 
    theme_bw() + theme(axis.text.x=element_blank()) + 
    scale_fill_manual(values=c(mycol1, mycol2, mycol3, mycol4)) + 
    xlab(xlabel) + ylab(ylabel) + guides(fill=FALSE) +
    scale_y_continuous(limits=c(lower, upper))# + facet_wrap(~filter)
  return(p)
}


panel_s2a <- cool_box_plot_utr('UTR10X', 'Technology', 'Fraction of bases 10x', 0, 1)
panel_s2b <- cool_box_plot_utr('UTR20X', 'Technology', 'Fraction of bases 20x', 0, 1)
panel_s2c <- cool_box_plot_utr('UTR30X', 'Technology', 'Fraction of bases 30x', 0, 1)
panel_s2d <- cool_box_plot_utr('UTR100X', 'Technology', 'Fraction of bases 100x',0, 1)
theme_set(theme_bw())
plot_grid(panel_s2a, panel_s2b, panel_s2c, panel_s2d, ncol=2, nrow=2,
          labels=c('A', 'B', 'C', 'D'))