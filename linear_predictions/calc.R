#!/usr/bin/env Rscript

library(relaimpo)

exones = read.table('all_exo_q10_mf.tsv', header=T, sep='\t')

for (i in c('agilent', 'illumina', 'roche', 'truseq', 'wgs')) {
	if (i != 'wgs') {
		fit_1 <- lm(exones[, paste0('MCOV_', i)] ~ exones$GC + exones$LEN + as.factor(exones[, paste0('INC_', i)]) + as.numeric(exones[, paste0('MF_', i)]))
	} else {
		fit_1 <- lm(exones[, 'MCOV_wgs'] ~ exones$GC + exones$LEN + as.numeric(exones[, paste0('MF_', i)]))
	}
	dput(fit_1, file = paste0('REGRESSION_', i))
	print(summary(fit_1))
	print(calc.relimp(fit_1, type=c("lmg", "first", "last")))
}
