#!/usr/bin/env Rscript

library(randomForest)

exones = read.table('all_exo_q10_mf.tsv', header=T, sep='\t')

for (i in c('agilent', 'illumina', 'roche', 'truseq', 'wgs')) {
	exones[, paste0('LCOV_', i)] = exones[, paste0('MCOV_', i)] < 0.1
}

for (i in c('agilent', 'illumina', 'roche', 'truseq', 'wgs')) {
	if (i != 'wgs') {
		fit_1 <- randomForest(as.factor(exones[, paste0('LCOV_', i)]) ~ exones$GC + exones$LEN + as.factor(exones[, paste0('INC_', i)]) + as.numeric(exones[, paste0('MF_', i)]), importance=T)
	} else {
		fit_1 <- randomForest(as.factor(exones[, 'LCOV_wgs']) ~ exones$GC + exones$LEN + as.numeric(exones[, paste0('MF_', i)]), importance=T)
	}
	dput(fit_1, file = paste0('FOREST_', i))
	print(summary(fit_1))
	print(importance(fit_1))
}
