#
#	gwasWeighted.R
#Tue Nov 17 15:52:43 CET 2020


packageDefinition = list(
	name = 'gwasWeighted',
	files = c('Rmeta.R', 'Rdata.R', 'Rsystem.R', 'Rfunctions.R', 'genomicRotation.R'),
	description = list(
		title = 'Weighted Bonferroni Correction in GWASs',
		author = 'Stefan B\uf6hringer <r-packages@s-boehringer.org>',
		description = 'Incorporate information into a GWAS analysis using related but not identical outcomes.',
		depends = c(),
		suggests = c(),
		license = 'LGPL',
		news = "0.1-0   Initial release"
	),
	git = list(
		readme = '## Installation\n```{r}\nlibrary(devtools);\ninstall_github("sboehringer/gwasWeighted")\n```\n',
		push = F,
		pushOnNewVersion = T,
		remote = 'https://github.com/sboehringer/gwasWeighted.git'
	)
);


extractWindow = function(dSource, chr, pos, windowSize) {
	window = dSource[dSource$chr == chr &
		dSource$pos >= pos - windowSize &
		dSource$pos < pos + windowSize, , ];
	return(window);
}

kernelStandard = function(window, center, windowSize) {
    weights = dnorm(window$pos, center, sd = windowSize); # we first centered at 0
    W = (window$P %*% weights / sum(weights))[1, 1]; #weighted sum
    return(W);
}

#Here's the section for standart weighted Bonferroni correction
#dTarget should be a dataframe with columns of: "marker", "P", "chr", "pos"
#dSources: list of data frames each with columns of: "W", "chr", "pos"

normMinLog = function(W)normalize(-log10(W))
normPower = function(W, power)normalize(W ** power)
combinerCorr = function(corrs)function(Ws)weight_by_corr(Ws, corrs)
# snps: can be data.frame from target study, needs: chr, pos
# combiner: min, combinerCorr(corr_coeff)
# weightTransform: identity, normMinLog
GWASweights = function(snps, dSources, windowSize = 0,
	kernel = kernelStandard, combiner = min, weightTransform = identity, power = 1, eps = 1e-5) {

	W = sapply(1:nrow(snps), function(i) { 
		snp = snps[i, ]; #snp = row of table
        
		Ws = sapply(dSources, function(dSource) {
			window = extractWindow(dSource, snp$chr, snp$pos, windowSize);
			if (length(window$chr) == 0) return(1);
			# kernel also computes weighted sum
			return(kernel(window, snp$pos, windowSize))
		});
		return(combiner(Ws));
	});
	# cave: must be compatible with combiner
	W = weightTransform(W);
	W = normPower(W, power);
	return(data.frame(snps, weight = W));
}

GWASweighted = function(dTarget, dSources, Wpower, ...) {
	Ws = GWASweights(dTarget, dSources, ...);
	P = bonferroniWeighted(dTarget$P, normPower(Ws$weight, Wpower));
	return(P);
}

weightCombiners = list(
	min = function(...)min,
	correlation = function(dTarget, dSources, precomputed) {
		correlations = if (notE(precomputed$correlations))
			precomputed$correlations else 
			corr_coeff(dTarget, dSources);
		return(function(Ws)weight_by_corr(Ws, correlations));
	}
);

runGWASweighted = function(dTarget, dSources, Wpower = 1, Kwidth = 5e4,
	weightCombiner = 'min', weightTransform = 'identity', precomputed = list()) {
	combiner = weightCombiners[[ weightCombiner ]](dTarget, dSource, precomputed);
	P = GWASweighted(dTarget, dSources,
		Wpower, windowSize = Kwidth, combiner = combiner,
		weightTransform = get(weightTransform)
	);
	return(P);
}
