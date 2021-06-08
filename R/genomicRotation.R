#
#	genomicRotation.R
#Mon Dec 17 15:45:47 CET 2018

#
#	<p> genomic rotation
#

# df: data.frame containing row-wise genomic information
# N: the shift number
genomicRotationIdcs = function(df, Nshift) {
	if (Nshift > nrow(df)) {
		#warning('Shift too large');
		Nshift = Nshift %% nrow(df);
	}
	#	print(Nshift);
	indeces = c(Seq(Nshift+1, nrow(df)), Seq(1, Nshift));
	return(indeces);
}
# columns: which columns to rotate
genomicRotation = function(df, Nshift, columns = 'P') {
	indeces = genomicRotationIdcs(df, Nshift);
	df[, columns] = df[indeces, columns, drop = F];
	df
}

# rotate by percentage
# perc: fraction of total SNP number by which to shift P-values (columns argument), rotation is performed
#	per chromosome
# orderBy: columns by which to order data frames first, needed to assure alignment of later aggregation
# columns: which columns to rotate
genomicRotationPerc = function(d, perc, columns = 'P', by = 'chr', orderBy = c('chr', 'pos')) {
	dO = d[order.df(d, orderBy), , drop = F];
	Nshift = round(perc * nrow(d));
	# rotation per chromosome
	idcsRot = by(dO, dO[[by]], function(d)genomicRotationIdcs(d, Nshift));
	# compute permutation
	Nrow = sapply(idcsRot, length);
	Ioff = pop(c(0, cumsum(Nrow)));
	idcsL = lapply(seq_along(Ioff), function(i)idcsRot[[i]] + Ioff[i]);
	idcs = unlist(lapply(seq_along(Ioff), function(i)idcsRot[[i]] + Ioff[i]));
	dO[, columns] = dO[idcs, columns, drop = F];
	return(dO);
}

# rotate several data.frames by percentage
genomicRotationsPerc = function(dfs, perc, columns = 'P', by = 'chr', orderBy = c('chr', 'pos')) {
	return(lapply(dfs, genomicRotationPerc, perc = perc, columns = columns, by = by, orderBy = orderBy));
}

# rotate several data.frames by different percentages, number of data.frames and percentages has to be the same
genomicRotationsPercs = function(dfs, percs, columns = 'P', by = 'chr', orderBy = c('chr', 'pos')) {
	return(lapply(seq_along(percs), function(i) {
		genomicRotationPerc(dfs[[i]], percs[i], by = by, orderBy = orderBy)
	}));
}

#
#	<p> data simulation
#

# Nstudies: total number of studies
# rotation: percentages of SNP positions between which to rotate
simulateNull = function(d0, Nstudies = 3, rotation = c(.004, 04), by = 'chr', orderBy = c('chr', 'pos')) {
	ds = rep(list(d0), Nstudies);
	rot = runif(Nstudies - 1, rotation[1], rotation[2]);
	dfs = genomicRotationsPercs(ds, c(0, rot), 'P', by = by, orderBy = orderBy);
	return(dfs);
}

modifyRegion = function(d0, region, rotation = c(0, 5), attenuation = .99, noise = .01, Nregion = 5e4) {
	Is = with(region, which(d0$chr == chr & d0$pos >= pos - Nregion & d0$pos <= pos + Nregion));
	d0R = d0[Is, , drop = F];
	Nsnps = nrow(d0R);
	rotA = runif(1, rotation[[1]]/Nsnps, rotation[[2]]/Nsnps);
	d0Rot = genomicRotationPerc(d0R, rotA, columns = 'P');
	Tmod = qchisq(d0Rot$P, 1, lower.tail = F) * attenuation * rlnorm(Nsnps, 0, noise);
	d0Rot$P = pchisq(Tmod, 1, lower.tail = F);
	return(list(Is = Is, data = d0Rot));
}

# Simulation idea: identify top-hits from target study, copy these signals to a number of discovery studies,
#	modify the signal by rotation, fixed attenuation, random pertubation
# rotationAlt: rotation of top-region after copying to discovery study, given in number of expected loci
# Ncopy: to how many discovery studies is a signal copied
# attenuation: fixed factor to attenuate signal from the target study on chisq scale
# noise: standard deviation of lognormal factor to scale signals before copying
simulateAltTopHits = function(d0, Nstudies = 3, rotation = c(.004, 04),
	Ptop = 2e-6, Nregion = 5e4,
	rotationAlt = c(1, 5), Ncopy = 1, attenuation = .95, noise = .05) {
	ds0 = simulateNull(d0, Nstudies, rotation);
	th = droplevels(snpsTopHitsFromTable(ds0[[1]], P ~ chr + pos, P = Ptop, Nregion = Nregion));
	for (i in Seq(1, nrow(th))) {
		for (j in 1:Ncopy) {
			d0Reg = modifyRegion(ds0[[1]], as.list(th[i, ]), rotationAlt, attenuation, noise, Nregion);
			ds0[[j + 1]][d0Reg$Is, ] = d0Reg$data;
		}
	}
	return(ds0);
}

#
#	<p> analysis
#

sumlog = function (p) {
	keep <- !is.na(p) & (p > 0) & (p <= 1);
	if (sum(keep) == 0) return(NA);
	if (sum(keep) == 1) return(list(p = p[keep], df = 2));
	lnp <- log(p[keep]);
	chisq <- (-2) * sum(lnp);
	df <- 2 * length(lnp);
	r = list(chisq = chisq, df = df, p = pchisq(chisq, df,  lower.tail = FALSE), validp = p[keep])
	return(r)
}

Sumlog = function(p)sumlog(p)$p

analyzeDataMeta = function(dfs, pars = list()) {
	Ps = do.call(cbind, list.kp(dfs, 'P'));
	Ina = which(apply(Ps, 1, function(p)any(is.na(p))));
	# <!> cave bias
	if (length(Ina) > 0) Ps = Ps[- Ina, , drop = F];
	Pmeta = apply(Ps, 1, Sumlog);
	r = list(Pmin = min(Pmeta, na.rm = T), TS = tailStrength(Pmeta), gini = Gini(Pmeta, na.rm = T));
	return(r);
}

analyzeDataWB = function(dfs, pars) {
	r = run_pipeline(dfs[[1]], dfs[-1], pars$power, pars$kernel, weightCombiner = 'min');
	return(r);
}

analyzeDataAll = function(dfs, pars) {
	r = list(analyzeDataMeta(dfs, pars), analyzeDataWB(dfs, pars));
	return(r);
}

#
#	<p> simulation runs
#

runSimulationSingle = function(d0, Nstudies = 3, rotation = c(.01, .05), ...) {
	ds = simulateNull(d0, Nstudies, rotation);
	r = analyzeDataAll(ds, ...);
	return(r);
}


runSimulationNull = function(d0, Nstudies = 3, Nrep = 10, rotation = c(.01, .05), ...) {
	r = lapply(Seq(1, Nrep), function(i) {
		runSimulationSingle(d0, Nstudies, rotation, ...)
	});
	return(r);
}

#Nstudies, rotation, rotationAlt, Ncopy, attenuation, noise
runSimulationAltSingle = function(d0, Nstudies = 3, rotation = c(.004, 04),
	Ptop = 2e-6, Nregion = 5e4,
	rotationAlt = c(1, 5), Ncopy = 1, attenuation = .95, noise = .05, ...) {

	ds = simulateAltTopHits(d0, Nstudies, rotation, Ptop, Nregion, rotationAlt, Ncopy, attenuation, noise);
	r = analyzeDataAll(ds, ...);
	return(r);
}

#
#	<p> data management
#

# in vector v, find index min j \in 1, ..., N so that v[1:j] contains at least U unique elements
uniqueIndex = function(v, U) {
	#Nu = sapply(seq_along(v), function(i)length(unique(data$chr[1:i])));
	# more efficient version
	u = c();
	for (i in seq_along(v)) {
		u = unique(c(u, v[i]));
		if (length(u) == U) return(i);
	}
	return(NA);
}

selectTopChrs = function(data, Nchrs = 3, fCols = P ~ chr + pos, P = 5e-4, Nregion = 5e4) {
	th = droplevels(snpsTopHitsFromTable(data, fCols, P = P, Nregion = Nregion));
	if (nrow(th) == 0) stop('no top hits identified');
	Ichr = uniqueIndex(th$chr, Nchrs);
	r = th[Seq(1, Ichr), , drop = F];
	return(r);
}

reduceDataSingle = function(d, Nchrs = 3, fCols = P ~ chr + pos, P = 5e-4, Nregion = 5e4, NselRange = NULL,
	chrs = NULL) {
	if (is.null(chrs)) chrs = selectTopChrs(d, Nchrs, fCols, P, Nregion);
	dSel = if (notE(NselRange)) {
		do.call(rbind, lapply(1:nrow(chrs), function(i) {
			subset(d, chr == chrs$chr[i] & pos >= chrs$pos[i] - NselRange & pos <= chrs$pos[i] + NselRange)
		}))
	} else merge(d, Df(chr = chrs$chr));
	Ina = which(is.na(dSel$P));
	if (length(Ina) > 0) dSel = dSel[-Ina, , drop = F];
	return(dSel);
}

# <!> names not homogenized
reduceData = function(ds, Nchrs = 3, fCols = P ~ chr + pos, P = 5e-4, Nregion = 5e4, NselRange = NULL) {
	chrs = selectTopChrs(ds[[1]], Nchrs, fCols, P, Nregion);
	dsSel = lapply(ds, function(d)reduceDataSingle(d, Nchrs, fCols, P, Nregion, NselRange, chrs));
	return(dsSel);
}

#
#	<p> parallel simulations null
#
dataColsDefault = c('chr', 'pos', 'marker', 'P');
simulateModelsWBNullScenarioSingle = function(Nstudies, pars, rotation, studyIndex, reduce) {
	LogS(2,
		'Nstudies: %{Nstudies}d, rotation: %{rot}s, studyIndex %{studyIndex}s',
		rot = join(rotation, ';')
	);
	path = parallelize_lookup(studyIndex);
	d = readTable(path, headerMap = list(posPhy = 'pos'))[dataColsDefault];
	d0 = do.call(reduceDataSingle, c(list(d), reduce));
	r = runSimulationSingle(d0, Nstudies, rotation, pars);
	return(r);
}
simulateModelsWBNullScenario = function(Nrep, Nstudies, pars, rotation, studyIndex, reduce) {
	r = Lapply(Seq(1, Nrep), function(i)
		simulateModelsWBNullScenarioSingle(Nstudies, pars, rotation, studyIndex, reduce)
	);
	return(r);
}


simulateModelsWBNull = function(modelList, parallel = F) {
	r = iterateModels(modelList, simulateModelsWBNullScenario, parallel = parallel);
	return(r);
}

#
#	<p> parallel simulations alt
#
simulateModelsWBAltScenarioSingle = function(Nstudies, pars, rotation, studyIndex, reduce,
	Ptop, Nregion, rotationAlt, Ncopy, attenuation, noise) {

	LogS(2,
		'Nstudies: %{Nstudies}d, rotation: %{rot}.2f, studyIndex %{studyIndex}d',
		rot = join('; ', rotation)
	);
	path = parallelize_lookup(studyIndex);
	d = readTable(path, headerMap = list(posPhy = 'pos'))[dataColsDefault];
	d0 = do.call(reduceDataSingle, c(list(d), reduce));
	r = runSimulationAltSingle(d0, Nstudies, rotation,
		Ptop, Nregion, rotationAlt, Ncopy, attenuation, noise, pars);
	return(r);
}
simulateModelsWBAltScenario = function(Nrep, Nstudies, pars, rotation, studyIndex, reduce,
	Ptop, Nregion, rotationAlt, Ncopy, attenuation, noise) {

	r = Lapply(Seq(1, Nrep), function(i)
		simulateModelsWBAltScenarioSingle(Nstudies, pars, rotation, studyIndex, reduce,
			Ptop, Nregion, rotationAlt, Ncopy, attenuation, noise)
	);
	return(r);
}


simulateModelsWBAlt = function(modelList, parallel = F) {
	r = iterateModels(modelList, simulateModelsWBAltScenario, parallel = parallel);
	return(r);
}

#
#	<p> summary
#
summaryMeta = function(r, PminData) {
	Pmin = list.kpu(r, '[[1]]$Pmin');
	Psupport = mean(Pmin < PminData, na.rm = T);
	return(Df(Analysis = 'meta', Psupport = Psupport));
}
summaryWB = function(r, PminData) {
	Pmin = filterList(list.kp(r, '[[2]]$[[1]]'), notE);
	PminA = aaply(laply(Pmin, as.matrix), c(2, 3), identity);
	Psupport = apply(PminA, 1:2, function(Pmin)mean(Pmin < PminData, na.rm = T));
	combs = merge.multi.list(as.list(dimnames(Psupport)));
	PsupportDf = apply(combs, 1, function(i)Psupport[i[1], i[2]]);
	combsChar = apply(combs, 1, function(i)Sprintf('E=%{E}d K=%{K}d', E = i[1], K = i[2]));
	return(Df(Analysis = combsChar, Psupport = PsupportDf));
}
summaryAll = function(r, PminData) {
	rbind(summaryMeta(r, PminData), summaryWB(r, PminData));
}
summaryPlot = function(sa, plotTitle = NULL) {
	pl = ggplot(data = sa, aes(x = Analysis, y = Psupport)) + geom_bar(stat="identity") +
		geom_hline(aes(yintercept = .5, color = 'red')) + xlab('Type of Analysis') +
		theme_bw() + theme(legend.position="none");
	if (notE(title)) pl = pl + ggtitle(plotTitle);
	return(pl);
}

studyData = function(modelList) {
	nlapply(modelList$studyIndex, function(n) {
		d = readTable(parallelize_lookup(n));
		return(d);
		#Pmin = min(d$P, na.rm = T);
		#return(list(P = Pmin));
	});
}

# ds = studyDescriptors(modelList)
summaryPlotsSave = function(r, ds,
	nameTemplate = 'results-sim/SimAlt-%{studyName}s-%{Ncopy}d-%{attenuation}d-%{noise}d.%{ext}s',
	plotType = 'pdf') {
	dTags = names(ds);
	plots = lapply(seq_along(r$results), function(i) {
		Istudy = r$models[i, 'studyIndex'];
		studyName = dTags[Istudy];
		summ = summaryAll(r$results[[i]], min(ds[[Istudy]]$P, na.rm = T));
		summL = list(summary = summ, minPdata = min(ds[[Istudy]]$P));
		pathSummary = with(as.list(r$models[i, ]), Sprintf(nameTemplate, ext = 'Rdata'));
		save(summL, file = pathSummary);
		summaryPlot(summ, plotTitle = Sprintf('Study %{studyName}s'))
	});
	ilapply(plots, function(pl, i) {
		studyName = dTags[r$models[i, 'studyIndex']];
		path = with(as.list(r$models[i, ]), Sprintf(nameTemplate, ext = plotType));
		plot_save(pl, plot_path = path)
	});
	return(plots);
}
