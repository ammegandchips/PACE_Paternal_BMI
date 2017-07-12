# Functions for the PACE paternal BMI project: Phase one

## Remove outliers using the Tukey (IQR*3) method
IQR.removal <- function(meth.matrix){
	require(matrixStats)
	rowIQR <- rowIQRs(meth.matrix, na.rm = T)
	row2575 <- rowQuantiles(meth.matrix, probs = c(0.25, 0.75), na.rm = T)
	maskL <- meth.matrix < row2575[,1] - 3 * rowIQR 
	maskU <- meth.matrix > row2575[,2] + 3 * rowIQR 
	meth.matrix[maskL] <- NA
	meth.matrix[maskU] <- NA
	meth.matrix
}

# EWAS function (including setting up the data, running sva, running the EWAS, extracting the relevant data)
ewas.function <- function(complete.pheno, meth.matrix, cell.counts, ewas.model, autosomal.sites){
	# prepare phenotype data for SVA
	print("preparing phenotype data for sva...")
	if(ewas.model=="minimal.pat" | ewas.model == "minimal.mat" | ewas.model=="cell.pat" | ewas.model == "cell.mat"| ewas.model == "full"){
	pheno <- droplevels(na.omit(complete.pheno[,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity","sex")]))
	pheno$sex <-NULL
		}else{
		if(ewas.model=="full.normal"){
		pheno <- droplevels(na.omit(complete.pheno[complete.pheno$mat.obesity==0,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity","sex")]))
		pheno$sex <-NULL
			}else{
			if(ewas.model=="full.boys"){
			pheno <- droplevels(na.omit(complete.pheno[complete.pheno$sex==0,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity")]))
				}else{
				if(ewas.model=="full.girls"){
				pheno <- droplevels(na.omit(complete.pheno[complete.pheno$sex==1,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity")]))
					}else{
					if(ewas.model=="full.normal.boys"){
					pheno <- droplevels(na.omit(complete.pheno[complete.pheno$sex==0 & complete.pheno$mat.obesity==0,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity")]))
						}else{
						if(ewas.model=="full.normal.girls"){
						pheno <- droplevels(na.omit(complete.pheno[complete.pheno$sex==1 & complete.pheno$mat.obesity==0,c("Sample_Name","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity")]))
		}}}}}}
	print("DONE!")
	# merge with cell counts
	print("merging with cell counts...")
	pheno <- na.omit(merge(pheno,cell.counts,by="Sample_Name"))
	print("DONE!!!")
	# match to methylation data
	print("matching methylation data to phenotype data...")
	meth.matrix<-meth.matrix[,match(pheno$Sample_Name,colnames(meth.matrix))]
	pheno<-pheno[match(colnames(meth.matrix),pheno$Sample_Name),]
	print("DONE!!")
	# generate surrogate variables to adjust for batch
	print("generating surrogate variables for batch...")
    require(sva)
    require(matrixStats)
    autosomal.sites = intersect(autosomal.sites, rownames(meth.matrix))
    meth.matrix.sva = meth.matrix[autosomal.sites,]
    var.idx = order(rowVars(meth.matrix.sva, na.rm=T), decreasing=T)[1:20000]
	k = which(is.na(meth.matrix.sva), arr.ind=TRUE)
	meth.matrix.sva[k] = rowMedians(meth.matrix.sva, na.rm=TRUE)[k[,1]]
	mod = model.matrix(reformulate(paste0("pheno$",colnames(pheno[-1]))))
	mod0 = mod[,-grep("pat.bmi.z|mat.bmi.z",colnames(mod))]
    sva.ret = sva(meth.matrix.sva, mod=mod, mod0=mod0, n.sv=10)
    SVs = as.data.frame(sva.ret$sv)
    colnames(SVs) <-paste0("SV",1:ncol(SVs))
	pheno <- cbind(pheno,SVs)
	print("DONE!!!!")
	# prepare phenotype data for EWAS	
	print("preparing data for EWAS...")
	if(ewas.model=="minimal.pat"){
	pheno <- pheno[,which(colnames(pheno) %in% c("pat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity",colnames(SVs)))]
		}else{
		if(ewas.model=="minimal.mat"){
		pheno <- pheno[,which(colnames(pheno) %in% c("mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity", colnames(SVs)))]
			}else{
			if(ewas.model=="cell.pat"){
			pheno <- pheno[,-which(colnames(pheno) %in% c("Sample_Name","mat.bmi.z"))]
				}else{
				if(ewas.model=="cell.mat"){
				pheno <- pheno[,-which(colnames(pheno) %in% c("Sample_Name","pat.bmi.z"))]
					}else{
					pheno <- pheno[,-which(colnames(pheno) %in% c("Sample_Name"))]
					}}}}
	print("DONE!!!!!")
	# run EWAS
	print("running EWAS...")
    require(limma)
	des = model.matrix(reformulate(paste0("pheno$",colnames(pheno))))
	fit = lmFit(meth.matrix, des)
	fit.ebayes = eBayes(fit)
	n = rowSums(!is.na(meth.matrix))
	se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,2])
    margin.error = (se * qt(0.975, df=fit.ebayes$df.total))
	res = data.frame(n=n,
	        coef=fit.ebayes$coefficient[,2],
	        ci.low.pat=fit.ebayes$coefficient[,2] - margin.error,
	        ci.high=fit.ebayes$coefficient[,2] + margin.error,
	        se=se,
			p=fit.ebayes$p.value[,2],
			fdr=p.adjust(fit.ebayes$p.value[,2],"fdr"),
			bonf=p.adjust(fit.ebayes$p.value[,2],"bonferroni"))
	if(grepl("full",ewas.model)){
		se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,3])
   		margin.error = (se * qt(0.975, df=fit.ebayes$df.total))
			res.mat = data.frame(n=n,
	        	coef=fit.ebayes$coefficient[,3],
	        	ci.low=fit.ebayes$coefficient[,3] - margin.error,
	        	ci.high=fit.ebayes$coefficient[,3] + margin.error,
	       		se=se,
				p=fit.ebayes$p.value[,3],
				fdr=p.adjust(fit.ebayes$p.value[,3],"fdr"),
				bonf=p.adjust(fit.ebayes$p.value[,3],"bonferroni"))
		res = cbind(res,res.mat)
		colnames(res) <- c(paste(colnames(res)[1:8],"pat",sep="."), paste(colnames(res)[9:16],"mat",sep="."))
		 			}
	print("DONE!!!!!! :)")
	list(model.name = ewas.model,
		model.variables = colnames(des),
		ewas.results = res)
	}


## Summarise the EWAS variables

summarise.ewas.variables <- function (complete.pheno, meth.matrix, cell.counts){
	# merge with cell counts
	pheno<- na.omit(merge(complete.pheno,cell.counts,by="Sample_Name"))
	# match to methylation data
	meth.matrix<-meth.matrix[,match(pheno$Sample_Name,colnames(meth.matrix))]
	pheno<-pheno[match(colnames(meth.matrix),pheno$Sample_Name),]
	# prepare phenotype data
	full.sample <- droplevels(na.omit(pheno[,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity","sex",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	normal.weight.mums <- droplevels(na.omit(pheno[pheno$mat.obesity==0,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity","sex",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	boys <- droplevels(na.omit(pheno[pheno$sex==0,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	girls <- droplevels(na.omit(pheno[pheno$sex==1,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	normal.weight.mums.of.boys <- droplevels(na.omit(pheno[pheno$sex==0 & pheno$mat.obesity==0,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	normal.weight.mums.of.girls <- droplevels(na.omit(pheno[pheno$sex==1 & pheno$mat.obesity==0,c("pat.bmi","mat.bmi","pat.bmi.z","mat.bmi.z","pat.age","pat.smoke","pat.ses","mat.age","mat.smoke","mat.ses","parity",colnames(cell.counts)[-which(colnames(cell.counts)=="Sample_Name")])]))
	# summarise data in "table 1" format
	require(tableone)
	full.sample.t1 <- CreateTableOne(data=full.sample,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity","sex"))
	normal.weight.mums.t1 <- CreateTableOne(data=normal.weight.mums,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity","sex"))
	boys.t1 <- CreateTableOne(data=boys,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity"))
	girls.t1 <- CreateTableOne(data=girls,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity"))
	normal.weight.mums.of.boys.t1 <- CreateTableOne(data=normal.weight.mums.of.boys,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity"))
	normal.weight.mums.of.girls.t1 <- CreateTableOne(data=normal.weight.mums.of.girls,factorVars=c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity"))
	# correlation between maternal and paternal BMI (also used to calculate adjustments for non-paternity up to 30% (DE1) as described in Lawlor 2008, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2265763/)
	cor.res <- cor.test(full.sample$pat.bmi, full.sample$mat.bmi, method="spearman")
	cor.res <- c(cor.res$estimate,  cor.res$p.value)
	pat.bmi.var <- var(full.sample$pat.bmi.z)
	mat.bmi.var <- var(full.sample$pat.bmi.z)
	pat.mat.bmi.covar <- cov(full.sample$pat.bmi.z, full.sample$mat.bmi.z, use="complete.obs")
	D <-sapply((1:30/100),function(P)((((1-P)*pat.bmi.var) / pat.mat.bmi.covar) * (((1-P)*pat.mat.bmi.covar) / mat.bmi.var )) * (((pat.bmi.var/(1*pat.mat.bmi.covar)) * (pat.mat.bmi.covar/mat.bmi.var))))
	DE1 <-D^(-1)
	cor.cov.de1 <- data.frame(stat=c("pat.mat.bmi.cor.rho","pat.mat.bmi.cor.p","pat.bmi.var","mat.bmi.var","pat.mat.bmi.covar",paste("D",1:30/100),paste("DE1",1:30/100)),
	values=c(cor.res,pat.bmi.var,mat.bmi.var,pat.mat.bmi.covar,D,DE1))
	#Association between paternal/maternal BMI and continuous covariates
	require(psych)
	assoc.cont<-corr.test(full.sample[,-which(colnames(full.sample) %in% c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity","sex"))],method="spearman")
	assoc.cont <- cbind(assoc.cont$r[,c("pat.bmi","mat.bmi")],
	assoc.cont$p[,c("pat.bmi","mat.bmi")])
	colnames(assoc.cont)<-c("pat.bmi.r","mat.bmi.r","pat.bmi.p","mat.bmi.p")
	#Association between paternal/maternal BMI and categorical covariates	
	pat.t.test<-apply(full.sample[,which(colnames(full.sample) %in% c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity","sex"))], 2, function(x) t.test(full.sample$pat.bmi~x))
	pat.t.test<-as.data.frame(t(rbind(rbind(sapply(pat.t.test,"[[","estimate"),sapply(pat.t.test,"[[","conf.int")),sapply(pat.t.test,"[[","p.value"))))
	colnames(pat.t.test)<-c("mean.in.group.0", "mean.in.group.1", "ci.lower.for.diff", "ci.upper.for.diff", "p.for.diff")
	pat.t.test$measure <-"pat.bmi"
	mat.t.test<-apply(full.sample[,which(colnames(full.sample) %in% c("mat.smoke","pat.smoke","mat.ses","pat.ses","parity","sex"))], 2, function(x) t.test(full.sample$mat.bmi~x))
	mat.t.test<-as.data.frame(t(rbind(rbind(sapply(mat.t.test,"[[","estimate"),sapply(mat.t.test,"[[","conf.int")),sapply(mat.t.test,"[[","p.value"))))
	colnames(mat.t.test)<-c("mean.in.group.0", "mean.in.group.1", "ci.lower.for.diff", "ci.upper.for.diff", "p.for.diff")
	mat.t.test$measure <-"mat.bmi"
	assoc.cat <-rbind(pat.t.test,mat.t.test)
	assoc.cat$difference <-assoc.cat$mean.in.group.1 - assoc.cat$mean.in.group.0
	# output
	list (full.sample.t1=full.sample.t1, 
	normal.weight.mums.t1=normal.weight.mums.t1, 
	boys.t1 =boys.t1, 
	girls.t1 =girls.t1, 
	normal.weight.mums.of.boys.t1 =normal.weight.mums.of.boys.t1, 
	normal.weight.mums.of.girls.t1 =normal.weight.mums.of.girls.t1, 
	cor.cov.de1 =cor.cov.de1,
	assoc.cont=assoc.cont, 
	assoc.cat=assoc.cat)
	}
	