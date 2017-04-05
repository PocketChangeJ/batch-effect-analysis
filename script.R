#The SVA package for removing batch effects and other unwanted variation in high-throughput
#experiments

#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
#browseVignettes("sva")

library(sva)
library(bladderbatch)
#The SVA package for removing batch effects and other unwanted variation in high-throughput
#experiments
library(pamr)
library(limma)

data(bladderdata)

#For the bladder cancer study, the variable of interest is cancer status. To begin we will assume no
#adjustment variables. The bladder data are stored in an expression set - a Bioconductor object used for
#storing gene expression data. The variables are stored in the phenotype data slot and can be obtained
#as follows:
 pheno = pData(bladderEset)
#The expression data can be obtained from the expression slot of the expression set.
 edata = exprs(bladderEset)
 
df=data.frame(pheno, t(edata))
attach(df)
normal_samples=sample[cancer=="Normal"]
cancer_samples=sample[-(cancer=="Normal")]
# boxplot(X1007_s_at[cancer=="Normal"] ~ batch[cancer=="Normal"])
detach(df)


# reshape to recreate the paper plots
m=dim(df)[2]
# df_long_normal_only=reshape(df[normal_samples,],idvar="sample",varying = list(5:m),direction = "long")
# names(df_long_normal_only)[6]="expr_value"

# save(df_long_normal_only, file="C:/Users/sadeh/Google Drive/2016 Spring/STAT646/batch effect/df_long_normal_only.RData")

load(file="C:/Users/sadeh/Google Drive/2016 Spring/STAT646/batch effect/df_long_normal_only.RData")


##############################
# since df_long is only normal samples of df transformed, no need to specify anymore:
# boxplot(df_long[normal_samples,5]~ sample[normal_samples])
par(mfrow=c(1,2))
boxplot(expr_value ~ sample , data=df_long_normal_only, 
        subset= batch==2,  
        col="green", xlab="normal samples: batch 2", ylab="expressions")

boxplot(expr_value ~ sample , data=df_long_normal_only, 
        subset= batch==3,
        col="red", xlab="normal samples: batch 3", ylab="expressions")


attach(df)
df[(batch==2) & (cancer=="Normal"),]$sample
# [1] 2 3 7 8
# 71020 71021 25 26

detach(df)

#spaghetti plot
par(mfrow=c(1,1))
ind=c(2,3,7,8,1,4,5,6)
ylim=range(df[,-(1:4)])
plot(1:8,df[ind,5+0], ylim=ylim, main="expression of first 100 genes among normal samples",
     xlab="1:4 batch 2, 5:8 batch 3 samples")
lines(1:8,df[ind,5+0], col=2)
for (i in 1:100){
  lines(1:8,df[ind,5+i], col=i)
}

ind=c(2,3,7,8,1,4,5,6)
ylim=range(df[,-(1:4)])
plot(1:8,df[ind,5+0], ylim=ylim, main="expression of first 10 genes among normal samples",
     xlab="1:4 batch 2, 5:8 batch 3 samples")
lines(1:8,df[ind,5+0], col=2)
for (i in 1:10){
  lines(1:8,df[ind,5+i], col=i)
}

# with complete 3 clusters
data=df[normal_samples,-(1:4)]

HC = hclust(dist(data))
par(mfrow=c(1,1))
plot(HC, main="complete linkage, normal samples")
# 20 21 25 26 class 1: batch 2
  
#with average 2 clusters: exactly 2 batches
HC = hclust(dist(data), method = "average")
par(mfrow=c(1,1))
plot(HC, main="average linkage, normal samples")

# #among features
# pc_of_features=prcomp(t(df[,-c(1:4)]))
# summary(pc_of_features)
# 
# #among samples
# pc=prcomp(data)
# summary(pc)

#Next we create the full model matrix - including both the adjustment variables and the variable of
#interest (cancer status). In this case we only have the variable of interest. Since cancer status has
#multiple levels, we treat it as a factor variable.
 mod = model.matrix(~as.factor(cancer), data=pheno)
#The null model contains only the adjustment variables. Since we are not adjusting for any other variables
#in this analysis, only an intercept is included in the model.
 mod0 = model.matrix(~1,data=pheno)
#Now that the model matrices have been created, we can apply the sva function to estimate batch and
#other artifacts.
 
#4. Applying the sva function to estimate batch and other artifacts 
#  The sva function performs two different steps. First it identities the number of latent factors that need
#  to be estimated. If the sva function is called without the n.sv argument specified, the number of
#  factors will be estimated for you. The number of factors can also be estimated using the num.sv.
 n.sv = num.sv(edata,mod,method="leek")
 n.sv
 
# Next we apply the sva function to estimate the surrogate variables:
 svobj = sva(edata,mod,mod0, n.sv=n.sv)
 
#  The sva function returns a list with four components, sv, pprob.gam, pprob.b, n.sv. 
#  sv is a matrix whose columns correspond to the estimated surrogate variables. They can be used in downstream
#  analyses as described below. 
#  pprob.gam is the posterior probability that each gene is associated with one or more latent variables [?]. 
#  pprob.b is the posterior probability that each gene is associated with the variables of interest [?]. 
#  n.sv is the number of surrogate variables estimated by the sva.

# 5.  Adjusting for surrogate variables using the f.pvalue function
# The f.pvalue function can be used to calculate parametric F-test p-values for each row of a data
# matrix. In the case of the bladder study, this would correspond to calculating a parametric F-test p-
# value for each of the 22,283 rows of the matrix. The F-test compares the models mod and mod0. They
# must be nested models, so all of the variables in mod0 must appear in mod. First we can calculate the
# F-test p-values for differential expression with respect to cancer status, without adjusting for surrogate
# variables, adjust them for multiple testing, and calculate the number that are significant with a Q-value
# less than 0.05.
 pValues = f.pvalue(edata,mod,mod0)
 qValues = p.adjust(pValues,method="BH")
 n.features=dim(edata)[1]
 sum(qValues<0.05)/n.features
  head(qValues)
 
#  Note that nearly 70% of the genes are strongly differentially expressed at an FDR of less than 5%
#  between groups. This number seems artificially high, even for a strong phenotype like cancer. Now
#  we can perform the same analysis, but adjusting for surrogate variables. The first step is to include
#  the surrogate variables in both the null and full models. The reason is that we want to adjust for the
#  surrogate variables, so we treat them as adjustment variables that must be included in both models.
#  Then P-values and Q-values can be computed as before.
 
 modSv = cbind(mod,svobj$sv)
 mod0Sv = cbind(mod0,svobj$sv)
 pValuesSv = f.pvalue(edata,modSv,mod0Sv)
 qValuesSv = p.adjust(pValuesSv,method="BH")
 
 sum(qValues<0.05)
 sum(qValuesSv<0.05)
 
 sum(qValues<0.05)/n.features
 sum(qValuesSv<0.05)/n.features
 
 ####they're almost same!!! 68% vs 67.7%  
 
 # Now these are the adjusted P-values and Q-values accounting for surrogate variables.
 
 
 
 # 6. Adjusting for surrogate variables using the limma package
 
#  The limma package is one of the most commonly used packages for differential expression analysis. The
#  sva package can easily be used in conjunction with the limma package to perform adjusted differential
#  expression analysis. The first step in this process is to fit the linear model with the surrogate variables
#  included.
  fit = lmFit(edata,modSv)
#  From here, you can use the limma functions to perform the usual analyses. As an example, suppose
#  we wanted to calculate differential expression with respect to cancer. To do that we first compute the
#  contrasts between the pairs of cancer/normal terms. We do not include the surrogate variables in the
#  contrasts, since they are only being used to adjust the analysis.
  contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
  fitContrasts = contrasts.fit(fit,contrast.matrix)
#  The next step is to calculate the test statistics using the eBayes function:
  eb = eBayes(fitContrasts)
  limma_out=topTableF(eb, adjust="BH", number=1000)
   


#    7. Applying the ComBat function to adjust for known batches
#    The ComBat function adjusts for known batches using an empirical Bayesian framework [1]. In order to
#    use the function, you must have a known batch variable in your dataset.
   
     batch = pheno$batch
#    Just as with sva, we then need to create a model matrix for the adjustment variables, including the
#    variable of interest. Note that you do not include batch in creating this model matrix - it will be included
#    later in the ComBat function. In this case there are no other adjustment variables so we simply fit an
#    intercept term.
     modcombat = model.matrix(~1, data=pheno)
#    Note that adjustment variables will be treated as given to the ComBat function. This means if you
#    are trying to adjust for a categorical variable with p different levels, you will need to give ComBat p-1
#    indicator variables for this covariate. We recommend using the model.matrix function to set these
#    up. For continuous adjustment variables, just give a vector in the containing the covariate values in a
#    single column of the model matrix.
#    We now apply the ComBat function to the data, using parametric empirical Bayesian adjustments.
     combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#    This returns an expression matrix, with the same dimensions as your original dataset. This new expres-
#      sion matrix has been adjusted for batch. Significance analysis can then be performed directly on the
#    adjusted data using the model matrix and null model matrix as described before:
     
     pValuesComBat = f.pvalue(combat_edata,mod,mod0)
     qValuesComBat = p.adjust(pValuesComBat,method="BH")
     
     combat_out=qValuesComBat[order(qValuesComBat, decreasing = TRUE)]
     
     # limma and combat mutual top genes
     sum(rownames(limma_out)[1:10]%in%names(combat_out[1:10]))
     sum(rownames(limma_out)[1:50]%in%names(combat_out[1:50]))
     sum(rownames(limma_out)[1:100]%in%names(combat_out[1:100]))
     
     sum(rownames(limma_out)[1:10]%in%names(combat_out[1:100]))
     sum(names(combat_out[1:10])%in%rownames(limma_out)[1:100])
     
     sum(rownames(limma_out)[1:100]%in%names(combat_out[1:1000]))
     sum(names(combat_out[1:100])%in%rownames(limma_out)[1:1000])
     
     sum(rownames(limma_out)[1:1000]%in%names(combat_out[1:1000]))
     
     sum(qValuesComBat<0.05)
     sum(qValuesComBat<0.05)/n.features
     
     data=data.frame(t(combat_edata[,normal_samples]))
     HC = hclust(dist(data), method = "average")
     par(mfrow=c(1,1))
     plot(HC, main="combat batch adjusted data: average linkage, normal samples")
     # [1] 2 3 7 8
     # 71020 71021 25 26
     plot(HC$height)

     #complete linkage
     HC = hclust(dist(data), method = "complete")
     par(mfrow=c(1,1))
     plot(HC, main="combat batch adjusted data:complete linkage, normal samples")
     # [1] 2 3 7 8
     # 71020 71021 25 26
     plot(HC$height)
     
#    These P-values and Q-values now account for the known batch effects included in the batch variable.
#    There are two additional options for the ComBat function. By default, it performs parametric empirical
#    Bayesian adjustments. If you would like to use nonparametric empirical Bayesian adjustments, use the
#    par.prior=FALSE option (this will take longer). Additionally, use the prior.plots=TRUE option to
#    give prior plots with black as a kernel estimate of the empirical batch eect density and red as the
#    parametric estimate. For example, you might chose to use the parametric Bayesian adjustments for
#    your data, but then can check the plots to ensure that the estimates were reasonable.

 #    Finally, we have now added the mean.only=TRUE option, that only adjusts the mean of the batch eects
#    across batches (default adjusts the mean and variance). This option is recommended for cases where
#    milder batch eects are expected (so no need to adjust the variance), or in cases where the variances
#    are expected to be dierent across batches due to the biology. For example, suppose a reseracher
#    wanted to project a knock-down genomic signature to be projected into the TCGA data. In this case,
#    the knockdowns samples may be very similar to each other (low variance) whereas the signature will
#    be at varying levels in the TCGA patient data. Thus the variances may be very dierent between the       
#    two batches (signature perturbation samples vs TCGA), so only adjusting the mean of the batch eect
#    across the samples might be desired in this case.

#      8. Removing known batch effects with a linear model
#      Direct adjustment for batch eects can also be performed using the f.pvalue function. In the bladder
#      cancer example, one of the known variables is a batch variable. This variable can be included as
#      an adjustment variable in both mod and mod0. Then the f.pvalue function can be used to detect
#      dierential expression. This approach is a simplied version of ComBat.
      modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
      mod0Batch = model.matrix(~as.factor(batch),data=pheno)
      pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
      qValuesBatch = p.adjust(pValuesBatch,method="BH")

      sum(qValuesBatch<0.05)
      sum(qValuesBatch<0.05)/n.features
