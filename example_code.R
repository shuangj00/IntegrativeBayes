# ***README***
# The following script is used to fit micribiome count data and covariate data to the 
# integrative Bayesian zero-inflated negative binomial hierarchical mixture model 
# proposed in the manuscript

# Before running the following code, please first load micribiome count data and covariate data.
# The necessary inputs should be 
# (1) a n-by-p count matrix Y, where n is the number of samples and p is the number 
# of taxa(feature)
# (2) a n-by-R covaritae matrix X, where R is the number of covariates
# (3) a n-dimensional vector z, indicating group allocation for n samples

# ========================================================================================
# ========================================================================================
# load functions & data matrices
# ========================================================================================
# ========================================================================================
Rcpp::sourceCpp('ZINBwCOV.cpp');
source('functions.R');
load("Example_data.Rdata");


# ========================================================================================
# ========================================================================================
# preprocessing
# ========================================================================================
# ========================================================================================
# keep the features that have at least 2 observed counts for both groups:
Y.input = Y.filter(Y.mat, zvec = z.vec, min.number = 2)[[2]]
# estimate the size factor s from the count matrix Y:
s.input = sizefactor.estimator(Y.mat)

# ========================================================================================
# ========================================================================================
# get true label for later visualization:
# ========================================================================================
# ========================================================================================
feature.remain = which(Y.filter(Y.mat, zvec = z.vec, min.number = 2)[[1]] == 1)
gamma.vec = gamma.vec[feature.remain]
delta.mat = delta.mat[,feature.remain]


# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm
# ========================================================================================
# ========================================================================================
S.iter = 10000
burn.in = 0.5
res = ZINBwCOV_main(Y_mat = Y.input,
                    z_vec = z.vec, 
                    s_vec = s.input,
                    X_mat = X.mat,
                    S = S.iter, burn_rate = burn.in)
# The MCMC outputs are stored in res
# $`mu0 est`: posterior mean(after burn-in) for the vector mu(0j)
# $`phi est`: posterior mean(after burn-in) for the dispersion parameter vector
# $`beta est`: posterior mean(after burn-in) for the Beta matrix
# $`gamma PPI`: PPI for all gamma(j) after burn-in
# $`delta PPI`: PPI for all delta(rj) fter burn-in
# $`R PPI`: PPI for all r(ij) after burn-in
# $`gamma sum`: sum of all gamma(j) for all iterations
# $`mukj full`: MCMC draws for mu(kj) after burn-in
# $`mu0 full`: MCMC draws for mu(0j) after burn-in
# $`beta full`: MCMC draws for beta(rj) after burn-in




# ========================================================================================
# ========================================================================================
# Visualiziting the results for two variable selection processes
# ========================================================================================
# ========================================================================================
## (1) Variable selection for discriminating features:
gamma_VS(res$`gamma PPI`, gamma.true = gamma.vec)
par(mar=c(5.1, 4.1, 4.1, 2.1))

## (2) Variable selection for significant feature-covariate association:
delta_ROC(as.vector(res$`delta PPI`), as.vector(abs(delta.mat)))
