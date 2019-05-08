#' Data simulation (null model & multiple SNPs)
#'
#' Simulate a dataset which includes genotypes, covariates and phenotypes.
#' @param N a numeric value: number of samples
#' @param nSNP a numeric value: number of SNPs
#' @param nCov a numeric value: number of covariates
#' @param maf a numeric value or a numeric vector of length nSNP: minor allele frequencies to simulate genotypes.
#' It should be between 0 and 0.5.
#' @param prev a numeric value: the expected proportion of cases among all samples.
#' It should be between 0 and 1.
#' @param TypeOfCov a character vector of length nCOV to specify the covariates types: "binary" ~ Bernoulli(0.5), "continuous" ~ Normal(0,1).
#' Default value is NULL, that is, 'rep(c("binary","continous"), length.out=nCov)'.
#' @param CoefOfCov a numeric vector of length nCOV to specify the covariates coefficients.
#' Default value is NULL, that is, 'rep(0.5, length.out=nCov)'.
#' @return an R list including the following elements
#' \item{Phen.mtx}{an R matrix (N * nCov+1) of phenotype and covariates}
#' \item{Geno.mtx}{an R matrix (N * nSNP) of genotypes}
#' @examples
#' # Specify all arguments
#' N = 10000;
#' nSNP = 100
#' nCov = 2
#' maf = 0.05;  # OR maf = runif(nSNP, 0, 0.5)
#' prev = 0.01;
#' # Data simulation process
#' Data.ls = data.simu.null(N, nSNP, nCov, maf, prev)
#' Data.ls$Phen.mtx[1:10,]
#' Data.ls$Geno.mtx[1:10,1:10]
#' @export
data.simu.null = function(N,
                          nSNP,
                          nCov,
                          maf,
                          prev,
                          TypeOfCov=NULL,
                          CoefOfCov=NULL)
{
  ############## Simulation of Genotype matrix
  if(max(maf)>0.5 | min(maf)<0) stop("Argument 'maf' should be between 0 and 0.5.")
  if(length(maf)==1){
    G=matrix(rbinom(N*nSNP,2,maf),N,nSNP)
  }else if(length(maf)==nSNP){
    G=sapply(maf,FUN=function(x){rbinom(N,2,x)})
  }else{
    stop("Argument 'maf' is to specify minor allele frequency for genotype simulation. It should be a numeric value or a numeric vector of length n.")
  }

  colnames(G) = paste0("rs",1:nSNP)

  ############## Simulation of Covariates matrix
  if(is.null(TypeOfCov)) TypeOfCov = rep(c("binary","continuous"), length.out=nCov)
  if(is.null(CoefOfCov)) CoefOfCov = rep(0.5, length.out=nCov)
  if(length(TypeOfCov)!=nCov | length(CoefOfCov)!=nCov) stop("'TypeOfCov' and 'CoefOfCov' should be vector of length 'nCov'.")

  Phen.mtx = c()
  for(i in 1:nCov){
    if(TypeOfCov[i]=="binary") Cov.tmp = rbinom(N,1,0.5)
    else if (TypeOfCov[i]=="continuous") Cov.tmp = rnorm(N)
    else stop("'TypeOfCov' should be 'binary' or 'continuous'.")
    Phen.mtx = cbind(Phen.mtx, Cov.tmp)
  }
  colnames(Phen.mtx) = paste0("Cov",1:nCov)

  ############## Simulation of Phenotype
  eta = apply(t(Phen.mtx) * CoefOfCov, 2, sum)
  alpha0 = uniroot(f=GetMean, c(-100, 100), prev=prev, eta=eta)  # calculate intercept term to keep prevalence as expected
  alpha0 = alpha0$root
  eta1 = eta+alpha0
  mu = exp(eta1)/(1+exp(eta1))
  y = rbinom(N,1,mu)
  real.prev = signif(mean(y==1),4)

  ##
  message(paste0("Expected prevalence is ",prev," and the real prevalence is ",real.prev,"."))
  ##

  Phen.mtx = data.frame(y=y, Phen.mtx)
  Geno.mtx = G
  Data.ls = list(Phen.mtx = Phen.mtx,
                 Geno.mtx = Geno.mtx)
  return(Data.ls)
}


# a lower function to calculate mean of Pr(y=1) for all subjects
GetMean=function(x, prev, eta){
  eta1 = eta+ x
  mu<-exp(eta1)/(1+exp(eta1))
  return(mean(mu)-prev)
}

### data simulation in which phenotype is affected by only G and covariates
data.simu.alt = function(N,
                         maf,
                         prev,
                         OR.G,
                         silent=T)
{
  G = rbinom(N,2,maf)

  E1 = rbinom(N,1,0.5)
  E2 = rnorm(N)

  betaG = log(OR.G)
  eta = E1 + E2 + betaG*G
  alpha0=uniroot(f=GetMean, c(-100, 100), prev=prev, eta=eta)  # calculate intercept term to keep prevalence as expected
  alpha0=alpha0$root
  eta1=eta+alpha0
  mu=exp(eta1)/(1+exp(eta1))
  y=rbinom(N,1,mu)
  real.prev=signif(mean(y==1),4)
  ##
  if(!silent)
    message(paste0("Pre-specified prevalence is ",prev," and the real prevalence is ",real.prev,"."))
  ##
  Data.mtx = data.frame(y=y, x1=E1, x2=E2, g=G)
  return(Data.mtx)
}


### data simulation in which phenotype is affected by only G and covariates
data.simu.alt.GE = function(N,
                           maf,
                           prev,
                           OR.G,
                           OR.GE,
                           silent=T)
{
  G = rbinom(N,2,maf)

  E1 = rbinom(N,1,0.5)
  E2 = rnorm(N)

  e1=rnorm(N)

  betaG = log(OR.G)
  betaGE = log(OR.GE)

  eta = E1 + E2 + betaG*G + betaGE*G*e1
  alpha0=uniroot(f=GetMean, c(-100, 100), prev=prev, eta=eta)  # calculate intercept term to keep prevalence as expected
  alpha0=alpha0$root
  eta1=eta+alpha0
  mu=exp(eta1)/(1+exp(eta1))
  y=rbinom(N,1,mu)
  real.prev=signif(mean(y==1),4)

  ##
  if(!silent)
    message(paste0("Pre-specified prevalence is ",prev," and the real prevalence is ",real.prev,"."))
  ##
  Data.mtx = data.frame(y=y, x1=E1, x2=E2, g=G, e1=e1, mu=mu)
  return(Data.mtx)
}

