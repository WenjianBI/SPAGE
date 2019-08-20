#' SaddlePoint Approximation implementation of GxE analysis
#'
#' Test for association between marginal GxE multiplicative interaction effect and dichotomous phenotypes.
#' @param obj.null output object of function SPAGE_Null_Model.
#' @param Envn.mtx a numeric environment matrix with each row as an individual and each column as an environmental factor.
#' Column names of environmental factors and row names of subject IDs are required.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#' Column names of genetic variations and row names of subject IDs are required. Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#' If the test statistic lies within the standard deviation cutoff of the mean, its p value is calculated based on normal distribution approximation, otherwise, its p value is calculated based on saddlepoint approximation.
#' @param impute.method a character string (default= "none") to specify the method to impute missing genotypes.
#' "bestguess" imputes missing genotypes as most likely values (0,1,2), "random" imputes missing genotypes by generating binomial(2,p) random variables (p is the MAF), and "fixed" imputes missing genotypes by assigning the mean genotype values (2p).
#' @param missing.cutoff a numeric value (default=0.15) to specify the cutoff of the missing rates of SNPs.
#' Any SNP with missing rates higher than the cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default=0) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param Firth.cutoff a numeric value (default=0, no Firth output) to specify the p-value cutoff for Firth's test.
#' Only when the SPA p-value less than the cutoff, Firth's test p-value is calculated.
#' @param BetaG.cutoff a numeric value (default=0.001) to specify the p-value cutoff for betaG estimation. See details for more information.
#' @param BetaG.SPA a logical value (default=F) to determine p.value.BetaG is calculated based on SPA (TRUE) or a normal distribution approximation (FALSE).
#' @param G.Model a character string (default="Add") to determine the genetic model. Options include "Add" (default, no change), "Dom" (g>=1: 1; g<1: 0) and "Rec" (g>1: 1; g<=1: 0). Be careful when dosage genotype data is used. We do not check MAF before transformation.
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rate}
#' \item{p.value.BetaG}{p value of the marginal genotypic effect based on normal distribution approximation (BetaG.SPA=F) or saddlepoint approximation (BetaG.SPA=T)}
#' \item{p.value.spa-xx}{p value of the marginal GxE effect based on saddlepoint approximation. xx is the name of the environmental factor}
#' \item{p.value.norm-xx}{p value of the marginal GxE effect based on the normal distribution approximation. xx is the name of the environmental factor}
#' \item{p.value.Firth-xx}{p value of the marginal GxE effect based on the Firth's test. xx is the name of the environmental factor. xx is the name of the environmental factor}
#' \item{Stat-xx}{test statistic of the marginal GxE effect. xx is the name of the environmental factor}
#' \item{Var-xx}{estimated variance of the marginal GxE effect. xx is the name of the environmental factor}
#' \item{z-xx}{z score of the marginal GxE effect. xx is the name of the environmental factor}
#' @details
#' Here we propose a scalable and accurate method, SPAGE (SaddlePoint Approximation implementation of G×E analysis), that is applicable for genome-wide scale phenome-wide G×E studies (PheWIS).
#' SPAGE fits a genotype-independent logistic model only once across the whole-genome analysis to reduce computation cost and uses a saddlepoint approximation (SPA) to calibrate the test statistics for analysis of phenotypes with unbalanced case-control ratios.
#' When genotypic effect is small or moderate (true for most of the variants), the method can control type I error rates well.
#' We first test for the marginal genotypic effect (normal approximation if Beta.SPA=F and SPA if Beta.SPA=T) and if the p value is less than the pre-given argument 'BetaG.cutoff', we will update the test statistic and p value.
#' @examples
#' # Specify all arguments
#' N = 1000
#' Data.ls = data.simu.null(N = N, nSNP = 10, nCov = 2, maf = 0.3, prev = 0.1)
#' subjectID = paste0("ID",1:N)
#' Phen.mtx = Data.ls$Phen.mtx
#' obj.null = SPAGE_Null_Model(y ~ Cov1 + Cov2, subjectID = subjectID, data = Phen.mtx, out_type = "D")
#' Envn.mtx = as.matrix(Phen.mtx)[,"Cov1",drop=FALSE]
#' Geno.mtx = Data.ls$Geno.mtx
#' rownames(Geno.mtx) = rownames(Envn.mtx) = subjectID
#' SPAGE(obj.null, Envn.mtx, Geno.mtx)
#' @export
#' @import SPAtest
SPAGE = function(obj.null,
                 Envn.mtx,
                 Geno.mtx,
                 Cutoff = 2,
                 impute.method = "none",
                 missing.cutoff = 0.15,
                 min.maf = 0,
                 Firth.cutoff = 0,
                 BetaG.cutoff = 0.001,
                 BetaG.SPA = F,
                 G.Model = "Add")
{
  ### check input arguments
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  Cutoff=Cutoff,
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  Firth.cutoff=Firth.cutoff,
                  BetaG.cutoff=BetaG.cutoff)

  check.input(obj.null, Envn.mtx, Geno.mtx, par.list)
  print("Warnings: please make sure subjects in Covariates, Genotype and Environmental factor dataset are of the same order.")
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output dataframe
  n.Envn = ncol(Envn.mtx)
  n.Geno = ncol(Geno.mtx)
  names.Envn = colnames(Envn.mtx)

  ### Start analysis
  print("Start Analyzing...")

  # Cycle for genotype matrix
  idx.total = idx.in.set = 1
  idx.in.set = 1;

  header = make.header(names.Envn)

  output.per.set = c()

  for(i in colnames(Geno.mtx)){

    g = Geno.mtx[,i]
    output = SPAGE.one.SNP(g,
                           obj.null,
                           Envn.mtx,
                           Cutoff,
                           impute.method,
                           missing.cutoff,
                           min.maf,
                           Firth.cutoff,
                           BetaG.cutoff,
                           BetaG.SPA,
                           G.Model)

    output.per.set = rbind(output.per.set, output)
  }

  colnames(output.per.set)=header;
  rownames(output.per.set)=colnames(Geno.mtx)
  return(output.per.set)
}

#' SaddlePoint Approximation implementation of GxE analysis (One-SNP-version)
#'
#' One-SNP-version SPAGE function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGE. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGE.
#' @examples
#' Data.ls = data.simu.null(N = 1000, nSNP = 10, nCov = 2, maf = 0.3, prev = 0.1)
#' Phen.mtx = Data.ls$Phen.mtx
#' obj.null = SPAGE_Null_Model(y~Cov1+Cov2, data=Phen.mtx, out_type="D")
#' Envn.mtx = Data.ls$Phen.mtx[,"Cov1",drop=FALSE]
#' ## Check help(SPAGE) to better understand the output?
#' SPAGE.one.SNP(Data.ls$Geno.mtx[,"rs1"], obj.null, Envn.mtx)
#' @export
#' @import SPAtest
SPAGE.one.SNP = function(g,
                         obj.null,
                         Envn.mtx,
                         Cutoff = 2,
                         impute.method = "none",
                         missing.cutoff = 0.15,
                         min.maf = 0,
                         Firth.cutoff = 0,
                         BetaG.cutoff = 0.001,
                         BetaG.SPA = F,
                         G.Model = "Add")
{
  if(G.Model=="Add"){}   # do nothing if G.Model is "Add"
  else{
    if(G.Model=="Dom") g=ifelse(g>=1,1,0)
    if(G.Model=="Rec") g=ifelse(g<=1,0,1)
    if(G.Model!="Dom" & G.Model!="Rec") stop("The argument of G.Model should be 'Add', 'Dom' or 'Rec'")
  }

  pval.cutoff.spa = (1-pnorm(Cutoff))*2  # transform standard deviation cutoff to p-value cutoff

  ### Conduct genotype imputation and calcuate summary statistics for genotype
  ug = impute.geno(g, impute.method)
  MAF = ug$MAF
  missing.rate = ug$missing.rate
  g = ug$g
  pos1 = which(g!=0)
  if(missing.rate > missing.cutoff | MAF < min.maf | length(pos1)<=1){

    output = c(MAF, missing.rate, NA)
    for(j in 1:ncol(Envn.mtx)){
      output = c(output, rep(NA,6))
    }
  }else{

    tilde.g = with(obj.null, c(g - XXWX_inv %*% (XW[,pos1] %*% g[pos1])))
    z.G = with(obj.null, sum(g[pos1] * y.mu[pos1]))
    var.G = with(obj.null, sum(tilde.g * W * tilde.g))
    stat.G = z.G^2/var.G
    pval.G = pchisq(stat.G, 1, lower.tail = F)
    if(pval.G < pval.cutoff.spa & BetaG.SPA){
      NAset = which(g==0)
      spa.G = try(get.spa.pvalue(tilde.g,NAset,obj.null$mu,obj.null$y), silent = T)
      if(class(pval.G)!="try-error") pval.G = spa.G$pval.spa
    }
    output = c(MAF, missing.rate, pval.G)

    if(pval.G > BetaG.cutoff){
      ### if not significant: we let beta.g = 0
      flag.mu.update = F
      mu = obj.null$mu
      W = obj.null$W
      mn.G = mn.GE = 0
    }else{
      flag.mu.update = T
      beta.g = try(score.beta.g(z.G, var.G, g, obj.null), silent = T)
      if(class(beta.g)=="try-error") beta.g = 0
      mu = update.mu(obj.null, beta.g, g)
      W = mu * (1-mu)
      var.G = sum(tilde.g * W * tilde.g)
      mn.G = sum(tilde.g * (mu-obj.null$mu))
    }

    # Cycle for environmental factor matrix
    for(j in colnames(Envn.mtx)){
      ### GxE term
      ge = rep(0,length(g))
      ge[pos1] = g[pos1] * Envn.mtx[pos1,j]
      tilde.ge = with(obj.null, c(ge - XXWX_inv %*% (XW[,pos1] %*% ge[pos1])))  # N-dimensional vector
      z.GE = with(obj.null, sum(ge[pos1] * y.mu[pos1]))

      if(flag.mu.update) mn.GE = sum(tilde.ge * (mu-obj.null$mu))
      cov.G.GE = sum(tilde.g * W * tilde.ge)
      var.GE = sum(tilde.ge * W * tilde.ge)
      var = var.GE-cov.G.GE^2/var.G

      z = z.GE - mn.GE - (cov.G.GE/var.G) * (z.G - mn.G)
      stat = z^2/var
      pval.norm = pchisq(stat, 1, lower.tail = F)

      ### Check normal approximation p-value to see if we conduct SPA test
      if(var==0){
        output = c(output, c(NA, NA, NA, z, var, stat))
        next
      }
      if(pval.norm < pval.cutoff.spa){
        G1 = tilde.ge-(cov.G.GE/var.G)*tilde.g
        NAset = which(g==0)
        spa.res = try(get.spa.pvalue(G1,NAset,mu,obj.null$y),silent = T)        # a character of "SPA" or "fastSPA"
        if(class(spa.res)=="try-error")
          pval.spa=NA
        else
          pval.spa = spa.res$pval.spa
        # is.converge.spa = spa.res$is.converge.spa
      }else{
        pval.spa = pval.norm
      }

      ### Firth test
      if(Firth.cutoff==0) pval.firth=NA
      else{
        if (pval.spa < Firth.cutoff){
          X.mtx = as.matrix(cbind(ge, g, obj.null$X))
          re.firth = SPAtest:::fast.logistf.fit(x = X.mtx,
                                                y = obj.null$y, firth = TRUE)
          beta = re.firth$beta[1]
          var.beta = re.firth$var[1, 1]
          stat.beta = beta^2/var.beta
          pval.firth = pchisq(stat.beta, 1, lower.tail = F)
        }else{
          pval.firth = pval.spa
        }
      }

      ### summarize all results
      output = c(output, c(pval.spa, pval.norm, pval.firth, z, var, stat))
    }
  }
  return(output)
}


check.input = function(obj.null,
                       Envn.mtx,
                       Geno.mtx,
                       par.list)
{
  if(class(obj.null)!="SPAGE_Null_Model") stop("Argument 'obj.null' should be output of function SPAGE_Null_Model.")
  if(!is.numeric(Envn.mtx)|!is.matrix(Envn.mtx)) stop("Input 'Envn.mtx' should be a numeric matrix.")
  if(is.null(colnames(Envn.mtx))) stop("Column names of 'Envn.mtx' should be given.")
  if(is.null(rownames(Envn.mtx))) stop("Row names of 'Envn.mtx' should be given.")
  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(any(rownames(Envn.mtx)!=rownames(Geno.mtx))) stop("Please check the sample order of 'Envn.mtx' and 'Geno.mtx'.")
  if(any(rownames(Envn.mtx)!=obj.null$subjectID)) stop("Please check the sample order of 'Envn.mtx' and 'Phen.mtx'.")
  if(!is.null(Geno.mtx)){
    if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")
    if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  }
  if(!is.numeric(par.list$Cutoff)|par.list$Cutoff<0) stop("Argument 'Cutoff' should be a numeric value greater than or equal to 0.")
  if(!is.element(par.list$impute.method,c("none","bestguess","random","fixed"))) stop("Argument 'impute.method' should be 'none', 'bestguess', 'random' or 'fixed'.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>1) stop("Argument 'min.maf' should be a numeric value between 0 and 1.")
  if(!is.numeric(par.list$Firth.cutoff)|par.list$Firth.cutoff<0|par.list$Firth.cutoff>1) stop("Argument 'Firth.cutoff' should be a numeric value between 0 and 1.")
  if(!is.numeric(par.list$BetaG.cutoff)|par.list$BetaG.cutoff<0|par.list$BetaG.cutoff>1) stop("Argument 'BetaG.cutoff' should be a numeric value between 0 and 1.")
}


make.header = function(names.Envn)
{
  header1 = c("MAF", "missing.rate", "p.value.BetaG")               # 3
  header2 = c("p.value.spa", "p.value.norm", "p.value.Firth",       # 3
              "Stat","Var","z")                                     # 3
  for(envn in names.Envn){
    header = c(header1, paste0(header2,"-",envn))
  }
  header = matrix(header,nrow=1)
  return(header)
}

impute.geno = function(g, impute.method)
{
  if(impute.method=="none"){
    N.na=0;
  }else{
    which.na = which(is.na(g))
    N.na = length(which.na)
  }

  if(N.na==0) MAF=mean(g)/2
  else MAF = mean(g[-which.na])/2

  if(MAF > 0.5){
    g = 2-g;
    MAF = 1-MAF
  }
  missing.rate = N.na/length(g)

  if(N.na>0){
    if(impute.method=="fixed") g[which.na] = MAF
    if(impute.method=="random") g[which.na] = rbinom(N.na, 2, MAF)
    if(impute.method=="bestguess") g[which.na] = round(2*MAF)
  }

  return(list(MAF=MAF,
              missing.rate=missing.rate,
              g=g))
}



####### calculate SPA test p-values based on SPAtest package (3.0.0)

get.spa.pvalue = function(G1,         # GxE interaction term
                          NAset,
                          mu,         # Updated mu vector based on genotype
                          y)
{

  res=NULL
  output = "P"
  nodes.fixed = NULL
  Cutoff = 0.1
  alpha = 5e-8

  q <- sum(G1 * y)
  g = G1

  if (length(NAset)/length(G1) < 0.5) {
    out <- SPAtest:::Saddle_Prob(q, mu = mu, g = g, Cutoff = Cutoff,
                                 alpha = alpha, output = output, nodes.fixed = nodes.fixed,
                                 nodes.init = nodes.init)
  }
  else {
    out <- SPAtest:::Saddle_Prob_fast(q, g = g, mu = mu, gNA = g[NAset],
                                      gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset],
                                      Cutoff = Cutoff, alpha = alpha, output = output,
                                      nodes.fixed = nodes.fixed, nodes.init = nodes.init)
  }

  pval.spa = out$p.value
  is.converge.spa = out$Is.converge
  res = list(pval.spa = pval.spa,
             is.converge.spa = is.converge.spa)
  return(res)
}

####### This is from function ScoreTest_SPA_wOR() from Rounak

score.beta.g = function(z.G,
                        var.G,
                        g,             # Updated GxE interaction term
                        obj.null)      # output of function fit.null.model()
{
  g = g-mean(g)

  y = obj.null$y
  mu = obj.null$mu
  W = obj.null$W

  Score = z.G
  VarF = var.G
  ncase = sum(y)
  ncontrol = length(y)-sum(y)

  mu1 = ncase/(ncase+ncontrol)
  VarW = mu1*(1-mu1)*sum(g^2)
  Scale = sqrt(VarF/VarW)
  g1 = g * Scale

  beta.g=score_solve(y, g1, Score)
  return(beta.g)
}

score_solve<-function(y, g1, qadj)
{
  g1 = as.vector(g1)
  gamma = c(0,0)
  rep = 1
  repeat{
    exp.eta1 = exp(gamma[1]+g1*gamma[2])
    mu1 = exp.eta1/(1+exp.eta1)
    W1 = mu1*(1-mu1)
    H22<-sum(g1^2*W1)
    H12<-sum(g1*W1)
    H11<-sum(W1)
    H<-solve(matrix(c(H11,H12,H12,H22),2,2))
    XW2 <- cbind(1,g1) * (W1)^0.5
    myQR <- qr(XW2)
    Q <- qr.Q(myQR)
    h <- (Q*Q) %*% rep(1, ncol(Q))
    s<-c(sum(y-mu1+h*(0.5-mu1)),qadj-sum(g1*mu1-g1*h*(0.5-mu1)))
    gammanew<-gamma+H%*%s
    if(sum(abs(gamma-gammanew))<10^-5 || rep>100)	break
    gamma=gammanew
    rep=rep+1
  }
  return(gammanew[2])
}

update.mu = function(obj.null, beta.g, tilde.g)
{
  mu0 = obj.null$mu

  X = obj.null$X

  eta1 = with(obj.null, c(X %*% coef + tilde.g * beta.g))
  exp.eta1 = exp(eta1)
  mu1 = exp.eta1/(1+exp.eta1)

  H = solve(t(X)%*%(mu1*(1-mu1)*X))

  rep = 1
  X.t = t(X)
  repeat{

    f = X.t %*% (mu1-mu0)
    d.coef = -1 * H %*% f
    d.eta = c(X %*% d.coef)
    eta1 = eta1+d.eta
    exp.eta1 = exp(eta1)
    mu1 = exp.eta1/(1+exp.eta1)
    if(mean(abs(d.eta)) < 10^-6 || rep>100)	break

    #### uncomment the following line for original NR algorithm
    # H = solve(t(X)%*%(mu1*(1-mu1)*X))
    rep = rep+1
  }

  return(mu1)
}

update.g = function(g, G.Model){
  if(G.Model=="Dom") return()
  if(G.Model=="Rec") return(ifelse(g<=1,0,1))
}


