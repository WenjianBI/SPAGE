#' Get parameters and residuals from the NULL model
#'
#' Fit a covariates-only model to compute model parameters and residuals for SPAGE.
#' @param formula an object of class “formula”: a symbolic description of the Covariates-Only model to be fitted.
#' @param subjectID a character vector to specify the subject IDs. When calculating SPAGE p-value, the subjects IDs will be used to match genotype and environmental exposure information.
#' @param data an optional data frame containing the variables in the model (default=NULL). If it is NULL, the variables are taken from 'environment(formula)'
#' @param out_type an indicator of the outcome type. "D" for the dichotomous outcome.
#' @return an R object to pass to function SPAGE for GxE analyses.
#' @examples
#' Data.ls = data.simu.null(N = 1000, nSNP = 10, nCov = 2, maf = 0.3, prev = 0.1)
#' Phen.mtx = Data.ls$Phen.mtx
#' obj.null = SPAGE_Null_Model(y~Cov1+Cov2, data=Phen.mtx, out_type="D")
#' @export
SPAGE_Null_Model = function(formula,
                            subjectID,
                            data=NULL,
                            out_type="D")
{
  check.outType(out_type)
  if(out_type=="D")
    res.glm = glm(formula, data, family=binomial)

  mu = res.glm$fitted.values  # the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.
  coef = res.glm$coefficients
  W = (1-mu)*mu
  X = model.matrix(res.glm)   # Design matrix
  WX = W*X
  XWX_inv = solve(t(X)%*%WX)           # inverse matrix of (X %*% W %*% X)

  XW = t(WX)
  y = res.glm$y
  XXWX_inv = X %*% XWX_inv

  obj.null = list(mu=mu, W=W, coef=coef, y=y, XW=XW, XXWX_inv = XXWX_inv, y.mu=y-mu, X=X, subjectID=subjectID)
  class(obj.null) = "SPAGE_Null_Model"

  return(obj.null)
}

#### Check OutType
check.outType = function(out_type)
{
  # if (out_type != "C" && out_type != "D") {
  #   stop("Invalid out_type!. Please use either \"C\" for the continous outcome or \"D\" for the dichotomous outcome.")
  # }
  if (out_type != "D") {
    stop("Invalid out_type!. Please use \"D\" for the dichotomous outcome.")
  }
}


