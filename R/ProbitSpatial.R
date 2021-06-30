# library(Rcpp)
# require(RcppEigen)
# require(Matrix)
# require(RANN)
# library(numDeriv)
# library(devtools)

#' Class of Spatial Probit Model.
#'
#' @slot beta numeric, the estimated parameters for the covariates.
#' @slot rho numeric, the estimated spatial autocorrelation parameter.
#' @slot lambda numeric, the estimated spatial error autocorrelation parameter.
#' @slot coeff numeric, all estimated parameters.
#' @slot loglik  numeric, the likelihood associated to the estimated model.
#' @slot formula \code{formula}.
#' @slot nobs numeric, number of observations.
#' @slot nvar numeric, number of covariates.
#' @slot y numeric, vector of observed dependent variable.
#' @slot X matrix, matrix of covariates.
#' @slot time numeric, estimation time.
#' @slot DGP character, DGP of the model (SAR, SEM or SARAR).
#' @slot method character, estimation method ("\code{conditional}" or 
#' 		"\code{full-lik}").
#' @slot varcov character, indicates the matrix used in the algorithm 
#' 		("\code{varcov}" or "\code{precision}").
#' @slot W SparseMatrix, the spatial weight matrix of y.
#' @slot M SparseMatrix, the spatial weight matrix of the disturbances.
#' @slot iW_CL numeric, the order of approximation used in the conditional 
#' 		method.
#' @slot iW_FL numeric, the order of approximation used inside the likelihood 
#'    function for the \code{full-lik} method.
#' @slot iW_FG numeric, the order of approximation used inside the gradient 
#'    functions for the \code{full-lik} method.
#' @slot reltol numeric, the relative convergence tolerance.
#' @slot prune numeric, the pruning for the gradient functions.
#' @slot env an \code{environment} containing information for use in later 
#'    function calls to save time.
#' @slot message a integer giving any additional information or NULL.
#'	
#' @export
setClass("ProbitSpatial", 
slots= list(beta		="numeric",
  rho		="numeric",
  lambda	="numeric",
  coeff		="numeric",
  loglik	="numeric",
  formula	="formula",
  nobs		="integer",
  nvar		="integer",
  y			="numeric",
  X			="matrix",
  time		="numeric",
  DGP		="character",
  method  	="character",
  varcov  	="character",
  W			="Matrix",
  M			="ANY",
  iW_CL		="numeric",
  iW_FL		="numeric",
  iW_FG		="numeric",
  reltol	="numeric",
  prune		="numeric",
  env		="environment",
  message	="numeric")         
)

#' Estimated coefficients of a spatial probit model.
#' 
#' Returns the coefficients estimated by a \code{ProbitSpatial} model.
#'
#' @param object an object of class \code{ProbitSpatial}.
#' @param ... ignored
#'
#' @return  It returns the value of the estimated parameters.
#' 
#' @method coef ProbitSpatial
coef.ProbitSpatial <- function(object,...) 
    {
	return(object@coeff)
    }

#' Effects of a spatial probit model.
#' 
#' Returns the marginal effects of a \code{ProbitSpatial} model.
#'
#' @usage effects_ProbitSpatial(object)
#' @param object an object of class \code{ProbitSpatial}.
#'
#' @return  It returns the marginal effects of the estimated 
#' 	\code{ProbitSpatial} model.
#'
#' @details The \code{effects} function has different outputs according to the 
#' 	DGP of the \code{ProbitSpatial} model:
#' \describe{
#'   \item{\code{"SAR"}}{The marginal effects of a spatial autoregressive 
#' 	model are more complicated than usual measurements of impacts for non 
#' 	spatial models. Here we follow LeSage and Pace and propose the following 
#' 	summaries for impact measures:}
#' 		\describe{
#'		\item{Average direct effects:}{the average over all the observations of 
#'       the effects of the change of an explanatory variable of a single 
#'       observation on the choice probability of that same observation.}
#'		\item{Average indirect effects:}{the average over all the observations 
#'       of the effect of a change on a explanatory variable on the choice 
#'       probability of the neighbouring observations.}
#'		\item{Average total effects:}{the sum of direct and indirect impacts.}
#' }
#'   \item{\code{"SEM"}}{marginal effects should be interpreted as if it were a 
#'		 standard probit model.}
#' }
#' @references
#' 	J. LeSage and R.K. Pace. \emph{Introduction to Spatial Econometrics}, CRC 
#' 	Press, chapter 10.1.6, 2009.
#'
#' @export 
effects_ProbitSpatial <- function(object) 
     {
idx1 <- match(c("(Intercept)","Intercept","(intercept)","intercept"),colnames(object@X))
if (any(!is.na(idx1))) {
	nvar <- object@nvar-1
    X <- object@X[,-idx1[!is.na(idx1)]]
    beta <- object@beta[-idx1[!is.na(idx1)]]
} else { 
	nvar <- object@nvar
    X <- object@X
    beta <- object@beta
    }
if (object@DGP == "SEM"){
    ME <- rep(0,nvar)
    dfit <- dnorm(fitted(object,type="link"))
    	for (i in 1:nvar){ME[i] <- mean(beta[i]*dfit)}
	ME <- data.frame(ME,row.names=names(beta))
	names(ME) <- c("average marginal effetcs")
} else {
    ME <- matrix(0,nvar,3)
    dfit <- dnorm(fitted(object,type="link"))
    D <- Matrix::diag(dfit)
    iW <- ApproxiW(object@W,object@rho,object@iW_CL)
    P <- D %*% iW 
    for (i in 1:nvar){
    	dPdi <- P %*% Matrix::diag(beta[i],object@nobs)
        #ME[i,3] <- mean(Matrix::rowSums(dPdi)) #Old version
		ME[i,3] <- (1-object@rho)^-1*mean(dfit)*beta[i] #New version 07/08/17
        ME[i,1] <- mean(Matrix::diag(dPdi))
        ME[i,2] <- ME[i,3]-ME[i,1]
        }
    ME <- data.frame(ME,row.names=names(beta))
    names(ME) <- c("direct","indirect","total")
    }
return(ME)
}

#' Extract spatial probit model fitted values.
#' 
#' Extract the fitted values of a \code{ProbitSpatial} model.
#'
#' @param object an object of class \code{ProbitSpatial}.
#' @param type the type of output: 
#' \describe{
#'   \item{\code{"link"}}{the value of the latent variable. Default.}
#'   \item{\code{"response"}}{probability.}
#'   \item{\code{"binary"}}{binary 0/1 output.}
#' }
#' @param cut the threshold probability for the \code{"binary"} type. 
#' 	Default is 0.5. 
#' @param ... ignored
#'
#' @return Returns the vector of fitted values of the \code{ProbitSpatial} model
#'
#' @method fitted ProbitSpatial

fitted.ProbitSpatial<-function(object,type=c("link","response","binary"),cut=0.5,...) 
{
type <- match.arg(type)
if (object@DGP %in% c("SAR","SARAR")){ 
	iW <- ApproxiW(object@W,object@rho,object@iW_CL)
    f <- iW %*% (object@X %*% object@beta)
    } else {
    f <- (object@X %*% object@beta) }

if (type=="link") {f<-as.numeric(f)}
if (type=="response") {f<-pnorm(as.numeric(f))}
if (type=="binary") {f<-(pnorm(as.numeric(f)) >= cut)*1}
return(f)
}

#' Extract names of ProbitSpatial class.
#'
#' Extract names of ProbitSpatial class.
#'
#' @param x an object of class \code{ProbitSpatial}.
#' @param ... ignored
#'
#' @return Returns the names of the \code{ProbitSpatial} object. 
#' 
#' @method names ProbitSpatial

names.ProbitSpatial <-  function(x,...){return(slotNames(x))}

#' Spatial probit model predictions.
#' 
#' Predicts of a \code{ProbitSpatial} model on a set \code{X} of covariates.
#' Works on both in-sample and out-of-sample using BLUP formula from Goulard et 
#' al. (2017)
#'
#' @param object an object of class \code{ProbitSpatial}.
#' @param X a matrix of explanatory variables. If oos=TRUE, it may contain more 
#'	observations than the dataset on which the model has been trained
#' @param type the type of output: 
#' \describe{
#'   \item{\code{"link"}}{the value of the latent variable. Default}
#'   \item{\code{"response"}}{probability.}
#'   \item{\code{"binary"}}{binary 0/1 output.}
#' }
#' @param cut the threshold probability for the \code{"binary"} type. 
#' 	Default is 0.5. 
#' @param oos logical. If TRUE, out-of-sample predictions are returned.
#' @param WSO W matrix containing weights of in-sample and 
#' out-of-sample data. Observations must be ordered in such a way that
#' the first elements belong to the in-sample data and the remaining ones
#' to the out-of-sample data.
#' @param ... ignored
#' 
#' @details If \code{oos=FALSE}, the function computes the predicted values for  #' the estimated model (same as \code{fitted}). Otherwise, it applies the BLUP #' formula of Goulard et al. (2017):
#' \deqn{\hat{y} = (\hat(y_S),\hat(y_O)),}
#' where the sub-indexes S and O refer, respectively, to the in-sample and
#' out-of-sample data. \eqn{\hat{y_S}} corresponds to fitted values, while
#' \eqn{\hat{y_O}} is computed as follows:
#' \deqn{\hat{y_O} = (I-\rho W)^{-1}(X\beta)-Q_{OO}^{-1}Q_{OS}(y_S-\hat{y_S}),}
#' where \eqn{Q} is the precision matrix of 
#' \eqn{\Sigma=\sigma^2((I-\rho W)'(I-\rho W))^{-1}.} and the sub-indexes OO and 
#' OS refer to the corresponding block matrices.
#'
#' @return Returns a vector of predicted values for the set \code{X} of 
#' 	covariates if \code{oos=FALSE} or the best linear unbiased predictors of the #' set \code{XOS} if \code{oos=TRUE}.
#'
#' @references
#' \describe{
#' \item{Goulard et al. (2017)}{M. Goulard, T. Laurent and C. Thomas-Agnan. 
#'	About predictions in spatial autoregressive models: optimal and almost
#'  optimal strategies. \emph{Spatial Economic Analysis} 12, 304-325, 2017.}}
#'
#' @method predict ProbitSpatial

predict.ProbitSpatial <- function(object,X,type=c("link","response","binary"),cut=0.5,oos=FALSE,WSO=NULL,...)
{
  if(oos==FALSE){
    if (ncol(X) != object@nvar) {stop("number of columns of X does not match 								 			 			
                                      the number of variables of the model")}	
    if (nrow(X) != object@nobs) {stop("number of observations of X does not
                                      match the one expected by the model. 			 
                                      Consider out-of-sample predictions.")}
    type <- match.arg(type)
    if (object@DGP == "SAR"){
      iW <- ApproxiW(object@W,object@rho,object@iW_CL)
      f <- iW %*% (X %*% object@beta)} else { f <- (X %*% object@beta)  }
    if (type=="link") {f<-as.numeric(f)}
    if (type=="response") {f<-pnorm(as.numeric(f))}
    if (type=="binary") {f<-(pnorm(as.numeric(f))>=cut)*1}
    return(f) } else {
      type <- match.arg(type)
      if (object@DGP == "SAR"){
        n<-dim(WSO)[1]
        ns<-object@nobs
        no<-n-ns
        XSO<-X
        Q<-Matrix::crossprod(Matrix::Diagonal(n,1)-object@rho * WSO)
        #t(diag(1,n)-object@rho*WSO)%*%(diag(1,n)-object@rho*WSO)
        QOO<-Q[(ns+1):n,(ns+1):n]
        QOS<-Q[(ns+1):n,1:ns]
        iW <- ApproxiW(WSO,object@rho,object@iW_CL)
        f <- iW %*% (XSO %*% object@beta)
        re<-Matrix::solve(QOO,QOS)%*%(object@y-f[1:ns])
        f[(ns+1):n]<-f[(ns+1):n]-re
      } else {f <- (X %*% object@beta)}
      if (type=="link") {f<-as.numeric(f)}
      if (type=="response") {f<-pnorm(as.numeric(f))}
      if (type=="binary") {f<-(pnorm(as.numeric(f))>=cut)*1}
      return(f)}
    }
    
#' Extract spatial probit model residuals.
#' 
#' Compute the residuals of an estimated \code{ProbitSpatial} model.
#' 
#' @param object an object of class \code{ProbitSpatial}.
#' @param ... ignored
#' 
#' @return Return a vector containing the generalised residuals of the 
#' 	\code{ProbitSpatial} model.
#'
#' @method residuals ProbitSpatial


residuals.ProbitSpatial <- function(object,...) 
    { object@y - fitted(object,type="response") }  
 
#' Fit a spatial probit model.
#' 
#' Approximate likelihood estimation of the probit model with spatial 
#' autoregressive (SAR), spatial error (SEM), spatial autoregressive with 
#' autoregressive disturbances (SARAR).
#' 
#' @usage ProbitSpatialFit(formula,data,W,
#'          DGP='SAR',method="conditional",varcov="varcov",
#'          M=NULL,control=list())
#' 
#' @param formula an object of class \code{formula}: a symbolic 
#' 	description of the model to be fitted.
#' @param data the data set containing the variables of the model.
#' @param W  the spatial weight matrix of class \code{"dgCMatrix"}.
#' @param DGP the data generating process of \code{data}: SAR, SEM, SARAR 
#' (Default is SAR).
#' @param method the optimisation method: \code{"conditional"} or 
#' 	\code{"full-lik"} (Defaul is \code{"conditional"}, see Details).
#' @param varcov the likelihood function is computed using the 
#' 	variance-covariance matrix (\code{"varcov"}) or the precision matrix 
#' 	(\code{"precision"})? Default is \code{"varcov"}.
#' @param M the second spatial weight matrix for SARAR models. Same class as W. 
#' @param control a list of control parameters. See Details.
#' 
#' @details 
#' 	The estimation is based on the approximate value of the true likelihood of 
#' 	spatial probit models. 
#' 	The DGP of the spatial autoregressive model (SAR) model is the following
#' 	\deqn{y = \rho Wy + X\beta + \epsilon,}
#' 	where the disturbances \eqn{\epsilon} are iid standard normally distributed, 
#' 	\eqn{W} is a sparse spatial weight matrix and \eqn{\rho} is the spatial lag 
#' 	parameter. The variance of the error term is equal 
#' 	to \eqn{\Sigma=\sigma^2((I_n-\rho W)^{-1}((I_n-\rho W)^{-1})^{t})}.
#' 	The DGP of the spatial error model (SEM) is as follows
#' 	\deqn{y = X\beta+u,}
#' 	\deqn{u = \rho W u + \epsilon,}
#' 	where the disturbances \eqn{\epsilon} are iid standard normally distributed, 
#' 	\eqn{W} is a sparse spatial weight matrix and \eqn{\rho} is the spatial 
#'	error  parameter. The variance of the error term  
#'	is equal to \eqn{\Sigma=\sigma^2((I_n-\rho W)^{-1}((I_n-\rho W 
#'	)^{-1})^{t})}.
#' 	The DGP of the spatial autoregressive model with autoregressive disturbances 
#' 	(SARAR) is as follows
#' 	\deqn{y = \rho Wy + X\beta + u,}
#' 	\deqn{u = \lambda M u + \epsilon,}
#' 	where the disturbances \eqn{\epsilon} are iid standard normally distributed, 
#' 	\eqn{W} and \eqn{M} are two sparse spatial weight matrix, while \eqn{\rho} 
#' 	 and \eqn{\lambda} are the spatial lag and spatial error parameters, 
#' 	 respectively. The variance of the error term  
#'	is equal to \eqn{\Sigma=\sigma^2((I_n-\rho W)^{-1}(I_n-\lambda 
#'	M)^{-1}((I_n-\lambda M)^{-1})^{t}((I_n-\rho W)^{-1})^{t})}.
#'
#' 	The approximation is inspired by the Mendell-Elston approximation 
#' 	of the multivariante normal probabilities (see References). It makes use of 
#' 	the Cholesky decomposition of the variance-covariance matrix \eqn{\Sigma}.
#' 
#' 	The \code{ProbitSpatialFit} command estimates the model by maximising the 
#' 	approximate log-likelihood. We propose two optimisation method:
#' 	\describe{
#'  \item{\code{"conditional"}:}{ it relies on a standard probit estimation  
#' 	which applies to the model estimated 
#' 	conditional on \eqn{\rho}.}
#'  \item{\code{"full-lik"}:}{ it minimises the full-log-likelihood using the 
#' 	analytical gradient functions (only available for SAR and SEM 
#' 	specification). The optimisation is performed by means of the 
#' 	\code{\link[stats]{optim}} function with \code{method = "BFGS"}.}
#' 	}
#' 	In both cases a \code{"conditional"} estimation is performed. If 
#' 	\code{method="conditional"}, then \code{ProbitSpatialFit} returns 
#' 	the results of this first estimation. In case \code{method="full-lik"},
#' 	the function tries to improve the log-likelihood by means of a further 
#' 	exploration around the value of the parameters found by the conditional 
#'	step.
#' 	The conditional step is usually very accurate and particularly fast. The 
#' 	second step is more time consuming and does not always improve the results
#' 	of the first step. We dissuade the user from using the full-likelihood 
#'	method 	for sample sizes bigger than ten thousands, since the computation of 
#'	the gradients is quite slow.  Simulation studies reported in Martinetti and 
#' 	Geniaux (2017) prove that the conditional estimation is highly reliable,
#'	even  if compared to the full-likelihood ones.
#' 
#' 	In order to reduce the computation time of the function
#'  \code{ProbitSpatialFit}, we propose a variant of the likelihood-function 
#' 	estimation that uses the inverse of the variance-covariance matrix (a.k.a. 
#' 	precision matrix). This variant applies to both the \code{"conditional"} and 
#' 	the \code{"full-lik"} methods and can be invoked by setting 
#' 	\code{varcov="precision"}. Simulation studies reported in Martinetti and 
#' 	Geniaux (2017) suggest that the accuracy of the results with the precision 
#' 	matrix are sometimes worst than the one with the true variance-covariance 
#' 	matrix, but the estimation time is considerably reduced.
#' 
#' 	The control argument is a list that can supply any of the following 
#' 	components:
#' 	\describe{
#'   \item{\code{iW_CL}}{the order of approximation of \eqn{(I_n-\rho W)^{-1}} 
#' 	used in the \code{"conditional"} method. Default is 6, while 0 means no 
#' 	approximation (it uses exact inversion of matrixes, not suitable for big 
#' 	sample sizes). See Martinetti and Geniaux (2017) for further references.}
#'   \item{\code{iW_FL}}{the order of approximation of \eqn{(I_n-\rho W)^{-1}} 
#' 	used in the computation of the likelihood function for the \code{"full-lik"} 
#' 	method. Default is 0, meaning no approximation.}
#'   \item{\code{iW_FG}}{the order of approximation of \eqn{(I_n-\rho W)^{-1}} 
#' 	used in the computation of the gradient functions for the \code{"full-lik"} 
#' 	method. Default is 0, meaning no approximation.}
#' 	\item{\code{reltol}}{relative convergence tolerance. It represents 	
#'	\code{tol} in \code{\link[stats]{optimize}} function  for 
#'	\code{method="conditional"} and  \code{reltol} in \code{\link[stats]{optim}} 
#'	function for \code{method="full-lik"}. Default is 1e-5.}
#'   \item{\code{prune}}{the pruning value used in the gradients. Default is 0,
#' 	meaning no pruning. Typacl values are around 1e-3 and 1e-6. They help 
#' 	reducing the estimation time of the gradient functions.}
#'   \item{\code{silent}}{Default is TRUE.}
#' }
#'
#' @return Return an object of class \code{ProbitSpatial}.
#' 
#'@references
#' \describe{
#' \item{Mendell and Elston (1974)}{N. Mendell and R. Elston. Multifactorial 
#' 	qualitative traits: genetic analysis and prediction of recurrence risks. 
#' \emph{Biometrics} 30, 41--57, 1974.}
#' 
#' \item{Martinetti and Geniaux (2017)}{D. Martinetti and G. Geniaux. 
#' 	Approximate likelihood estimation of spatial probit models. \emph{Regional 
#' 	Science and Urban Economics} 64, 30-45, 2017.}}
#'
#' @examples
#' \donttest{
#' n <- 1000
#' nneigh <- 3
#' rho <- 0.5
#' beta <- c(4,-2,1)
#' W <- generate_W(n,nneigh,seed=123)
#' X <- cbind(1,rnorm(n,2,2),rnorm(n,0,1))
#' colnames(X) <- c("intercept","X1","X2")
#' y <- sim_binomial_probit(W=W,X=X,beta=beta,rho=rho,model="SAR")
#' d <- as.data.frame(cbind(y,X))
#' mod <- ProbitSpatialFit(y~X1+X2,d,W,
#'        DGP='SAR',method="conditional",varcov="varcov")
#'}
#' @export

ProbitSpatialFit<-function(formula,data,W,DGP='SAR',method="conditional",varcov="varcov",M=NULL,control=list()){
  con <- list(iW_CL=6,iW_FL=0,iW_FG=0,reltol=1e-05,prune=1e-4,silent=TRUE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.extract(mf, "response")
  Y <- as.numeric(Y>0)
  X <- model.matrix(mt, mf)  
 idx1<-match(c("(Intercept)","Intercept","(intercept)","intercept"),colnames(X))
  if (any(!is.na(idx1))) {  colnames(X)[idx1[!is.na(idx1)]] <- "(Intercept)"}
  myenv <- new.env()
  try(rm(myenv),silent =TRUE)
  myenv <- new.env()
  myenv[["appiWCL"]] <- con$iW_CL
  myenv[["appiWFL"]] <- con$iW_FL
  myenv[["appiNFG"]] <- con$iW_FG
  myenv[["eps"]] <- con$prune
  if(is.null(W) | any(abs(Matrix::rowSums(W)-1)>1e-12)) stop('W must be a valid row normalized spatial weight matrix')
  if(class(W)=="matrix") W <- Matrix::Matrix(W)
  myenv[["WW"]] <- W
  myenv[["MM"]] <- M
  myenv[["ind"]] <- X
  myenv[["de"]] <- Y
  myenv[["reltol"]] <- con$reltol
  message <- 0
  
  ##### Conditional method
  if(varcov=='precision') method_sigma='UP' else method_sigma='UC'
  if (DGP=="SAR" ){DGP_1='SAR'}
  if (DGP=="SEM" ){DGP_1='SEM'}
  if (DGP=="SARAR"){DGP_1='SARAR'}
  lik <- get(paste('conditional',DGP_1,method_sigma,sep='_'))
  llik <- get(paste('lik',DGP_1,method_sigma,sep='_'))
  init_cond <- Sys.time()
  out <- suppressWarnings(lik(myenv))
  fint_cond <- Sys.time() 
  tim_cond <- as.numeric(difftime(fint_cond,init_cond,units='secs'))
  mycoef_cond <- unlist(out$par)
  myenv$l_cond <- out$l
	if (DGP_1 == 'SARAR') {second_place=c("rho","lambda")
	} else 	{second_place='rho'}
  names(mycoef_cond) <- c(colnames(X),second_place)
  retourcond=1
  
  ##### Full-likelihood method (if requested)
  if(method == "full-lik" && DGP %in% c("SAR","SEM"))
  { retourcond=0
    if (con$prune==0) method_grad <- "FG" else method_grad <- "AG"
    ggrad <- get(paste('grad',DGP,method_sigma,method_grad,sep='_'))
    init_FL <- Sys.time()
    out = optim(mycoef_cond,llik,ggrad,myenv,method = "BFGS",control = list(reltol = myenv$reltol))
    fint_FL <- Sys.time()
    tim_FL <- as.numeric(difftime(fint_FL,init_FL,units='secs'))
    if(!is.list(out)) {
      retourcond=1
      cat('Convergence failed with Full Maximum Likelihood: try to increase iW_FL and/or decrease prune (in case of approximate gradients) and/or stick to the results of conditional likelihood estimation')
      message=1
    } else {
      mycoef_FL=unlist(out$par)
      names(mycoef_FL) <- c(colnames(X),'rho')
      myenv$l_FL=out$value
    }
  }
  
  ##### Standard deviation for betas estimates conditional on rho
  ##### CHECK THIS PART FOR SARAR SEM 
  if (retourcond==1){mycoef <- mycoef_cond
	  lhat<-myenv$l_cond} else {
	  mycoef <- mycoef_FL
	  lhat <- myenv$l_FL}
  k <- ncol(X)
  beta <- mycoef[1:k]
  rho <- mycoef[k+1]
  if (DGP=="SARAR") lambda <- mycoef[k+2]
  if (con$silent==FALSE){
    # fit contains xstar= iW %*% X and v = sqrt(Matrix::diag(iW))
    iW <- ApproxiW(W,rho,con$iW_CL)
	if (DGP %in% c("SAR","SARAR")) {
		xstar <- iW %*% X
		} else {xstar <- X}
	if (DGP %in% c("SAR","SEM")) {
	    v <- sqrt(Matrix::diag(Matrix::tcrossprod(iW)))
		} else {
		iM <- ApproxiW(M,lambda,con$iW_CL)
		v <- sqrt(Matrix::diag(Matrix::tcrossprod(iW %*% iM)))
		}
    xb <- as.numeric(xstar %*% beta)/v
    p <- pnorm(xb)
    p[which(p==1)] <- 0.9999
    p[which(p==0)] <- 0.0001
    g <- (dnorm(xb)^2)/(p * (1 - p))
    gmat <- as.matrix((sqrt(g)/v) * xstar)
    vmat1 <- Matrix::solve(Matrix::crossprod(gmat))
    semat1 <- sqrt(Matrix::diag(vmat1))
    
    # Likelihood-ratio test for rho parameter
	Beta0 <- coef(glm.fit(X, Y,family=binomial(link='probit')))
    if (DGP=="SARAR"){
	LR_rho<-2*(llik(c(Beta0,0,lambda),myenv)-lhat)
	LR_lambda<-2*(llik(c(Beta0,rho,0),myenv)-lhat)
	LR <- c(LR_rho,LR_lambda)
	} else { 
	LR <- 2*(llik(c(Beta0,0),myenv)-lhat)}
    cat("St. dev. of beta conditional on rho and Lik-ratio of rho", "\n")
    
    if (DGP=="SARAR") {fin_res<-c(beta,rho,lambda)} else {fin_res<-c(beta,rho)}
    outmat <- cbind(fin_res, c(semat1,LR*NA), c(beta/semat1,LR), 
                    c(2 * (1 - pnorm(abs(beta)/semat1)),pchisq(LR, 1, lower.tail = FALSE)))
    colnames(outmat) <- c("Estimate", "Std. Error", "z-value", 
                          "Pr(>|z|)")
	if (DGP_1 == 'SARAR') {
	second_place=c("rho","lambda")
	} else 	{second_place='rho'}
    rownames(outmat) <- c(colnames(X),second_place)
    # print(outmat)
    }
  
  out <- new("ProbitSpatial", 
	beta	= mycoef[1:k],
    rho		= mycoef[1+k],
    lambda 	= ifelse(DGP=="SARAR",mycoef[2+k],0),
    coeff	= mycoef,
    loglik	= ifelse(retourcond==1,myenv$l_cond,myenv$l_FL),
    formula	= formula,
    nobs	= ncol(W),
    nvar	= k,
    y		= Y,
    X		= X,
    time	= ifelse(retourcond==1,tim_cond,tim_FL),
    DGP		= DGP,	
    method  = method,
    varcov  = varcov,
    W		= W,
    M 		= M,
    iW_CL	= con$iW_CL,
    iW_FL	= con$iW_FL,
    iW_FG	= con$iW_FG,
    reltol	= con$reltol,
    prune	= con$prune,
    env		= myenv,
    message	= message)
  return(out)
}  

#' Probit with Spatial Dependence, SAR, SEM, and SARAR Models.
#'
#' \code{ProbitSpatial} package allows to fit spatial autoregressive (SAR) and 
#'	spatial error (SEM) probit models. It also provides functions to simulated 
#'	spatial binary data, an empirical data set and different methods for the 
#'	diagnostic of the estimated model.
#'	
#'	The main function of this package is \code{ProbitSpatialFit}. It allows to 
#'	fit both SAR and SEM models for big datasets in a reasonable time. The 
#'	function is based on the maximisation of the approximate likelihood 
#'	function. The approximation is inspired by the Mendell and Elston algorithm 
#'	for computing multivariate normal probabilities and take advantage of the 
#'	sparsity of the spatial weight matrix. Two methods are available for the 
#'	estimation of the model parameter: the first one is known as conditional 
#'	method (see Case (1992)) and performs relatively well in terms of accuracy 
#'	of the estimated parameters and is very rapid. The second method, that 
#'	minimises the full-log-likelihood, is slower but it should be more accurate. 
#'	Monte Carlo experiments on simulated data reported in Martinetti and Geniaux 
#'	(2017) showed that the full-log-likelihood approach is not always 
#'	overperforming the conditional method in terms of accuracy. At the present 
#'	stage, our suggestion is to use the conditional method for a first 
#'	estimation and only attempt the full-likelihood approach in a second moment, 
#'	when the dataset size is not bigger than a few thousands. 
#'
#'	Another feature of the \code{ProbitSpatialFit} function is the possibility 
#'	to fit the model using the precision matrix instead of the 
#'	variance-covariance matrix, since it is usually sparser and hence allows 
#'	faster computations (see LeSage and Pace (2009)).
#'
#'	The output of \code{ProbitSpatialFit} function is an object of class 
#'	\code{ProbitSpatial}, for which the methods  
#'	\code{residuals}, \code{fitted}, \code{effects}, \code{predict} and 
#'	\code{coef} are available.
#' 
#'	The package also contains the function \code{sim_binomial_probit} that 
#'	allows to simulate data samples of both SAR and SEM models. It can be used 
#'	to replicate the Monte Carlo experiments reported in Martinetti and Geniaux 
#'	(2017) as well as the experiment of Calabrese and Elkink (2014). 
#'	An empirical data set \code{\link{Katrina}} on the reopening decisions of 
#'	firms in the aftermath of the Katrina Hurricane in New Orleans is also 
#'	available (LeSage et al.(2011)).
#'  
#'  Other packages in CRAN repository on the same subject are 
#'	\code{McSpatial} (McMillen (2013)) and 
#'	\code{spatialprobit} (Wilhelm and Godinho de Matos 
#'	(2013)).
#'	
#'	The core functions of the present package have been coded using the 
#'	\code{Rcpp} and \code{RcppEigen} libraries (Bates and Eddelbuettel (2013)), 
#'	that allow direct interchange of rich R objects between R and C++.
#'  	
#'	@author Davide Martinetti \email{davide.martinetti@@inra.fr} and 
#' 	Ghislain Geniaux  \email{ghislain.geniaux@@inra.fr}
#'	
#'	@references  
#'	\describe{
#'	\item{Bates and Eddelbuettel (2013)}{D. Bates and D. Eddelbuettel. Fast and 
#'	elegant numerical linear algebra using the RcppEigen package. \emph{Journal 
#'	of Statistical Software} 52, 1--24, 2013.}
#' \item{Case (1992)}{A. C. Case. Neighborhood Influence and Technological 
#'	Change. \emph{Regional Science and Urban Economics} 22, 491--508, 1992.}
#' \item{Calabrese and Elkink (2014)}{R. Calabrese and J.A. Elkink. Estimators 
#'	of binary spatial autoregressive models: a Monte Carlo study. \emph{Journal 
#'	of Regional Science} 54, 664--687, 2014.}
#' \item{LeSage and Pace (2009)}{J. LeSage and R.K. Pace. \emph{Introduction to 
#'	Spatial Econometrics}, CRC Press, chapter 10.1.6, 2009.}
#' \item{LeSage et al. (2011)}{P. LeSage, R. K. Pace, N. Lam, R. Campanella and 
#'	X. Liu. New Orleans business recovery in the aftermath of Hurricane 
#'	Katrina. \emph{Journal of the Royal Statistical Society A} 174, 1007--1027, 
#'	2011.}
#' \item{Martinetti and Geniaux (2017)}{D. Martinetti and G. Geniaux. 
#' 	Approximate likelihood estimation of spatial probit models. \emph{Regional 
#' 	Science and Urban Economics} 64, 30-45, 2017.}
#' \item{McMillen (2013)}{D. McMillen. McSpatial: Nonparametric spatial data 
#'	analysis. R package version 2.0, 
#'	\url{http://CRAN.R-project.org/package=McSpatial}, 2013.}
#' \item{Mendell and Elston (1974)}{N. Mendell and R. Elston. Multifactorial 
#' 	qualitative traits: genetic analysis and prediction of recurrence risks. 
#'  \emph{Biometrics} 30, 41--57, 1974.}
#' \item{Wilhelm and Godinho de Matos (2013)}{S. Wilhelm and M. Godinho de 
#'	Matos. Estimating Spatial Probit Models in R. \emph{The R Journal} 5, 
#'	130--143, 2013.}
#' }
"_PACKAGE"
 
#' Conditional SAR UC. 
#' 
#' Performs conditional estimation of SAR model with variance-covariance matrix.
#' @usage conditional_SAR_UC(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export

conditional_SAR_UC <- function(myenv)
{
 logl <- function(rho) { 
 -lik_SAR_UC_conditional(rho,myenv)$l}
 rho = optimize(logl, lower = -1, upper = 1, tol= myenv$reltol)
 out <- list(lik_SAR_UC_conditional(rho$minimum,myenv)$beta,rho$minimum)
 names(out) <- c("beta", "rho")
 return(list(par=out,l=rho$objective))}
 
#' Conditional SAR UP. 
#'
#' Performs conditional estimation of SAR model with precision matrix.
#' @usage conditional_SAR_UP(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export
 conditional_SAR_UP <- function(myenv)
{ 
 logl <- function(rho) { 
 -lik_SAR_UP_conditional(rho,myenv)$l}
 rho = optimize(logl, lower = -1, upper = 1, tol= myenv$reltol)
 out <- list(lik_SAR_UP_conditional(rho$minimum,myenv)$beta,rho$minimum)
 names(out) <- c("beta", "rho")
 return(list(par=out,l=rho$objective))}
 
#' Conditional SEM UC. 
#' 
#' Performs conditional estimation of SEM model with variance-covariance matrix.
#' @usage conditional_SEM_UC(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export
 conditional_SEM_UC <- function(myenv)
{  logl <- function(rho) { 
 -lik_SEM_UC_conditional(rho,myenv)$l}
 rho = optimize(logl, lower = -1, upper = 1, tol= myenv$reltol)
 out <- list(lik_SEM_UC_conditional(rho$minimum,myenv)$beta,rho$minimum)
 names(out) <- c("beta", "rho")
 return(list(par=out,l=rho$objective))}
 
#' Conditional SEM UP. 
#' 
#' Performs conditional estimation of SEM model with precision matrix.
#' @usage conditional_SEM_UP(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export
conditional_SEM_UP <- function(myenv)
{  logl <- function(rho) { 
 -lik_SEM_UP_conditional(rho,myenv)$l}
 rho = optimize(logl, lower = -1, upper = 1, tol= myenv$reltol)
 out <- list(lik_SEM_UP_conditional(rho$minimum,myenv)$beta,rho$minimum)
 names(out) <- c("beta", "rho")
 return(list(par=out,l=rho$objective))}
 
#' Conditional SARAR UC. 
#' 
#' Performs conditional estimation of SARAR model with variance-covariance 
#' matrix.
#' @usage conditional_SARAR_UC(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export 
   conditional_SARAR_UC <- function(myenv)
  {
    logl <- function(par) { -lik_SARAR_UC_conditional(par,myenv)$l}
    inipar<-c(0,0)
    Ui<-matrix(c(1,0,-1,0,0,1,0,-1),4,2)
    bvec <- c(-0.99, -0.99, -0.99,-0.99)
    res<-constrOptim(inipar,logl,NULL,ui=Ui,ci=bvec, control=list(reltol=1e-4))
   out<-list(par=c(lik_SARAR_UC_conditional(res$par,myenv)$beta,
	   res$par[1],res$par[2]),res$value)
    names(out) <- c("par","l")
    return(out)
  }
  
#' Conditional SARAR UP. 
#' 
#' Performs conditional estimation of SARAR model with precision matrix.
#' @usage conditional_SARAR_UP(myenv)
#' 
#' @param myenv an \code{environment}.
#' @return the log-likelihood and the estimated parameters.
#' 
#' @details We discourage the direct use of this function.
#' 
#' @export 
   conditional_SARAR_UP <- function(myenv)
  {
    logl <- function(par) { -lik_SARAR_UP_conditional(par,myenv)$l}
    inipar<-c(0,0)
    Ui<-matrix(c(1,0,-1,0,0,1,0,-1),4,2)
    bvec <- c(-0.99, -0.99, -0.99,-0.99)
    res<-constrOptim(inipar,logl,NULL,ui=Ui,ci=bvec, control=list(reltol=1e-4))
   out<-list(par=c(lik_SARAR_UP_conditional(res$par,myenv)$beta,
	   res$par[1],res$par[2]),res$value)
    names(out) <- c("par","l")
    return(out)
  }  
 

#' Generate a random spatial weight matrix. 
#'
#' Generate a spatial weight matrix of given size and number of nearest 
#'	neighbors from randomly-located observations on the unit square.
#' 
#' @usage generate_W(n, nneigh, seed=123)
#' 
#' @param n the size of the matrix.
#' @param nneigh the number of nearest neighbors.
#' @param seed an integer to set the seed for the random generated
#' 	locations.
#' 
#' @return a matrix of class \code{dgCMatrix} (sparse matrix). 
#' 
#' @details The output matrix has zero diagonal and it is row-standardised. 
#' 	The \code{n} observations are allocated randomly in the unit square.  
#' 	For each observation, the \code{nneigh} closests observations w.r.t. the 	
#' 	Euclidean distance are assigned with a weight equal to 1/\code{nneigh}.
#' 
#' @seealso \code{\link{sim_binomial_probit}}.
#'
#' @examples
#' W <- generate_W(100,4,seed=12)
#' 
#'
#' @export

generate_W <-function(n,nneigh,seed=123)
{
  set.seed(seed)
  coord <- cbind(runif(n,0,1),runif(n,0,1))
  k=nneigh+1
  nb1 <- RANN::nn2(as.matrix(coord), k=k ,treetype = c("bd"))
  W <- Matrix::sparseMatrix(i=rep(seq_along(rep(k,n)),rep(k,n)),j=t(nb1$nn.idx),x=1/k) 
  diag(W) <- 0
  if(class(W)=='matrix') W <- Matrix::Matrix(W)
  ret <- Matrix::summary(as(W, "dgCMatrix") )
  ti <- tapply(ret$x,ret$i,function(x) sum(x,na.rm=TRUE))
  ret$x <- as.numeric(ret$x/ti[match(ret$i,as.numeric(names(ti)))])
  W <- Matrix::sparseMatrix(i = ret$i, j =  ret$j,x= ret$x,dims=dim(W))
  Matrix::drop0(W)
  W
}


#' Simulate the dependent variable of a SAR/SEM/SARAR model.
#' 
#' The function \code{sim_binomial_probit} is used to generate the dependent 
#' 	variable of a spatial binomial probit model, where all the data and 
#' 	parameters of the model can be modified by the user. 
#' 
#' @usage sim_binomial_probit(W,X,beta,rho,model="SAR",M=NULL,lambda=NULL,
#' sigma2=1,ord_iW=6,seed=123)
#' 
#' @param W the spatial weight matrix (works for \code{"SAR"} and 
#' 	\code{"SEM"} models).
#' @param X the matrix of covariates.
#' @param beta the value of the covariates parameters.
#' @param rho the value of the spatial dependence parameter (works for 
#' 	\code{"SAR"} and \code{"SEM"} models).
#' @param model the type of model, between \code{"SAR"}, \code{"SEM"}, 
#' 	\code{"SARAR"} (Default is \code{"SAR"}).
#' @param M the second spatial weight matrix (only if \code{model} is 
#' 	\code{"SARAR"}).
#' @param lambda the value of the spatial dependence parameter (only if 
#' 	\code{model} is \code{"SARAR"}).
#' @param sigma2  the variance of the error term (Defaul is 1).
#' @param ord_iW the order of approximation of the matrix 
#' 	\eqn{(I_n-\rho W)^{-1}}.
#' @param seed to set the random generator seed of the error term.
#' 
#' @return a vector of zeros and ones 
#' 
#' @details The \code{sim_binomial_probit} generates a vector of dependent 
#' 	variables for a spatial probit model. It allows to simulate the following 
#' 	DGPs (Data Generating Process):
#' 	SAR
#' 	\deqn{z = (I_n-\rho W)^{-1}(X\beta+\epsilon)	}
#' 	SEM
#' 	\deqn{z = X\beta+(I_n-\rho W)^{-1}\epsilon	}
#' 	SARAR
#' 	\deqn{z = (I_n-\rho W)^{-1}(X\beta+(I_n-\lambda M)^{-1}\epsilon)	}
#' 	where \eqn{\epsilon} are independent and normally distributed with mean zero 
#' 	and variance \code{sigma2} (default is 1).
#' 
#' 	The matrix \code{X} of covariates, the corresponding parameters \code{beta}, 
#' 	the spatial weight matrix \code{W} and the corresponding spatial dependence 
#' 	parameter \code{rho} need to be passed by the user. Eventually, the same 
#' 	applies for \code{lambda} and \code{M} for the SARAR model.
#' 	
#' 	The matrix \eqn{(I_n-\rho W)^{-1}} is computed using the 
#' 	\code{ApproxiW} function, that can either invert \eqn{(I_n-\rho W)} 
#' 	exactely, if \code{order_iW=0} (not suitable for \code{n} bigger than 1000),  
#' 	or using the Taylor approximation 
#' 	\deqn{(I_n-\rho W)^{-1}= I_n+\rho W+\rho^2 W^2+\ldots 	}
#' 	of order \code{order_iW} (default is approximation of order 6).
#' 
#' @seealso \code{\link{generate_W}}, \code{\link{ProbitSpatialFit}}.
#'
#' @examples
#' n <- 500
#' nneigh <- 3
#' rho <- 0.5
#' beta <- c(4,-2,1)
#' W <- generate_W(n,nneigh)
#' X <- cbind(1,rnorm(n,2,2),rnorm(n,0,1))
#' #SAR
#' y <- sim_binomial_probit(W,X,beta,rho,model="SAR") #SAR model
#' #SEM
#' y <- sim_binomial_probit(W,X,beta,rho,model="SEM") #SEM model
#' #SARAR
#' M <- generate_W(n,nneigh,seed=1)
#' lambda <- -0.5
#' y <- sim_binomial_probit(W,X,beta,rho,model="SARAR",M=M,lambda=lambda) 
#'
#'
#' @export

sim_binomial_probit<-function(W,X,beta,rho,model="SAR",M=NULL,lambda=NULL,
					sigma2=1,ord_iW=6,seed=123){
  set.seed(seed)
  k <- length(beta)
  n <- dim(W)[1]
  iW <- ApproxiW(W, rho, ord_iW) 
  if (!is.null(M)){  iM <- ApproxiW(M, lambda, ord_iW) }
  if (dim(X)[1] != dim(W)[1]) (cat("X and W need to have the same number of rows"))
  if (dim(X)[2] != length(beta))(cat("the number of columns of X has to equal to the length of beta"))
  e_sim <- rnorm(n,0,sigma2) 
  if (model=="SAR"){z <- iW%*%(X %*% beta + e_sim)}
  if (model=="SEM"){z <- (X %*% beta + iW%*%e_sim)}
  if (model=="SARAR"){z <- iW%*%(X %*% beta + iM%*%e_sim)}
  Y_sim <- as.double(z >= 0)
  return(Y_sim)
}


#' Spatial probit model summaries.
#' 
#' Print the results of a \code{ProbitSpatial} model.
#'
#' @param object an object of class \code{ProbitSpatial}.
#' @param covar should the statistics be computed with the matrix of 	
#'	variance of the parametes or not. Default is FALSE, hence Likelihood-ratio 
#' 	statistics are printed.
#' @param ... further arguments
#'
#' @return This functions does not return any value.
#'
#' @details The \code{summary} function prints
#' \describe{
#'   \item{Model}{Featurs on the model and dataset.}
#'   \item{Time}{Estimation time.}
#'   \item{Statistics}{Standard errors of the estimated parameters. If 
#' \code{covar=TRUE}, it uses the matrix of variance of the parameters, else the 
#' likelihood ratio test.}
#'   \item{Accuracy}{Confusion Matrix and accuracy of the estimated model.}
#' }
#'
#' @method summary ProbitSpatial
#' @export

summary.ProbitSpatial <- function(object,covar=FALSE,...) {
 cat("-- Univariate conditional estimation of spatial probit --\n\n")
 cat("Sample size = ", object@nobs,"\n")
 cat("Number of covariates = ",object@nvar,"\n")
 cat("DGP = ", object@DGP,"\n")
 cat("estimation method = ", object@method,"\n")
 vc<-ifelse(object@varcov=="varcov", "Var-Covar Matrix", "Precision Matrix")
 cat("Variance covariance = ", vc,"\n" )
 if (object@iW_CL>0) {cat("order of approx. of iW in the conditional step = ", object@iW_CL,"\n")}
 if (object@method=="full-lik"){ 
  if (object@iW_FL>0) {cat("order of approximation of iW in the likelihood function = ", object@iW_FL,"\n")}
  if (object@iW_FG>0) {cat("order of approximation of iW in the gradient function = ", object@iW_FG,"\n")}
  if (object@prune>0) {cat("pruning in the gradient functions = ", object@prune,"\n")}
  }
 cat("Execution time = ", object@time,"\n\n")
                        cat("-----------------------------------------------\n\n")

 mod_covar <- ifelse(object@varcov=='varcov',"UC","UP")
 if (covar == TRUE){
 mycoef <- object@coeff
 lik <- get(paste('lik',object@DGP,mod_covar, sep='_'))
 H <- numDeriv::hessian(lik,x=mycoef,env=object@env)
 se <- sqrt(diag(abs(solve(H))))
 outmat <- cbind(mycoef, se, mycoef/se, 2 * (1 - pnorm(abs(mycoef)/se)))
 colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>z)")
 rownames(outmat) <- names(mycoef)
 cat("Unconditional standard errors with variance-covariance matrix\n\n")
 } else { lik<-get(paste('conditional',object@DGP,mod_covar, sep='_'))
 llik<-get(paste('lik',object@DGP, mod_covar,sep='_'))
 env1 <- object@env
 Beta0 <-   suppressWarnings(coef(glm.fit(object@X,object@y,intercept=FALSE,family=binomial(link='probit'))))
 if (object@DGP=="SARAR"){	
 	LR_rho=-2*(llik(c(Beta0,0,object@lambda),env1)-object@loglik)
  	LR_lambda=-2*(llik(c(Beta0,object@rho,0),env1)-object@loglik)
  	LR <- c(LR_rho,LR_lambda)
   } else { LR=-2*(llik(c(Beta0,0),env1)-object@loglik) }
 
 LR_beta <- c()
 for(i in 1:object@nvar){
 	XX <- env1$ind
    env1$ind <- as.matrix(XX[,-i])
    lc=lik(env1)
    env1$ind <- as.matrix(XX)
    LR_beta <- c(LR_beta,-2*(lc$l-object@loglik))}		
  LRtheta <- abs(as.numeric(c(LR_beta,LR)))
  outmat <- cbind(object@coeff, LRtheta, pchisq(LRtheta, 1, lower.tail = FALSE))
  colnames(outmat) <- c("Estimate", "LR test", "Pr(>z)")
  if(object@DGP=="SARAR"){rownames(outmat)<-c(colnames(object@X),'rho','lambda') 
  } else { rownames(outmat) <- c(colnames(object@X),'rho') }

  cat("Unconditional standard errors with likelihood-ratio test\n")
  }
  print(outmat)
 
  cat("\n-----------------------------------------------\n")
  f<-fitted(object,type="binary")
  y<-object@y
  TP <- sum(f==1 & y==1)
  TN <- sum(f==0 & y==0)
  FP <- sum(f==1 & y==0)
  FN <- sum(f==0 & y==1)
  conf_matrix<-matrix(c(TP,FN,FP,TN),2,2)
  colnames(conf_matrix) <- c("pred 1","pred 0")
  rownames(conf_matrix) <- c("true 1","true 0")
  cat("Confusion Matrix:\n")
  print(conf_matrix)
  cat("Accuracy:\t", (TP+TN)/(TP+TN+FP+FN), "\n")
  cat("Sensitivity:\t", (TP)/(TP+FN), "\t Specificity:\t",(TN)/(FP+TN),"\n")
  cat("Pos Pred Value:\t", (TP)/(TP+FP), "\t Neg Pred Value:",(TN)/(TN+FN),"\n")
 }

#'	New Orleans business recovery in the aftermath of Hurricane Katrina.
#'	
#'	This dataset has been used in the LeSage et al. (2011) paper entitled "New 
#'	Orleans business recovery in the aftermath of Hurricane Katrina" to study 
#'	the decisions of shop owners to reopen business after Hurricane Katrina. The 
#'	dataset contains 673 observations on 3 streets in New Orleans and can be 
#'	used to estimate the spatial probit models and to replicate the findings in 
#'	the paper.
#'	
#' @usage data(Katrina)
#' @docType data
#' @format
#' Katrina is a data frame with 673 observations on the following 15 
#' variables:
#' \describe{
#'   \item{\code{code}}{a numeric vector}
#'   \item{\code{long}}{longitude coordinate of store}
#'   \item{\code{lat}}{latitude coordinate of store}
#'   \item{\code{street1}}{a numeric vector}
#'   \item{\code{medinc}}{median income}
#'   \item{\code{perinc}}{a numeric vector}
#'   \item{\code{elevation}}{a numeric vector}
#'   \item{\code{flood}}{flood depth (measured in feet)}
#'   \item{\code{owntype}}{type of store ownership: "sole proprietorship" vs. 
#'		"local chain" vs. "national chain"}
#'   \item{\code{sesstatus}}{socio-economic status of clientele (1-5): 1-2 = low #'		status customers, 3 = middle, 4-5 = high status customers}
#'   \item{\code{sizeemp}}{"small size" vs. "medium size" vs. "large size" 
#'   	firms}
#'   \item{\code{openstatus1}}{a numeric vector}
#'   \item{\code{openstatus2}}{a numeric vector}
#'   \item{\code{days}}{days to reopen business}
#'   \item{\code{street}}{1=Magazine Street, 2=Carrollton Avenue, 3=St. Claude 
#'   	Avenue}
#' }
#' 
#' Katrina is a data frame with 673 observations on the following 13 variables.
#' \describe{
#'   \item{\code{long}}{longitude coordinate of store}
#'   \item{\code{lat}}{latitude coordinate of store}
#'   \item{\code{flood_depth}}{flood depth (measured in feet)}
#'   \item{\code{log_medinc}}{log median income}
#'   \item{\code{small_size}}{binary variable for "small size" firms}
#'   \item{\code{large_size}}{binary variable for "large size" firms}
#'   \item{\code{low_status_customers}}{binary variable for low socio-economic 
#'		status of clientele}
#'   \item{\code{high_status_customers}}{binary variable for high socio-economic 
#'   	status of clientele}
#'   \item{\code{owntype_sole_proprietor}}{a binary variable indicating "sole 
#'		proprietor" ownership type}
#'   \item{\code{owntype_national_chain}}{a binary variable indicating 
#'		"national_chain" ownership type}
#'   \item{\code{y1}}{reopening status in the very short period 0-3 months; 
#'   	1=reopened, 0=not reopened}
#'   \item{\code{y2}}{reopening status in the period 0-6 months; 1=reopened, 
#'		0=not reopened}
#'   \item{\code{y3}}{reopening status in the period 0-12 months; 1=reopened, 
#'		0=not reopened}
#' }
#'
#' @details
#'
#'	The Katrina dataset contains the data found on the website before some 
#'	of the variables are recoded. For example, the socio-economic status of 
#'	clientele is coded as 1-5 in the raw data, but only 3 levels will be used in 
#'	estimation: 1-2 = low status customers, 3 = middle, 4-5 = high status 
#'	customers. Hence, with "middle" as the reference category, Katrina contains 
#'	2 dummy variables for low status customers and high status customers.
#'
#'	The dataset Katrina is the result of these recoding operations and can be 
#'	directly used for model estimation.
#'
#' @note
#'	When definining the reopening status variables y1 (0-3 months), y2 (0-6 
#'	months), and y3 (0-12 months) from the days variable, the Matlab code 
#'	ignores the seven cases where days=90. To be consistent with the number of 
#'	cases in the paper, we define y1,y2,y3 in the same way: y1=sum(days < 90), 
#'	y2=sum(days < 180 & days != 90), y3=sum(days < 365 & days != 90). So this is 
#'	not a bug, its a feature.
#'
#' @source
#'	The raw data was obtained from the Royal Statistical Society dataset website 
#'	\url{www.blackwellpublishing.com/rss/Volumes/Av174p4.htm} and brought 
#'	to RData format by Wilhelm and Godinho de Matos (2013). 
#'  
#' @examples
#' \dontrun{
#' 	data(Katrina)
#' 	attach(Katrina)
#' 	table(y1) # 300 of the 673 firms reopened during 0-3 months horizon, p.1016
#' 	table(y2) # 425 of the 673 firms reopened during 0-6 months horizon, p.1016
#' 	table(y3) # 478 of the 673 firms reopened during 0-12 months horizon, p.1016
#' 	detach(Katrina)
#'
#' 
#' 	# replicate LeSage et al. (2011), Table 3, p.1017
#' 	require(spdep)
#'  
#' 	# (a) 0-3 months time horizon
#' 	# LeSage et al. (2011) use k=11 nearest neighbors in this case
#' 	nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11))
#' 	listw <- nb2listw(nb, style="W")
#' 	W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
#' 
#' 	fit1_cond <- ProbitSpatialFit(y1 ~ flood_depth + log_medinc + small_size + 
#' 		large_size +low_status_customers +  high_status_customers + 
#' 		owntype_sole_proprietor + owntype_national_chain, 
#' 		W=W1, data=Katrina, DGP='SAR', method="conditional", varcov="varcov")
#' 	summary(fit1_cond)
#'
#' 	fit1_FL <- ProbitSpatialFit(y1 ~ flood_depth + log_medinc + small_size + 
#' 		large_size +low_status_customers +  high_status_customers + 
#' 		owntype_sole_proprietor + owntype_national_chain, 
#' 		W=W1, data=Katrina, DGP='SAR', method="full-lik", varcov="varcov")
#' 	summary(fit1_FL)
#'
#' 	fit1_cond_10nn <- ProbitSpatialFit(y1 ~ flood_depth+ log_medinc+ small_size+
#' 		large_size +low_status_customers +  high_status_customers + 
#' 		owntype_sole_proprietor + owntype_national_chain, 
#' 		W=W1, data=Katrina, DGP='SAR', method="conditional", varcov="varcov",
#' 		control=list(iW_CL=10))
#' 	summary(fit1_cond_10nn)
#'
#'	# (b) 0-6 months time horizon
#'	# LeSage et al. (2011) use k=15 nearest neighbors
#'	nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=15))
#'	listw <- nb2listw(nb, style="W")
#'	W2 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
#'	
#'	fit2_cond <- ProbitSpatialFit(y2 ~ flood_depth + log_medinc + small_size + 
#'		large_size + low_status_customers + high_status_customers + 
#'		owntype_sole_proprietor + owntype_national_chain, 
#'		W=W2, data=Katrina, DGP="SAR", method="full-lik", varcov="varcov")
#'	summary(fit2_cond)  
#'	
#'	fit2_FL <- ProbitSpatialFit(y2 ~ flood_depth + log_medinc + small_size + 
#'		large_size + low_status_customers + high_status_customers + 
#'		owntype_sole_proprietor + owntype_national_chain, 
#'		W=W2, data=Katrina, DGP="SAR", method="full-lik", varcov="varcov")
#'	summary(fit2_FL)  
#'	
#'	# (c) 0-12 months time horizon
#'	# LeSage et al. (2011) use k=15 nearest neighbors as in 0-6 months
#'	W3 <- W2
#'	fit3_cond <- ProbitSpatialFit(y3 ~ flood_depth + log_medinc + small_size + 	
#'		large_size + low_status_customers + high_status_customers + 
#'		owntype_sole_proprietor + owntype_national_chain, 
#'		W=W3, data=Katrina, DGP="SAR", method="conditional", varcov="varcov")
#'	summary(fit3_cond)
#'	
#'	fit3_FL <- ProbitSpatialFit(y3 ~ flood_depth + log_medinc + small_size + 
#'		large_size + low_status_customers + high_status_customers + 
#'		owntype_sole_proprietor + owntype_national_chain, 
#'		W=W3, data=Katrina, DGP="SAR", method="full-lik", varcov="varcov")
#'	summary(fit3_FL)
#'	
#'	# replicate LeSage et al. (2011), Table 4, p.1018
#'	# SAR probit model effects estimates for the 0-3-month time horizon
#'	effects(fit1_cond)  
#'	
#'	# replicate LeSage et al. (2011), Table 5, p.1019
#'	# SAR probit model effects estimates for the 0-6-month time horizon
#'	effects(fit2_cond)
#'	
#'	# replicate LeSage et al. (2011), Table 6, p.1020
#'	# SAR probit model effects estimates for the 0-12-month time horizon
#'	effects(fit3_cond)
#'	}
#'	
#'	@references
#'	\describe{
#'  \item{LeSage et al. (2011)}{P. LeSage, R. K. Pace, N. Lam, R. Campanella and 
#'	X. Liu. New Orleans business recovery in the aftermath of Hurricane 
#'	Katrina. \emph{Journal of the Royal Statistical Society A}, 174, 1007--1027, 
#'	2011.}
#'  \item{Wilhelm and Godinho de Matos (2013)}{S. Wilhelm and M. Godinho de 
#'	Matos. Estimating Spatial Probit Models in R. \emph{The R Journal} 5, 
#'	130--143, 2013.}
#'	}
"Katrina"



#' Extract from ProbitSpatial class.
#' 
#' extract a slot from \code{ProbitSpatial} class object.
#' 
#' @param x an object of class \code{ProbitSpatial}.
#' @param name of the slot.
#'
#' @return The content of the slot.
#'
#' @method $ ProbitSpatial
setMethod("$",
    signature(x = "ProbitSpatial"),
    function (x, name) 
    {
        slot(x,name)
    }
)
