
// includes from the plugin
#include <RcppEigen.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP ApproxiW( SEXP WW, SEXP la, SEXP order) ;
SEXP lik_SAR_UC( SEXP th, SEXP env) ;
SEXP grad_SAR_UC_FG( SEXP th, SEXP env) ;
SEXP grad_SAR_UC_AG( SEXP th, SEXP env) ;
SEXP lik_SAR_UC_conditional( SEXP th, SEXP env) ;
SEXP lik_SAR_UP( SEXP th, SEXP env) ;
SEXP grad_SAR_UP_FG( SEXP th, SEXP env) ;
SEXP grad_SAR_UP_AG( SEXP th, SEXP env) ;
SEXP lik_SAR_UP_conditional( SEXP th, SEXP env) ;
SEXP lik_SEM_UC( SEXP th, SEXP env) ;
SEXP grad_SEM_UC_FG( SEXP th, SEXP env) ;
SEXP grad_SEM_UC_AG( SEXP th, SEXP env) ;
SEXP lik_SEM_UC_conditional( SEXP th, SEXP env) ;
SEXP lik_SEM_UP( SEXP th, SEXP env) ;
SEXP grad_SEM_UP_FG( SEXP th, SEXP env) ;
SEXP grad_SEM_UP_AG( SEXP th, SEXP env) ;
SEXP lik_SEM_UP_conditional( SEXP th, SEXP env) ;
}

// definition

SEXP ApproxiW( SEXP WW, SEXP la, SEXP order ){
BEGIN_RCPP
using Eigen::SparseMatrix;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 int ord = Rcpp::as<int>(order) ;
 double lambda = Rcpp::as<double>(la) ;
 const SparseMatrix<double> W(as<SparseMatrix<double> >(WW));
 int n = W.rows();
 SparseMatrix<double> A = W;
 double b = lambda;
 typedef Eigen::Triplet<double> T;
 std::vector<T> tripletList;
 tripletList.reserve(n);
 for(int j = 0; j < n; ++j)
{
 tripletList.push_back(T(j,j,1));
}
 SparseMatrix<double> iW(n,n);
 iW.setFromTriplets(tripletList.begin(), tripletList.end());
 iW = iW+lambda*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 b = lambda*b;
 iW = iW + b*A;
}
 //SparseMatrix<double> iWW = iW.transpose()*iW;
 //return List::create(Named("iW")=iW,Named("iWW")=iWW);
 return wrap(iW);
END_RCPP
}


SEXP lik_SAR_UC( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWFL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 double  ep = e["eps"]; 
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);

if  ( abs(rho)>1 ){
    return wrap("Error");
    Rcpp::stop("Parameter rho out of bounds"); }
 
 //Compute the inverse of I-rho*W. If "ord"==0, use exact method, otherwise use Taylor serie appro of order "ord"
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> iW(n,n);
 if (ord==0){
  SparseMatrix<double> IW = I-rho*W;
  Eigen::SparseLU<SparseMatrix<double> > solver;
  solver.compute(IW);
  iW = solver.solve(I);
 iW.prune(ep,1);
  } else {
  SparseMatrix<double> A = W;
  double rr = rho; 
  iW = I + rho*W;
  for(int j = 2; j < ord; ++j)
  {
    A = A * W;
    rr = rho*rr;
    iW = iW + rr*A;
  }
 }
 
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 //Set Xstar<-(iW%*%indep)/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
 if ( (rs.array() < 0).any() ){
    return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = iW*indep;
 Xstar = Xstar.array()/(rs.replicate(1,k-1)).array();
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Sigma=Sigma.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Sigma
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Sigma);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixL();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
if  ( (CC.diagonal().array() == 0).any() ){
    return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n);
 VectorXd hata = VectorXd::Zero(n);
 VectorXd hatb = VectorXd::Zero(n);
 VectorXd VecU = VectorXd::Zero(n);
 VectorXd Vecg = VectorXd::Zero(n);
 double logP =0;
 double U = 0;
 int i = 0;
 hata[i] = (lo[i])/CC.coeff(i,i);
 hatb[i] = (up[i])/CC.coeff(i,i);
 U = R::pnorm(hatb[i],0,1,1,0)- R::pnorm(hata[i],0,1,1,0);
 if (U == 0){
    return wrap("Error");}
 VecU[i] = U; 
 mu[i] = (R::dnorm(hata[i],0,1,0)-R::dnorm(hatb[i],0,1,0))/U;
 logP += log(U);
 
 for (int j = 1; j < n; ++j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 Vecg[j] = g;
 hata[j] = (lo[j]-g)/CC.coeff(j,j);
 hatb[j] = (up[j]-g)/CC.coeff(j,j);
 U = R::pnorm(hatb[j],0,1,1,0)- R::pnorm(hata[j],0,1,1,0);
if (U == 0){
    return wrap("Error");}
 VecU[j] = U; 
 mu[j] = (R::dnorm(hata[j],0,1,0)-R::dnorm(hatb[j], 0,1,0))/U;
 logP += log(U);
 }
 e.assign("eiW",iW);
 e.assign("eSigma",Sigma);
 e.assign("eXstar",Xstar);
 e.assign("exb",xb);
 e.assign("eCC",CC);
 e.assign("eVecU",VecU);
 e.assign("eVecg",Vecg);
 e.assign("ehata",hata);
 e.assign("ehatb",hatb);
 e.assign("emu",mu);
 e.assign("elo",lo);
 e.assign("eup",up);
 e.assign("elogP",logP);
 e.assign("efirstorder",firstorder);
 e.assign("eAMDord",AMDord);
 
 return wrap(-logP);
END_RCPP
}


SEXP grad_SAR_UC_FG( SEXP th, SEXP env ){
BEGIN_RCPP
 
 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = indep.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict lower part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyLower>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 }  
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 SparseMatrix<double> M = Cm1*dSigma*(Cm1.transpose());
 MatrixXd dense_diag_M = (0.5*M.diagonal()).asDiagonal();
 SparseMatrix<double> diag_M = dense_diag_M.sparseView();
 SparseMatrix<double> SLM = M.triangularView<StrictlyLower>();
 SparseMatrix<double,Eigen::RowMajor> dC= CC*(SLM+diag_M);
 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigma.diagonal().array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num1=(buf1*iW*indep).array() * (Sigmadiagsqrt.replicate(1,k-1)).array();
 VectorXd dsqrtdSigma = dSigmadiag.array() / (2*Sigmadiagsqrt).array();
 MatrixXd num2=(iW*indep).array()*(dsqrtdSigma.replicate(1,k-1)).array();
 MatrixXd num = num1-num2;
 num = num.array()/(Sigmadiag.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb= MatrixXd::Zero(n,k-1); 
 VectorXd dUdr= VectorXd::Zero(n);
 VectorXd dLUdr= VectorXd::Zero(n);
 VectorXd dadr= VectorXd::Zero(n);
 VectorXd dbdr= VectorXd::Zero(n);
 VectorXd dhatadr= VectorXd::Zero(n);
 VectorXd dhatbdr= VectorXd::Zero(n);
 VectorXd dmudr= VectorXd::Zero(n);
 VectorXd dgdr = VectorXd::Zero(n);
 
 int i = 0;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*CC.coeff(i,i)-lo(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*CC.coeff(i,i)-up(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 
 for (int i = 1; i < n; ++i) {
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double gdr = 0;
 double cii = CC.coeff(i,i);
 double cii2 = cii*cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 gdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 gdr += it2.value()*mu(it2.index());}
 gdr = gdr - dC.coeff(i,i)*mu(i);
 
 dgdr(i) = (gdr*cii - Vecg(i)*dC.coeff(i,i))/cii2;
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = (zz-g)/cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*cii-lo(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = (zz-g)/cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*cii-up(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 return wrap(grad);
END_RCPP
}


SEXP grad_SAR_UC_AG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 double  ep = e["eps"]; 
 
 CC.prune(ep,1);
 Sigma.prune(ep,1);
 iW.prune(ep,1);
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = indep.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict lower part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyLower>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 } 
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 dSigma.prune(ep,1);
 Cm1.prune(ep,1);
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 SparseMatrix<double> M = Cm1*dSigma*(Cm1.transpose());
 MatrixXd dense_diag_M = (0.5*M.diagonal()).asDiagonal();
 SparseMatrix<double> diag_M = dense_diag_M.sparseView();
 SparseMatrix<double> SLM = M.triangularView<StrictlyLower>();
 SparseMatrix<double,Eigen::RowMajor> dC= CC*(SLM+diag_M);
 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigma.diagonal().array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num1=(buf1*iW*indep).array() * (Sigmadiagsqrt.replicate(1,k-1)).array();
 VectorXd dsqrtdSigma = dSigmadiag.array() / (2*Sigmadiagsqrt).array();
 MatrixXd num2=(iW*indep).array()*(dsqrtdSigma.replicate(1,k-1)).array();
 MatrixXd num = num1-num2;
 num = num.array()/(Sigmadiag.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb= MatrixXd::Zero(n,k-1); 
 VectorXd dUdr= VectorXd::Zero(n);
 VectorXd dLUdr= VectorXd::Zero(n);
 VectorXd dadr= VectorXd::Zero(n);
 VectorXd dbdr= VectorXd::Zero(n);
 VectorXd dhatadr= VectorXd::Zero(n);
 VectorXd dhatbdr= VectorXd::Zero(n);
 VectorXd dmudr= VectorXd::Zero(n);
 VectorXd dgdr = VectorXd::Zero(n);
 
 int i = 0;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*CC.coeff(i,i)-lo(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*CC.coeff(i,i)-up(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 
 for (int i = 1; i < n; ++i) {
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double gdr = 0;
 double cii = CC.coeff(i,i);
 double cii2 = cii*cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 gdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 gdr += it2.value()*mu(it2.index());}
 gdr = gdr - dC.coeff(i,i)*mu(i);
 
 dgdr(i) = (gdr*cii - Vecg(i)*dC.coeff(i,i))/cii2;
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = (zz-g)/cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*cii-lo(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = (zz-g)/cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*cii-up(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 return wrap(grad);
END_RCPP
}


SEXP lik_SAR_UC_conditional( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 double theta = Rcpp::as<double>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWCL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 int n = W.rows();
 int k = indep.cols();
 double rho = theta;
 SparseMatrix<double> A = W;
 double rr = theta; 
 SparseMatrix<double> iW(n,n);
 iW.setIdentity();
 //Compute (I-rho*W)^-1
 iW = iW+rho*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 rr = rho*rr;
 iW = iW + rr*A;
}
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 Environment Mat("package:Matrix"); 
 
 //Set Xstar<-(iW%*%indep)/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
//if ( (rs.array() < 0).any() ){ return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = iW*indep;
 Xstar = Xstar.array()/(rs.replicate(1,k)).array();
 
 Environment speed("package:speedglm");
 Function sglm = speed["speedglm.wfit"];
 Environment stats("package:stats");
 Function binomial = stats["binomial"];
 
 Rcpp::List fit = sglm(_["y"] =dep,_["X"] =Xstar,_["intercept"] =0,_["family"] =binomial("probit"));
 SEXP fitcoef = fit[0];
 VectorXd beta = as<VectorXd>(fitcoef);
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Sigma=Sigma.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Sigma
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Sigma);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixL();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
 //if  ( (CC.diagonal().array() == 0).any() ){ return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n); 
 double logP = 0;
 double a = 0;
 double b = 0;
 double U = 0;
 int i = 0;
 a = (lo[i])/CC.coeff(i,i);
 b = (up[i])/CC.coeff(i,i);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
 //if (U == 0){return wrap("Error");}
 mu[i] = (R::dnorm(a,0,1,0)-R::dnorm(b,0,1,0))/U;
 logP += log(U);
 
 for (int j = 1; j < n; ++j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 a = (lo[j]-g)/CC.coeff(j,j);
 b = (up[j]-g)/CC.coeff(j,j);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
  //if (U == 0){ return wrap("Error");}
 mu[j] = (R::dnorm(a,0,1,0)-R::dnorm(b, 0,1,0))/U;
 logP += log(U);
 }
 
 
 return List::create(Named("l") = logP,Named("beta") = beta);
END_RCPP
}


SEXP lik_SAR_UP( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env);  
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWFL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 double  ep = e["eps"]; 
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 if  ( abs(rho)>1 ){
    return wrap("Error");
    Rcpp::stop("Parameter rho out of bounds"); }
 
 //Compute the inverse of I-rho*W. If "ord"==0, use exact method, otherwise use Taylor serie appro of order "ord"
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> iW(n,n);
 if (ord==0){
 SparseMatrix<double> IW = I-rho*W;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(IW);
 iW = solver.solve(I);
 iW.prune(ep,1);
 } else {
 SparseMatrix<double> A = W;
 double rr = rho; 
 iW = I + rho*W;
 for(int j = 2; j < ord; ++j)
  {
    A = A * W;
    rr = rho*rr;
    iW = iW + rr*A;
  }
 }
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 //Set Xstar<-(iW%*%indep)/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
 if ( (rs.array() < 0).any() ){
    return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = iW*indep;
 Xstar = Xstar.array()/(rs.replicate(1,k-1)).array();
 VectorXd xb = Xstar*beta;
 
 // Compute the precision matrix Omega (the inverse of Sigma (iWW))
 SparseMatrix<double> Omega = I-rho*W;
 Omega = Omega.transpose()*Omega;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength, _["decreasing"]=1);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Omega=Omega.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Precision matrix Omega
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Omega);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixU();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
if  ( (CC.diagonal().array() == 0).any() ){
    return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n);
 VectorXd hata = VectorXd::Zero(n);
 VectorXd hatb = VectorXd::Zero(n);
 VectorXd VecU = VectorXd::Zero(n); 
 VectorXd Vecg = VectorXd::Zero(n); 
 double logP = 0;
 double U = 0;
 int i = n-1;
 hata[i] = (lo[i])*CC.coeff(i,i);
 hatb[i] = (up[i])*CC.coeff(i,i);
 U = R::pnorm(hatb[i],0,1,1,0)- R::pnorm(hata[i],0,1,1,0);
 if (U == 0){
    return wrap("Error");}
 VecU[i] = U;
 if (U == 0){
    return wrap("Error");}
 mu[i] = (R::dnorm(hata[i],0,1,0)-R::dnorm(hatb[i],0,1,0))/(U*CC.coeff(i,i));
 logP += log(U);
 
 for (int j = n-2; j >= 0; --j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 Vecg[j] = g;
 hata[j] = (lo[j]*CC.coeff(j,j))+g;
 hatb[j] = (up[j]*CC.coeff(j,j))+g;
 U = R::pnorm(hatb[j],0,1,1,0)- R::pnorm(hata[j],0,1,1,0);
 if (U == 0){
    return wrap("Error");}
 VecU[j] = U;
 mu[j] = (((R::dnorm(hata[j],0,1,0)-R::dnorm(hatb[j], 0,1,0))/U)-g)/CC.coeff(j,j);
 logP += log(U);
 }
 e.assign("eiW",iW);
 e.assign("eSigma",Sigma);
 e.assign("eOmega",Omega);
 e.assign("eXstar",Xstar);
 e.assign("exb",xb);
 e.assign("eCC",CC);
 e.assign("eVecU",VecU);
 e.assign("eVecg",Vecg);
 e.assign("ehata",hata);
 e.assign("ehatb",hatb);
 e.assign("emu",mu);
 e.assign("elo",lo);
 e.assign("eup",up);
 e.assign("elogP",logP);
 e.assign("efirstorder",firstorder);
 e.assign("eAMDord",AMDord);
 
 return wrap(-logP);
END_RCPP
}


SEXP grad_SAR_UP_FG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using Eigen::Lower;
 using Eigen::StrictlyUpper;
 using Eigen::Upper;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 SparseMatrix<double> Omega = e["eOmega"]; 
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p.inverse());
 Omega = Omega.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict upper part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyUpper>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 } 
 
 //Compute derivative of Omega in rho
 // d(Omega) = d( (I-rho*W)^t * (I-rho*W)) = -( W^t * (I-rho*W) + (W^t * (I-rho*W))^t )
 SparseMatrix<double> bufOmega = -(W.transpose()*(I-rho*W));
 SparseMatrix<double> dOmega = bufOmega + SparseMatrix<double>(bufOmega.transpose());
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 // dC = f(Cm1*dSigma*(Cm1^t)) * C^t
 // where f(A) takes half of the diagonal of A and the strict upper part of A 
 SparseMatrix<double> R = 0.5*(dOmega.cwiseProduct(I));
 R = R + SparseMatrix<double>(dOmega.triangularView<StrictlyLower>());  
 SparseMatrix<double> Z = Cm1*(R*Cm1.transpose());
 SparseMatrix<double> M = Z.triangularView<StrictlyUpper>() + SparseMatrix<double>(Z.triangularView<StrictlyLower>().transpose()) + SparseMatrix<double>((Z.cwiseProduct(I)));
 SparseMatrix<double,Eigen::RowMajor> dC= M*(CC.transpose());
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigma.diagonal().array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num1=(buf1*iW*indep).array() * (Sigmadiagsqrt.replicate(1,k-1)).array();
 VectorXd dsqrtdSigma = dSigmadiag.array() / (2*Sigmadiagsqrt).array();
 MatrixXd num2=(iW*indep).array()*(dsqrtdSigma.replicate(1,k-1)).array();
 MatrixXd num = num1-num2;
 num = num.array()/(Sigmadiag.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 //VectorXd dgdr = VectorXd::Zero(n); 
 
 int i = n-1;
 double Cii = CC.coeff(i,i);
 double Cii2 = Cii*Cii;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)*Cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*Cii*dadb.row(i)+phi_a(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii+lo(i)*dC.coeff(i,i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*dhatadr(i)*VecU(i)*Cii+phi_a(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/Cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = (hatb(i)*phi_b(i)*VecU(i)*Cii*dbdb.row(i)+phi_b(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii+up(i)*dC.coeff(i,i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)*Cii+phi_b(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 } 
 
 for (int i = n-2; i >= 0; --i) { 
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double dgdr = 0;
 double dmdr = 0;
 double mr = 0;
 Cii = CC.coeff(i,i);
 Cii2 = Cii*Cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 dgdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 dgdr += it2.value()*mu(it2.index());}
 dgdr = dgdr - dC.coeff(i,i)*mu(i);
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = zz*Cii+g;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -( (hata(i)*phi_a(i)*dadb.row(i)*VecU(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii + lo(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 mr = phi_a(i)/VecU(i) - Vecg(i);
 dmdr = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i)) -dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2;
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = zz*Cii+g;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = ((hatb(i)*phi_b(i)*dbdb.row(i)*VecU(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii + up(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 mr = -phi_b(i)/VecU(i) - Vecg(i);
 dmdr = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i)) - dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2; ;
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP grad_SAR_UP_AG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using Eigen::Lower;
 using Eigen::StrictlyUpper;
 using Eigen::Upper;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 SparseMatrix<double> Omega = e["eOmega"]; 
 
 double  ep = e["eps"]; 
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p.inverse());
 Omega = Omega.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 CC.prune(ep,1);
 Sigma.prune(ep,1);
 iW.prune(ep,1);
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict upper part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyUpper>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 } 
 
 //Compute derivative of Omega in rho
 // d(Omega) = d( (I-rho*W)^t * (I-rho*W)) = -( W^t * (I-rho*W) + (W^t * (I-rho*W))^t )
 SparseMatrix<double> bufOmega = -(W.transpose()*(I-rho*W));
 SparseMatrix<double> dOmega = bufOmega + SparseMatrix<double>(bufOmega.transpose());
 
 Cm1.prune(ep,1);
 dOmega.prune(ep,1);
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 // dC = f(Cm1*dSigma*(Cm1^t)) * C^t
 // where f(A) takes half of the diagonal of A and the strict upper part of A 
 SparseMatrix<double> R = 0.5*(dOmega.cwiseProduct(I));
 R = R + SparseMatrix<double>(dOmega.triangularView<StrictlyLower>());  
 SparseMatrix<double> Z = Cm1*(R*Cm1.transpose());
 SparseMatrix<double> M = Z.triangularView<StrictlyUpper>() + SparseMatrix<double>(Z.triangularView<StrictlyLower>().transpose()) + SparseMatrix<double>((Z.cwiseProduct(I)));
 SparseMatrix<double,Eigen::RowMajor> dC= M*(CC.transpose());
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigma.diagonal().array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num1=(buf1*iW*indep).array() * (Sigmadiagsqrt.replicate(1,k-1)).array();
 VectorXd dsqrtdSigma = dSigmadiag.array() / (2*Sigmadiagsqrt).array();
 MatrixXd num2=(iW*indep).array()*(dsqrtdSigma.replicate(1,k-1)).array();
 MatrixXd num = num1-num2;
 num = num.array()/(Sigmadiag.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 //VectorXd dgdr = VectorXd::Zero(n); 
 
 int i = n-1;
 double Cii = CC.coeff(i,i);
 double Cii2 = Cii*Cii;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)*Cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*Cii*dadb.row(i)+phi_a(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii+lo(i)*dC.coeff(i,i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*dhatadr(i)*VecU(i)*Cii+phi_a(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/Cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = (hatb(i)*phi_b(i)*VecU(i)*Cii*dbdb.row(i)+phi_b(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii+up(i)*dC.coeff(i,i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)*Cii+phi_b(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 } 
 
 for (int i = n-2; i >= 0; --i) { 
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double dgdr = 0;
 double dmdr = 0;
 double mr = 0;
 Cii = CC.coeff(i,i);
 Cii2 = Cii*Cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 dgdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 dgdr += it2.value()*mu(it2.index());}
 dgdr = dgdr - dC.coeff(i,i)*mu(i);
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = zz*Cii+g;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -( (hata(i)*phi_a(i)*dadb.row(i)*VecU(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii + lo(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 mr = phi_a(i)/VecU(i) - Vecg(i);
 dmdr = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i)) -dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2;
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = zz*Cii+g;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = ((hatb(i)*phi_b(i)*dbdb.row(i)*VecU(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii + up(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 mr = -phi_b(i)/VecU(i) - Vecg(i);
 dmdr = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i)) - dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2; ;
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP lik_SAR_UP_conditional( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 double theta = Rcpp::as<double>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWCL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 int n = W.rows();
 int k = indep.cols();
 double rho = theta;
 
 SparseMatrix<double> A = W;
 double rr = rho; 
 SparseMatrix<double> iW(n,n);
 iW.setIdentity();
 //Compute (I-rho*W)^-1
 iW = iW+rho*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 rr = rho*rr;
 iW = iW + rr*A;
}
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 // Compute the precision matrix Omega (the inverse of Sigma (iWW))
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> Omega = I-rho*W;
 Omega = Omega.transpose()*Omega;
 
 //Set Xstar<-(iW%*%indep)/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
//if ( (rs.array() < 0).any() ){   return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = iW*indep;
 Xstar = Xstar.array()/(rs.replicate(1,k)).array();
 
 Environment speed("package:speedglm");
 Function sglm = speed["speedglm.wfit"];
 Environment stats("package:stats");
 Function binomial = stats["binomial"];
 
 Rcpp::List fit = sglm(_["y"] =dep,_["X"] =Xstar,_["intercept"] =0,_["family"] =binomial("probit"));
 SEXP fitcoef = fit[0];
 VectorXd beta = as<VectorXd>(fitcoef);
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength, _["decreasing"]=1);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Omega=Omega.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Precision matrix Omega
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Omega);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixU();
 VectorXi AMDord = Ch.permutationPinv().indices();
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
//if  ( (CC.diagonal().array() == 0).any() ){return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n); 
 double logP = 0;
 double a = 0;
 double b = 0;
 double U = 0;
 int i = n-1;
 a = (lo[i])*CC.coeff(i,i);
 b = (up[i])*CC.coeff(i,i);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[i] = (R::dnorm(a,0,1,0)-R::dnorm(b,0,1,0))/(U*CC.coeff(i,i));
 logP += log(U);
 
 for (int j = n-2; j >= 0; --j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 a = (lo[j]*CC.coeff(j,j))+g;
 b = (up[j]*CC.coeff(j,j))+g;
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[j] = (((R::dnorm(a,0,1,0)-R::dnorm(b, 0,1,0))/U)-g)/CC.coeff(j,j);
 logP += log(U);
 }
 
 return List::create(Named("l") = logP, Named("beta")=beta);
END_RCPP
}


SEXP lik_SEM_UC( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Eigen::SparseVector;
 using Rcpp::Environment; 
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWFL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 double  ep = e["eps"]; 
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 if  ( abs(rho)>1 ){
    return wrap("Error");
    Rcpp::stop("Parameter rho out of bounds"); }
 
 //Compute the inverse of I-rho*W. If "ord"==0, use exact method, otherwise use Taylor serie appro of order "ord"
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> iW(n,n);
 if (ord==0){
 SparseMatrix<double> IW = I-rho*W;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(IW);
 iW = solver.solve(I);
 iW.prune(ep,1);
 } else {
 SparseMatrix<double> A = W;
 double rr = rho; 
 iW = I + rho*W;
 for(int j = 2; j < ord; ++j)
    {
    A = A * W;
    rr = rho*rr;
    iW = iW + rr*A;
    }
 }
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 //Set Xstar<-indep/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
if ( (rs.array() < 0).any() ){
    return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = indep;
 Xstar = Xstar.array()/(rs.replicate(1,k-1)).array();
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Sigma=Sigma.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Sigma
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Sigma);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixL();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
if  ( (CC.diagonal().array() == 0).any() ){
    return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n);
 VectorXd hata = VectorXd::Zero(n);
 VectorXd hatb = VectorXd::Zero(n);
 VectorXd VecU = VectorXd::Zero(n); 
 VectorXd Vecg = VectorXd::Zero(n); 
 double logP = 0;
 double U = 0;
 int i = 0;
 hata[i] = (lo[i])/CC.coeff(i,i);
 hatb[i] = (up[i])/CC.coeff(i,i);
 U = R::pnorm(hatb[i],0,1,1,0)- R::pnorm(hata[i],0,1,1,0);
 if (U == 0){
    return wrap("Error");}
 VecU[i] = U; 
 mu[i] = (R::dnorm(hata[i],0,1,0)-R::dnorm(hatb[i],0,1,0))/U;
 logP += log(U);
 
 for (int j = 1; j < n; ++j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 Vecg[j] = g;
 hata[j] = (lo[j]-g)/CC.coeff(j,j);
 hatb[j] = (up[j]-g)/CC.coeff(j,j);
 U = R::pnorm(hatb[j],0,1,1,0)- R::pnorm(hata[j],0,1,1,0);
if (U == 0){
    return wrap("Error");}
 VecU[j] = U; 
 mu[j] = (R::dnorm(hata[j],0,1,0)-R::dnorm(hatb[j], 0,1,0))/U;
 logP += log(U);
 }
 e.assign("eiW",iW);
 e.assign("eSigma",Sigma);
 e.assign("eXstar",Xstar);
 e.assign("exb",xb);
 e.assign("eCC",CC);
 e.assign("eVecU",VecU);
 e.assign("eVecg",Vecg);
 e.assign("ehata",hata);
 e.assign("ehatb",hatb);
 e.assign("emu",mu);
 e.assign("elo",lo);
 e.assign("eup",up);
 e.assign("elogP",logP);
 e.assign("efirstorder",firstorder);
 e.assign("eAMDord",AMDord);
 
 return wrap(-logP);
END_RCPP
}


SEXP grad_SEM_UC_FG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Eigen::SparseVector;
 using Rcpp::Environment; 
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict lower part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyLower>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 }  
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf = iW*W*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 SparseMatrix<double> M = Cm1*dSigma*(Cm1.transpose());
 MatrixXd dense_diag_M = (0.5*M.diagonal()).asDiagonal();
 SparseMatrix<double> diag_M = dense_diag_M.sparseView();
 SparseMatrix<double> SLM = M.triangularView<StrictlyLower>();
 SparseMatrix<double,Eigen::RowMajor> dC= CC*(SLM+diag_M);
 
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal(); 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigmadiag.array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num = indep.array()*((dSigmadiag.replicate(1,k-1)).array()); 
 VectorXd den = 2*(Sigmadiag.array()*Sigmadiagsqrt.array()); 
 num = -num.array()/(den.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 VectorXd dgdr = VectorXd::Zero(n);
 
 int i = 0;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*CC.coeff(i,i)-lo(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*CC.coeff(i,i)-up(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 
 for (int i = 1; i < n; ++i) {
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double gdr = 0;
 double cii = CC.coeff(i,i);
 double cii2 = cii*cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 gdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 gdr += it2.value()*mu(it2.index());}
 gdr = gdr - dC.coeff(i,i)*mu(i);
 
 dgdr(i) = (gdr*cii - Vecg(i)*dC.coeff(i,i))/cii2;
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = (zz-g)/cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*cii-lo(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = (zz-g)/cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*cii-up(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 //SEXP sumLU_b = csum(dUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 //double grad_r =  dUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP grad_SEM_UC_AG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Eigen::SparseVector;
 using Rcpp::Environment;  
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 double  ep = e["eps"];
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 CC.prune(ep,1);
 Sigma.prune(ep,1);
 iW.prune(ep,1);
 
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict lower part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyLower>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 } 
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf = iW*W*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 dSigma.prune(ep,1);
 Cm1.prune(ep,1);
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 SparseMatrix<double> M = Cm1*dSigma*(Cm1.transpose());
 MatrixXd dense_diag_M = (0.5*M.diagonal()).asDiagonal();
 SparseMatrix<double> diag_M = dense_diag_M.sparseView();
 SparseMatrix<double> SLM = M.triangularView<StrictlyLower>();
 SparseMatrix<double,Eigen::RowMajor> dC= CC*(SLM+diag_M);
 
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal(); 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigmadiag.array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num = indep.array()*((dSigmadiag.replicate(1,k-1)).array()); 
 VectorXd den = 2*(Sigmadiag.array()*Sigmadiagsqrt.array()); 
 num = -num.array()/(den.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 VectorXd dgdr = VectorXd::Zero(n);
 
 int i = 0;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*CC.coeff(i,i)-lo(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/CC.coeff(i,i);
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*CC.coeff(i,i)-up(i)*dC.coeff(i,i))/(CC.coeff(i,i)*CC.coeff(i,i));
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 
 for (int i = 1; i < n; ++i) {
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double gdr = 0;
 double cii = CC.coeff(i,i);
 double cii2 = cii*cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 gdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 gdr += it2.value()*mu(it2.index());}
 gdr = gdr - dC.coeff(i,i)*mu(i);
 
 dgdr(i) = (gdr*cii - Vecg(i)*dC.coeff(i,i))/cii2;
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = (zz-g)/cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*dadb.row(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = (dadr(i)*cii-lo(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i));
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = (zz-g)/cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) =(hatb(i)*phi_b(i)*VecU(i)*dbdb.row(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i));
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = (dbdr(i)*cii-up(i)*dC.coeff(i,i))/cii2 - dgdr(i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*VecU(i)*dhatbdr(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i));
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 //SEXP sumLU_b = csum(dUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 //double grad_r =  dUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP lik_SEM_UC_conditional( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 double theta = Rcpp::as<double>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWCL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 int n = W.rows();
 int k = indep.cols();
 double rho = theta;
 
 SparseMatrix<double> A = W;
 double rr = theta; 
 SparseMatrix<double> iW(n,n);
 iW.setIdentity();
 //Compute (I-rho*W)^-1
 iW = iW+rho*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 rr = rho*rr;
 iW = iW + rr*A;
}
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 //Set Xstar<-indep/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
//if ( (rs.array() < 0).any() ){return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = indep;
 Xstar = Xstar.array()/(rs.replicate(1,k)).array();
 
 Environment speed("package:speedglm");
 Function sglm = speed["speedglm.wfit"];
 Environment stats("package:stats");
 Function binomial = stats["binomial"];
 
 Rcpp::List fit = sglm(_["y"] =dep,_["X"] =Xstar,_["intercept"] =0,_["family"] =binomial("probit"));
 SEXP fitcoef = fit[0];
 VectorXd beta = as<VectorXd>(fitcoef);
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Sigma=Sigma.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Sigma
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Sigma);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixL();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
//if  ( (CC.diagonal().array() == 0).any() ){return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n);
 double logP = 0;
 double a = 0;
 double b = 0;
 double U = 0;
 int i = 0;
 a = (lo[i])/CC.coeff(i,i);
 b = (up[i])/CC.coeff(i,i);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[i] = (R::dnorm(a,0,1,0)-R::dnorm(b,0,1,0))/U;
 logP += log(U);
 
 for (int j = 1; j < n; ++j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 a = (lo[j]-g)/CC.coeff(j,j);
 b = (up[j]-g)/CC.coeff(j,j);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[j] = (R::dnorm(a,0,1,0)-R::dnorm(b, 0,1,0))/U;
 logP += log(U);
 }
 
 return List::create(Named("l") = logP,Named("beta") = beta);
END_RCPP
}


SEXP lik_SEM_UP( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Rcpp::List; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWFL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 double  ep = e["eps"]; 
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 if  ( abs(rho)>1 ){
    return wrap("Error");
    Rcpp::stop("Parameter rho out of bounds"); }
 
 //Compute the inverse of I-rho*W. If "ord"==0, use exact method, otherwise use Taylor serie appro of order "ord"
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> iW(n,n);
 if (ord==0){
 SparseMatrix<double> IW = I-rho*W;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(IW);
 iW = solver.solve(I);
 iW.prune(ep,1);
 } else {
 SparseMatrix<double> A = W;
 double rr = rho; 
 iW = I + rho*W;
 for(int j = 2; j < ord; ++j)
    {
     A = A * W;
    rr = rho*rr;
    iW = iW + rr*A;
    }
 }
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 //Set Xstar<-indep/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
if ( (rs.array() < 0).any() ){
    return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = indep;
 Xstar = Xstar.array()/(rs.replicate(1,k-1)).array();
 VectorXd xb = Xstar*beta;
 
 // Compute the precision matrix Omega (the inverse of Sigma (iWW))
 SparseMatrix<double> Omega = I-rho*W;
 Omega = Omega.transpose()*Omega;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength, _["decreasing"]=1);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Omega=Omega.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Precision matrix Omega
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Omega);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixU();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
if  ( (CC.diagonal().array() == 0).any() ){
    return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n); 
 VectorXd hata = VectorXd::Zero(n); 
 VectorXd hatb = VectorXd::Zero(n); 
 VectorXd VecU = VectorXd::Zero(n); 
 VectorXd Vecg = VectorXd::Zero(n); 
 double logP = 0;
 double U = 0;
 int i = n-1;
 hata[i] = (lo[i])*CC.coeff(i,i);
 hatb[i] = (up[i])*CC.coeff(i,i);
 U = R::pnorm(hatb[i],0,1,1,0)- R::pnorm(hata[i],0,1,1,0);
 if (U == 0){
    return wrap("Error");}
 VecU[i] = U;
 mu[i] = (R::dnorm(hata[i],0,1,0)-R::dnorm(hatb[i],0,1,0))/(U*CC.coeff(i,i));
 logP += log(U);
 
 for (int j = n-2; j >= 0; --j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 Vecg[j] = g;
 hata[j] = (lo[j]*CC.coeff(j,j))+g;
 hatb[j] = (up[j]*CC.coeff(j,j))+g;
 U = R::pnorm(hatb[j],0,1,1,0)- R::pnorm(hata[j],0,1,1,0);
if (U == 0){
    return wrap("Error");}
 VecU[j] = U;
 mu[j] = (((R::dnorm(hata[j],0,1,0)-R::dnorm(hatb[j], 0,1,0))/U)-g)/CC.coeff(j,j);
 logP += log(U);
 }
 
 e.assign("eiW",iW);
 e.assign("eSigma",Sigma);
 e.assign("eOmega",Omega);
 e.assign("eXstar",Xstar);
 e.assign("exb",xb);
 e.assign("eCC",CC);
 e.assign("eVecU",VecU);
 e.assign("eVecg",Vecg);
 e.assign("ehata",hata);
 e.assign("ehatb",hatb);
 e.assign("emu",mu);
 e.assign("elo",lo);
 e.assign("eup",up);
 e.assign("elogP",logP);
 e.assign("efirstorder",firstorder);
 e.assign("eAMDord",AMDord);
 
 return wrap(-logP);
END_RCPP
}


SEXP grad_SEM_UP_FG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using Eigen::Lower;
 using Eigen::StrictlyUpper;
 using Eigen::Upper;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Eigen::SparseVector;
 using Rcpp::Environment; 
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 SparseMatrix<double> Omega = e["eOmega"]; 
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p.inverse());
 Omega = Omega.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict upper part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyUpper>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 } 
 
 //Compute derivative of Omega in rho
 // d(Omega) = d( (I-rho*W)^t * (I-rho*W)) = -( W^t * (I-rho*W) + (W^t * (I-rho*W))^t )
 SparseMatrix<double> bufOmega = -(W.transpose()*(I-rho*W));
 SparseMatrix<double> dOmega = bufOmega + SparseMatrix<double>(bufOmega.transpose());
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 // dC = f(Cm1*dSigma*(Cm1^t)) * C^t
 // where f(A) takes half of the diagonal of A and the strict upper part of A 
 SparseMatrix<double> R = 0.5*(dOmega.cwiseProduct(I));
 R = R + SparseMatrix<double>(dOmega.triangularView<StrictlyLower>());  
 SparseMatrix<double> Z = Cm1*(R*Cm1.transpose());
 SparseMatrix<double> M = Z.triangularView<StrictlyUpper>() + SparseMatrix<double>(Z.triangularView<StrictlyLower>().transpose()) + SparseMatrix<double>((Z.cwiseProduct(I)));
 SparseMatrix<double,Eigen::RowMajor> dC= M*(CC.transpose());
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal(); 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigmadiag.array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num = indep.array()*((dSigmadiag.replicate(1,k-1)).array()); 
 VectorXd den = 2*(Sigmadiag.array()*Sigmadiagsqrt.array()); 
 num = -num.array()/(den.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 //VectorXd dgdr = VectorXd::Zero(n); 
 
 int i = n-1;
 double Cii = CC.coeff(i,i);
 double Cii2 = Cii*Cii;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)*Cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*Cii*dadb.row(i)+phi_a(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii+lo(i)*dC.coeff(i,i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*dhatadr(i)*VecU(i)*Cii+phi_a(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/Cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = (hatb(i)*phi_b(i)*VecU(i)*Cii*dbdb.row(i)+phi_b(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii+up(i)*dC.coeff(i,i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)*Cii+phi_b(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 } 
 
 for (int i = n-2; i >= 0; --i) { 
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double dgdr = 0;
 double dmdr = 0;
 double mr = 0;
 Cii = CC.coeff(i,i);
 Cii2 = Cii*Cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 dgdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 dgdr += it2.value()*mu(it2.index());}
 dgdr = dgdr - dC.coeff(i,i)*mu(i);
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = zz*Cii+g;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -( (hata(i)*phi_a(i)*dadb.row(i)*VecU(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii + lo(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 mr = phi_a(i)/VecU(i) - Vecg(i);
 dmdr = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i)) -dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2;
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = zz*Cii+g;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = ((hatb(i)*phi_b(i)*dbdb.row(i)*VecU(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii + up(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 mr = -phi_b(i)/VecU(i) - Vecg(i);
 dmdr = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i)) - dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2; ;
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP grad_SEM_UP_AG( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using Eigen::Lower;
 using Eigen::StrictlyUpper;
 using Eigen::Upper;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 VectorXd theta = Rcpp::as<VectorXd>(th);
 Environment e(env);  
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiNFG"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 SparseMatrix<double> iW = e["eiW"]; 
 SparseMatrix<double> Sigma = e["eSigma"]; 
 MatrixXd Xstar = e["eXstar"];
 VectorXd xb = e["exb"];
 SparseMatrix<double,Eigen::RowMajor> CC= e["eCC"];
 VectorXd VecU = e["eVecU"];
 VectorXd Vecg = e["eVecg"];
 NumericVector hata = e["ehata"];
 NumericVector hatb = e["ehatb"];
 VectorXd mu = e["emu"];
 VectorXd lo = e["elo"];
 VectorXd up = e["eup"];
 VectorXi firstorder = e["efirstorder"];
 VectorXi AMDord = e["eAMDord"];
 double  logP = e["elogP"]; 
 
 SparseMatrix<double> Omega = e["eOmega"]; 
 
 double  ep = e["eps"]; 
 
 CC.prune(ep,1);
 Sigma.prune(ep,1);
 iW.prune(ep,1);
 
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 PermutationMatrix<Dynamic,Dynamic> p = p1*p2;
 
 //Permutate rows and columns of matrixes
 W=W.twistedBy(p.inverse());
 iW=iW.twistedBy(p.inverse());
 Sigma = Sigma.twistedBy(p.inverse());
 Omega = Omega.twistedBy(p2.inverse());
 //Permutate rows  of matrixes or vectors
 dep = p.inverse()*dep;
 indep = p.inverse()*indep;
 xb = p.inverse()*xb;
 Xstar = p.inverse()*Xstar;
 
 int n = W.rows();
 int k = theta.size();
 double rho = theta(k-1);
 VectorXd beta = theta.head(k-1);
 
 //Compute C^-1
 // C = D(I+N), where D is diag(C) and N such that D*N=M (M =strict upper part of C)
 // C^-1 = (I+N)^-1*D^-1, where (I+N)^-1=sum_{i=0}^n (-1)^i*N^i
 SparseMatrix<double,Eigen::RowMajor> ICm1(n,n);
 ICm1.setIdentity();
 SparseMatrix<double> Dm1 = (CC.cwiseProduct(ICm1)).cwiseInverse();  
 SparseMatrix<double> N= Dm1*(CC.triangularView<StrictlyUpper>()); 
 
 //Compute (I+N)^-1. If "ord"==0, then it uses exact solve method. Else, use taylor appro of order "ord"
 SparseMatrix<double> I(n,n); 
 I.setIdentity();
 SparseMatrix<double> Cm1(n,n); 
 if (ord==0){
 SparseMatrix<double> iN = I+N;
 Eigen::SparseLU<SparseMatrix<double> > solver;
 solver.compute(iN);
 Cm1 = solver.solve(Dm1); 
 } else {
 SparseMatrix<double> Npower = N;
 SparseMatrix<double> iNm1_sparse = I-N;
 int m1 = -1;
 for(int j = 2; j < ord; ++j)
{
 Npower = Npower * N;
 m1= -1*m1;
 iNm1_sparse = iNm1_sparse + m1*Npower;
}
 Cm1 = iNm1_sparse*Dm1;  
 }
 
 //Compute derivative of Omega in rho
 // d(Omega) = d( (I-rho*W)^t * (I-rho*W)) = -( W^t * (I-rho*W) + (W^t * (I-rho*W))^t )
 SparseMatrix<double> bufOmega = -(W.transpose()*(I-rho*W));
 SparseMatrix<double> dOmega = bufOmega + SparseMatrix<double>(bufOmega.transpose());
 
 Cm1.prune(ep,1);
 dOmega.prune(ep,1); 
 
 //Compute derivative of C in rho <- dC
 //Use BAYESIAN FILTERING AND SMOOTHING - Simo Sarkka - Theorem A.1
 // dC = f(Cm1*dSigma*(Cm1^t)) * C^t
 // where f(A) takes half of the diagonal of A and the strict upper part of A 
 SparseMatrix<double> R = 0.5*(dOmega.cwiseProduct(I));
 R = R + SparseMatrix<double>(dOmega.triangularView<StrictlyLower>());  
 SparseMatrix<double> Z = Cm1*(R*Cm1.transpose());
 SparseMatrix<double> M = Z.triangularView<StrictlyUpper>() + SparseMatrix<double>(Z.triangularView<StrictlyLower>().transpose()) + SparseMatrix<double>((Z.cwiseProduct(I)));
 SparseMatrix<double,Eigen::RowMajor> dC= M*(CC.transpose());
 
 //Compute derivative of Sigma in rho <- dSigma
 // dSigma = iW*W*iW*iWt+iW*iWt*Wt*iWt
 SparseMatrix<double> buf1 = iW*W;
 SparseMatrix<double> buf = buf1*Sigma;
 SparseMatrix<double> dSigma = buf+SparseMatrix<double>(buf.transpose());
 
 //Diagonal of Sigma 
 VectorXd Sigmadiag = Sigma.diagonal(); 
 //Diagonal of Sigma sqrt
 VectorXd Sigmadiagsqrt = Sigmadiag.array().sqrt();
 //Diagonal of dSigma
 VectorXd dSigmadiag = dSigma.diagonal();
 
 //Compute dxb/dr
 MatrixXd num = indep.array()*((dSigmadiag.replicate(1,k-1)).array()); 
 VectorXd den = 2*(Sigmadiag.array()*Sigmadiagsqrt.array()); 
 num = -num.array()/(den.replicate(1,k-1)).array();
 VectorXd dxbdr = num*beta;
 
 NumericVector phi_a = dnorm(hata,0,1,0);
 NumericVector phi_b = dnorm(hatb,0,1,0);
 NumericVector Phi_a = pnorm(hata,0,1,1,0);
 NumericVector Phi_b = pnorm(hatb,0,1,1,0);
 
 MatrixXd dUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dLUdb = MatrixXd::Zero(n,k-1);
 MatrixXd dmudb = MatrixXd::Zero(n,k-1);
 MatrixXd dadb = MatrixXd::Zero(n,k-1);
 MatrixXd dbdb = MatrixXd::Zero(n,k-1);
 
 VectorXd dUdr = VectorXd::Zero(n);
 VectorXd dLUdr = VectorXd::Zero(n);
 VectorXd dadr = VectorXd::Zero(n);
 VectorXd dbdr = VectorXd::Zero(n);
 VectorXd dhatadr = VectorXd::Zero(n);
 VectorXd dhatbdr = VectorXd::Zero(n);
 VectorXd dmudr = VectorXd::Zero(n);
 //VectorXd dgdr = VectorXd::Zero(n); 
 
 int i = n-1;
 double Cii = CC.coeff(i,i);
 double Cii2 = Cii*Cii;
 if (dep(i) == 0){
 //Part for derivative on beta
 dadb.row(i) = Xstar.row(i)*Cii;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -(hata(i)*phi_a(i)*VecU(i)*Cii*dadb.row(i)+phi_a(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Part for derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii+lo(i)*dC.coeff(i,i);
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 dmudr(i) = -(hata(i)*phi_a(i)*dhatadr(i)*VecU(i)*Cii+phi_a(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 }
 else{
 //Derivative on beta
 dbdb.row(i) = Xstar.row(i)/Cii;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = (hatb(i)*phi_b(i)*VecU(i)*Cii*dbdb.row(i)+phi_b(i)*Cii*dUdb.row(i))/(VecU(i)*VecU(i)*Cii2);
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii+up(i)*dC.coeff(i,i);
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 dmudr(i) = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)*Cii+phi_b(i)*(dUdr(i)*Cii+VecU(i)*dC.coeff(i,i)))/(VecU(i)*VecU(i)*Cii2);
 } 
 
 for (int i = n-2; i >= 0; --i) { 
 VectorXd g = VectorXd::Zero(k-1);
 VectorXd zz = Xstar.row(i);
 double dgdr = 0;
 double dmdr = 0;
 double mr = 0;
 Cii = CC.coeff(i,i);
 Cii2 = Cii*Cii;
 
 SparseVector<double> vec = CC.row(i);
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += (it.value()*dmudb.row(it.index()));
 dgdr +=  it.value()*dmudr(it.index());}
 
 SparseVector<double> rowdc = dC.row(i);
 for (SparseVector<double>::InnerIterator it2(rowdc); it2; ++it2){
 dgdr += it2.value()*mu(it2.index());}
 dgdr = dgdr - dC.coeff(i,i)*mu(i);
 
 if (dep(i) == 0){
 //Derivative on beta
 dadb.row(i) = zz*Cii+g;
 dUdb.row(i) = -phi_a(i)*dadb.row(i);
 dLUdb.row(i) = dUdb.row(i)/(1-Phi_a(i));
 dmudb.row(i) = -( (hata(i)*phi_a(i)*dadb.row(i)*VecU(i)+phi_a(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dadr(i) = dxbdr(i);
 dhatadr(i) = dadr(i)*Cii + lo(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = -phi_a(i)*dhatadr(i); 
 dLUdr(i) = dUdr(i)/(1-Phi_a(i)); 
 mr = phi_a(i)/VecU(i) - Vecg(i);
 dmdr = -(hata(i)*phi_a(i)*VecU(i)*dhatadr(i)+phi_a(i)*dUdr(i))/(VecU(i)*VecU(i)) -dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2;
 }
 else{
 //Derivative on beta 
 dbdb.row(i) = zz*Cii+g;
 dUdb.row(i) = phi_b(i)*dbdb.row(i);
 dLUdb.row(i)= dUdb.row(i)/Phi_b(i);
 dmudb.row(i) = ((hatb(i)*phi_b(i)*dbdb.row(i)*VecU(i)+phi_b(i)*dUdb.row(i))/(VecU(i)*VecU(i)) -g)/Cii;
 //Derivative on rho
 dbdr(i) = dxbdr(i);
 dhatbdr(i) = dbdr(i)*Cii + up(i)*dC.coeff(i,i) + dgdr;
 dUdr(i) = phi_b(i)*dhatbdr(i); 
 dLUdr(i) = dUdr(i)/Phi_b(i); 
 mr = -phi_b(i)/VecU(i) - Vecg(i);
 dmdr = (hatb(i)*phi_b(i)*dhatbdr(i)*VecU(i)+phi_b(i)*dUdr(i))/(VecU(i)*VecU(i)) - dgdr;
 dmudr(i) = (dmdr*Cii-mr*dC.coeff(i,i))/Cii2; ;
 } 
 }
 Environment Mat("package:Matrix");
 Environment base("package:base");
 Function csum = Mat["colSums"];
 SEXP sumLU_b = csum(dLUdb);
 VectorXd grad_b = as<VectorXd>(sumLU_b);
 
 double grad_r =  dLUdr.sum();
 Function comb = base["c"];
 SEXP gradd= comb(grad_b,grad_r);
 VectorXd grad = as<VectorXd>(gradd);  
 grad = grad/logP;
 
 return wrap(grad);
END_RCPP
}


SEXP lik_SEM_UP_conditional( SEXP th, SEXP env ){
BEGIN_RCPP

 using Eigen::PermutationMatrix;
 using Eigen::Dynamic;
 using Eigen::SparseMatrix;
 using Eigen::StrictlyLower;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using Eigen::VectorXi;
 using Eigen::Map; 
 using Rcpp::Environment; 
 using Eigen::SparseVector;
 typedef Eigen::Triplet<double> T;
 
 double theta = Rcpp::as<double>(th);
 Environment e(env); 
 SparseMatrix<double> W = e["WW"];
 int ord = e["appiWCL"];
 VectorXd dep = e["de"];
 MatrixXd indep = e["ind"];
 
 int n = W.rows();
 int k = indep.cols();
 double rho = theta;
 
 SparseMatrix<double> A = W;
 double rr = rho; 
 SparseMatrix<double> iW(n,n);
 iW.setIdentity();
 //Compute (I-rho*W)^-1
 iW = iW+rho*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 rr = rho*rr;
 iW = iW + rr*A;
}
 //Compute (I-rho*W)^-1  * t( (I-rho*W)^-1 ) i.e. Sigma
 SparseMatrix<double> Sigma = iW*iW.transpose();
 
 // Compute the precision matrix Omega (the inverse of Sigma (iWW))
 SparseMatrix<double> I(n,n);
 I.setIdentity();
 SparseMatrix<double> Omega = I-rho*W;
 Omega = Omega.transpose()*Omega;
 
 //Set Xstar<-(iW%*%indep)/sqrt(diag(Sigma))
 VectorXd rs = Sigma.diagonal();
//if ( (rs.array() < 0).any() ){return wrap("Error");}
 rs = rs.array().sqrt();
 MatrixXd Xstar = indep;
 Xstar = Xstar.array()/(rs.replicate(1,k)).array();
 
 Environment speed("package:speedglm");
 Function sglm = speed["speedglm.wfit"];
 Environment stats("package:stats");
 Function binomial = stats["binomial"];
 
 Rcpp::List fit = sglm(_["y"] =dep,_["X"] =Xstar,_["intercept"] =0,_["family"] =binomial("probit"));
 SEXP fitcoef = fit[0];
 VectorXd beta = as<VectorXd>(fitcoef);
 VectorXd xb = Xstar*beta;
 
 VectorXd lo(n),up(n),intLength(n);
 for (int i = 0 ; i<n ; i++){
 if (dep[i] == 0){
 lo[i] = xb[i];
 up[i] = R_PosInf;  
 intLength[i] = 1-R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 else{
 lo[i] = R_NegInf;
 up[i] = xb[i];  
 intLength[i] = R::pnorm(xb[i]/rs[i],0,1,1,0);
 }
 }
 
 Environment base("package:base");
 Function ordF = base["order"];
 SEXP fo = ordF(intLength, _["decreasing"]=1);
 VectorXi firstorder = as<VectorXi>(fo);
 firstorder = firstorder.array()-1;
 
 //Reordering of lo, up and Sigma via Permutation matrix
 PermutationMatrix<Dynamic,Dynamic> p1(firstorder);
 Omega=Omega.twistedBy(p1.inverse());
 lo = p1.inverse()*lo;
 up = p1.inverse()*up;
 
 //Computing Cholesky Decomposition of Precision matrix Omega
 using Eigen::LLT;
 typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;
 const SpChol Ch(Omega);
 SparseMatrix<double,Eigen::RowMajor> CC = Ch.matrixU();
 VectorXi AMDord = Ch.permutationPinv().indices();
 
 //Reorder lo and up according to AMDord
 PermutationMatrix<Dynamic,Dynamic> p2(AMDord);
 lo = p2.inverse()*lo;
 up = p2.inverse()*up;
//if  ( (CC.diagonal().array() == 0).any() ){return wrap("Error");}
 
 // Computing the loglikelihood with Genz Trinh method of uniconditional approximation
 using Eigen::SparseVector;
 VectorXd mu = VectorXd::Zero(n); 
 double logP = 0;
 double a = 0;
 double b = 0;
 double U = 0;
 int i = n-1;
 a = (lo[i])*CC.coeff(i,i);
 b = (up[i])*CC.coeff(i,i);
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[i] = (R::dnorm(a,0,1,0)-R::dnorm(b,0,1,0))/(U*CC.coeff(i,i));
 logP += log(U);
 
 for (int j = n-2; j >= 0; --j) {
 SparseVector<double> vec = CC.row(j);
 double g = 0;
 for (SparseVector<double>::InnerIterator it(vec); it; ++it){
 g += mu[it.index()]*it.value();}
 a = (lo[j]*CC.coeff(j,j))+g;
 b = (up[j]*CC.coeff(j,j))+g;
 U = R::pnorm(b,0,1,1,0)- R::pnorm(a,0,1,1,0);
//if (U == 0){return wrap("Error");}
 mu[j] = (((R::dnorm(a,0,1,0)-R::dnorm(b, 0,1,0))/U)-g)/CC.coeff(j,j);
 logP += log(U);
 }
 
 return List::create(Named("l") = logP, Named("beta")=beta);
END_RCPP
}



