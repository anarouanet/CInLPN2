// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <stdlib.h> /* srand, rand */
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
//===========================================================================================================
//function print matrix
int f_mat_print(arma::mat& B)
{
  int J = B.n_cols;
  int I = B.n_rows;
  printf("\n \n===========================\n affichage de la matrice\n============================\n \n");
  printf(" dimension ::::%d * %d\n\n", I,J);
  for(int i =0; i<I; i++)
  {
    printf("| ");
    for(int j =0; j<J; j++)
    {
      printf("%.4f ",B(i,j));
    }
    printf("| \n");
  }
  printf("\n \n===========================\n fin de l'affichage\n============================\n \n");
  return 1;
}

//===========================================================================================================
// function to test if a matrix is inversible
arma::mat f_inv_mat(arma::mat& B)
{
  int g = 0;
  mat C;
  try{

    C = inv_sympd(B);
    //C = inv(B);

  }
  catch(const std::runtime_error& e){
    g=1;
    printf("Matrix B is not inversible : %d \n", g);
    C = zeros(B.n_cols,B.n_cols);
  }
  //   f_mat_print(C);
  if(g==1){
    printf("Matrix B is not inversible : %d \n", g);
  }

  return (C);
}


//======================================================================
//' Function that vectorises a matrix by rows
//' 
//' @param M a matrix
//' @export
// [[Rcpp::export]]
arma::vec vectorise(arma::mat& M){
  int n_c = M.n_cols;
  int n_r = M.n_rows;
  vec v = zeros(n_r*n_c);
  int pp=0;
  for(int i =0; i<n_r; i++){
    for(int k =0; k<n_c; k++){
      v(pp) = M(i,k);
      pp++;
    }
  }
  return(v);
}

//===========================================================================================================
arma::vec InnerProd(arma::vec v1, arma::vec v2){
  int n = v1.size();
  vec ip = zeros(n);
  for(int i=0; i < n; i++){
    ip[i] = v1[i]*v2[i];
  }
  return (ip);
}


//===========================================================================================================
//' Function that creates a K-block diagonal matrix with non diagonal element fixed to 0
//' 
//' @param Kvector a vector of length K
//' 
//' @return a diagonal matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat KmatDiag(arma::vec& Kvector){
  // Kvector : K-vector of paramters used to construct the matrix Kmat
  for(int i=0; i < (int)Kvector.size(); i++){
    Kvector[i] = pow(Kvector[i],2);
  }
  mat Kmat = diagmat(Kvector);
  return (Kmat);
}


//'  Function that computes a symetric D matric from it Cholesky L
//'  
//' @param q an integer
//' @param qvector a vector of length q
//' 
//' @return a symetric positive define matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat DparChol(int q, arma::vec& qvector){
  // qvector : vecteur de param?tres de taille K pour la construction de la matrice
  mat L = zeros<mat>(q,q);
  int pc=0 ;//  pour index? un param?tre dans le vecteur paraChol
  for( int j=0;j< q; ++j){
    for(int i=j; i<q; ++i){
      L(i,j) = qvector[pc];
      pc+=1;
    }
  }
  mat D = L*L.t();
  return (D);
}

//===========================================================================================================
/* ***********************************
Function aijt
Description:
aijt a function that computes the coefficient a_ij of the transition matrix A at time t from
model.matrix and vector of parameters
*/

double aijt(int t, arma::vec alpha_ijl,  arma::mat modA_mat){
  //modA_mat: model.matrix
  // remember that, length(alpha_ijl) = ncol(modA_mat)
  double a_ij_t = as_scalar(modA_mat(span(t,t),span(0,modA_mat.n_cols-1))*alpha_ijl);
  return(a_ij_t);
}



//'  Function that computes all coefficients of the transition matrix at time t
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t an integer indicating the time at which coefficients are computed 
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a vector of elements of the transition matrix
//' @export
//' 
// [[Rcpp::export]]
arma::vec vecaijt( int K, int t, arma::vec& vec_alpha_ij, arma::mat& modA_mat){
  // K = number of outcomes
  // L = taille de la matrice modA_mat.
  int L = modA_mat.n_cols;
  arma::vec vec_a_ij_t = zeros<vec>(K*K);
  int pp = 0; // index of loop
  for( int k = 0 ; k < K*K; k++){
    vec_a_ij_t(k) = aijt(t, vec_alpha_ij(span(pp,pp+L-1)), modA_mat);
    pp+=L;
  }
  return(vec_a_ij_t);

}

//===========================================================================================================

//'  Function that  constructs the transition matrix at time t
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t an integer indicating the time at which coefficients are computed 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return Id + DeltaT*A where A is the transition matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat ConstrA(int K, int t, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat){

  vec a_ij_t = vecaijt(K, t, vec_alpha_ij, modA_mat);
  mat A = zeros<mat>(K,K);
  vec vect = ones(K);
  mat Id = diagmat(vect);
  int p=0;

  for(int i = 0;i< K;i++){
    for(int j = 0;j< K;j++){
      A(i,j) = a_ij_t(p);
      p++;
    }
  }
  return(Id + DeltaT*A); // Atild = ((1/DeltaT)*Id + A)
}

//=======================================================================================================

//'  Function that  creates a matrix K,K*(max(tau_i)-1) containing  sub-matrices (A_t)_(t=0,tau_i-1)
//'  
//' @param K an integer representing the size of K*K matrix
//' @param tau_i vector of integers indicating times at which coefficients are computed 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transition matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix containing matrix of form Id + DeltaT*A
//' @export
//' 
// [[Rcpp::export]]
arma::mat GmatA0totaui(int K, arma::vec& vec_alpha_ij, arma::vec& tau_i, double DeltaT,  arma::mat modA_mat) {
  int siz_tau_i = tau_i.size();

  mat G_mat_A_0_to_tau_i = zeros<mat>(K,K*siz_tau_i);
  for (int i=0; i< siz_tau_i; i++) {
    G_mat_A_0_to_tau_i(span(0,K-1),span(K*i,K*i+K-1)) = ConstrA(K, tau_i[i], DeltaT, vec_alpha_ij, modA_mat);
  }
  return (G_mat_A_0_to_tau_i);
}

//=======================================================================================================

//'  Function that  computes the product of A(t) for t1 to t2
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t1 indicates the started discretized time
//' @param t2 indicates the ended discretized time
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat ProdA(int K, int t2, int t1, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat){
  // t2 = end time
  // t1 = start time
  // t1 <= t2
  mat pA = ConstrA(K,t1, DeltaT, vec_alpha_ij, modA_mat);
  for(int i = (t1+1); i <= t2; i++){
    pA = pA*ConstrA(K, i, DeltaT, vec_alpha_ij, modA_mat);
  }
  return(pA);
}


//===========================================================================================
//'  Function that  creates a big matrix containing  Prod(A_t)t=t_ini,tau.
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t_ini indicates the started discretized time
//' @param tau vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat GmatprodAstotau( int K, arma::vec& vec_alpha_ij, arma::vec& tau,
                           int t_ini, double DeltaT, arma::mat modA_mat) {
  // t_ini = t initial
  int nb_t = tau.size()-t_ini;
  int ii = 0; //index de boucle
  mat G_mat_prod_A_s_to_tau = zeros<mat>(K,K*nb_t);
  for (int i=t_ini; i< (int)tau.size(); i++) {
    G_mat_prod_A_s_to_tau(span(0,(K-1)),span(ii,(ii+K-1))) = ProdA(K, tau[i], t_ini, DeltaT, vec_alpha_ij, modA_mat);
    ii +=K;
  }
  return (G_mat_prod_A_s_to_tau);
}

//===========================================================================================
//'  Function that  creates a big matrix ts_G_mat_prod_A_0_to_tau containing  Prod(A_t)t=0,tau.
//'  
//' @param K an integer representing the size of K*K matrix
//' @param tau vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat tsGmatprodA0totau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, double DeltaT, arma::mat modA_mat) {
  //Creation of matrix G_mat_prod_A_0_to_tau that contains all products  A(j) from t_j a Tmax: t_j \in 0, Tmax
  int T = tau.size();
  mat ts_G_mat_prod_A_0_to_tau = zeros<mat>(K*T,K*T);
  for(int s=0; s<T; s++){
    ts_G_mat_prod_A_0_to_tau(span(K*s,K*(s+1)-1),span(K*s,K*T-1)) = GmatprodAstotau( K, vec_alpha_ij, tau, s, DeltaT, modA_mat);
  }
  return(ts_G_mat_prod_A_0_to_tau);
}

//===========================================================================================================
/* ***********************************
Function matHit
Description:
matHit a function that construct the observation matrix H_i(t) at time t.
Dependances: none
*/

arma::mat matHit(arma::vec X_i_t){
  // X_i_t : vector of observation
  int K = (int)X_i_t.size();
  int l=0; //loop variable
  int p=0; //loop variable

  for( int i=0; i< K; i++){
    if(!isnan(X_i_t[i])){
      p++;
    }
  }
  mat H_i_t = zeros(p,K); //initialisation of matrix H_i(t)
  // filling ofa matrix H_i(t)
  for( int i=0; i< K; i++){
    if(!isnan(X_i_t(i))){
      H_i_t(l,i) =1;
      l++;
    }
  }
  if(p==0){
    mat H_i_t = zeros(K,K);
  }
  return(H_i_t);
}


// function that provides the vector of times (from tau) that correspond to observed outcomes values (X_ik)
arma::vec matTik(arma::vec X_ik, arma::vec tau){
  // X_ik : vector of observation

  int T = (int)X_ik.size();
  int p=0; //loop variable

  for( int i=0; i< T; i++){
    if(!isnan(X_ik(i))){
      p++;
    }
  }

  vec T_i_k = zeros<vec>(p); //initialisation of matrix H_i(t)
  int ii=0;
  for( int i=0; i< T; i++){
    if(!isnan(X_ik(i))){
      T_i_k(ii)=tau(i);
      ii++;
    }
  }
  
  return(T_i_k);
}
//===========================================================================================
//' After vectorising the vzector Yi, this function returns
//' a vector indicating missing values : 1 = observed value, 0 = missing value
//' 
//' @param Yi a matrix with possibly NAs
//' 
//' @return a vector of elements 0,1
//' @export
//' 
// [[Rcpp::export]]
arma::vec compoYiNA(arma::mat& Yi){
  int K = Yi.n_cols;
  int ni = Yi.n_rows;
  vec compo_Yi_NA = ones(K*ni);
  int p=0;
  for( int i=0; i< ni; i++){
    for( int j=0; j< K; j++){
      if(isnan(Yi(i,j))){
        compo_Yi_NA(p) =0;
      }
      p++;
    }
  }
  return(compo_Yi_NA);
}

//===========================================================================================
//' Function that returns Yi (a vector) without NAs values
//' 
//' @param Yi a vector with possibly NAs
//' 
//' @return a vector 
//' @export
//' 
// [[Rcpp::export]]
arma::vec YiwoNA(arma::vec Yi){
  int Ti = Yi.size();
  mat NAs = Yi.elem(find_nonfinite(Yi)); // repere NAs values in the vector Yi
  vec Yi_wo_NA = zeros(Ti-NAs.n_rows*NAs.n_cols);
  int p=0;
  for( int i=0; i< Ti; i++){
    if(!isnan(Yi(i))){
      Yi_wo_NA(p) =Yi(i);
      p++;
    }

  }
  return(Yi_wo_NA);
}
//===========================================================================================================

/* ***********************************
Function matNui
Description:
matNui a function that construct the matrix nu_t_j, the expectation of processes at time t_j
Dependance: none
*/

//===========================================================================================
//'  Function that constructs the matrix nu_t_j, the expectation of the processes at time t_j
//'  
//' @param nD an integer indicating the number of processes
//' @param tau_i vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matrix containing at time t
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat matNui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                 arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i){
  // matNu_i : matrix of size xi.n_rows*nD containing expectation of processes
  int n_cols_xi = xi.n_cols;
  int mi=tau_i.size(); // number of observations
  int T = max(tau_i)+1;
  mat matNu_i = zeros(mi,nD);
  mat Mu_t = zeros(nD,1);
  int i = 0;
  for(int t=0; t< T; t++){
    if(t==0){
      Mu_t = x0i*alpha_mu0;
    }
    else{
      Mu_t = DeltaT*xi(span(t*nD,(t+1)*nD-1), span(0,n_cols_xi-1))*alpha_mu
      + G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),nD*(t-1)+nD-1))*Mu_t;
    }
    if(t ==(int)tau_i(i)){
      matNu_i(span(i,i), span(0,nD-1)) = Mu_t.t();
      i++;
    }
  }
  return (matNu_i);
}

//===========================================================================================
//'  Function that constructs the matrix matYtild_ui (Lambda_i|ui), the expectation of the processes at time t_j given ui
 //'  
 //' @param nD an integer indicating the number of processes
 //' @param tau_i vector of integers indicating times 
 //' @param DeltaT double that indicates the discretization step
 //' @param x0i model.matrix for baseline's submodel
 //' @param xi model.matrix for change's submodel
 //' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
 //' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
 //' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
 //' matrix containing at time t
 //' @param zi model.matrix for change's submodel random effects
 //' @param ui matrix containing the random intercepts on the baseline level
 //' @param vi matrix containing the random intercepts on the change model
 //' 
 //' @return a matrix
 //' @export
 //' 
 // [[Rcpp::export]]
 arma::mat matYtild_ui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                  arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi, arma::colvec& ui, arma::colvec& vi){
   // matNu_i : matrix of size xi.n_rows*nD containing expectation of processes
   int n_cols_xi = xi.n_cols;
   int n_cols_zi = zi.n_cols;
   
   int mi=tau_i.size(); // number of observations
   int T = max(tau_i)+1;
   mat matNu_i = zeros(mi,nD);
   mat Mu_t = zeros(nD,1);
   int i = 0;
   for(int t=0; t< T; t++){
     if(t==0){
       Mu_t = x0i*alpha_mu0 + ui;
     }
     else{
       Mu_t = DeltaT*(xi(span(t*nD,(t+1)*nD-1), span(0,n_cols_xi-1))*alpha_mu + zi(span(t*nD,(t+1)*nD-1), span(0,n_cols_zi-1))*vi)
       + G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),nD*(t-1)+nD-1))*Mu_t;
     }
     if(t ==(int)tau_i(i)){
       matNu_i(span(i,i), span(0,nD-1)) = Mu_t.t();
       i++;
     }
   }
   return (matNu_i);
 }

/* ***********************************
Function f_Yi_r_NA_by0
Description:
f_Yi_r_NA_by0 a function that replace NA by 0.0. Just for computatinal need
Dependences : none
*/

//===========================================================================================
//' Function that replaces NAs by 0.0 just for computatinal need
//' 
//' @param Yi a matrix with possibly NAs
//' 
//' @return a matrix 
//' @export
//' 
// [[Rcpp::export]]
arma::mat f_Yi_r_NA_by0(arma::mat& Yi){
  int K = Yi.n_cols;
  int ni = Yi.n_rows;
  vec compo_Yi_NA = ones(K*ni);
  int p=0;
  for( int i=0; i< ni; i++){
    for( int j=0; j< K; j++){
      if(isnan(Yi(i,j))){
        Yi(i,j) = 0.0;
      }
      p++;
    }
  }
  return(Yi);
}

//===========================================================================================================
/* ***********************************
Function YiNui
Description:
YiNui a function that compute the difference est une fonction (mat_Yi - mat_Nu_i)
delete missing values (NA) and return a vector
Dependances: matNui
*/

//===========================================================================================
//' Function that computes the difference (mat_Yi - mat_Nu_i), deletes missing values (NAs) and 
//' returns a vector. mat_Yi is the outcome and mat_Nu_i is the expectation
//'  
//' @param nD an integer indicating the number of processes
//' @param matrixP a matrix that matches markers to latent processes
//' @param tau a vector of integers indicating times 
//' @param tau_i a vector of integers indicating times for individual i
//' @param DeltaT double that indicates the discretization step
//' @param Yi a matrix of the outcomes
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matric containing at time t
//' 
//' @return a vector
//' @export
//' 
// [[Rcpp::export]]
arma::vec YiNui(int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i){
  // Yi : matrice of observation of subject i
  mat Nu_cp = matNui(nD, tau, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
  mat Nu_cp_i = zeros(Yi.n_rows,nD);
  
  for(int i=0; i<(int)tau_i.size(); i++){
    Nu_cp_i.row(i) = Nu_cp.row(tau_i(i));
  }

  // Nu_cp : Expectation of overall times
  mat M =  Yi - Nu_cp_i*matrixP.t();

  return(YiwoNA(vectorise(M)));
}

//===========================================================================================================
/* ***********************************
Function normalCDF
Description:
 normalCDF a function that computes the cumulative distribution function of a standard normal
*/

double normalCDF(double value, bool lower_tail = true)
{
  //double M_SQRT1_2=sqrt(0.5);
  double proba = 0.5 * erfc(-value * M_SQRT1_2);

  if(!lower_tail)
    proba = 1- proba;         
  
  return proba;
}

//===========================================================================================================
/* ***********************************
Function f_marker
Description:
 f_marker a function that computes f(Yi|Lambda_i) with cumulative probit function
delete missing values (NA) and returns a double
Dependances: matNui
*/

bool Isnotnan(double i) {
  return !isnan(i);
}

//===========================================================================================
//' Function that computes the difference (mat_Yi - mat_Nu_i), delates missing values (NAs) and 
//' returns a vector. mat_Yi is the outcomes and mat_Nu_i is the expectation
//'  
//' @param Lambdai a matrix of dimension nT x nD containing the sampled lambda_i
//' @param nD an integer indicating the number of processes
//' @param matrixP a matrix that matches markers to latent processes
//' @param tau a vector of integers indicating times 
//' @param tau_i a vector of integers indicating times for individual i
//' @param DeltaT double that indicates the discretization step
//' @param Ytildi vector of individual transformed outcomes
//' @param YtildPrimi vector of individual transformed derivatives
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param paraSig variances of marker-specific measurement errors
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matric containing at time t
//' @param paraEtha2 transformation parameters
//' @param if_link: link function indicator, 0 if linear, 1 if splines, 2 if thresholds
//' @param zitr: minY and maxY of observed ordinal Y
//' @param ide indicator if the values between zitr(0) and zitr(1) are observed in Y 
//' @param paras_k: number of parameters for link function for each marker k
//' @param K2_lambda_t: vector indicating to which latent process corresponds each value of Lambdai
//' @param K2_lambda: vector indicating to which latent process corresponds each marker
//' 
//' @return a double
//' @export
//' 

// [[Rcpp::export]]
double f_marker(arma::mat& Lambdai, int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Ytildi, arma::mat& YtildPrimi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::vec& paraSig, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::colvec& paraEtha2, arma::vec& if_link, arma::colvec& zitr, 
                arma::mat& ide, arma::vec& paras_k, arma::vec& K2_lambda_t, arma::vec& K2_lambda){

  int K = matrixP.n_rows;
  double logvrais = 0.0;

  int param = 0;
  int K_t = 0;

  for(int k=0; k<K; k++){
    // Retrieve part of Ytildi and Lambdai that corresponds to marker k
    vec Ytildik=Ytildi.col(k);

    // if(!std::all_of(Ytildi.col(k).begin(), Ytildi.col(k).end(), Isnotnan) ){
    //   int ind=0;
    //   
    //   for(int b = 0 ; b < Ytildi.n_rows; b++){
    //     if(!isnan(Ytildi(b,k))){
    //       Ytildik(ind)=Ytildi(b,k);
    //       ind++;
    //     }
    //   }
    //   Ytildik =Ytildik.head(ind);
    // }else{
    //   Ytildik=Ytildi.col(k);
    // }

    // vec K2_lambda_t[j]: which latent process linked to the jth observation in lambda
    // K2_lambda[k]: which latent process linked to the kth markers
    int nT = K2_lambda_t.size() ;
    vec Lambdai_k(Ytildik.size());
    int ind=0;
    for(int b = 0 ; b < nT; b++){
      if(K2_lambda_t(b)==K2_lambda(k)){
        Lambdai_k(ind)=Lambdai(b,k);
        ind++;
      }
    }


     if (if_link[k]==2){// thresholds
       
       double inf;
       double sup;
       bool lower=true;
       double vrais=1;
       
       vec params_thresholds(paras_k[k]);
       params_thresholds[0] = paraEtha2[param];
       
       for(int r=1; r<paras_k[k]; r++){
         params_thresholds[r] = params_thresholds[r-1] + pow(paraEtha2[param + r],2);
       }
      
       for(int j=0; j<(int)tau_i.size(); j++){

         if (Ytildik(j,k)==zitr[K_t*2]){ // if Ytildi == minY
           double gamma0 =  params_thresholds[0] - Lambdai_k[j]; 
           
           vrais  *= normalCDF(gamma0, lower);

         }else{
           sup = params_thresholds[0];
           inf = sup;

           int pp=0;
           for(int p=0; p<(int)(zitr[K_t*2 + 1]-zitr[K_t*2]-1); p++){
             
             if(ide(p)==1){
               sup = params_thresholds[pp+1];

               if(Ytildik(j,k) == (zitr[K_t*2]+p)){
                 double gamma1 = sup - Lambdai_k[j];
                 double gamma0 = inf - Lambdai_k[j];
                 vrais  *= (normalCDF(gamma1) - normalCDF(gamma0));
                 
                 //cout << " p "<<p << " pp "<<pp << endl;
                 //cout << " sup "<<sup << " gamma1 "<<gamma1 << " normalCDF(gamma1) "<<normalCDF(gamma1)<<endl;
                 //cout << " inf "<<inf << " gamma0 "<<gamma0 << " normalCDF(gamma0) "<<normalCDF(gamma0)<<endl;
                 //cout << " vrais "<<vrais << " CDF1-CDF0 "<<(normalCDF(gamma1) - normalCDF(gamma0)) << endl;

                 // fout << "Ytildij: "<< Ytildi(j,k) << " zitr[K_t*2]+p: "<<zitr[K_t*2]+p 
                 //      << " sup "<<sup << " inf "<<inf << " param "<< param
                 //      << " params_thresholds[p+1] "<< params_thresholds[p+1] << " vrais "<<vrais
                 //      << " temp "<<temp<< " gamma1 "<<gamma1<< " gamma0 "<<gamma0
                 //      <<" Ngamma1 "<<normalCDF(gamma1)<< " Ngamma0 "<<normalCDF(gamma0)<<endl;
               }
               inf = sup;
               pp += 1;
             }
           }
           
           if (Ytildik(j,k)==zitr[K_t*2 + 1]){ // if Ytildi == maxY
             double gamma0 =  sup - Lambdai_k[j]; 
             vrais  *= (1-normalCDF(gamma0));
           }
         }
       }
      K_t +=1;
     logvrais += log(vrais);

     }else if (if_link[k]==0){// linear

       // ##### computering of the likelihood ##########################
       double log_det= Lambdai_k.size()*log(paraSig[k]);
       //double log_Jac_Phi=log(params_lin[0])*Lambdai_k.size();
       double log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi.col(k)))));
       vec y=Ytildik-Lambdai_k;
       
//       vrais = -0.5*(sum(Lambdai_k.size())*log(2*M_PI) + log_det + as_scalar(y.t()*y/paraSig[k])) ;
       logvrais += -0.5*(Lambdai_k.size()*log(2*M_PI) + log_det + as_scalar(y.t()*y/paraSig[k])) + log_Jac_Phi;

       //dmvnorm(Ytilde, mean = Ytild[i+1,], sigma = sigmae2*diag(length(mu))) / a^(length(mu))
     }
     param += paras_k[k];
   }

   // int mk=ncol(matrixP); //number of ordinal markers
  // for(int i=0; i<(int)tau_i.size(); i++){
  //   for(int j=0; j<mk; i++){
  // 
  //   }
  // }
  // -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i))
  // 
  // 
  // // Yi : matrice of observation of subject i
  // mat Nu_cp = matNui(nD, tau, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
  // mat Nu_cp_i = zeros(Yi.n_rows,nD);
  // for(int i=0; i<(int)tau_i.size(); i++){
  //   Nu_cp_i.row(i) = Nu_cp.row(tau_i(i));
  // }
  // // Nu_cp : Expectation of overall times
  // mat M =  Yi - Nu_cp_i*matrixP.t();

  return(logvrais);
}
//===========================================================================================
//'  Function that constructs the matrix lambda_t_j, the value of the processes at time t_j, given the random effects
//'  
//' @param nD an integer indicating the number of processes
//' @param tau_i vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel (beta)
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel (gamma)
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matrix containing at time t
//' @param randomeffects vector of individual random effects
//' @param zi model.matrix for random change's submodel
//' @param ordered indicator if tau_i is ordered or not (usually ordered, except for GK nodes)
//' 
//' @return a vector
//' @export
//' 
// [[Rcpp::export]]
arma::vec matNui_ui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                 arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::vec& randomeffects, arma::mat& zi,
                 bool ordered){
  // matNu_i : matrix of size xi.n_rows*nD containing expectation of processes
  
  int n_cols_xi = xi.n_cols;
  int n_cols_zi = zi.n_cols;
  int mi=tau_i.size(); // number of observations
  int T = max(tau_i)+1;
  mat matNu_i = zeros(mi,nD);
  mat Mu_t = zeros(nD,1);

  colvec ui=randomeffects.subvec( 0, nD-1 );//randomeffects(linspace(0, nD-1, nD));
  colvec vi;
  if(randomeffects.size()>nD){
    vi=randomeffects.subvec( nD, randomeffects.size()-1 );
  }else{
    vi = zeros(1);
  }

  //Verify Xi, Zi with randomeffects.size()>1
  int i = 0;
  for(int t=0; t< T; t++){
    if(t==0){
      Mu_t = x0i*alpha_mu0+ ui;
    }else{
      if(randomeffects.size()>nD){
        if(1>2&t==(T-1)) {
          std::cout << " mattRE "<<t << " "<<T<<endl;
          //std::cout << " xi "<< xi <<endl;
          std::cout << " xi "<<xi.n_rows<< " "<<n_cols_xi<<endl;
          std::cout << " t*nD "<<t*nD << " (t+1)*nD-1 "<<(t+1)*nD-1<<endl;
          std::cout << " n_cols_xi-1 "<<n_cols_xi-1<<endl;
          std::cout << " tau_i-1 "<<tau_i.t()<<endl;
          
        }
        Mu_t = DeltaT*(xi(span(t*nD,(t+1)*nD-1), span(0,n_cols_xi-1))*alpha_mu +
          zi(span(t*nD,(t+1)*nD-1), span(0,n_cols_zi-1))*vi) +
          G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),nD*(t-1)+nD-1))*Mu_t;
      }else{
        Mu_t = DeltaT*(xi(span(t*nD,(t+1)*nD-1), span(0,n_cols_xi-1))*alpha_mu) +
          G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),nD*(t-1)+nD-1))*Mu_t;
      }

    }

    if(ordered){
      if(t ==(int)tau_i(i)){
        matNu_i(span(i,i), span(0,nD-1)) = Mu_t.t();
        i++;
      }
    }else{
      for(int ii=0; ii< tau_i.size(); ii++){
        if(t ==(int)tau_i(ii)){
          matNu_i(span(ii,ii), span(0,nD-1)) = Mu_t.t();
        }
      }
    }
  }
  return (YiwoNA(vectorise(matNu_i)));
}

//===========================================================================================================
double f_marker_ui(arma::vec& Lambdai, int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Ytildi, arma::mat& YtildPrimi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::vec& paraSig, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::colvec& paraEtha2, arma::vec& if_link, arma::colvec& zitr, 
                arma::mat& ide, arma::vec& paras_k, arma::vec& K2_lambda_t, arma::vec& K2_lambda){
  
  int K = matrixP.n_rows;
  double logvrais = 0.0;

  int param = 0;
  int K_t = 0;

  for(int k=0; k<K; k++){
    // Retrieve part of Ytildi and Lambdai that corresponds to marker k
    vec Ytildik=Ytildi.col(k);
    
    if(!std::all_of(Ytildi.col(k).begin(), Ytildi.col(k).end(), Isnotnan) ){
      int ind=0;
      
      for(int b = 0 ; b < Ytildi.n_rows; b++){
        if(!isnan(Ytildi(b,k))){
          Ytildik(ind)=Ytildi(b,k);
          ind++;
        }
      }
      Ytildik =Ytildik.head(ind);
    }else{
      Ytildik=Ytildi.col(k);
    }
    
    // K2_lambda_t: which latent process linked to the observation in lambda
    // K2_lambda: which latent process linked to the K markers
    
    vec Lambdai_k(Ytildik.size());
    int ind=0;
    for(int b = 0 ; b < K2_lambda_t.size(); b++){
      if(K2_lambda_t(b)==K2_lambda(k)){
        Lambdai_k(ind)=Lambdai(b);
        ind++;
      }
    }
    
    
    if (if_link[k]==2){// thresholds
      
      double inf;
      double sup;
      bool lower=true;
      double vrais=1;
      
      vec params_thresholds(paras_k[k]);
      params_thresholds[0] = paraEtha2[param];
      
      for(int r=1; r<paras_k[k]; r++){
        params_thresholds[r] = params_thresholds[r-1] + pow(paraEtha2[param + r],2);
      }
      
      for(int j=0; j<(int)tau_i.size(); j++){
        
        if (Ytildik(j,k)==zitr[K_t*2]){ // if Ytildi == minY
          double gamma0 =  params_thresholds[0] - Lambdai_k[j]; 
          
          vrais  *= normalCDF(gamma0, lower);
          
        }else{
          sup = params_thresholds[0];
          inf = sup;
          
          int pp=0;
          for(int p=0; p<(int)(zitr[K_t*2 + 1]-zitr[K_t*2]-1); p++){
            
            if(ide(p)==1){
              sup = params_thresholds[pp+1];
              
              if(Ytildik(j,k) == (zitr[K_t*2]+p)){
                double gamma1 = sup - Lambdai_k[j];
                double gamma0 = inf - Lambdai_k[j];
                vrais  *= (normalCDF(gamma1) - normalCDF(gamma0));
                
                //cout << " p "<<p << " pp "<<pp << endl;
                //cout << " sup "<<sup << " gamma1 "<<gamma1 << " normalCDF(gamma1) "<<normalCDF(gamma1)<<endl;
                //cout << " inf "<<inf << " gamma0 "<<gamma0 << " normalCDF(gamma0) "<<normalCDF(gamma0)<<endl;
                //cout << " vrais "<<vrais << " CDF1-CDF0 "<<(normalCDF(gamma1) - normalCDF(gamma0)) << endl;
                
                // fout << "Ytildij: "<< Ytildi(j,k) << " zitr[K_t*2]+p: "<<zitr[K_t*2]+p 
                //      << " sup "<<sup << " inf "<<inf << " param "<< param
                //      << " params_thresholds[p+1] "<< params_thresholds[p+1] << " vrais "<<vrais
                //      << " temp "<<temp<< " gamma1 "<<gamma1<< " gamma0 "<<gamma0
                //      <<" Ngamma1 "<<normalCDF(gamma1)<< " Ngamma0 "<<normalCDF(gamma0)<<endl;
              }
              inf = sup;
              pp += 1;
            }
          }
          
          if (Ytildik(j,k)==zitr[K_t*2 + 1]){ // if Ytildi == maxY
            double gamma0 =  sup - Lambdai_k[j]; 
            vrais  *= (1-normalCDF(gamma0));
          }
        }
      }
      K_t +=1;
      logvrais += log(vrais);
      
    }else if (if_link[k]==0){// linear
      
      // vec params_lin(paras_k[k]);
      // for(int r=0; r<paras_k[k]; r++){
      //   params_lin[r] = paraEtha2[param + r];
      // }
      
      // ##### computering of the likelihood ##########################
      double log_det= Lambdai_k.size()*log(paraSig[k]);
      //double log_Jac_Phi=log(params_lin[0])*Lambdai_k.size();
      double log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi.col(k)))));
      vec y=Ytildik-Lambdai_k;
      
      //       vrais = -0.5*(sum(Lambdai_k.size())*log(2*M_PI) + log_det + as_scalar(y.t()*y/paraSig[k])) ;
      logvrais += -0.5*(Lambdai_k.size()*log(2*M_PI) + log_det + as_scalar(y.t()*y/paraSig[k])) + log_Jac_Phi;
      
      //dmvnorm(Ytilde, mean = Ytild[i+1,], sigma = sigmae2*diag(length(mu))) / a^(length(mu))
    }
    param += paras_k[k];
  }
  
  // int mk=ncol(matrixP); //number of ordinal markers
  // for(int i=0; i<(int)tau_i.size(); i++){
  //   for(int j=0; j<mk; i++){
  // 
  //   }
  // }
  // -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i))
  // 
  // 
  // // Yi : matrice of observation of subject i
  // mat Nu_cp = matNui(nD, tau, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
  // mat Nu_cp_i = zeros(Yi.n_rows,nD);
  // for(int i=0; i<(int)tau_i.size(); i++){
  //   Nu_cp_i.row(i) = Nu_cp.row(tau_i(i));
  // }
  // // Nu_cp : Expectation of overall times
  // mat M =  Yi - Nu_cp_i*matrixP.t();
  
  return(logvrais);
}

// /*============================================================
// Computes the hazard risk or survival if surv == true
// ==============================================================*/
//' gammaX vector of linear predictors for 1 and 2 transitions (including association on random effects if assoc <=2)
//' surv Computation of survival funtion (1) or risk function (0)
//' trans index for computation of survival function on all transitions (-1), on first transition(0), or second transition (1) for nE=2
//' when 
double fct_risq_base(double t,  int status, arma::vec& param_basehaz, int basehaz, arma::vec& knots_surv, 
                     int nE, arma::vec& gammaX, bool surv, int trans){
  double out=0;
  
  if(basehaz==0){
    
    if(nE == 1){
      if(surv){ //Computation of survival
        //out = exp(-pow(t/param_basehaz(1),param_basehaz(0))* exp(gammaX(0)));
        out = exp(-pow(t*param_basehaz(0),param_basehaz(1))* exp(gammaX(0)));
      }else{
        // b0'=b1, b1'=1/b0
        //out = param_basehaz(0)/param_basehaz(1)*pow(t/param_basehaz(1),(param_basehaz(0)-1))* exp(gammaX(0));
        out = param_basehaz(0)*param_basehaz(1)*pow(t*param_basehaz(0),(param_basehaz(1)-1))* exp(gammaX(0));
      }
      
    }else{ //nE==2
      if(surv){
        if(trans==-1){
          //out = exp(-pow(t/param_basehaz(1),param_basehaz(0))* exp(gammaX(0)))* exp(-pow(t/param_basehaz(3),param_basehaz(2))* exp(gammaX(1)));
          out = exp(-pow(t*param_basehaz(0),param_basehaz(1))* exp(gammaX(0)))* exp(-pow(t*param_basehaz(2),param_basehaz(3))* exp(gammaX(1)));
        }else if(trans==0){
          //out = exp(-pow(t/param_basehaz(1),param_basehaz(0))* exp(gammaX(0)));
          out = exp(-pow(t*param_basehaz(0),param_basehaz(1))* exp(gammaX(0)));
        }else if(trans==1){
          //out = exp(-pow(t/param_basehaz(3),param_basehaz(2))* exp(gammaX(1)));
          out = exp(-pow(t*param_basehaz(2),param_basehaz(3))* exp(gammaX(1)));
        }

      }else{ //hazard risk 
        if(status == 0){
          out = 1;
        }else if(status == 1){
          //out = param_basehaz(0)/param_basehaz(1)*pow(t/param_basehaz(1),(param_basehaz(0)-1))* exp(gammaX(0));
          out = param_basehaz(0)*param_basehaz(1)*pow(t*param_basehaz(0),(param_basehaz(1)-1))* exp(gammaX(0));
        }else if(status == 2){
          //out = param_basehaz(2)/param_basehaz(3)*pow(t/param_basehaz(3),(param_basehaz(2)-1))* exp(gammaX(1));
          out = param_basehaz(2)*param_basehaz(3)*pow(t*param_basehaz(2),(param_basehaz(3)-1))* exp(gammaX(1));
        }
      }
    }

    
  }else{ // Splines
    
  }

  
  return(out);
}

// /*============================================================
// Computes the hazard risk or survival function for a vector or timepoints 
// ==============================================================*/
//' ptGK_delta vector of projections of GK nodes onto grid of delta
//' ptGK vector of individual GK nodes
//' alpha vector of association parameters
//' delta_i event status 
//' survfunc indicator if output is survival function or hazard risk
//' trans index for computation of survival/risk function on all transitions (-1), on first transition(0), or second transition (1)

vec fct_pred_curlev_slope(arma::vec& ptGK_delta, arma::vec& ptGK, arma::colvec& xti1, arma::colvec& xti2, 
                          arma::mat& xti1_intY, arma::mat& xti2_intY, arma::vec& ui_r, int delta_i, arma::colvec& param_surv, arma::colvec& param_surv_intY, int assoc, 
                          int nD, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, arma::colvec& alpha_mu,
                          arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi, arma::vec& param_basehaz, int basehaz, arma::vec& knots_surv, 
                          arma::vec& gamma_X, int nE, bool survfunc, int trans){
  
  vec out = zeros(ptGK_delta.size());
  vec curlev = zeros(ptGK_delta.size());
  vec curslope = zeros(ptGK_delta.size());

  int nA = 1; // length of vector of association parameters
  if(assoc == 2 || assoc == 5)
    nA++;
  nA *= nD;
  vec alpha = zeros(nA*nE) ; // association params for event 1, then for event 2

  for( int i=0; i<nE; i++){
    for(int j=0; j<nA; j++){
      alpha(i*nA+j)= param_surv(xti1.size() + xti2.size()*i + i*nA+j);
    }
  }

  // prediction Y(t)
  if(assoc == 3 || assoc == 5){ 
    curlev = matNui_ui(nD, ptGK_delta, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi, false);
  }

  // std::cout << " ptGK_delta "<< ptGK_delta.t()
  //           << " alpha "<<alpha.t()
  //           << " curlev "<<curlev.t()
  //           << " ui_r "<<ui_r.t()
  //           << " DeltaT "<<DeltaT<<endl;
  
  if(assoc == 4 || assoc == 5){
    std::cout << " to develop !"<<endl;
  }
  
  // Computation of the survival/risk function
  if(survfunc){ // survival function, never used so far
    std::cout << " verif use of fct_pred_curlev_slope for survival computation!! ";
    mat cumrisq(ptGK_delta.size(), nE);
    
    for(int j=0; j<nE; j++){
      for( int i=0; i<ptGK.size(); i++) // Cumulative risk with regression parameters
        cumrisq(i,j) = -log(fct_risq_base(ptGK(i), 0, param_basehaz, basehaz, knots_surv, nE, gamma_X, true, j));
    
    }//nE
    
    for( int i=0; i<ptGK.size(); i++){
      for( int j=0; j<nE; j++){
        double temp = alpha(j)*curlev(i);

        if(param_surv_intY.size()>0){
          std::cout << " ajout interactions  ";
          if(j==0){
            for(int jj=0; jj<xti1_intY.size(); jj++){
              //temp += param_surv_intY(jj)*curlev(i)*xti1_intY(jj);
            }
          }
        }
        
        out(i) *= exp(-cumrisq(i,j)*exp(temp));//Changer alpha avec nD
      }
    }
    


  }else{ // Computation of the risk

    if(trans==-1){

      mat risq = zeros(ptGK_delta.size(),nE);
      mat alphaY = zeros(ptGK.size(), nE);
      
      for( int j=0; j<nE; j++){
        // alpha * Yt
        for( int i=0; i<ptGK.size(); i++){
          for( int nd=0; nd<nD; nd++)
            alphaY(i, j) += alpha(j*(nD)+nd)*curlev(i*nD+nd);
        }

        for( int i=0; i<ptGK.size(); i++){// Instantaneous risk with regression parameters
          risq(i,j) = fct_risq_base(ptGK(i), j+1, param_basehaz, basehaz, knots_surv, nE, gamma_X,false, j);
          
          double temp = alphaY(i,j);

          if(param_surv_intY.size()>0){
            if(j==0){
              for(int jj=0; jj<xti1_intY.size(); jj++){
                for( int nd=0; nd<nD; nd++){
                  temp += param_surv_intY(jj*nD+nd)*curlev(i*nD+nd)*xti1_intY(0,jj);
                }
              }
            }
            if(j==1){
              for(int jj=0; jj<xti2_intY.size(); jj++){
                for( int nd=0; nd<nD; nd++){
                  temp += param_surv_intY(xti1_intY.n_cols*nD+jj*nD+nd)*curlev(i*nD+nd)*xti2_intY(0,jj);
                }
              }
            }
          }

          out(i) += risq(i,j)*exp(temp); // curlev(i) = exp(alpha Yt)
        } 
      }

    }else {

      // alpha * Yt
      vec alphaY= zeros(ptGK.size());
      for( int i=0; i<ptGK.size(); i++){
        for( int nd=0; nd<nD; nd++){
          alphaY(i) += alpha(trans*nD+nd)*curlev(i*nD+nd);// trans=0 or 1
        }
      }

      vec risq = zeros(ptGK_delta.size());
      for( int i=0; i<ptGK.size(); i++) // Instantaneous risk with regression parameters for transition trans
        risq(i) = fct_risq_base(ptGK(i), delta_i, param_basehaz, basehaz, knots_surv, nE, gamma_X, false, trans); //trans = delta_i-1

      for( int i=0; i<ptGK.size(); i++) {
        double temp=alphaY(i);
        if(param_surv_intY.size()>0){
          if(trans==0){
            for(int jj=0; jj<xti1_intY.size(); jj++){
              for( int nd=0; nd<nD; nd++){
                temp += param_surv_intY(jj*nD+nd)*curlev(i*nD+nd)*xti1_intY(0,jj);
              }
            }
          }
          if(trans==1){
            for( int nd=0; nd<nD; nd++){
              for(int jj=0; jj<xti2_intY.size(); jj++){
                         // << "param corr "<< param_surv_intY(xti1_intY.n_cols*nd+jj)<<endl;//xti1_intY.n_cols*nD+jj*nD+nd
                temp += param_surv_intY(xti1_intY.n_cols*nD+jj*nD+nd)*curlev(i*nD+nd)*xti2_intY(0,jj);
              }
            }
          }
        }
        
        out(i) += risq(i)*exp(temp);
      }
    }
  }
  return(out);
}

//===========================================================================================================

//===========================================================================================================

// /*============================================================
// Computes the survival function using Gauss Konrod quadrature
// ==============================================================*/
// combines a 7-point Gauss rule with a 15-point Kronrod rule (Kahaner, Moler & Nash 1989, §5.5).
// Gauss points are incorporated into the Kronrod points
double fct_surv_Konrod(double t_i, arma::colvec& xti1, arma::colvec& xti2, arma::mat& xti1_intY, arma::mat& xti2_intY, arma::vec& ui_r, arma::vec& param_basehaz, int basehaz, arma::colvec& param_surv, arma::colvec& param_surv_intY, arma::vec& knots_surv, int assoc,
                       int nD, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi,
                       int nE, arma::vec& gamma_X){

  // Konrod weights
  vec wgk(8);
  wgk(0)=0.022935322010529224963732008058970;
  wgk(1)=0.063092092629978553290700663189204;
  wgk(2)=0.104790010322250183839876322541518;
  wgk(3)=0.140653259715525918745189590510238;
  wgk(4)=0.169004726639267902826583426598550;
  wgk(5)=0.190350578064785409913256402421014;
  wgk(6)=0.204432940075298892414161999234649;
  wgk(7)=0.209482141084727828012999174891714;
  
  vec wgk_15(15);
  wgk_15(0)=wgk(0)  ;
  wgk_15(1)=wgk(0)  ;
  wgk_15(2)=wgk(1);
  wgk_15(3)=wgk(1);
  wgk_15(4)=wgk(2);
  wgk_15(5)=wgk(2);
  wgk_15(6)=wgk(3);
  wgk_15(7)=wgk(3);
  wgk_15(8)=wgk(4);
  wgk_15(9)=wgk(4);
  wgk_15(10)=wgk(5);
  wgk_15(11)=wgk(5);
  wgk_15(12)=wgk(6);
  wgk_15(13)=wgk(6);
  wgk_15(14)=wgk(7);
  
  // Konrod nodes
  vec ptGK(15);
  ptGK(0) = 0.991455371120812639206854697526329;
  ptGK(1) = -0.991455371120812639206854697526329;
  ptGK(2) = 0.949107912342758524526189684047851;
  ptGK(3) = -0.949107912342758524526189684047851;
  ptGK(4) = 0.864864423359769072789712788640926;
  ptGK(5) = -0.864864423359769072789712788640926;
  ptGK(6) = 0.741531185599394439863864773280788;
  ptGK(7) = -0.741531185599394439863864773280788;
  ptGK(8) = 0.586087235467691130294144838258730;
  ptGK(9) = -0.586087235467691130294144838258730;
  ptGK(10) = 0.405845151377397166906606412076961;
  ptGK(11) = -0.405845151377397166906606412076961;
  ptGK(12) = 0.207784955007898467600689403773245;
  ptGK(13) = -0.207784955007898467600689403773245;
  ptGK(14) = 0.000000000000000000000000000000000;
  
  vec risq_GK_event = zeros(15); // baseline hazard
  vec pred_GK_event = zeros(15); // prediction Y(t)

  double surv=0;
  vec fct_pred_surv = zeros(15);
  vec ptGK_ti = zeros(15); // GK nodes that correspond to [0, t_i] 
  vec deltaT_ptGK_ti = zeros(15); // projection of GK nodes onto the deltaT grid

  for( int i=0; i<15; i++){
    // projection of GK nodes onto [0, ti]
    ptGK_ti(i) = ptGK(i)*(t_i-0)/2+(t_i+0)/2;
    // Correspondance ptGK_ti and tau_i for computation of current Y
    deltaT_ptGK_ti(i) = round(ptGK_ti(i)/ (double) DeltaT); // search for last tau_i before ptGK_ti(i)
  }

  risq_GK_event= fct_pred_curlev_slope(deltaT_ptGK_ti, ptGK_ti, xti1, xti2, xti1_intY, xti2_intY, ui_r, 1, param_surv, param_surv_intY, assoc, 
                                       nD, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, param_basehaz, basehaz, knots_surv, 
                                       gamma_X, nE, false,-1); //trans=-1

  double cumrisk=0;
  for( int i=0; i<15; i++){
    cumrisk += wgk_15(i)*risq_GK_event(i);
  }

  cumrisk *= t_i/2;
  surv=exp(-cumrisk);

  return(surv); 
}

//===========================================================================================================

// /*============================================================
// Computes the individual likelihood of time-to-event, conditionally on the random effects
// ==============================================================*/
//
//' ui_r vector of individual random effects
//' t_i individual time-to-event
//' delta_i individual status of event
//' xti1 vector of individual covariates for first event
//' xti2 vector of individual covariates for competing event
//' param_surv regression parameters
//' param_basehaz parameters for baseline hazard function
//' basehaz type of baseline hazard function
//' knots_surv vector of knots if basehaz == Splines
//' assoc function of the random effects that captures association 
//' //'    (0: random intercept, 1: random slope, 2: random intercept and slope, 3: current value, 4: current slope, 5: current value and slope)
//' truncation boolean, indicating if left truncation or not
//' 
arma::vec f_survival_ui(arma::vec& ui_r, double t_0i, double t_i,int delta_i, arma::colvec& xti1, arma::colvec& xti2, 
                     arma::mat& xti1_intY, arma::mat& xti2_intY,  arma::colvec& param_surv, arma::colvec& param_surv_intY,
                     arma::vec& param_basehaz, int basehaz, arma::vec& knots_surv, int assoc, bool truncation,
                     int nD, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                     arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi,
                     int nE){

  mat surv (1,1);
  surv(0,0) = 0;
  mat haz(1,1);
  haz(0,0) = 1;
  vec fti = ones(2);
  vec  gamma_X(nE);
  mat gammaX = zeros(nE,1);

  int nA = 1;
  if(assoc == 2 || assoc == 5) //random intercept + slope
    nA ++;
  nA *= nD;

  if(xti1.size()>0)
    gammaX(span(0,0), span(0,0)) = xti1.t()*param_surv(span(0, xti1.size()-1));

  if(nE==2 && xti2.size()>0 )
    gammaX(span(1,1), span(0,0)) = xti2.t()*param_surv(span(xti1.size()+ nA, xti1.size()+ nA+ xti2.size()-1));
    //gammaX(span(1,1), span(0,0)) = xti2.t()*param_surv(span(xti1.size() + nA, xti1.size() + nA + xti2.size()-1));

  if(assoc <= 2){// random intercept (0), random slope (1) or both (2)

    vec p_surv(nA*nE);
    for(int j=0; j<nE; j++){
      p_surv(span(j*nA,j*nA+nA-1))=param_surv(span(j*(nA+xti1.size()+xti2.size())+(j+1)%2*xti1.size(),j*(nA+xti1.size()+xti2.size())+(j+1)%2*xti1.size()+nA-1));
    }

    
    //int tp = xti1.size();
    for( int j=0; j<nE; j++){
      if(assoc == 0){
        vec ui_surv(nD);
        int nRE = ui_r.size()/nD;
        for(int d=0; d<nD; d++)
          ui_surv(d)=ui_r(d*nRE+d);
        
        gammaX(span(j,j), span(0,0)) += ui_surv*p_surv;
        
      }else if(assoc == 1){
        vec ui_surv(nD);
        int nRE = ui_r.size()/nD;
        for(int d=0; d<nD; d++)
          ui_surv(d)=ui_r(d*nRE+1);

        gammaX(span(j,j), span(0,0)) += ui_r(1)*param_surv(xti1.size()+ xti2.size()+j);
        
      }else if(assoc == 2){

        gammaX(span(j,j), span(0,0)) += ui_r.t()*p_surv(span(j*nA,j*nA+nA-1));
      }
    }

    for(int j=0; j<nE; j++)
      gamma_X(j) = gammaX(j,0);
    surv(0,0) = fct_risq_base(t_i, 0, param_basehaz, basehaz, knots_surv, nE, gamma_X, true, -1);

    if(delta_i>0){
      haz(0,0) = fct_risq_base(t_i, delta_i, param_basehaz, basehaz, knots_surv, nE, gamma_X, false, -1);
    }else{
      haz(0,0)=1;
    }
    
    fti(0) = surv(0,0)*haz(0,0);
    
  }else{

    for( int j=0; j<nE; j++)
      gamma_X(j) = gammaX(j,0);

    vec event(1);
    event(0)=t_i;
    vec deltaT_ptGK_ti(1);
    deltaT_ptGK_ti(0) = round(t_i/ (double) DeltaT);
    vec hazard(1);
    hazard.fill(1);

    if(delta_i!=0)
      hazard = fct_pred_curlev_slope(deltaT_ptGK_ti, event, xti1, xti2, xti1_intY, xti2_intY, ui_r, delta_i, param_surv, param_surv_intY, assoc, 
                                    nD, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, param_basehaz, basehaz, knots_surv, 
                                    gamma_X, nE, false, delta_i-1); // trans = delta_i-1
    //double test = fct_risq_base(t_i, 0, param_basehaz, basehaz, knots_surv, nE, gamma_X, true, -1);
    double surv = 1;

    surv = fct_surv_Konrod(t_i, xti1, xti2, xti1_intY, xti2_intY, ui_r, param_basehaz, basehaz, param_surv, param_surv_intY, knots_surv, assoc,
                           nD, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, nE, gamma_X);
    fti(0) = surv*hazard(0,0);

    double surv0 = 1;
    if(truncation)
      surv0 = fct_surv_Konrod(t_0i, xti1, xti2, xti1_intY, xti2_intY, ui_r, param_basehaz, basehaz, param_surv, param_surv_intY, knots_surv, assoc,
                           nD, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, nE, gamma_X);
    fti(1) = surv0;

    // std::cout << " gamma_X "<<gamma_X<<endl
    //           << " event(0) "<<event(0)<<endl
    //           << " deltaT_ptGK_ti(0) "<<deltaT_ptGK_ti(0)<<endl
    //           << " delta_i "<<delta_i<<endl
    //           << " hazard "<< hazard<<endl
    //           << " surv "<<surv<<endl
    //           << " surv0 "<<surv0<<endl;
    
  }
  
  
  
 return(fti); 
}

//===========================================================================================================

// /*============================================================
// generate one random gaussian variable with mean and stddev
// ==============================================================*/
// double rand_normal(double mean, double stddev)
// {//Box muller method
//   static double n2 = 0.0;
//   static int n2_cached = 0;
//   if (!n2_cached)
//   {
//     double x, y, r;
//     do
//     {
//       x = 2.0*rand()/RAND_MAX - 1;
//       y = 2.0*rand()/RAND_MAX - 1;
//
//       r = x*x + y*y;
//     }
//     while (r == 0.0 || r > 1.0);
//     {
//       double d = sqrt(-2.0*log(r)/r);
//       double n1 = x*d;
//       n2 = y*d;
//       double result = n1*stddev + mean;
//       n2_cached = 1;
//       return result;
//     }
//   }
//   else
//   {
//     n2_cached = 0;
//     return n2*stddev + mean;
//   }
// }
//
// /*============================================================
// Generate multivariate gaussian vector with mean m and SD = LLt
// ==============================================================*/
// // [[Rcpp::export]]
// arma::mat mvnorm(int seed, arma::vec& m, arma::mat& SD){
//   mat L = chol(SD).t();
//   // mat L = SD.t();
//   int d = m.size();
//   vec xx = zeros(d); // for indpt standard gaussian vector
//   for( int i=0; i<d; i++){
//     srand (seed*33*i);
//     xx[i] = rand_normal(0.0,1.0);
//   }
//   return(m + L*xx);
// }

/*============================================================
 Generate multivariate gaussian vector with mean m and SD = LLt
==============================================================*/
arma::mat mvnorm(int seed, arma::vec m, arma::mat SD){
  Rcpp::Environment base("package:DynNet");
  Rcpp::Function g = base["f_mvrnorm"];
  vec x = as<arma::vec>(wrap(g(Rcpp::_["seed"] = seed, Rcpp::_["m"] = m, Rcpp::_["sd"] = SD)));
  return(x);
}

// /*================================================================
// Function MC : Compute the prediction in real scale using Monte Carlo approch
//  Dependences : none
// ==================================================================*/
// // [[Rcpp::export]]
// arma::vec MC(int K, int nr, arma::vec& mu, arma:: mat& SD, List& knots, arma::vec& ParaTransformY, int degree){
//   Rcpp::Environment base("package:DynNet");
//   Rcpp::Function f = base["R_MC"];
//   arma::vec yi = as<arma::vec>(wrap(f(K, nr, mu, SD, knots, ParaTransformY, degree)));
//   return(yi);
// }


//===========================================================================================
//' Function that transforms a vector to a matrix
//' 
//' @param y a vector 
//' @param K an integer indicating the number of columns of the returned matrix
//' @param m_i an integer indicating the number of rows of the returned matrix
//' 
//' @return a matrix 
//' @export
//' 
// [[Rcpp::export]]
arma::mat VecToMat(arma::vec& y, int K, int m_i){

  mat mat_y = zeros(m_i,K);
  for(int p = 0; p < m_i; p++){
    mat_y.row(p) = y(span(p*K, ((p+1)*K-1))).t();
  }

  return(mat_y);
}
