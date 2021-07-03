#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <eigen3/Eigen/Dense>
#include <stdexcept>
#include <math.h>
#include <vector>
#include <algorithm>
#include "genericfun.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

//using Eigen::MatrixXd;
/*
 individual contribution to the log-likelihood
 */

double Loglikei(int K, int nD, arma::mat matrixP, int m_i, arma::vec tau, arma::vec tau_i, arma::mat Ytildi, arma::mat YtildPrimi,
                arma::mat x0i, arma::mat z0i, arma::mat xi, arma::mat zi, arma::colvec alpha_mu0,
                arma::colvec alpha_mu, arma::mat matDw,  arma::mat matDw_u, arma::mat matDu,
                arma::mat matB, arma::mat Sig, arma::mat G_mat_A_0_to_tau_i,
                arma::mat G_mat_prod_A_0_to_tau,  double DeltaT){

  double loglik_i = 0.e0;
  // ###### compute  Yi - E(Yi) ##### deleting missing values #####
  vec Ytildi_nu_i = YiNui(nD, matrixP, tau, tau_i, DeltaT, Ytildi, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
  
  //  ##### computing of the  matrix  matVY_i #####
  int sizeYi = Ytildi_nu_i.size();
  mat matVY_i = zeros(sizeYi,sizeYi); // initialization of the  matVY_i to zero
  vec k_i = zeros<vec>(m_i);// vector of number of observed marker at each observation time
  int q = zi.n_cols; // number of random effets on the slope
  
  // found the max of the vector tau_i
  double maxTau_i = 0.0;
  for (unsigned i = 0; i < tau_i.size(); i++){
    if (tau_i[i] > maxTau_i)
      maxTau_i = tau_i[i];
  }
  mat GrdZi = zeros<mat>((maxTau_i+1)*nD, q);
  int p_j=0; // loop variable
  int p_k =0; // loop variable
  vec vect = ones<vec>(nD);
  double log_Jac_Phi=0.e0;
  
  // ##### computering of GrdZi ####################################
  //  GrdZi : this matrix contains all model.matrix at each time of tau_i
  for(int t = 0; t<= maxTau_i; t++){
    if(t==0){
      GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = 0.e0*zi(span(t*nD,(t+1)*nD-1), span(0,q-1));
    }
    else{
      GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = (DeltaT*zi(span(t*nD,(t+1)*nD-1), span(0,q-1)) +
        G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),t*nD-1))*GrdZi(span((t-1)*nD,t*nD-1), span(0,q-1)));
    }
  }
  
  // ######################################################################################################
  for( int j =0 ; j < m_i; j++){
    p_k = p_j;
    // ###### computing of the matrix H_i_t_j  ##############
    mat matH_i_t_j = matHit(Ytildi.row(j).t()); // Yi.row(j) = vector of observations at time j;
    for( int k =j ; k < m_i; k++){
      // ###### computering of the matrix H_i_t_k ##############
      mat matH_i_t_k = matH_i_t_j;
      if(k!=j) {
        matH_i_t_k = matHit(Ytildi.row(k).t()); // Yi.row(k) = vector of observations at time k;
      }
      if(sum(sum(matH_i_t_k)) == 0){ // #### this mean that there is no observation at this time
        k_i(k) = 0;
      }
      else{
        k_i(k) = matH_i_t_k.n_rows;
      }
      
      // ###########################################################
      if( (k_i(j) != 0) & (k_i(k) != 0)){
        //
        // ###### Computering of the componante relate to the RE ########
        // start: computering of phi_0_j_0 et de phi_0_k_0###############
        mat phi_0_j_0 = diagmat(vect);
        mat phi_0_k_0 = diagmat(vect);
        if(tau_i(j)>0){
          phi_0_j_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(j)-1),nD*(tau_i(j)-1)+nD-1));
        }
        if(tau_i(k)>0){
          phi_0_k_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(k)-1),nD*(tau_i(k)-1)+nD-1));
        }
        
        matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
          ((phi_0_j_0*z0i)*matDw* (z0i*phi_0_k_0).t() +
          (phi_0_j_0*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
          (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*phi_0_k_0).t() +
          (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t()
          )*matrixP.t())*matH_i_t_k.t();
        
        //// For simplicty we define GrdZi(span(0*K,(0+1)*K-1), span(0,q-1))) = 0 : as if we spefify
        // model for slope from t_0 but with covariate set to 0
        
        if( k == j ){
          // ##### Fill the diagonal of the matrix VY_i by adding  \Sigma  #########################
          matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*Sig*matH_i_t_j.t();
        }
        
        // ###### VY_i is a symetric matrix; so we fill lower triangular matri by transpose the upper triangular part #########
        if(k != j){
          matVY_i(span((p_k), (p_k+k_i(k)-1)), span(p_j,(p_j+k_i(j)-1))) = matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))).t();
        }
      }
      p_k += k_i(k); // incrementation of p_k
    }
    p_j += k_i(j);// incrementation of p_j
  }
  // ######################### computing logJac #################################
  log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi))));
  double abs_det_matVY_i = abs(det(matVY_i));
  
  // ##### computering of the likelihood ##########################
  loglik_i = -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i)) + log_Jac_Phi;

  return(loglik_i);
  // return(log_Jac_Phi);
}




bool compFun0(int i) {
  return i > 0;
}


double normalCDF(double value)
{
  return 0.5 * erfc(-value * M_SQRT1_2);
}

/*
 individual contribution to the log-likelihood
 Generate latent process for antithetic MC: \int f(Y|Lambda) f(Lambda) dLambda = 1/N sum^N/2 [f(Y|Lambda_n)+f(Y|\tilde{Lambda}_n)]
 Covariance covariance matrix of the vector Lambda_i(t_01, t_02, t11, t12) with t_{occasion,process}
 */

double Loglikei_GLM(int K, int nD, arma::mat& matrixP, int m_i, arma::vec& tau, arma::vec tau_i, arma::mat Ytildi, arma::mat YtildPrimi,
                    arma::mat x0i, arma::mat z0i, arma::mat xi, arma::mat zi, arma::colvec& alpha_mu0,
                    arma::colvec& alpha_mu, arma::mat& matDw,  arma::mat& matDw_u, arma::mat& matDu,
                    arma::mat& matB, arma::mat& Sig, arma::mat& G_mat_A_0_to_tau_i,
                    arma::mat& G_mat_prod_A_0_to_tau,  double DeltaT, arma::vec& ParamTransformY, arma::vec& df, arma::vec& if_link,  
                    arma::vec& zitr,  arma::mat& ide, arma::vec& paras_k,
                    double t_0i, double t_i, int delta_i, arma::vec& xti1, arma::vec& xti2, int basehaz, arma::vec& knots_surv, bool survival, arma::vec& param_surv, arma::vec& param_basehaz, int assoc, bool truncation, 
                    arma::mat& seq_i,  int type_int, arma::vec& ind_seq_i, int MCnr, int sub, int nE){

  vec Ytildi_nu_i;
  mat matVY_i;
  vec k_i = zeros<vec>(m_i);// vector of number of observed marker at each observation time
  vec PNu_cp_i;
  mat matV_i;
  mat sigMSM;
  int check=0; // 1 both, 2 close likelihood only
  //cout << " max(if_link) "<<max(if_link)<< " survival "<< survival << " check "<< check <<endl;

  if((max(if_link) < 2 && !survival) || check==1){//|| check==1){
    // ###### compute  Yi - E(Yi) ##### deleting missing values #####
    Ytildi_nu_i = YiNui(nD, matrixP, tau, tau_i, DeltaT, Ytildi, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);

    //  ##### Compute mean PNu_cp_i = E(P.Lambda_i) #####
    mat Nu_cp = matNui(nD, tau, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
    mat Nu_cp_i = zeros(Ytildi.n_rows,nD);
    for(int i=0; i<(int)tau_i.size(); i++){
      Nu_cp_i.row(i) = Nu_cp.row(tau_i(i));
    }
    // Nu_cp : Expectation of overall times
    PNu_cp_i =YiwoNA(vectorise(Nu_cp_i*matrixP.t()));
    
    //  ##### computing of the  matrix  matVY_i #####
    int sizeYi = Ytildi_nu_i.size();
    matVY_i = zeros(sizeYi,sizeYi); // initialization of the  matVY_i to zero
    matV_i = zeros(m_i*nD, m_i*nD); // Computation of the covariance matrix
    
    //mat matVY_icheck = zeros(sizeYi,sizeYi); // initialization of the  matVY_i to zero
    int q = zi.n_cols; // number of random effets on the slope
    
    // found the max of the vector tau_i
    double maxTau_i = 0.0;
    for (unsigned i = 0; i < tau_i.size(); i++){
      if (tau_i[i] > maxTau_i)
        maxTau_i = tau_i[i];
    }
    mat GrdZi = zeros<mat>((maxTau_i+1)*nD, q);
    int p_j=0; // loop variable
    int p_k =0; // loop variable
    vec vect = ones<vec>(nD);
    int add = 0;
    vec MSigmaM=zeros<vec>(sizeYi);
    int aa = 0;
    
    // ##### computering of GrdZi ####################################
    //  GrdZi : this matrix contains all model.matrix at each time of tau_i
    for(int t = 0; t<= maxTau_i; t++){
      if(t==0){
        GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = 0.e0*zi(span(t*nD,(t+1)*nD-1), span(0,q-1));
      }
      else{
        GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = (DeltaT*zi(span(t*nD,(t+1)*nD-1), span(0,q-1)) +
          G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),t*nD-1))*GrdZi(span((t-1)*nD,t*nD-1), span(0,q-1)));
      }
    }
    
    for( int j =0 ; j < m_i; j++){
      for( int k =j ; k < m_i; k++){
        mat phi_0_j_0 = diagmat(vect);
        mat phi_0_k_0 = diagmat(vect);
        if(tau_i(j)>0){
          phi_0_j_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(j)-1),nD*(tau_i(j)-1)+nD-1));
        }
        if(tau_i(k)>0){
          phi_0_k_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(k)-1),nD*(tau_i(k)-1)+nD-1));
        }
        
        if(nD>1)
          add = nD-1;
        
        matV_i(span(j*nD,j*nD+add),span(k*nD,k*nD+add)) =(phi_0_j_0*z0i)*matDw* (z0i*phi_0_k_0).t() +
          (phi_0_j_0*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
          (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*phi_0_k_0).t() +
          (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t()
          ;
        
        if(k != j){
          matV_i(span(k*nD,k*nD+add), span(j*nD,j*nD+add)) = matV_i(span(j*nD,j*nD+add),span(k*nD,k*nD+add)).t();
        }
      }
    }
    
    // ######################################################################################################
    for( int j =0 ; j < m_i; j++){
      
      p_k = p_j;
      // ###### computing of the matrix H_i_t_j  ##############
      mat matH_i_t_j = matHit(Ytildi.row(j).t()); // Yi.row(j) = vector of observations at time j;
      
      for( int k =j ; k < m_i; k++){
        
        // ###### computering of the matrix H_i_t_k ##############
        mat matH_i_t_k = matH_i_t_j;
        if(k!=j) {
          matH_i_t_k = matHit(Ytildi.row(k).t()); // Yi.row(k) = vector of observations at time k;
        }
        if(sum(sum(matH_i_t_k)) == 0){ // #### this mean that there is no observation at this time
          k_i(k) = 0;
        }
        else{
          k_i(k) = matH_i_t_k.n_rows;
        }
        
        // ###########################################################
        if( (k_i(j) != 0) & (k_i(k) != 0)){
          
          //
          // ###### Computering of the componant relating to the RE ########
          // start: computering of phi_0_j_0 et de phi_0_k_0###############
          mat phi_0_j_0 = diagmat(vect);
          mat phi_0_k_0 = diagmat(vect);
          if(tau_i(j)>0){
            phi_0_j_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(j)-1),nD*(tau_i(j)-1)+nD-1));
          }
          if(tau_i(k)>0){
            phi_0_k_0 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(k)-1),nD*(tau_i(k)-1)+nD-1));
          }
          
          matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
            ((phi_0_j_0*z0i)*matDw* (z0i*phi_0_k_0).t() +
            (phi_0_j_0*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
            (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*phi_0_k_0).t() +
            (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t()
            )*matrixP.t())*matH_i_t_k.t();
          
          //cout << "matVY_icheck\n "<< matVY_icheck(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1)));
          //cout << "matV_i2 "<< matH_i_t_j*( matrixP*
          //  (matV_i(span(j*nD,j*nD+add),span(k*nD,k*nD+add)))*matrixP.t())*matH_i_t_k.t();
          
          //matVY_icheck(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
          //  (matV_i(span(j*nD,j*nD+add),span(k*nD,k*nD+add)))*matrixP.t())*matH_i_t_k.t();
          
          //// For simplicty we define GrdZi(span(0*K,(0+1)*K-1), span(0,q-1))) = 0 : as if we specify
          // model for slope from t_0 but with covariate set to 0
          
          if( k == j ){
            // ##### Fill the diagonal of the matrix VY_i by adding  \Sigma  #########################
            mat MSM=matH_i_t_j*Sig*matH_i_t_j.t();
            matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += MSM;
            //matVY_icheck(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*Sig*matH_i_t_j.t();
            MSigmaM(span(aa, aa +  matH_i_t_j.n_rows -1))= MSM.diag();
            aa = aa + matH_i_t_j.n_rows;
          }
          
          // ###### VY_i is a symetric matrix; so we fill lower triangular matri by transpose the upper triangular part #########
          if(k != j){
            matVY_i(span((p_k), (p_k+k_i(k)-1)), span(p_j,(p_j+k_i(j)-1))) = matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))).t();
            //matVY_icheck(span((p_k), (p_k+k_i(k)-1)), span(p_j,(p_j+k_i(j)-1))) = matVY_icheck(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))).t();
          }
        }
        p_k += k_i(k); // incrementation of p_k
      }
      p_j += k_i(j);// incrementation of p_j
    }
    
    vec ViY = vectorise(matVY_i);
    //vec ViY_check = vectorise(matVY_icheck);
    sigMSM =diagmat(MSigmaM);
    
  }else{
    int p_j=0; // loop variable
    int p_k =0; // loop variable

    for( int j =0 ; j < m_i; j++){
      p_k = p_j;
      // ###### computing of the matrix H_i_t_j  ##############
      mat matH_i_t_j = matHit(Ytildi.row(j).t()); // Yi.row(j) = vector of observations at time j;
      
      for( int k =j ; k < m_i; k++){
        // ###### computering of the matrix H_i_t_k ##############
        mat matH_i_t_k = matH_i_t_j;
        if(k!=j) {
          matH_i_t_k = matHit(Ytildi.row(k).t()); // Yi.row(k) = vector of observations at time k;
        }
        if(sum(sum(matH_i_t_k)) == 0){ // #### this mean that there is no observation at this time
          k_i(k) = 0;
        }
        else{
          k_i(k) = matH_i_t_k.n_rows;
        }
        p_k += k_i(k); // incrementation of p_k
      }
      p_j += k_i(j);// incrementation of p_j
    }
  }
  
  //double ya=0;
  //for( int j =0 ; j < ViY.size(); j++){
  //  ya += abs(ViY(j)-ViY_check(j));
  //}
  
  double lvrais = 0;
  double loglik_i=0;
  double loglik_i2=0;
  double loglik_i0 = 0;
  double log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi))));

  if((max(if_link) == 1 && !survival)|| check==1 ){
    // To check the linear closed form of likelihood
    double abs_det_matVY_i = abs(det(matVY_i));
    
    // ##### computering of the likelihood ########################## all linear
    loglik_i = -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i)) + log_Jac_Phi;
    loglik_i0 = -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i)) ;//+ log_Jac_Phi;

    //cout << " loglik_i "<<loglik_i<< " k_i "<<k_i.t();
    //mat checkB = matVY_i-sigMSM;
    //mat Vvi=sigMSM;
    loglik_i2 = -0.5*(sum(k_i)*log(2*M_PI) + log(det(sigMSM)) + as_scalar(Ytildi_nu_i.t()*inv_sympd(sigMSM)*Ytildi_nu_i));// + log_Jac_Phi;
    //cout << " loglik_i2 "<<loglik_i2<< " log_Jac_Phi "<< log_Jac_Phi <<endl;
    lvrais = loglik_i;
  }

  if(check==1 || max(if_link)>1 || survival){

    // MatrixXd precMat = matV_i;
    // LLT<MatrixXd> lltOfA(matV_i); // compute the Cholesky decomposition of A
    // MatrixXd L = lltOfA.matrixL();
    // double logDetPrecMat=  2*L.diagonal().array().log().sum();
    //fout2 << " logDetPrecMat " << logDetPrecMat<<endl;
    
    //precMat = L.inverse().transpose()*L.inverse();
    
    // Computation choleski
    //cout << " Add cholesky for V_Lambda_i !"<<endl;
    //mat Chol_matV_i = mvrnormArma(matV_i);//chol(matV_i); //chol.t()*chol = mat
    vec K2_lambda = zeros<vec>(K); // which latent process linked to the K markers
    
    for(int j =0 ; j < K; j++){
      int jjj=0;  
      while(matrixP(j,jjj)== 0){
        jjj ++;
      }
      K2_lambda[j]=jjj;
    }
    
    int sizeLambda=0;
    vec K2_lambda_t = zeros<vec>(sum(k_i)); // which latent process linked to the observation in lambda
    vec ind_lambda = zeros<vec>(sum(k_i));  // which indices of PNu_cp_i to include into lambda
    int ind=0;
    //sum(k_i) number of observations for subject i

    for(int j =0 ; j < Ytildi.n_rows; j++){
      vec markers = zeros<vec>(K);
      for(int b = 0 ; b < Ytildi.n_cols; b++){
        if(!isnan(Ytildi(j,b))){
          if(markers[b]==0){
            ind_lambda[sizeLambda]=ind;
            K2_lambda_t[sizeLambda]= K2_lambda[b];
            sizeLambda ++;
          }else{
            markers[b]=1;
          }
        }
        ind ++;
      }
    }
    
    vec paraSig(Ytildi.n_cols);
    
    int integral_ui = 1;
    lvrais = 0;
    
    if(integral_ui ==1 || survival){
      int nq = matDw_u.n_cols + matDw.n_cols ;
      mat var_RE = zeros<mat>(nq, nq);
      mat ui;
      //Verify it works with nD>1, may have to introduce nD somewhere
      //  matDw       matDw_u
      //  matDw_u.t() matDu
      for(int j =0 ; j < nq; j++){
        for(int jj =0 ; jj < nq; jj++){
          if(j < matDw.n_cols & jj < matDw.n_cols){
            var_RE(j, jj) = matDw(j, jj);
            
          }else if(j >= matDw.n_cols & jj < matDw.n_cols){
            var_RE(j, jj) = matDw_u(jj, j-matDw.n_cols);
            
          }else if(j < matDw.n_cols & jj >= matDw.n_cols){
            var_RE(j, jj) = matDw_u(j, jj-matDw.n_cols);
            
          }else if(j >= matDw.n_cols & jj >= matDw.n_cols){
            var_RE(j, jj) = matDu(j-matDw.n_cols, jj-matDw.n_cols);
            
          }
        }
      }

      mat chol_var_RE = chol(var_RE).t();

      if(type_int == -1){// MC
        mat uii = chol_var_RE * randn< Mat<double> >(chol_var_RE.n_rows, MCnr);
        ui = uii.t();
      }else if(type_int==0){ // AMC
        cout << " to develop !"<<endl;
      }else {//QMC
        ui = seq_i * chol_var_RE.t();
        //mat uii = chol_var_RE * randn< Mat<double> >(chol_var_RE.n_rows, MCnr);
        //ui = uii.t();
      }              

      //cout << " Vi "<< ZI*var_RE*ZI.t()<<endl<<endl;

      bool aMC=true;

      if(aMC){
        double out1;
        double out2;
        double vrais =0;
 
        //double log_Jac_Phi=0;
        double vraisr_surv=1;
        double vrais_surv_check=1;
        double vrais_survtot =0;
        double vraisY_tot =0;
        double min_lvraisr=0;
        double max_lvraisr=0;
        int pb_QMC = 0;

        //log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi))));
        
        for(int nr=0; nr < MCnr; nr++){
          
          int k_t=0;
          vec ui_r = ui.row(nr).t();
          if(nD>1)
            cout << " choisir les bons ui_r"<<endl;
          
          //ui_r.fill(0);
          //cout << " ui_r "<<ui_r.t();
          
          vec Lambda_nr;
          //vec Lambda_nr = matNui_ui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi, true);
     
          double lvraisr=0;
          int kk = 0;
          for (int k = 0 ; k < K; k++){
            vec ParaTransformYk = ParamTransformY(span(kk, (kk+df[k]-1)));
            vec tau_ik = matTik(Ytildi.col(k), tau_i);
            Lambda_nr = matNui_ui(nD, tau_ik, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi, true);
            vec Ytildik = YiwoNA(vectorise(Ytildi.col(k)));
            int nik = YiwoNA(Ytildik).size();
            mat Sig_k = eye(nik,nik)*Sig(k,k);
            
            if(type_int == -1){ //-1 MC 0 AMC 
              cout << " develop likelihood computation with integral ui for MC or AMC "<<endl;


                
                //Computation \Lambda_nr
                //vec Lambda_nr = lambda_ui(ui_r, Xi, beta, A, delta, gamma);
                //double out2 = f_marker(Lambda_nrMC, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParaTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
                //lvrais += out2;
              
              //lvrais /= MCnr;
              
            }else if(type_int > 0){//QMC`
              if(if_link(k)<2){// If linear or splines
                //vec Lambda_nr = matNui_ui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi);

                //vec Ytildi_nu_i_ui = vectorise(Ytildi)-Lambda_nr;
                vec Ytildi_nu_i_uik = Ytildik-Lambda_nr;
                out2 = -0.5*(nik*log(2*M_PI) + log(det(Sig_k)) + as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(Sig_k)*Ytildi_nu_i_uik));

                //cout << nik<<  " out2 "<< out2 <<" Sig_k "<< Sig_k(0,0) << " Ytildik" << Ytildik.t()<<" Ytildi_nu_i_uik "<< Ytildi_nu_i_uik.t();
              
                // if(nr==0){
                //   cout << " sum(k_i) "<< sum(k_i) 
                //        <<  " log(2*M_PI) "<< log(2*M_PI)
                //        <<" log(det(sigMSM)) "<< log(det(sigMSM)) 
                //        << " produit "<< as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(sigMSM)*Ytildi_nu_i_uik)<<endl;
                //   cout << " Ytildi_nu_i_uik "<< Ytildi_nu_i_uik.t()<<endl
                //        << " inv_sympd(sigMSM) "<< inv_sympd(sigMSM) <<endl<<endl;
                // }
                if(out2<=min_lvraisr){
                  min_lvraisr = out2;
                }
                
                if(out2>=max_lvraisr){
                  max_lvraisr = out2;
                }
                
                if(exp(out2)==0)
                  pb_QMC++;
                lvraisr += (out2);

              }else if(if_link(k)==2){ // thresholds

                double phi1;
                double phi2;
                double lvraisk=0;
                for (int j = 0 ; j < nik; j++){//NA!!
                  if(k_i(j)>0){ // change: if k_i(j,k)==1
                    
                    // double binf = ParaTransformYk(0);
                    // double bsup = ParaTransformYk(0);
                    // 
                    // for (int m = 0 ; m < (zitr(2*k_t+1)-zitr(2*k_t)+1); m++){
                    //   if(Ytildi(j, k)<(zitr(2*k_t) + m)){
                    //      binf = binf +  pow(ParaTransformYk_b(m), 2);
                    //   }
                    // }

                    int mm=0; //index in ParaTransformYk
                    double inf=-10;
                    double sup=-10;
                    
                    if(Ytildi(j, k)==(zitr(2*k_t))){
                      double value = (ParaTransformYk(0)-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm)-Lambda_nr(j))/abs(Sig(k,k));
                      phi1 = normalCDF(value);
                      phi2 = 0;
                    }else{
                      
                      inf=ParaTransformYk(0);
                      sup=ParaTransformYk(0);
                      for (int m = 0 ; m < (df[k]-1); m++){
                        inf = sup;
                        sup += pow(ParaTransformYk(m+1),2);
                        
                        if(Ytildi(j, k)==(zitr(2*k_t) + m+1)){
                          double value = (sup-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm)-Lambda_nr(j))/abs(Sig(k,k));
                          phi1 = normalCDF(value);
                          value = (inf-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm-1)-Lambda_nr(j))/abs(Sig(k,k));
                          phi2 = normalCDF(value);
                        }
                      }
                      if(Ytildi(j, k)==(zitr(2*k_t+1))){
                        phi1 = 1;
                        double value = (inf-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm-1)-Lambda_nr(j))/abs(Sig(k,k));
                        phi2 = normalCDF(value);
                      }
                    }
                    
                    if(phi1< phi2){
                      cout << " j "<< j //<< " PT "<< ParaTransformYk.t()<<endl
                           << " exp " <<(exp(phi1) - exp(phi2))
                           << " phi1 "<< phi1
                           << " phi2 "<< phi2 <<endl; 
                    }

                    lvraisk += log(phi1-phi2);

                    // if(phi1==phi2 || isinf(lvraisk) )
                    //   cout << " nr "<< nr <<" k " << k <<"j"<<j<<" Ytildi(j, k) "<<Ytildi(j, k)<< " k_i "<<k_i.t()
                    //        << " zitr(2*k_t) "<< zitr(2*k_t) << " zitr(2*k_t+1) "<< zitr(2*k_t+1) <<" inf "<<inf << " sup "<<sup
                    //        << " Lambda_nr(j) "<< Lambda_nr(j)<< " phi1 "<<phi1 << " phi2 "<<phi2<< " log(phi1-phi2) "<< log(phi1-phi2) <<" lvraisk "<<lvraisk << " lvraisr " << lvraisr<< " lvrais " << lvrais<< " vrais " << vrais<<endl;
                      
                    if(2<1){
                      double inf;
                      double sup;
                      for (int m = 0 ; m <= (zitr(2*k_t+1)-zitr(2*k_t)); m++){ // m values between Ymin and Ymax
                        cout << " m2 "<<m<<endl;
                        //if(ide(kk + max(mm-1, 0))==1){ // check ide contains ide of k=0, k=1, ...
                        if(Ytildi(j, k)==(zitr(2*k_t) + m)){
                          
                          if(m==0){
                            double value = (sup-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm)-Lambda_nr(j))/abs(Sig(k,k));
                            phi1 = normalCDF(value);
                            phi2 = 0;
                          }else if(m==(zitr(2*k_t+1)-zitr(2*k_t))){
                            phi1 = 1;
                            double value = (inf-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm-1)-Lambda_nr(j))/abs(Sig(k,k));
                            phi2 = normalCDF(value);
                          }else{
                            double value = (sup-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm)-Lambda_nr(j))/abs(Sig(k,k));
                            phi1 = normalCDF(value);
                            value = (inf-Lambda_nr(j))/abs(Sig(k,k));//(ParaTransformYk(mm-1)-Lambda_nr(j))/abs(Sig(k,k));
                            phi2 = normalCDF(value);
                          }
                          
                          if(phi1< phi2){
                            cout << " j "<< j //<< " PT "<< ParaTransformYk.t()<<endl
                                 << " exp " <<(exp(phi1) - exp(phi2))
                                 << " phi1 "<< phi1
                                 << " phi2 "<< phi2 <<endl; 
                          }
                          
                          //vrais *= (exp(phi1) - exp(phi2));
                          // cout << " j "<< j << " phi1 "<< exp(phi1) - exp(phi2) <<  " phi1 "<<exp(phi1) << " phi2 "<< exp(phi2) <<endl
                          //     <<" ParaTransformYk(m+1) " << ParaTransformYk(m+1) <<" ParaTransformYk(m) " << ParaTransformYk(m) << "Sig(k,k) "<< Sig(k,k) <<endl;
                          
                          //cout << "m "<<m<<" ParaTransformY(ind_m) "<<ParaTransformY(ind_m) << " lambda "<<Lambda_nr(j)<<endl;
                          lvraisk += log(phi1-phi2);
                          
                          if(nr<0){
                            cout << nr << " m "<<m<<" lvraisr " << lvraisr << " phi1 "<< phi1 << " phi2 "<< phi2<<endl; 
                          }
                          //if(nr==10)
                          //  cout << nr << " ParaTransformY "<< ParaTransformY.t()<<endl; 
                          
                        }//if Y=m
                        mm++;
                      } // for m in Mk
                    }//2<1
 
                  }// if Yj observed
                } // for j
                k_t ++;
                lvraisr += lvraisk;
              }// if threshold
            } // if QMC
            //if(nr < 10){
            //  cout << " out2 "<< out2 <<endl;
            //        << " pi " << sum(k_i)*log(2*M_PI) 
            //        << " logdet " << log(det(sigMSM)) 
            //        << " exp "<<  as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(sigMSM)*Ytildi_nu_i_uik)
            //        << " ui_r "<< ui_r.t()
            //       << " vectorise(Ytildi) "<< vectorise(Ytildi).t()
            //        << " Lambda_nr.t() "<< Lambda_nr.t()
            //       << "Ytildi_nu_i_uik.t()  "<< Ytildi_nu_i_uik.t()<<endl;
            //}
            //double out2 = f_marker(Lambda_nr, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParaTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
            kk += df[k];
          }//k


          if(survival){
            vraisr_surv = f_survival_ui(ui_r, t_0i, t_i, delta_i, xti1, xti2, param_surv, param_basehaz, basehaz, knots_surv, assoc, truncation,
                                        nD, tau, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, nE);
          }

          vrais_survtot += vraisr_surv;
          vraisY_tot += exp(lvraisr);
          //vrais += exp(lvraisr)*vraisr_surv;
          vrais += exp(lvraisr+log(vraisr_surv));

          if(2<1){//verification survival likelihood if all regression parameters  =0 and baseline = Weibull
            double s1=exp(-pow(t_i/param_basehaz(1),param_basehaz(0))-pow(t_i/param_basehaz(3),param_basehaz(2)));
            double lambdat=1;
            if(delta_i==1)
              lambdat = param_basehaz(0)/param_basehaz(1)*pow(t_i/param_basehaz(1),param_basehaz(0)-1);
            if(delta_i==2)
              lambdat = param_basehaz(2)/param_basehaz(3)*pow(t_i/param_basehaz(3),param_basehaz(2)-1);
            double s0=exp(-pow(t_0i/param_basehaz(1),param_basehaz(0))-pow(t_0i/param_basehaz(3),param_basehaz(2)));
            vrais_surv_check=s1/s0*lambdat;
            //cout<< " vrais_surv_check "<<vrais_surv_check<< " s1 "<<s1 << " s0 "<<s0 << " lambdat "<< lambdat<<endl;
          }
        }//nr
        
        if(check){
          double minY = pow(10,10);
          for(int j = 0 ; j < Ytildi.size(); j++){
            if(Ytildi(j)<minY)
              minY = Ytildi(j);
          }
          //if(abs(loglik_i- log(vraisY_tot/MCnr) )>1)
            cout << " diffY "<<loglik_i-log_Jac_Phi- log(vraisY_tot/MCnr)<< " MCnr "<<MCnr<<" pb_QMC" <<pb_QMC<< " minY "<< minY << " vrais / MCnr "<<vrais / MCnr;
            cout << " loglik_i0 "<< loglik_i0<< " loglik_i2 "<<loglik_i2<<" log(vraisY_tot/MCnr) "<< log(vraisY_tot/MCnr)<< " log_Jac_Phi "<<log_Jac_Phi << " vrais "<<vrais<<endl;
            
            //cout << " diffT "<<log(vrais_surv_check)-log(vrais_survtot/MCnr) <<" vrais_t_check "<<vrais_surv_check
                 //<< " vrais_survtot/MCnr "<<vrais_survtot/MCnr<< // << " vraisTi "<< vrais_survtot/MCnr << endl;
          
            //cout << " loglik_i "<<loglik_i<< " log(vraisY_tot/MCnr) "<< log(vraisY_tot/MCnr)<< " diffT "<<log(vrais_surv_check)-log(vrais_survtot/MCnr) <<" pb_QMC" <<pb_QMC<< " minY "<< minY << " vrais / MCnr "<<endl;// << " vraisTi "<< vrais_survtot/MCnr << endl;
          //cout << " max_lvraisr "<< max_lvraisr << " min_lvraisr "<< min_lvraisr << " pb_QMC" <<pb_QMC<< " log_Jac_Phi "<< log_Jac_Phi <<endl;
        }

        
        vrais /= MCnr;
        lvrais += log(vrais) + log_Jac_Phi;
        
        if(pb_QMC > MCnr*0.15*K)
          lvrais=-pow(10,10);
      }
      

      
    }else{ //Integral on lambda ------------------
      ind_lambda = ind_lambda.head(sizeLambda);
      int ind2=0;
      ind=0;
      vec mean_lambdai(sizeLambda);
      
      if(sizeLambda != PNu_cp_i.size()){ // if one Y obs = nan or there are doublons in lambda_k(t)
        
        int ind3=0;
        for(int b = 0 ; b < PNu_cp_i.size(); b++){
          if(binary_search(ind_lambda.begin(), ind_lambda.end(), b)){
            mean_lambdai[ind3]= PNu_cp_i[b];
            ind3++;
          }else{
            matV_i.shed_col(b);
            matV_i.shed_row(b);
          }
        }
      }else{
        mean_lambdai=PNu_cp_i;
      }
      
      //LLT<MatrixXd> lltOfA(precMat); // compute the Cholesky decomposition of A
      //MatrixXd L = lltOfA.matrixL();
      //logDetPrecMat=  2*L.diagonal().array().log().sum();
      //precMat = L.inverse().transpose()*L.inverse();
      // cout << " det(matV_i) "<< det(matV_i) <<endl;
      // mat U;
      // vec s;
      // mat V;
      // 
      // svd(U,s,V,matV_i);
      // 
      // mat X= U*diagmat(s)*V.t();
      // 
      // cout << " matV_i-X "<< matV_i-X <<endl;
      
      //cout << " chol(matV_i) "<< chol(matV_i) <<endl;
      mat D;
      
      const bool chol_status = op_chol::apply_direct(D, matV_i, 1);  // '1' means "lower triangular"
      
      if(chol_status == false)
      {
        // C is not symmetric positive definite, so find approximate square root of C
        
        colvec eigval;  // NOTE: eT is constrained to be real (ie. float or double) in fn_mvnrnd.hpp
        mat eigvec;
        
        const bool eig_status = eig_sym_helper(eigval, eigvec, matV_i, 'd', "mvnrnd()");
        
        if(eig_status == false)  { return false; }
        
        double*   eigval_mem    = eigval.memptr();
        const uword eigval_n_elem = eigval.n_elem;
        
        // since we're doing an approximation, tolerate tiny negative eigenvalues
        
        const double tol = double(-100) * Datum<double>::eps * norm(matV_i, "fro");
        
        if(arma_isfinite(tol) == false)  { return false; }
        
        for(uword i=0; i<eigval_n_elem; ++i){
          const double val = eigval_mem[i];
          
          if( (val < tol) || (arma_isfinite(val) == false) )  { return false; }
        }
        
        for(uword i=0; i<eigval_n_elem; ++i)  { if(eigval_mem[i] < double(0))  { eigval_mem[i] = double(0); } }
        
        mat DD = eigvec * diagmat(sqrt(eigval));
        D.steal_mem(DD); //matV_i = D*D.t();
      }
      
      
      mat MatVi2=matV_i+0.0001*eye(matV_i.n_rows, matV_i.n_rows);
      mat invmatV_i= inv_sympd(MatVi2);
      mat D_bis = D+0.0001*eye(matV_i.n_rows, matV_i.n_rows);
      
      double log_det=log(det(MatVi2));
      //double log_det_D=2*log(det(D_bis));
      //cout << " log_det "<< log_det << " vs "<< log(det(matV_i))<< " vs " <<log_det_D << endl;
      // double logDetPrecMat=  2*L.diagonal().array().log().sum();
      //fout2 << " logDetPrecMat " << logDetPrecMat<<endl;
      
      //precMat = L.inverse().transpose()*L.inverse();
      //mat out = D * randn< Mat<double> >(D.n_rows, MCnr);
      
      mat Lambda;
      
      if(type_int<=0){// MC (-1) or aMC (0) (QMC = 1)
        Lambda = D * randn< Mat<double> >(D.n_rows, MCnr);
      }else{ // QMC
        //Lambda =  (D* seq_i.t()).t();
        Lambda =  seq_i;
      }
      
      bool aMC=true;
      double lvrais2 =0;
      
      if(aMC){
        int MCnr_d;
        
        if(type_int == -1){ //MC
          
          for(int nr=0; nr < 0; nr++){
            vec zero_lambda = zeros<vec>(mean_lambdai.size());
            //vec Lambda_nrMC = mvnrnd(mean_lambdai, matV_i);//chol(matV_i); //chol.t()*chol = mat
            vec Lambda_nrMC(Lambda.n_rows);// changer ! matrice nT x nD
            for(int i=0; i < Lambda.n_rows; i++){
              Lambda_nrMC(i) = mean_lambdai(i) + Lambda(i,nr);
            }
            
            //Lambda_nr = mvnrnd(mean_lambdai, matV_i);
            double out2 = f_marker(Lambda_nrMC, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParamTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
            lvrais += out2;
          }
          lvrais /= MCnr;
        }else if(type_int > 0 ){ //QMC
          MCnr_d = MCnr;
          double cstant = 0;
          vec lvrais_nr(MCnr);
          // Approximation of the likelihood by antithetic MC
          //Number of pairs of draws for antithetic MC
          //cout << " det(matV_i) "<<det(matV_i)<<endl;
          //vec s = svd( matV_i );
          //  cout << " svd "<<s.t()<<endl;
          //type_int = 0 if antithetic MC, 1 if halton, 2 if sobol
          for(int nr=0; nr < MCnr; nr++){
            vec zero_lambda = zeros<vec>(mean_lambdai.size());
            
            //vec Lambda_nrMC = mvnrnd(mean_lambdai, matV_i);//chol(matV_i); //chol.t()*chol = mat
            vec Lambda_nrMC(Lambda.n_cols);
            vec resid_nrMC(Lambda.n_cols);
            for(int i=0; i < Lambda.n_cols; i++){
              Lambda_nrMC(i) = mean_lambdai(i) + Lambda(nr,i);
              resid_nrMC(i) = Lambda(nr,i);
            }
            
            //Lambda_nr = mvnrnd(mean_lambdai, matV_i);
            double logout2 = f_marker(Lambda_nrMC, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParamTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
            
            
            double test = -0.5*(Lambda_nrMC.size()*log(2*M_PI) + log_det + as_scalar(resid_nrMC.t()*invmatV_i*resid_nrMC)) ;
            lvrais_nr(nr)=logout2;
            //if(type)
            // Sampling of Lambda_i
            //arma::vec sampled = mvnorm((nr+1)*666, zeros(matV_i.n_cols), eye(matV_i.n_cols, matV_i.n_cols));
            //arma::vec Lambda_nr = PNu_cp_i + Chol_matV_i*sampled;
            //vec Lambda_nr0 = mvnrnd(zero_lambda, matV_i);//chol(matV_i); //chol.t()*chol = mat
            //for(int i=0; i < Lambda.n_cols; i++){
            //  Lambda_nrMC(i) = mean_lambdai(i) - Lambda(nr,i);
            //}      //Lambda_nr = mvnrnd(mean_lambdai, matV_i);
            
            //double out = f_marker(Lambda_nrMC, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParaTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
            //vrais += out;
          }
          
          //cstant=max(lvrais_nr)-700;//(min(lvrais_nr)+max(lvrais_nr))/2;
          cstant=0;
          double vrais=0;
          //cstant=0;
          int n0=0;
          for(int nr=0; nr < MCnr_d; nr++){
            vrais += exp(lvrais_nr(nr)-cstant);
            if(exp(lvrais_nr(nr)-cstant)<pow(10,-5)){
              n0 +=1;
            }
          }
          
          cout << " n0 "<< n0<< " cstant "<< cstant <<
            " vrais "<< vrais <<
              " log(vrais) "<< log(vrais) << endl;
          
          lvrais =log(vrais)+cstant;
          lvrais2 =lvrais - log(MCnr_d-n0);
          lvrais -=log(MCnr/2);
        }
        
      }
    }// --- end if integration on lambda
  }
  return(lvrais);
}



/*
 individual contribution to the log-likelihood for the latent processes
 */

double Loglikei_latent(arma::colvec& Latent_i, int& K, int& nD, arma::mat& matrixP, int& m_i, arma::vec& tau, arma::vec& tau_i, arma::mat& Ytildi, arma::mat& YtildPrimi,
                       arma::mat& x0i, arma::mat& z0i, arma::mat& xi, arma::mat& zi, arma::colvec& alpha_mu0,
                       arma::colvec& alpha_mu, arma::mat& matDw,  arma::mat& matDw_u, arma::mat& matDu,
                       arma::mat& matB, arma::mat& Sig, arma::mat& G_mat_A_0_to_tau_i,
                       arma::mat& G_mat_prod_A_0_to_tau,  double& DeltaT){
  
  
  double loglik=0.e0;
  mat Vi; 
  //= Var_latent_i(K, nD, matrixP, m_i, tau, tau_i, Ytildi, YtildPrimi,
  //              x0i, z0i, xi, zi, alpha_mu0,
  //              alpha_mu, matDw,  matDw_u, matDu,
  //              matB, Sig, G_mat_A_0_to_tau_i,
  //              G_mat_prod_A_0_to_tau, DeltaT);
  double abs_det_matV_i = abs(det(Vi));
  
  vec mu_i = matNui(nD, tau_i, DeltaT, x0i, alpha_mu0,
                    xi, alpha_mu,G_mat_A_0_to_tau_i);
  
  loglik = 0;//-0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matV_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i));
  
  return(loglik);
}

bool compFun2(int i) {
  return i <2;
}
bool compFun1(int i) {
  return i == 2;
}


//===================================================================================================
/* ******************************************************
 Function f_Loglik: log-likelihood of the observed data
 */
//===========================================================================================
//' Function that computes the log-likelihood of the observed data
//'  
//' @param K an integer indicating the number of markers
//' @param nD an integer indicating the number of latent processes
//' @param mapping indicates which outcome measured which latent process, it is a mapping table between
//' outcomes and latents processes
//' @param paraOpt initial values for model parameters
//' @param paraFixe values associated to the index of parameters to be constrained
//' @param posfix position of parameters to be constrained
//' @param m_is vector of numbers of visit occasions for individuals
//' @param Mod_MatrixY model.matrix from markers transformation submodels
//' @param Mod_MatrixYprim model.matrix from the derivates of markers transformation submodels
//' @param df vector of numbers of parameters for each transformation model
//' @param nb_paraD number of paramerters of the variance-covariance matrix of random effects
//' @param x0 model.matrix for baseline's fixed submodel
//' @param x model.matrix for change's fixed submodel
//' @param z0 model.matrix for baseline's random effects submodel
//' @param z model.matrix for change's random effects submodel
//' @param q0 a vector of number of random effects on each initial latent process level
//' @param q a vector of number of random effects on each change latent process over time
//' @param if_link indicates if non linear link is used to transform an outcome
//' @param tau a vector of integers indicating times (including maximum time)
//' @param tau_is a vector of integers indicating times for individuals
//' @param modA_mat model.matrix for elements of the transistion matrix
//' @param DeltaT double that indicates the discretization step  
//' @param data_surv matrix of Tentry, Event, StatusEvent and covariates for survival models
//' @param basehaz baseline hasard function type
//' @param knots_surv knots for splines modelling the baseline hazard function
//' @param np_surv number of parameters in the survival sub-models !! change if nE>1 !!
//' @param survival boolean indicating if survival analysis
//' @param assoc type of association between outcomes and times-to-events
//' @param truncation boolean indicating if left trucation 
//' @param nE number of events
//' @param Xsurv1 design matrix for first event
//' @param Xsurv2 design matrix for second event
//' @param zitr min and max of ordinal outcomes
//' @param ide vector of observed values for ordinal outcomes
//' @param modA_mat_predGK_t design matrix for computing predictions of Y on [0;ti] in Gauss Konrod for all subject
//' @param modA_mat_predGK_t0 design matrix for computing predictions of Y on [0;t0i] in Gauss Konrod for all subject
//' @param pt_GK_t Gauss-Konrod nodes for integration on [0;ti] for all subject
//' @param pt_GK_t0 Gauss-Konrod nodes for integration on [0;t0i] for all subject
//' @return double 
//' @export
//' 
// [[Rcpp::export]]
double Loglik(int K, int nD, arma::vec& mapping, arma::vec& paraOpt, arma::vec& paraFixe, arma::vec& posfix, 
              arma::vec& paras_k, arma::mat& sequence, int type_int, arma::vec& ind_seq_i, int MCnr,  arma::vec& nmes, arma::vec& m_is,
              arma::mat& Mod_MatrixY, arma::mat& Mod_MatrixYprim, arma::vec& df, arma::mat& x,
              arma::mat& z, arma::vec& q, int nb_paraD, arma::mat& x0, arma::mat& z0,
              arma::vec& q0, 
              arma::mat& data_surv, int basehaz, arma::vec& knots_surv, arma::vec& np_surv, bool survival, int assoc, bool truncation,
              int nE, arma::mat& Xsurv1, arma::mat& Xsurv2,
              arma::vec& if_link, arma::vec& zitr, arma::vec& ide, 
              arma::vec& tau, arma::vec& tau_is, arma::mat& modA_mat, double DeltaT){

  double loglik = 0.e0;
  double loglik2 = 0.e0;
  arma::mat matrixP = zeros(K,nD);
  for(int k = 0; k<K; k++){
    matrixP(k,(mapping(k)-1)) = 1.e0;
  }
  int m = (int)tau.size();
  int N = (int)m_is.size();
  int p=0; //loop to read observations Y
  int ncol_x = (int)x.n_cols; // number of parameters for the mean slope (DeltaX)
  int ncol_x0 = (int)x0.n_cols; // number of parameters for mean of processes at baseline (X0)
  int ncol_z = (int)z.n_cols; // number of parameters for randoms effects on the slope
  int ncol_z0 = (int)z0.n_cols; // number of parameters for randoms effects on baseline processes values
  int L = (int)modA_mat.n_cols; // number of parameters for the transition matrix A
  
  // Identification of constraint parameters and non-constraint parameters
  int Nb_para = (int)posfix.size();
  vec paras = zeros(Nb_para);
  int Opt=0;
  int Fixe=0;
  for(int i=0; i < Nb_para ; i++){
    if(posfix[i]==0){
      paras[i] = paraOpt[Opt];
      Opt +=1;
    }
    else{
      paras[i] = paraFixe[Fixe];
      Fixe +=1;
    }
  }
  
  //Identification of groups of parameters
  int ipara =0;
  colvec alpha_mu0 = paras(span(ipara,ipara+ncol_x0-1));
  ipara += ncol_x0;
  colvec alpha_mu = paras(span(ipara,ipara+ncol_x-1));
  ipara += ncol_x;
  colvec alpha_D = paras(span(ipara,ipara + nb_paraD-1));
  ipara += nb_paraD;
  vec vec_alpha_ij = paras(span(ipara,ipara+L*nD*nD-1));
  ipara += L*nD*nD;
  vec paraB = zeros(nD);
  vec paraSig = paras(span(ipara,ipara+K-1));
  ipara += K;
  int nbParaTransformY = Mod_MatrixY.n_cols;
  colvec ParaTransformY = paras(span(ipara,ipara+nbParaTransformY-1));
  ipara += nbParaTransformY;
  
  int nq_s=0; // number parameters for baseline function
  colvec param_surv = zeros<vec>(sum(np_surv));
  colvec param_basehaz;
  if(survival){
    if(basehaz==0){
      nq_s = 2;
    }else{
      nq_s = knots_surv.size();
    }
    
    param_surv = zeros<vec>(sum(np_surv));
    param_basehaz = zeros<vec>(nq_s*nE);
    
    for(int k=0; k<nE;k++){
      param_basehaz(span(k*nq_s, (nq_s-1)*(1-k) + (nq_s*nE-1)*k)) = paras(span(ipara,ipara + nq_s -1));
      param_surv(span(k*np_surv(0), (np_surv(0)-1)*(1-k) + (sum(np_surv)-1)*k)) = paras(span(ipara + nq_s,ipara + nq_s + np_surv(k)-1));
      ipara += (np_surv(k) + nq_s);
    }
  }

  int nb_RE = sum(sum(q0)+sum(q));
  Mat<double> matD = DparChol(nb_RE, alpha_D);
  int n_cols_matD = matD.n_cols;
  Mat<double> matDw = matD(span(0,nD-1),span(0,nD-1));
  Mat<double> matDw_u = matD(span(0,nD-1),span(nD,n_cols_matD-1));
  Mat<double> matDu = (matD(span(nD,n_cols_matD-1),span(nD,n_cols_matD-1)));
  Mat<double> matB = KmatDiag(paraB); // structured variance-covariance matrice
  Mat<double> Sig = KmatDiag(paraSig); // noise

  // Computering of Ytild and YtildPrim, the transformed marker and it derivate from
  // Y, model.matrix and transformation parameters
  Mat<double> Ytild = zeros(Mod_MatrixY.n_rows,K);
  Mat<double> YtildPrim = zeros(Mod_MatrixYprim.n_rows,K);
  
  // re-writing transformed parameters
  // if_link = 1 : mean non linear link

  int kk = 0;
  int kkp = 0;
  colvec ParamTransformY(sum(df)); // transformed parameters

  for(int k=0; k<K;k++){
    if(if_link[k] == 1){ //splines
      ParamTransformY(span(kk, kk+df[k]-1)) = exp(ParaTransformY(span(kk, (kk+df[k]-1))));
      ParamTransformY[kk] = log(ParamTransformY[kk]);
      Ytild.col(k) = Mod_MatrixY.cols(kk, (kk+df[k]-1))*ParamTransformY(span(kk, kk+df[k]-1));
      YtildPrim.col(k) = Mod_MatrixYprim.cols(kkp, (kkp+df[k]-2))*ParamTransformY(span(kk+1, kk+df[k]-1));
      kk += df[k];
      kkp += df[k]-1;
    }else if(if_link[k] == 2){ //thresholds
      
      Ytild.col(k) = Mod_MatrixY.col(kk);
      ParamTransformY(kk) = ParaTransformY(kk);
      for(int j= 1; j < df[k]; j++ ){
        //ParamTransformY(kk+j)= ParamTransformY(kk+j-1) + ParaTransformY(kk+j)*ParaTransformY(kk+j);
        ParamTransformY(kk+j)= ParaTransformY(kk+j);
      }

      YtildPrim.col(k).fill(1);
      
      kk += df[k];
      kkp += df[k]-1;
    }else{ // linear
      ParamTransformY(span(kk, (kk+df[k]-1))) = ParaTransformY(span(kk, (kk+df[k]-1)));
      ParamTransformY[kk] = - ParamTransformY[kk]/ParamTransformY[kk+1];
      ParamTransformY[kk+1] = 1.e0/ParamTransformY[kk+1];
      
      Ytild.col(k) = Mod_MatrixY.cols(kk, (kk+df[k]-1))*ParamTransformY(span(kk, (kk+df[k]-1)));
      //cout << " Y "<<Mod_MatrixY.cols(kk, (kk+df[k]-1))<<endl;

      YtildPrim.col(k) = ParamTransformY[kk+1]*Mod_MatrixYprim.cols(kkp, (kkp+df[k]-2));
      kk += df[k];
      kkp += df[k]-1;
    }
  }
  
  //Computering of log-likelihood as sum of individuals contributions
  double loglik0=0;
  for(int n= 0; n < N ; n++ ){

    //if(n%200==0)
      //cout << "indiv "<< n <<endl;
    // printf("\n %d \n",(n+1));
    //Creation of matrix G_mat_prod_A_0_to_tau that contains all products  A(j) from t_i a Tmax: t_i \in 0, Tmax
      mat G_mat_prod_A_0_to_tau = GmatprodAstotau(nD, vec_alpha_ij, tau, 0, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));
      mat G_mat_A_0_to_tau_i = GmatA0totaui(nD, vec_alpha_ij, tau, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));
    //if( std::all_of(if_link.begin(), if_link.end(), compFun2) ){
    

    //mat G_mat_A_0_to_t_i = GmatA0totaui(nD, vec_alpha_ij, deltaT_ptGK_ti(span(n*15, (n+1)*15-1), DeltaT, modA_mat_predGK_ti(span(n*15,((n+1)*15-1)), span(0,(L-1))));
    //mat G_mat_A_0_to_t_0i = GmatA0totaui(nD, vec_alpha_ij, deltaT_ptGK_t0i(span(n*15, (n+1)*15-1), DeltaT, modA_mat_predGK_t0i(span(n*15,((n+1)*15-1)), span(0,(L-1))));
    
    double t_i=0.0; 
    double t_0i=0.0; 
    double delta_i=0.0;
    vec xti1 = zeros<vec>(Xsurv1.n_cols); 
    vec xti2 = zeros<vec>(Xsurv2.n_cols); 

    if(survival){
      t_0i = data_surv(n,0);
      t_i = data_surv(n,1);
      delta_i = data_surv(n, 2);
      
      int n_assoc = 1;
      if(assoc == 2 || assoc == 5){
        n_assoc ++;
      }
      
      for(int k= 0; k < Xsurv1.n_cols; k++ ){
        xti1(k) = Xsurv1(n, k);
      }
      
      if(nE==2){
        for(int k= 0; k < Xsurv2.n_cols; k++ ){
          xti2(k) = Xsurv2(n, k);
        }
      }

      //modA_mat_predGK_ti = modA_mat_predGK_t(span(n*15, (n+1)*15-1), span(0, modA_mat_predGK_t.n_cols-1));//zeros(15, modA_mat_predGK_t.n_cols);
      //modA_mat_predGK_t0i = modA_mat_predGK_t0(span(n*15, (n+1)*15-1), span(0, modA_mat_predGK_t0.n_cols-1));//zeros(15, modA_mat_predGK_t0.n_cols);
    }

    //  double out0 = Loglikei(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
    //                      YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
    //                      z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
    //                      z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
    //                      matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT);
    // loglik0 += out0;
    //}
    // }else if( std::all_of(if_link.begin(), if_link.end(), compFun1) ){
    //  std::cout << "All the elements are equal to 2.\n";
    //cout << " n "<<n;
    double out1 =0;
    out1= Loglikei_GLM(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
                                 YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
                                 z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
                                 z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
                                 matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT, ParamTransformY, df, if_link, zitr, ide, paras_k,
                                 t_0i, t_i, delta_i, xti1, xti2, basehaz, knots_surv, survival, param_surv, param_basehaz, assoc, truncation,
                                 sequence, type_int, ind_seq_i, MCnr, n, nE);
      
      loglik += out1;

    
    

    //cout << " n "<<n << " out1 "<<out1<<endl;
    

    // double out2 =  Loglikei(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
    //                     YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
    //                     z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
    //                     z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
    //                     matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT);
    // loglik2 += out2;
    
    
    p += m_is[n];
  }
  //cout << "   "<<loglik0 << " "<< loglik<<endl;
  
  // int n = 0;
  // mat G_mat_prod_A_0_to_tau = tsGmatprodA0totau(K, vec_alpha_ij, tau, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));
  // mat G_mat_A_0_to_tau_i = GmatA0totaui(K, vec_alpha_ij, tau, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));
  // double  loglik = Loglikei(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
  //                         YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
  //                         z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
  //                         z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
  //                         matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT);

  return(loglik);
  // return(YtildPrim); 
}
// 
// arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma, int chol);
// RcppExport SEXP ConservativeEstimates_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP cholSEXP) {
//   BEGIN_RCPP
//   Rcpp::RObject __result;
//   Rcpp::RNGScope __rngScope;
//   Rcpp::traits::input_parameter< int >::type n(nSEXP);
//   Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
//   Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
//   Rcpp::traits::input_parameter< int >::type chol(cholSEXP);
//   __result = Rcpp::wrap(mvrnormArma(n, mu, sigma, chol));
//   return __result;
//   END_RCPP
// }
// // trmvrnorm_rej_cpp
// arma::mat trmvrnorm_rej_cpp(int n, arma::vec mu, arma::mat sigma, arma::vec lower, arma::vec upper, int verb);
// RcppExport SEXP ConservativeEstimates_trmvrnorm_rej_cpp(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP verbSEXP) {
//   BEGIN_RCPP
//   Rcpp::RObject __result;
//   Rcpp::RNGScope __rngScope;
//   Rcpp::traits::input_parameter< int >::type n(nSEXP);
//   Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
//   Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
//   Rcpp::traits::input_parameter< arma::vec >::type lower(lowerSEXP);
//   Rcpp::traits::input_parameter< arma::vec >::type upper(upperSEXP);
//   Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
//   __result = Rcpp::wrap(trmvrnorm_rej_cpp(n, mu, sigma, lower, upper, verb));
//   return __result;
//   END_RCPP
// }
