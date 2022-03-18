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

double Loglikei_GLM(int K, int nD, arma::mat& matrixP, arma::vec& mapping, int m_i, arma::vec& tau, arma::vec tau_i, arma::mat Ytildi, arma::mat YtildPrimi,
                    arma::mat x0i, arma::mat z0i, arma::mat xi, arma::mat zi, arma::colvec& alpha_mu0,
                    arma::colvec& alpha_mu, arma::mat& matDw,  arma::mat& matDw_u, arma::mat& matDu,
                    arma::mat& matB, arma::mat& Sig, arma::mat& G_mat_A_0_to_tau_i,
                    arma::mat& G_mat_prod_A_0_to_tau,  double DeltaT, arma::vec& ParamTransformY, arma::vec& df, arma::vec& if_link,  
                    arma::vec& zitr,  arma::mat& ide, arma::vec& paras_k,
                    double t_0i, double t_i, int delta_i, arma::vec& xti1, arma::vec& xti2, arma::mat& xti1_intY, arma::mat& xti2_intY, int basehaz, arma::vec& knots_surv, bool survival, arma::vec& param_surv, arma::vec&  param_surv_intY, arma::vec& param_basehaz, int assoc, bool truncation, 
                    arma::mat& seq_i,  int type_int, arma::vec& ind_seq_i, int MCnr, int sub, int nE, int add_diag_varcov, arma::vec& q_nD){
  
  vec Ytildi_nu_i;
  mat matVY_i;
  vec k_i = zeros<vec>(m_i);// vector of number of observed marker at each observation time
  vec PNu_cp_i;
  mat sigMSM;
  double expotrick=0; // CHANGE
  if(nD==1)
    expotrick=0;
  
  int check=2; // 1 both QMC and closed likelihood for comparison of individual likelihood (type_int, MCnr must be defined!) 
               // 2 close likelihood if links = linear/splines (except if MCnr>0), QMC if links = thresholds or if survival = T
               // 3 MC integration (even if links = linear/splines and survival = F)
  int printa=0; //1 to print some intermediate values
  int q = sum(q_nD);
  
  if((max(if_link) < 2 && !survival)&& check <3 || check==1){

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
    //mat matVY_icheck = zeros(sizeYi,sizeYi); // initialization of the  matVY_i to zero
    //int q = zi.n_cols; // number of random effets on the slope
    
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
    vec MSigmaM=zeros<vec>(sizeYi);
    int aa = 0;
    
    // ##### computering of GrdZi ####################################
    //  GrdZi : this matrix contains all model.matrix at each time of tau_i
    if(q>0){
      for(int t = 0; t<= maxTau_i; t++){
        
        if(t==0){
          GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = 0.e0*zi(span(t*nD,(t+1)*nD-1), span(0,q-1));
        }
        else{
          GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = (DeltaT*zi(span(t*nD,(t+1)*nD-1), span(0,q-1)) +
            G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),t*nD-1))*GrdZi(span((t-1)*nD,t*nD-1), span(0,q-1)));
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
          // matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
          //   ((phi_0_j_0*z0i)*matDw* (z0i*phi_0_k_0).t() +
          //   (phi_0_j_0*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
          //   (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*phi_0_k_0).t() +
          //   (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t()
          //   )*matrixP.t())*matH_i_t_k.t();
          
          matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
            ((phi_0_j_0*z0i)*matDw* (z0i*phi_0_k_0).t())*matrixP.t())*matH_i_t_k.t();
          
          
          if(q>0){
            matVY_i(span(p_j,(p_j+k_i(j)-1)), span((p_k), (p_k+k_i(k)-1))) += matH_i_t_j*( matrixP*
              ((phi_0_j_0*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
              (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*phi_0_k_0).t() +
              (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t()
              )*matrixP.t())*matH_i_t_k.t();
          }
          //  matDw       matDw_u
          //  matDw_u.t() matDu
          
          
          //// For simplicty we define GrdZi(span(0*K,(0+1)*K-1), span(0,q-1))) = 0 : as if we specify
          // model for slope from t_0 but with covariate set to 0
          
          if( k == j){
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
    
  }

  if(max(if_link) == 2 || survival || MCnr>0){ //if check==1, already done
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
  
  double lvrais = 0;
  double loglik_i=0;
  double loglik_i2=0;
  double loglik_i0 = 0;
  double log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi))));

  if((max(if_link) < 2 && !survival)&& check <3 || check==1){

    // To check the linear closed form of likelihood
    double abs_det_matVY_i = abs(det(matVY_i));
    
    // ##### computering of the likelihood ########################## all linear
    loglik_i = -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i)) + log_Jac_Phi;

    //loglik_i0 = -0.5*(sum(k_i)*log(2*M_PI) + log(abs_det_matVY_i) + as_scalar(Ytildi_nu_i.t()*inv_sympd(matVY_i)*Ytildi_nu_i)) ;//+ log_Jac_Phi;

    //loglik_i2 = -0.5*(sum(k_i)*log(2*M_PI) + log(det(sigMSM)) + as_scalar(Ytildi_nu_i.t()*inv_sympd(sigMSM)*Ytildi_nu_i)) + log_Jac_Phi;
    lvrais = loglik_i;
    
  }
  
  if(max(if_link) == 2 || survival || check==1 || MCnr>0){
    
    vec K2_lambda = zeros<vec>(K); // which latent process linked to the K markers
    K2_lambda = mapping-1;
    // for(int j =0 ; j < K; j++){
    //   int jjj=0;  
    //   while(matrixP(j,jjj)== 0){
    //     jjj ++;
    //   }
    //   K2_lambda[j]=jjj;
    // }

    int sizeLambda=0;
    vec K2_lambda_t = zeros<vec>(sum(k_i)); // which latent process linked to the observation in lambda
    int ind=0;
    //sum(k_i) number of observations for subject i
    
    for(int j =0 ; j < Ytildi.n_rows; j++){
      for(int b = 0 ; b < Ytildi.n_cols; b++){
        if(!isnan(Ytildi(j,b))){
          K2_lambda_t[sizeLambda]= K2_lambda[b];
          sizeLambda ++;
        }
        ind ++;
      }
    }

    lvrais = 0;
    int nq = matDw.n_cols + q ;
    mat var_RE = zeros<mat>(nq, nq);
    mat ui;

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

    if(det(var_RE)<0){
      for(int j =0 ; j < var_RE.n_cols; j++){
        var_RE(j, j) += add_diag_varcov;
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
      //ui = chol_var_RE * seq_i.t() ;
      //mat uii = chol_var_RE * randn< Mat<double> >(chol_var_RE.n_rows, MCnr);
      //ui = uii.t();
    }              

    bool aMC=true;
    
    if(aMC){
      double out2;
      double vrais =0;
      
      //double log_Jac_Phi=0;
      vec vraisr_surv=ones(2);
      double surv0=0;
      double vrais_no_surv = 0;
      double vrais_surv_check=1;
      double vrais_survtot =0;
      double vraisY_tot =0;
      double min_lvraisr=0;
      double max_lvraisr=-pow(10,9);
      int nullm=0;
      int nullp=0;
      
      //log_Jac_Phi = sum(log(YiwoNA(vectorise(YtildPrimi))));
      
      for(int nr=0; nr < MCnr; nr++){
        double lvraisr=0;
        int k_t=0;
        vec ui_r = ui.row(nr).t();
        
        // ###inverse
        if(1>2){
          vec ui_r_temp = ui_r;
          ui_r(0) = ui_r_temp(1);
          ui_r(1) = ui_r_temp(0);
          ui_r(2) = ui_r_temp(3);
          ui_r(3) = ui_r_temp(2);
        }

        //vec Lambda_nr = matNui_ui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi, true);
        
        // if(nD>1){
        //   ui_r = zeros<vec>(nb_RE);
        //   vec ui_rall = ui.row(nr).t();
        //   int d_act = K2_lambda(k);// latent process corresponding to marker k
        //   
        //   int t_re=0;
        //   for(int i=0; i < nb_RE; i++){
        //     ui_r(i) = ui_rall(t_re + d_act);
        //     t_re += nD;
        //   }
        // }else{
        //   ui_r = ui.row(nr).t();
        // }
        //ui_r.fill(0);
        vec Lambda_nr = matNui_ui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, ui_r, zi, true);

        int kk = 0;
        for (int k = 0 ; k < K; k++){
          
          vec ParaTransformYk = ParamTransformY(span(kk, (kk+df[k]-1)));
          vec tau_ik = matTik(Ytildi.col(k), tau_i);
          vec Ytildik = YiwoNA(vectorise(Ytildi.col(k)));
          
          int nik = YiwoNA(Ytildik).size();
          nik = Ytildik.size();
          vec Lambda_nrk(nik);
          mat Sig_k = eye(nik,nik)*Sig(k,k);
          
          for (int j = 0 ; j < nik; j++){
            int jj=0;
            while(tau_i(jj)<tau_ik(j)){
              jj++;
            }
            if(tau_i(jj)==tau_ik(j)){
              Lambda_nrk(j)=Lambda_nr(jj*nD + K2_lambda[k]);
            }else{
              cout  << " problem definition tau_ik"<<endl;
              cout << " k "<<k<<" j "<<j<< " tau_i "<<tau_i.t()
                   <<" jj "<<jj<< " tau_ik "<<tau_ik.t()<<endl;
            }
          }
          
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
              vec Ytildi_nu_i_uik = Ytildik-Lambda_nrk;
              out2 = -0.5*(nik*log(2*M_PI) + log(det(Sig_k)) + as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(Sig_k)*Ytildi_nu_i_uik));

              if(printa==1&& nr<0){
                cout <<" k :"<<k<< " nr "<<nr<< endl<<endl<<"lvraisr" << lvraisr<<" out2 "<<out2<< " nik "<<nik <<  " log(det(Sig_k)) "<<log(det(Sig_k))<< " scalar "<< as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(Sig_k)*Ytildi_nu_i_uik) << " log_Jac_Phi "<<log_Jac_Phi<<endl;
                cout <<  " ui_r " <<ui_r.t();
                cout <<  " Ytildik " <<Ytildik.t();
                cout <<  " Lambda_nrk " <<Lambda_nrk.t();
                cout <<  " Ytildi_nu_i_uik " <<Ytildi_nu_i_uik.t()<<endl<<endl;
                //   << " matVY_i "<<matVY_i;
              }
              // if(nr==0){
              //   cout << " sum(k_i) "<< sum(k_i) 
              //        <<  " log(2*M_PI) "<< log(2*M_PI)
              //        <<" log(det(sigMSM)) "<< log(det(sigMSM)) 
              //        << " produit "<< as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(sigMSM)*Ytildi_nu_i_uik)<<endl;
              //   cout << " Ytildi_nu_i_uik "<< Ytildi_nu_i_uik.t()<<endl
              //        << " inv_sympd(sigMSM) "<< inv_sympd(sigMSM) <<endl<<endl;
              // }

              lvraisr += out2;
              
            }else if(if_link(k)==2){ // thresholds
              double phi1;
              double phi2;
              double vraisk=1;
              for (int j = 0 ; j < nik; j++){//NA!!
                
                if(k_i(j)>0){ // change: if k_i(j,k)==1
                  
                  // double binf = ParaTransformYk(0);
                  // double bsup = ParaTransformYk(0);
                  // 
                  // for (int m = 0 ; m < (zitr(2*k_t+1)-zitr(2*k_t)+1); m++){
                  //   if(Ytildik(j)<(zitr(2*k_t) + m)){
                  //      binf = binf +  pow(ParaTransformYk_b(m), 2);
                  //   }
                  // }
                  
                  int mm=0; //index in ParaTransformYk
                  double inf=-10;
                  double sup=-10;
                  
                  if(Ytildik(j)==(zitr(2*k_t))){
                    double value = (ParaTransformYk(0)-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5));//(ParaTransformYk(mm)-Lambda_nrk(j))/abs(Sig(k,k));
                    phi1 = normalCDF(value);
                    phi2 = 0;

                    if(printa)
                      cout << nr  << " j "<< j<< " m "<<zitr(2*k_t) << " value "<<value << " ui_r "<<ui_r.t()
                           << "PT "<<ParaTransformYk(0)<< " Lambda_nrk(j) "<<Lambda_nrk(j)<< " phi1 "<<phi1<<" phi "<<phi1-phi2<<endl;
                  }else{
                    
                    inf=ParaTransformYk(0);
                    sup=ParaTransformYk(0);
                    for (int m = 0 ; m < (df[k]-1); m++){
                      sup += pow(ParaTransformYk(m+1),2);
                      if(Ytildik(j)==(zitr(2*k_t) + m+1)){
                        double value = (sup-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5));//(ParaTransformYk(mm)-Lambda_nrk(j))/abs(Sig(k,k));
                        phi1 = normalCDF(value);
                        value = (inf-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5));//(ParaTransformYk(mm-1)-Lambda_nrk(j))/abs(Sig(k,k));
                        phi2 = normalCDF(value);
                        
                        if(printa)
                          cout << nr << " j "<< j<< " m "<<zitr(2*k_t) + m+1 << " inf "<< inf << " sup "<<sup<< " ui_r "<<ui_r.t()<<" value1 "<< (sup-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5))
                               << "value2 "<<(inf-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5))
                               << " phi1 "<<phi1<< " phi2 "<<phi2<<" phi "<<phi1-phi2<<endl;
                      }
                      inf = sup;
                    }
                    if(Ytildik(j)==(zitr(2*k_t+1))){
                      phi1 = 1;
                      double value = (inf-Lambda_nrk(j))/abs(pow(Sig(k,k),0.5));//(ParaTransformYk(mm-1)-Lambda_nrk(j))/abs(Sig(k,k));
                      phi2 = normalCDF(value);
                    
                      if(printa)
                        cout << nr << " j "<< j<< " m "<<zitr(2*k_t+1)<< " inf "<< inf 
                             << "value "<<value<< " ui_r "<<ui_r.t()
                             << " phi2 "<<phi2<< " phi "<<phi1-phi2<<endl;
                    }
                  }
                  
                  if(phi1< phi2){
                    cout << " phi1< phi2: j "<< j //<< " PT "<< ParaTransformYk.t()<<endl
                         << " exp " <<(exp(phi1) - exp(phi2))
                         << " phi1 "<< phi1
                         << " phi2 "<< phi2 
                         << " j "<< j << " k "<< k
                        << " Ytildi "<<Ytildik(j)<< " zitr "<<zitr.t(); 
                    cout << " Ytildi "<<Ytildi.t();
                  }
                  
                  vraisk *= (phi1-phi2);

                  if(phi1==phi2 && printa==1)
                    cout << nr << " j "<<j<<" Ytildik(j) "<<Ytildik(j)<< " phi "<<(phi1-phi2)<< " vraisk " <<vraisk<<endl;
                  
                  // if(phi1==phi2 || isinf(lvraisk) )
                  //   cout << " nr "<< nr <<" k " << k <<"j"<<j<<" Ytildik(j) "<<Ytildik(j)<< " k_i "<<k_i.t()
                  //        << " zitr(2*k_t) "<< zitr(2*k_t) << " zitr(2*k_t+1) "<< zitr(2*k_t+1) <<" inf "<<inf << " sup "<<sup
                  //        << " Lambda_nrk(j) "<< Lambda_nrk(j)<< " phi1 "<<phi1 << " phi2 "<<phi2<< " log(phi1-phi2) "<< log(phi1-phi2) <<" lvraisk "<<lvraisk << " lvraisr " << lvraisr<< " lvrais " << lvrais<< " vrais " << vrais<<endl;
                  
                }// if Yj observed
              } // for j
              k_t ++;
              lvraisr += log(vraisk);
            }// if threshold
          } // if QMC
          //if(nr < 10){
          //  cout << " out2 "<< out2 <<endl;
          //        << " pi " << sum(k_i)*log(2*M_PI) 
          //        << " logdet " << log(det(sigMSM)) 
          //        << " exp "<<  as_scalar(Ytildi_nu_i_uik.t()*inv_sympd(sigMSM)*Ytildi_nu_i_uik)
          //        << " ui_r "<< ui_r.t()
          //       << " vectorise(Ytildi) "<< vectorise(Ytildi).t()
          //        << " Lambda_nrk.t() "<< Lambda_nrk.t()
          //       << "Ytildi_nu_i_uik.t()  "<< Ytildi_nu_i_uik.t()<<endl;
          //}
          //double out2 = f_marker(Lambda_nrk, nD, matrixP, tau, tau_i, DeltaT, Ytildi, YtildPrimi, x0i, alpha_mu0, xi, paraSig, alpha_mu, G_mat_A_0_to_tau_i, ParaTransformY, if_link, zitr, ide, paras_k, K2_lambda_t, K2_lambda);
          kk += df[k];
        }//k
        
        if(survival){
          vraisr_surv = f_survival_ui(ui_r, t_0i, t_i, delta_i, xti1, xti2, xti1_intY, xti2_intY, param_surv, param_surv_intY, param_basehaz, basehaz, knots_surv, assoc, truncation,
                                      nD, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i, zi, nE);
        }

        vrais_survtot += vraisr_surv(0);
        vraisY_tot += exp(lvraisr);
        //vrais += exp(lvraisr)*vraisr_surv;

        vrais += exp(lvraisr+log(vraisr_surv(0))-expotrick);
        vrais_no_surv += exp(lvraisr);
        //cout << " vrais "<<vrais<< " lvrais "<<exp(lvraisr)<< " lvrais_surv0: "<<log(vraisr_surv(0))<<endl;
        //vrais += vraisr_surv(0);//CHANGER
        //vrais += exp(lvraisr);
        surv0 += vraisr_surv(1);
        
        if(min_lvraisr>=(lvraisr-expotrick))
          min_lvraisr = lvraisr-expotrick;
        if(max_lvraisr<=(lvraisr-expotrick))
          max_lvraisr=lvraisr-expotrick;
        if((lvraisr-expotrick) < -7*pow(10,2))
          nullm++;
        if((lvraisr-expotrick) > 7*pow(10,2))
          nullp++;
      }//nr

      if(check==1){
        double minY = pow(10,10);
        for(int j = 0 ; j < Ytildi.size(); j++){
          if(Ytildi(j)<minY)
            minY = Ytildi(j);
        }
        //if(abs(loglik_i- log(vraisY_tot/MCnr) )>1)
        //cout << " diffT "<<log(vrais_surv_check)-log(vrais_survtot/MCnr) <<" vrais_t_check "<<vrais_surv_check
        //<< " vrais_survtot/MCnr "<<vrais_survtot/MCnr<< // << " vraisTi "<< vrais_survtot/MCnr << endl;
        
        //cout << " loglik_i "<<loglik_i<< " log(vraisY_tot/MCnr) "<< log(vraisY_tot/MCnr)<< " diffT "<<log(vrais_surv_check)-log(vrais_survtot/MCnr) <<" pb_QMC" <<pb_QMC<< " minY "<< minY << " vrais / MCnr "<<endl;// << " vraisTi "<< vrais_survtot/MCnr << endl;
        //cout << " max_lvraisr "<< max_lvraisr << " min_lvraisr "<< min_lvraisr << " pb_QMC" <<pb_QMC<< " log_Jac_Phi "<< log_Jac_Phi <<endl;
      }
      if(1>2){//verification survival likelihood if all regression parameters  =0 and baseline = Weibull
        double s1=exp(-pow(t_i/param_basehaz(1),param_basehaz(0))-pow(t_i/param_basehaz(3),param_basehaz(2)));
        double lambdat=1;
        if(delta_i==1)
          lambdat = param_basehaz(0)/param_basehaz(1)*pow(t_i/param_basehaz(1),param_basehaz(0)-1);
        if(delta_i==2)
          lambdat = param_basehaz(2)/param_basehaz(3)*pow(t_i/param_basehaz(3),param_basehaz(2)-1);
        double s0=exp(-pow(t_0i/param_basehaz(1),param_basehaz(0))-pow(t_0i/param_basehaz(3),param_basehaz(2)));
        vrais_surv_check=s1/s0*lambdat;
        cout<< " vrais_surv_check "<<vrais_surv_check<< " s1 "<<s1 << " s0 "<<s0 << " lambdat "<< lambdat << " t_i "<<t_i<<endl
            << " log(vrais_survtot/MCnr) "<<log(vrais_survtot/MCnr)
            << " param_basehaz "<<param_basehaz.t() ;
      }
      surv0 /= MCnr;
      //vrais /= MCnr;
      vrais_no_surv /= MCnr;
      lvrais = -log(MCnr) + expotrick + log(vrais) + log_Jac_Phi - log(surv0); 
    
       if(printa==1 && check==1){
        cout << " diffY "<<loglik_i- lvrais<< " loglik_i "<< loglik_i << " lvrais "<<lvrais<< " log(surv0) "<<log(surv0) <<endl;
        //<< " MCnr "<<MCnr<< " minY "<< minY << " vrais / MCnr "<<vrais / MCnr;
        //cout << " loglik_i "<< loglik_i-log_Jac_Phi<< " loglik_i2 "<<loglik_i2<<" log(vraisY_tot/MCnr) "<< log(vraisY_tot/MCnr)<< " log_Jac_Phi "<<log_Jac_Phi << " vrais "<<vrais<<endl<<endl<<endl;
      }
    }
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
//' @param cholesky logical indicating if the variance covariance matrix is parameterized using the cholesky (TRUE, by default) or the correlation (FALSE)
//' @param q a vector of number of random effects on each change latent process over time
//' @param if_link indicates if non linear link is used to transform an outcome
//' @param tau a vector of integers indicating times (including maximum time)
//' @param tau_is a vector of integers indicating times for individuals
//' @param modA_mat model.matrix for elements of the transistion matrix
//' @param DeltaT double that indicates the discretization step  
//' @param data_surv matrix of Tentry, Event, StatusEvent and covariates for survival models
//' @param data_surv_intY matrix of covariates in interaction with Y in survival models
//' @param nYsurv number of covariates in interaction with Y in survival models for all nEvents
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
              arma::vec& q0, bool cholesky,
              arma::mat& data_surv, arma::mat& data_surv_intY, arma::vec& nYsurv,
              int basehaz, arma::vec& knots_surv, arma::vec& np_surv, bool survival, int assoc, bool truncation,
              int nE, arma::mat& Xsurv1, arma::mat& Xsurv2,
              arma::vec& if_link, arma::vec& zitr, arma::vec& ide, 
              arma::vec& tau, arma::vec& tau_is, arma::mat& modA_mat, double DeltaT){

  double loglik = 0.e0;
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
  
  int n_assoc = 1;
  if(survival){
    if(assoc == 2 || assoc == 5)
      n_assoc ++;
  }
  
  colvec param_surv_intY = zeros<vec>(sum(nYsurv)*nD*n_assoc);
  colvec param_basehaz;

  if(survival){
    
    if(basehaz==0){
      nq_s = 2;
    }else{
      nq_s = knots_surv.size();
    }

    vec np_surv0 = np_surv - nYsurv*nD*n_assoc;
    param_surv = zeros<vec>(sum(np_surv0));
    param_basehaz = zeros<vec>(nq_s*nE);

    for(int k=0; k<nE;k++){
      param_basehaz(span(k*nq_s, (nq_s-1)*(1-k) + (nq_s*nE-1)*k)) = paras(span(ipara,ipara + nq_s -1));
      ipara += nq_s;
      param_surv(span(k*np_surv0(0), (np_surv0(0)-1)*(1-k) + (sum(np_surv0)-1)*k)) = paras(span(ipara,ipara + np_surv0(k)-1));

      if(sum(nYsurv)>0)
        param_surv_intY(span(k*(np_surv(0)-np_surv0(0)), ((np_surv(0)-np_surv0(0))-1)*(1-k) + (sum(np_surv-np_surv0)-1)*k)) = paras(span(ipara+np_surv0(k),ipara + np_surv(k)-1));

      ipara += np_surv(k);
    }
  }

  int nb_RE = sum(sum(q0)+sum(q));
  Mat<double> matD=zeros(nb_RE,nb_RE);
  // alpha_D contains initial parameters (corr)
  mat prmea;

  if(cholesky==false){
    int ii=0;
    prmea = zeros(nb_RE, nb_RE); //triangular inf
    for(int i=0; i<nb_RE;i++){
      for(int j=i; j<nb_RE;j++){
        prmea(j,i)=alpha_D(ii);
        ii++;
      }
    }
    
    for(int i=0; i<nb_RE;i++){
      for(int j=0; j<=i;j++){
        if(i==j){
          matD(i,j)=prmea(i,j)*prmea(i,j);
        }else{
          matD(i,j)=(exp(prmea(i,j))-1)/(exp(prmea(i,j))+1);
          matD(j,i)=matD(i,j);
        }
      }
    }
    
    for(int i=0; i<nb_RE;i++){
      for(int j=0; j<nb_RE;j++){
        if(i != j ){
          matD(i,j)= matD(i,j)*pow(matD(i,i),0.5)*pow(matD(j,j),0.5);
        }
      }
    }
    
    // prmea = zeros(nb_RE, nb_RE);
    // int ii=0;
    // for(int i=0; i<nb_RE;i++){
    //   for(int j=i; j<nb_RE;j++){
    //     prmea(j,i)=alpha_D(ii);
    //     ii++;
    //   }
    // }
    // cout << " prmea0 "<<prmea;
    // 
    // mat DI=zeros(nb_RE, nb_RE);
    // DI.diag() = prmea.diag();
    // prmea = prmea + prmea.t() - DI;
    // cout << " prmea "<<prmea;
    // 
    // colvec sea = abs(prmea.diag());
    // mat corr = (exp(prmea)-1)/(exp(prmea)+1);
    // for(int i=0; i<nb_RE;i++)
    //   corr(i,i)=1;
    // 
    // matD = corr;
    // for(int i=0; i<nb_RE;i++)
    //   matD.col(i)=matD.col(i)*sea(i);
    // for(int i=0; i<nb_RE;i++)
    //   matD.row(i)=matD.row(i)*sea(i);
    // //cout << det(matD)<<" matD "<<matD;
    // cout << " matD "<<matD;
    
  }else{
    matD = DparChol(nb_RE, alpha_D);
  }
  
  int n_cols_matD = matD.n_cols;
  Mat<double> matDw = matD(span(0,nD-1),span(0,nD-1));
  Mat<double> matDw_u;
  Mat<double> matDu;
  
  if(n_cols_matD>nD){
    matDw_u = matD(span(0,nD-1),span(nD,n_cols_matD-1));
    matDu = (matD(span(nD,n_cols_matD-1),span(nD,n_cols_matD-1)));
  }else{
    matDw_u = zeros(1,1);
    matDu = zeros(1,1);
  }

  Mat<double> matB = KmatDiag(paraB); // structured variance-covariance matrice

  int add_diag_varcov = 0; // something to add on the diagonal of var_RE  to make it positive definite for cholesky computation in MC
  
  if(max(if_link)>1 || survival){//check condition
    int nq = matDw_u.n_cols + matDw.n_cols ;
    mat ui;
    //  matDw       matDw_u
    //  matDw_u.t() matDu
    // Verify var_RE is positive definite to compute the cholesky.
    double det_var_RE = det(matD);
    // Add something on the diagonal if varcov matrix is not positive definite
    if(det_var_RE<0){
      
      cout << " var_RE \n"<<matD;
      cout << " prmea \n"<<prmea;
      cout << " alpha_D "<<alpha_D.t();
      
      int p=4;
      while(det_var_RE<0 & p > 0){
        double add = pow(10,-p);
        for(int i=0; i < nq; i++){
          matD(i,i) += add;
        }
        det_var_RE = det(matD);
        if(det_var_RE<0)
          p-=1;
      }
      cout <<" proposed varcov not positive definite: "<<pow(10,-p)<< " added on diagonal"<<endl;
      add_diag_varcov = p;
      //  matDw       matDw_u
      //  matDw_u.t() matDu
      //Add on matDw and  matDu
      // for(int i=0; i < matDw.n_cols; i++){
      //   matDw(i,i) += pow(10,-p);
      // }
      // for(int i=0; i < matDu.n_cols; i++){
      //   matDu(i,i) += pow(10,-p);
      // }
    }
  }
  
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
      YtildPrim.col(k) = ParamTransformY[kk+1]*Mod_MatrixYprim.cols(kkp, (kkp+df[k]-2));
      //cout << " paramTransf "<<ParaTransformY(span(kk, (kk+df[k]-1))).t()<<endl;
      //cout << " Yprim "<<1<<endl;
      //cout << " YtildPrim "<<0.9174<<endl;
      kk += df[k];
      kkp += df[k]-1;
    }
  }
  if(survival){
    for(int j= 0; j < param_basehaz.size() ; j++ )
      param_basehaz(j) = pow(param_basehaz(j),2);
  }
  
  
  //Computering of log-likelihood as sum of individuals contributions
  for(int n= 0; n < N ; n++){ //nsubjects
    //if(n%200==0)
    
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
    mat xti1_intY(1, nYsurv(0)); 
    mat xti2_intY(1, nYsurv(1)); 

    if(survival){
      t_0i = data_surv(n,0);
      t_i = data_surv(n,1);
      delta_i = data_surv(n, 2);
      
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
      
      if(sum(nYsurv)>0){
        xti1_intY = data_surv_intY(span(n,n), span(0,nYsurv(0)-1));
        if(nE==2){
          xti2_intY = data_surv_intY(span(n,n), span(nYsurv(0), sum(nYsurv)-1));
        }
      }
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

    
    double out1 =0;
    out1= Loglikei_GLM(K, nD, matrixP, mapping, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
                       YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
                       z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
                       z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
                       matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT, ParamTransformY, df, if_link, zitr, ide, paras_k,
                       t_0i, t_i, delta_i, xti1, xti2, xti1_intY, xti2_intY, basehaz, knots_surv, survival, param_surv, param_surv_intY, param_basehaz, assoc, truncation,
                       sequence, type_int, ind_seq_i, MCnr, n, nE, add_diag_varcov, q);      
  
    loglik += out1;
    
    // double out2 =  Loglikei(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
    //                     YtildPrim(span(p,(p+m_is(n)-1)), span(0,(K-1))), x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))),
    //                     z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))), x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))),
    //                     z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),alpha_mu0, alpha_mu, matDw, matDw_u, matDu,
    //                     matB, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau,  DeltaT);
    // loglik2 += out2;
    //cout << "n:"<< n<< " loglik "<<  out1<<endl;
    
    p += m_is[n];
  }
  if(loglik <= -pow(10,9))
    loglik=-pow(10,9);
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
