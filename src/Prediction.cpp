#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <vector>
#include "genericfun.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

/*
Individual prediction
*/
arma::mat predi(int K, int nD, arma::mat matrixP, int m_i, arma::vec tau, arma::vec tau_i, arma::mat Ytild_i,
                arma::mat x0i, arma::mat z0i, arma::mat xi, arma::mat zi, arma::colvec alpha_mu0,
                arma::colvec alpha_mu, arma::mat matDw, arma::mat matDw_u, arma::mat matDu,
                arma::mat Sig, arma::mat& G_mat_A_0_to_tau_i, arma::mat& G_mat_prod_A_0_to_tau, double DeltaT,
                arma::mat& GrilleY, arma:: mat& GrilleYtild, arma::vec ParaTransformY, arma::vec if_link, arma::vec df, arma::vec minY, arma::vec maxY,
                List& knots, arma::vec degree, int MCnr, double eps){

  // to call R fonction from C++ code :
  Rcpp::Environment base("package:CInLPN2");
  Rcpp::Function f = base["f_trSpline"]; // make visible R function f_trSpline

  Rcpp::Function g = base["f_mvrnorm"]; // make visible R function f_mvrnorm


  //Ytild_i := observations  Xi := latent processes
  //m_i = number of observation;
  mat pred_MYtildFull = zeros(m_i,K); // marginal predictions
  mat pred_SSYtildFull = zeros(m_i,K); // subject-specific predictions
  mat Res_pred_MYtild = zeros(m_i,K); // marginal residuals
  mat Respred_SSYtild = zeros(m_i,K); // subject-specific residuals
  //--------------------
  // loop variables
  int p_Xj=0;
  int p_Xk =0;
  int p_Yj=0;
  int p_Yk =0;
  //--------------------
  int q = zi.n_cols; // number of random effects on the slopes
  // found the max of the vector tau_i
  double maxTau = 0.0;
  for (unsigned i = 0; i < tau_i.size(); i++){
    if (tau_i[i] > maxTau)
      maxTau = tau_i[i];
  }
  mat GrdZi = zeros((maxTau+1)*nD, q);
  mat pred = zeros(m_i,(7*K)); // matrix of all predictions (marginal and subject-specific)
  // lay of pred row : 1- pred_MYtildFull  2-ytild_i, 3- pred_MYtildFull,
  // 4- Res_pred_MYtild, 5- pred_SSYtildFull 6- pred_SSYtildFull, 7- Respred_SSYtild
  vec vect = ones(nD);







  //==============================================================================================================
  // ###### Compute of marginal and conditional distribution of transformed onservations ########
  // ##### here X define the netweork of latent processes
  //==============================================================================================================

  // marginal distribution  : Expectation (Nu_Ytild_i) and variance-covariance matrix (matVYtild_iFull)

  mat mat_Nu_i = matNui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);

  pred_MYtildFull = mat_Nu_i*matrixP.t(); // marginal prediction in transformed scale


  vec Ytild_i_nu_i = YiNui(nD, matrixP, tau, tau_i,  DeltaT, Ytild_i, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i); // (Ytild_i - mat_Nu_i)

  vec compo_Yi_NA = compoYiNA(Ytild_i); // begining of  VYtild_i computering===================
  int d = compo_Yi_NA.size();
  int m = sum(compo_Yi_NA);

  // ==== Computering of GrdZi ==============================================================
  for(int t = 0; t<= maxTau; t++){
    if(t==0){
      GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = 0*zi(span(t*nD,(t+1)*nD-1), span(0,q-1));
    }
    else{
      GrdZi(span(t*nD,(t+1)*nD-1), span(0,q-1)) = DeltaT*(zi(span(t*nD,(t+1)*nD-1), span(0,q-1)) +
        G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),t*nD-1))*GrdZi(span((t-1)*nD,t*nD-1), span(0,q-1)));
    }
  }
  mat matVX_i = zeros(nD*m_i,nD*m_i); // initialisation of matrix matVX_i ?
  mat VYtild_iFull = zeros(K*m_i,K*m_i);
  mat VXYtild_iFull = zeros(nD*m_i,K*m_i);
  for( int j =0 ; j < m_i; j++){
    p_Xk = p_Xj;
    p_Yk = p_Yj;
    for( int k =j ; k < m_i; k++){
      mat prodA_0_to_t_j_1 = diagmat(vect);
      mat prodA_0_to_t_k_1 = diagmat(vect);
      if(tau_i(j)>0){
        prodA_0_to_t_j_1 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(j)-1),nD*(tau_i(j)-1)+nD-1));
      }
      if(tau_i(k)>0){
        prodA_0_to_t_k_1 = G_mat_prod_A_0_to_tau(span(0,nD-1),span(nD*(tau_i(k)-1),nD*(tau_i(k)-1)+nD-1));
      }
      matVX_i(span(p_Xj,(p_Xj+nD-1)), span((p_Xk), (p_Xk+nD-1))) +=
        pow(DeltaT, (tau_i(j)+tau_i(k)))*(prodA_0_to_t_j_1*z0i)*matDw*(z0i*prodA_0_to_t_k_1).t() +
        pow(DeltaT, tau_i(j))*(prodA_0_to_t_j_1*z0i)*matDw_u*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t() +
        pow(DeltaT, tau_i(k))*(GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDw_u.t()*(z0i*prodA_0_to_t_k_1).t() +
        (GrdZi(span(tau_i(j)*nD,(tau_i(j)+1)*nD-1), span(0,q-1)))*matDu*(GrdZi(span(tau_i(k)*nD,(tau_i(k)+1)*nD-1), span(0,q-1))).t();

      VYtild_iFull(span(p_Yj,(p_Yj+K-1)), span((p_Yk), (p_Yk+K-1)))  = matrixP*matVX_i(span(p_Xj,(p_Xj+nD-1)), span((p_Xk), (p_Xk+nD-1)))*matrixP.t();

      VXYtild_iFull(span(p_Xj,(p_Xj+nD-1)), span((p_Yk), (p_Yk+K-1))) = matVX_i(span(p_Xj,(p_Xj+nD-1)), span((p_Xk), (p_Xk+nD-1)))*matrixP.t();

      // ###### VY_i is a symetric matrix; so we fill lower triangular matri by transpose the upper triangular part #########
      if(k != j){
        matVX_i(span(p_Xk, (p_Xk+nD-1)), span(p_Xj,(p_Xj+nD-1))) = matVX_i(span(p_Xj,(p_Xj+nD-1)), span(p_Xk, (p_Xk+nD-1))).t();
        VYtild_iFull(span(p_Yk,(p_Yk+K-1)), span((p_Yj), (p_Yj+K-1))) = VYtild_iFull(span(p_Yj,(p_Yj+K-1)), span((p_Yk), (p_Yk+K-1)));
      }

      p_Xk += nD; // incrementing p_k
      p_Yk += K; // incrementing p_k
    }
    p_Xj += nD;// incrementating p_j
    p_Yj += K;// incrementating p_j

  }

  //==========================================================================================================================
  // Marginal distribution of Ytild_i : Expectation (Nu_Ytild_i)
  // variance-covariance matrix (VYtild_iFull)

  vec Nu_Ytild_i = vectorise(mat_Nu_i*matrixP.t());
  vec vc = ones(m_i);
  mat A = diagmat(vc);
  mat Sig_i = kron(A,Sig); //
  VYtild_iFull += Sig_i;

  //====================================================================================================================
  // Conditional distribution  Ytild_i : Expectation (Nu_YtildCond_i) and variance-covariance matrix Sig_i
  // In case that K = D, then P = Id_K (identity matrix


  mat VYtild_i = zeros(m,m);

  p_Yj=0; // initializing loop variable
  p_Yk =0; // initializing loop variable
  for(int j =0; j < d; j++){
    p_Yk = p_Yj;
    for(int k =j; k< d; k++){
      if(compo_Yi_NA(j)*compo_Yi_NA(k)==1){
        VYtild_i(p_Yj,p_Yk) = VYtild_iFull(j,k);
        if(p_Yj!= p_Yk){
          VYtild_i(p_Yk,p_Yj) = VYtild_i(p_Yj,p_Yk);
        }
        p_Yk++;
      }
    }
    if(p_Yj!= p_Yk){p_Yj++;}
  } // end computering  VYtild_i ========================================================


  mat VXYtild_i = zeros(nD*m_i,m); //begining of computering VXYtild_i=============
  p_Yj=0; // initialisation variable de boucle
  for(int j =0; j < nD*m_i; j++){
    p_Yk =0; //initialisation variable de boucle
    for(int k =0; k< d; k++){
      if(compo_Yi_NA(k)==1){
        VXYtild_i(p_Yj,p_Yk) = VXYtild_iFull(j,k);
        p_Yk++;
      }
    }
    p_Yj++;
  }// end of computering VXYtild_i ========================================================

  vec Xi_hat = vectorise(mat_Nu_i) + VXYtild_i*inv_sympd(VYtild_i)*Ytild_i_nu_i;
  mat mat_Xi_hat = zeros(m_i,nD);

  for(int t=0; t< m_i; t++){
    mat_Xi_hat.row(t) = Xi_hat(span(t*nD,((t+1)*nD-1))).t();
  }
  mat mat_Nu_YtildCond_i = mat_Xi_hat*matrixP.t();
  pred_SSYtildFull = mat_Nu_YtildCond_i; // subject-specific pedictions in transformed scale

  vec Nu_YtildCond_i = vectorise(mat_Nu_YtildCond_i);


  //===================================================================================
  ///* marginal and subject-specific prediction prediction in real scales

  vec predMyi = zeros(m_i*K); //  Marginal Prediction
  vec predSSyi = zeros(m_i*K); //  Subject-specific Prediction

  for(int nr=0; nr < MCnr; nr++){

    // Marginal pred
    vec ytildM_i = mvnorm((nr+1)*666, Nu_Ytild_i,VYtild_iFull);
    mat mat_ytildM_i = VecToMat(ytildM_i, K, m_i);

    // SS pred
    vec ytildSS_i = mvnorm((nr+1)*666, Nu_YtildCond_i,Sig_i);
    mat mat_ytildSS_i = VecToMat(ytildSS_i, K, m_i);

    //Computering inverse of ytild_i by the link function H_k
    mat mat_yMn1 = zeros(m_i,K); // yn1 = H^-1(mat_ytildM_i)
    mat mat_ySSn1 = zeros(m_i,K); // yn1 = H^-1(mat_ytildSS_i)


    // Computering of(H_k)^-1(ytildSS_ik) = y_ikfor linear transformations =======================
    //===================================================================================================
    int kk = 0;
    for(int k=0; k<K; k++){
      if(if_link[k] == 0){
        vec ParamTransformYk = ParaTransformY(span(kk, (kk+df[k]-1)));
        mat_yMn1.col(k) = mat_ytildM_i.col(k)*ParamTransformYk[1] + ParamTransformYk[0]*ones(m_i);
        mat_ySSn1.col(k) = mat_ytildSS_i.col(k)*ParamTransformYk[1] + ParamTransformYk[0]*ones(m_i);
      }
      kk += df[k];
    }//===========end inverse linear link function==============================================


    // Iterative computering of (H_k)^-1(ytild_ik) = y_ik for I-splines basis link function
    // using Newton-Raphson algorithm
    //===================================================================================================
    // initializing the solution from grid values
    mat mat_yMn0 = mat_yMn1;
    mat mat_ySSn0 = mat_ySSn1;

    //
    for(int k=0; k<K; k++){
      if(if_link[k] == 1){
        for(int j=0; j< m_i;j++){
          // for marginal pred
          int lM = 0;
          while((mat_ytildSS_i(j,k) > GrilleYtild(lM,k)) & (lM< (m_i-1))){
            lM += 1;
          }
          mat_yMn0(j,k) = GrilleY(lM,k);

          // for subject-specific pred
          int lSS = 0;
          while((as_scalar(mat_ytildSS_i(j,k)) > as_scalar(GrilleYtild(lSS,k))) & (lSS< (m_i-1))){
            lSS += 1;
            //printf("l and k : %d, \t %d \n", l,k);
          }
          mat_ySSn0(j,k) = GrilleY(lSS,k);
        }
      }
    }
    vec yMn0 = vectorise(mat_yMn0);
    vec ySSn0 = vectorise(mat_ySSn0); // end initializing yMn0 and ySSn0


    double eps_nM =1.0;
    double eps_nSS =1.0; //convergence criteria for computering inverse by NR

    int itr =0;
    while((eps_nM>eps) & (eps_nSS>eps) & (itr < 10000)){
      // computering of H(yn) and HPrim(yn)
      // from vectoriel form to matrix form of yn

      // mat mat_yn0 = VecToMat(yn0, K, m_i);
      kk = 0; //remise a zero du compteur
      for(int k=0; k<K; k++){
        if(if_link[k] == 1){
          colvec ParamTransformYk = exp(ParaTransformY(span(kk, (kk+df[k]-1))));
          ParamTransformYk[0] = log(ParamTransformYk[0]);
          //Call of R function that computes the transformation by I-splines basis and  derivates of yn
          vec knots_k = as<vec>(knots[k]);
          vec mat_yMn0_colk = mat_yMn0.col(k);
          vec mat_ySSn0_colk = mat_ySSn0.col(k);
          // for marginal pred
          mat res_tr = as<arma::mat>(wrap(f(Rcpp::_["y"] = mat_yMn0_colk, Rcpp::_["minY"] = as_scalar(minY[k]),
                                            Rcpp::_["maxY"] = maxY[k],  Rcpp::_["knots"] = knots_k,
                                            Rcpp::_["degree"] = as_scalar(degree[k]), Rcpp::_["paras"] = ParamTransformYk)));

          mat_yMn1.col(k) = mat_yMn0.col(k) - InnerProd(res_tr.col(1), (res_tr.col(0)-mat_ytildM_i.col(k)));

          //for subject specific pred
          res_tr = as<arma::mat>(wrap(f(Rcpp::_["y"] = mat_ySSn0_colk, Rcpp::_["minY"] = as_scalar(minY[k]),
                                        Rcpp::_["maxY"] = maxY[k],  Rcpp::_["knots"] = knots_k,
                                        Rcpp::_["degree"] = as_scalar(degree[k]), Rcpp::_["paras"] = ParamTransformYk)));
          mat_ySSn1.col(k) = mat_ySSn0.col(k) - InnerProd(res_tr.col(1), (res_tr.col(0)-mat_ytildSS_i.col(k)));
        }
        kk += df[k];
      }
      vec yMn1 = vectorise(mat_yMn1);
      vec ySSn1 = vectorise(mat_ySSn1);
      eps_nM = as_scalar((yMn1-yMn0).t()*(yMn1-yMn0));
      eps_nSS = as_scalar((ySSn1-ySSn0).t()*(ySSn1-ySSn0));
      yMn0 = yMn1;
      ySSn0 = ySSn1;
      itr +=1; // incrementing iteration
    }
    predMyi = predMyi+yMn0;
    predSSyi = predSSyi+ySSn0;
    if(itr ==10000){ printf("Warnings!!! \n Maximun of iterations (10000) reached without convergence during inverse computing with Newton Raphson algorithm \n");}
  }
  predMyi = predMyi/MCnr;
  predSSyi = predSSyi/MCnr;//============= end invertion by NR======================================================

  mat pred_MYFull = VecToMat(predMyi,K, m_i);
  mat pred_SSYFull = VecToMat(predSSyi,K, m_i);

  // //=============================================================================================
  // //marginal and subject-specific residuals in transformed scale
  Res_pred_MYtild = Ytild_i - pred_MYtildFull; // marginal residuals
  Respred_SSYtild = Ytild_i - pred_SSYtildFull; // subject-specific residuals
  for( int k = 0; k< d; k++){
    if(compo_Yi_NA[k]==0){
      Res_pred_MYtild[k] = NAN;
      Respred_SSYtild[k] = NAN;
    }
  }

  ///=================================================================================
  // formated output for the predictions
  int ii = 0; //loop variable
  for(int k=0; k<K; k++){
    pred.col(ii) = pred_MYFull.col(k);
    pred.col(ii+1) = Ytild_i.col(k);
    pred.col(ii+2) = pred_MYtildFull.col(k);
    pred.col(ii+3) = Res_pred_MYtild.col(k);

    pred.col(ii+4) = pred_SSYFull.col(k);
    pred.col(ii+5) = pred_SSYtildFull.col(k);
    pred.col(ii+6) = Respred_SSYtild.col(k);


    ii = ii + 7;
  }

  return (pred);
  // return (pred_SSYFull);
}


//==============================================================================================================================
/* ******************************************************
Predictions for overall individuals
*/
//===========================================================================================
//' Function that computes the predictions (marginal and subject-specific) for individuals
//'  
//' @param K an integer indicating the number of markers
//' @param nD an integer indicating the number of latent processes
//' @param mapping indicates which outcome measured which latent process, it is a mapping table between
//' outcomes and latents processes
//' @param paras values of model parameters
//' @param m_is vector of numbers of visit occasions for individuals
//' @param Mod_MatrixY model.matrix from markers transformation submodels
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
//' @param MCnr an integer that indicates the number of sample for MC method  
//' @param minY a vector of minima of outcomes
//' @param maxY a vector of maxima of outcomes
//' @param knots indicates position of knots used to transform outcomes
//' @param degree indicates degree of the basis of splines
//' @param epsPred convergence criteria for prediction using MC method
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat pred(int K, int nD, arma::vec& mapping, arma::vec& paras, arma::vec& m_is,
               arma::mat& Mod_MatrixY, arma::vec df, arma::mat& x, arma::mat& z, arma::vec& q, bool cholesky,
               int nb_paraD, arma::mat& x0, arma::mat& z0, arma::vec& q0, arma::vec if_link, arma::vec tau,
               arma::vec& tau_is, arma::mat& modA_mat, double DeltaT, int MCnr, arma::vec minY, arma::vec maxY,
               List& knots, arma::vec degree, double epsPred){

  // appel de fonctions externe R
  Rcpp::Environment base("package:CInLPN2");
  Rcpp::Function f = base["f_trSpline"];


  //printf("Begining of predictions n \n");
  mat pred_Y = zeros(sum(m_is),7*K);
  arma::mat matrixP = zeros(K,nD);
  for(int k = 0; k<K; k++){
    matrixP(k,(mapping(k)-1)) = 1.e0;
  }

  int m = tau.size();
  int N = m_is.size();
  int p=0; //loop variables
  int ncol_x = (int)x.n_cols; // number of parameters for the mean slope (DeltaX)
  int ncol_x0 = (int)x0.n_cols; // number of parameters for mean of processes at baseline (X0)
  int ncol_z = (int)z.n_cols; // number of parameters for randoms effects on the slope
  int ncol_z0 = (int)z0.n_cols; // number of parameters for randoms effects on baseline processes values
  int L = modA_mat.n_cols; // number of parameters for the transition matrix A

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
  vec paraSig = paras(span(ipara,ipara+K-1));
  ipara += K;
  int nbParaTransformY = Mod_MatrixY.n_cols;
  colvec ParaTransformY = paras(span(ipara,ipara+nbParaTransformY-1));
  ipara += nbParaTransformY;

  int nb_RE = sum(sum(q0)+sum(q));
  Mat<double> matD;
  // alpha_D contains initial parameters (corr)
  if(cholesky==false){
    mat prmea = zeros(nb_RE, nb_RE);
    int ii=0;
    for(int i=0; i<nb_RE;i++){
      for(int j=0; j<=i;j++){
        prmea(i,j)=alpha_D(ii);
        ii++;
      }
    }
    
    mat DI=zeros(nb_RE, nb_RE);
    DI.diag() = prmea.diag();
    prmea = prmea + prmea.t() - DI;
    
    colvec sea = abs(prmea.diag());
    mat corr = (exp(prmea)-1)/(exp(prmea)+1);
    for(int i=0; i<nb_RE;i++)
      corr(i,i)=1;
    
    matD = corr;
    for(int i=0; i<nb_RE;i++)
      matD.col(i)=matD.col(i)*sea(i);
    for(int i=0; i<nb_RE;i++)
      matD.row(i)=matD.row(i)*sea(i);
    
  }else{
    matD = DparChol(nb_RE, alpha_D);
  }

  int n_cols_matD = matD.n_cols;
  mat matDw = matD(span(0,nD-1),span(0,nD-1));
  mat matDw_u = DeltaT*matD(span(0,nD-1),span(nD,n_cols_matD-1));
  mat matDu = DeltaT*matD(span(nD,n_cols_matD-1),span(nD,n_cols_matD-1))*DeltaT;
  mat Sig = KmatDiag(paraSig); // noice


  Mat<double> Ytild = zeros(Mod_MatrixY.n_rows,K);
  Mat<double> GrilleYtild = zeros(N,K);
  Mat<double> GrilleY = zeros(N,K);
  int kk = 0;
  int kkp = 0;
  for(int k=0; k<K;k++){
    for(int l=0; l<N; l++){
      GrilleY(l,k) = minY[k] + (maxY[k]-minY[k])*(l/N);
    }

    if(if_link[k] == 1){
      colvec ParamTransformYk = exp(ParaTransformY(span(kk, (kk+df[k]-1))));
      ParamTransformYk[0] = log(ParamTransformYk[0]);
      Ytild.col(k) = Mod_MatrixY.cols(kk, (kk+df[k]-1))*ParamTransformYk;
      vec knots_k = as<vec>(knots[k]);
      vec GrilleY_colk = GrilleY.col(k);
      mat res_tr = as<arma::mat>(wrap(f(Rcpp::_["y"] = GrilleY_colk, Rcpp::_["minY"] = as_scalar(minY[k]),
                                        Rcpp::_["maxY"] = maxY[k], Rcpp::_["knots"] = knots_k,
                                        Rcpp::_["degree"] = as_scalar(degree[k]), Rcpp::_["paras"] = ParamTransformYk)));
      GrilleYtild.col(k) = res_tr.col(0);
      kk += df[k];
      kkp += df[k]-1;
    }
    else{
      colvec ParamTransformYk = ParaTransformY(span(kk, (kk+df[k]-1)));
      ParamTransformYk[0] = - ParamTransformYk[0]/ParamTransformYk[1];
      ParamTransformYk[1] = 1.e0/ParamTransformYk[1];

      Ytild.col(k) = Mod_MatrixY.cols(kk, (kk+df[k]-1))*ParamTransformYk;
      GrilleYtild.col(k) = (GrilleY.col(k)-ParamTransformYk[0])/ParamTransformYk[1];
      kk += df[k];
      kkp += df[k]-1;
    }

  }

  // predictions=============================================
  for(int n= 0; n < N; n++ ){
    //Creation of matrix G_mat_prod_A_0_to_tau that contains all products  A(j) from t_j a Tmax: t_j \in 0, Tmax
    mat G_mat_prod_A_0_to_tau = GmatprodAstotau(nD, vec_alpha_ij, tau, 0, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));
    mat G_mat_A_0_to_tau_i = GmatA0totaui(nD, vec_alpha_ij, tau, DeltaT, modA_mat(span(n*m,((n+1)*m-1)), span(0,(L-1))));

    pred_Y.rows(p,(p+m_is[n]-1)) = predi(K, nD, matrixP, m_is(n), tau, tau_is(span(p,(p+m_is(n)-1))), Ytild(span(p,(p+m_is(n)-1)), span(0,(K-1))),
                x0(span(n*nD,(n+1)*nD-1), span(0,(ncol_x0-1))), z0(span(n*nD,(n+1)*nD-1), span(0,(ncol_z0-1))),
                x(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_x-1))), z(span(n*nD*m,((n+1)*nD*m-1)), span(0,(ncol_z-1))),
                alpha_mu0, alpha_mu, matDw, matDw_u, matDu, Sig, G_mat_A_0_to_tau_i, G_mat_prod_A_0_to_tau, DeltaT,
                GrilleY, GrilleYtild, ParaTransformY, if_link, df, minY, maxY, knots, degree, MCnr, epsPred);

    p += m_is[n];
  }
  
  return(pred_Y);
}
