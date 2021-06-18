arma::vec InnerProd(arma::vec v1, arma::vec v2);
arma::mat KmatDiag(arma::vec& Kvector);
arma::mat DparChol(int q, arma::vec& qvector);
double aijt(int t, arma::vec alpha_ijl,  arma::mat modA_mat);
arma::vec vecaijt( int K, int t, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
arma::mat ConstrA(int K, int t, double DeltaT, arma::vec vec_alpha_ij, arma::mat modA_mat);
arma::mat ProdA(int K, int t2, int t1, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
arma::mat GmatA0totaui( int K, arma::vec& vec_alpha_ij, arma::vec& tau_i, double DeltaT, arma::mat modA_mat);
arma::mat GmatprodAstotau( int K, arma::vec& vec_alpha_ij, arma::vec& tau, int t_ini, double DeltaT, arma::mat modA_mat);
arma::mat tsGmatprodA0totau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, double DeltaT, arma::mat modA_mat);
arma::mat matHit(arma::vec X_i_t);
arma::vec compoYiNA(arma::mat& Yi);
arma::vec YiwoNA(arma::vec Yi);
arma::mat matNui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                 arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);
arma::vec matNui_ui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                    arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::vec& randomeffects, arma::mat& zi, bool ordered);
arma::vec YiNui(int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);
double normalCDF(double value, bool lower_tail);
double f_marker(arma::mat& Lambdai, int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& YtildPrimi, arma::mat& x0i, arma::colvec& alpha_mu0,
             arma::mat& xi, arma::vec& paraSig, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::colvec& paraEtha2, arma::vec& if_link, arma::colvec& zitr, 
             arma::mat& ide, arma::vec& paras_k, arma::vec& K2_lambda_t, arma::vec& K2_lambda);
arma::mat transformY(arma::mat& Y, arma::colvec& paraEtha0, arma::colvec& paraEtha1);
arma::vec vectorise(arma::mat& M);
arma::mat mvnorm(int seed, arma::vec m, arma::mat SD);
arma::mat VecToMat(arma::vec& y, int K, int m_i);
arma::mat f_inv_mat(arma::mat& B);
int f_mat_print( arma::mat& B);
double f_survival_ui(arma::vec& ui_r, double t_0i, double t_i, int delta_i, arma::colvec&xti1, arma::colvec&xti2, arma::vec& param_surv, arma::vec& param_basehaz, int basehaz, arma::vec& knots_surv, int assoc, bool truncation, int nD, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                     arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi, int nE);
double fct_surv_Konrod(double t_i, double t_0i, arma::colvec&xti1, arma::colvec&xti2, arma::vec& ui_r, int delta_i, arma::vec& param_basehaz, int basehaz, arma::colvec& param_surv, 
                       arma::vec& knots_surv, int assoc, bool truncation, int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, 
                       arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi, int nE, arma::vec& gamma_X);
arma::vec fct_pred_curlev_slope(arma::vec& ptGK_delta, arma::vec& ptGK, arma::colvec&xti1, arma::colvec&xti2, arma::vec& ui_r, int delta_i, arma::colvec& param_surv, int assoc, int nD, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i, arma::mat& zi, arma::vec& param_basehaz, int basehaz, arma::vec& knots_surv, arma::vec& gamma_X, int nE, bool survfunc);