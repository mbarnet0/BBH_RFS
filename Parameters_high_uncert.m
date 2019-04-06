%%% This file combine all parameters so far for HJB, ODE, SPE and SCC
clear all
close all
clc

filename = [pwd, '/climate_pars_high_uncert_BUMP'];

%% Uncertainty Pars
% theta = 0;
% kappa = 5000;
theta = 100;
kappa = 0.1;
if kappa>10
    RobDrift = 0;
else
    RobDrift = 1./(2.*kappa);
end

%% Emission Parameters

%%% Castruccio:
M = csvread('par_values.csv');
par_rho = M(:,1);
par_phi = M(:,2);
par_beta0 = M(:,3);
par_beta1 = M(:,4);
par_sigma_z_2 = M(:,5);

phi = mean(par_phi);
var_phi = var(par_phi);
rho = mean(par_rho);
var_rho = var(par_rho);
beta0 = 0;
par_beta1rho = par_beta1.*par_rho;
beta1rho_ave = mean(par_beta1rho);
beta1rho_min = min(par_beta1rho);
beta1rho_max = max(par_beta1rho);
var_beta1rho = var(par_beta1rho);
cov_beta1rho_rho = cov(par_beta1rho,par_rho);
cov_beta1rho_rho = cov_beta1rho_rho(1,2);
sigma_z = mean(sqrt(par_sigma_z_2));

psi_PS = 0.07;
var_psi_PS = 0.005;
par_psi_PS = psi_PS.*ones(size(par_rho));

psi_hotel = 0.0064;
var_psi_hotel = 0.0005;
par_psi_hotel = psi_hotel.*ones(size(par_rho));

%%% Matthews
%%% Use McD file for updated lambda parameters
McD = csvread('TCRE_MacDougallEtAl2017_update.csv');
McD = McD./1000;
lambda_McD = mean(McD(:,1));
var_lambda_McD = var(McD(:,1));
lambda_max_McD = max(McD(:,1));
skew_lambda_McD = skewness(McD(:,1)).*std(McD(:,1)).^3;
kur_lambda_McD = kurtosis(McD(:,1)).*std(McD(:,1)).^4;
par_lambda_McD = McD(:,1);

%%% Old version of lambda parameters
N = csvread('matthews_param_estimates_update.csv');
lambda_CCR = mean(N(:,3));
var_lambda_CCR = var(N(:,3));
lambda_max_CCR = max(N(:,3));
skew_lambda_CCR = skewness(N(:,3)).*std(N(:,3)).^3;
kur_lambda_CCR = kurtosis(N(:,3)).*std(N(:,3)).^4;
par_lambda_CCR = N(:,3);

lambda_CCR_ppm = mean(N(:,4));
var_lambda_CCR_ppm = var(N(:,4));
lambda_max_CCR_ppm = max(N(:,4));
skew_lambda_CCR_ppm = skewness(N(:,4)).*std(N(:,4)).^3;
kur_lambda_CCR_ppm = kurtosis(N(:,4)).*std(N(:,4)).^4;
par_lambda_CCR_ppm = N(:,4);

lambda_TCR = mean(N(:,7));
var_lambda_TCR = var(N(:,7));
lambda_max_TCR = max(N(:,7));
skew_lambda_TCR = skewness(N(:,7)).*std(N(:,7)).^3;
kur_lambda_TCR = kurtosis(N(:,7)).*std(N(:,7)).^4;
par_lambda_TCR = N(:,7);

lambda_TCR_ppm = mean(N(:,8));
var_lambda_TCR_ppm = var(N(:,8));
lambda_max_TCR_ppm = max(N(:,8));
skew_lambda_TCR_ppm = skewness(N(:,8)).*std(N(:,8)).^3;
kur_lambda_TCR_ppm = kurtosis(N(:,8)).*std(N(:,8)).^4;
par_lambda_TCR_ppm = N(:,8);

lambda_TCR_hotel = mean(N(:,5));
% lambda_TCR_hotel = min(N(:,5));
% lambda_TCR_hotel = max(N(:,5));
var_lambda_TCR_hotel = var(N(:,5));
lambda_max_TCR_hotel = max(N(:,5));
skew_lambda_TCR_hotel = skewness(N(:,5)).*std(N(:,5)).^3;
kur_lambda_TCR_hotel = kurtosis(N(:,5)).*std(N(:,5)).^4;
par_lambda_TCR_hotel = N(:,5);

% %%% default using TCR values for now
% lambda_PS = lambda_TCR;
% var_lambda_PS = var_lambda_TCR;
% par_lambda_PS = par_lambda_TCR;
% 
% lambda_hotel = lambda_TCR_hotel;
% var_lambda_hotel = var_lambda_TCR_hotel;
% par_lambda_hotel = par_lambda_TCR_hotel;

%% Economic/capital
% % % % delta = 0.02;
delta = 0.025;
% delta = 0.01;
% % % % alphaO = 0.05;
% alphaO = 0.05;
alphaO= 0.04;
xi_m = alphaO;
xi_m_3state = 1;
xi_o = alphaO;
xi_g = 1-alphaO;
sigma_g = 0.02;
sigma_o = 0.02;
% sigma_d = 0.005;
% sigma_d = 0.01;
sigma_d = 0.005;
sigma_m = 0.01;
sigma_n = 0.01;
A_O = 0.1;
A_G = 0.1;

alpha = alphaO;

%%% capital
%%% quadratic adj cost
delt = 0.02;
nu2 = 10;
% EW
nu1Oew = 1-nu2.*A_O;
nu1Gew = 1-nu2.*A_G;
% BWY
nu1Obwy = 1+nu2.*(delt-A_O);
nu1Gbwy = 1+nu2.*(delt-A_G);
% EW
mu_oew = A_O.*nu1Oew + nu2.*0.5.*(A_O.^2)-delt;
mu_gew = A_G.*nu1Gew + nu2.*0.5.*(A_G.^2)-delt;
% BWY
mu_obwy = A_O.*nu1Obwy + nu2.*0.5.*(A_O.^2-delt.^2)-delt;
mu_gbwy = A_G.*nu1Gbwy + nu2.*0.5.*(A_G.^2-delt.^2)-delt;

nu1 = nu1Obwy;
mu_o = mu_obwy;
mu_g = mu_gbwy;

%%% Log adj cost
Gamma = 0.2;
Theta = 0.2;
Alpha = -0.045;
% % % Alpha = -0.06;
% Alpha = -0.07133;

%%% for Expl model
% Gamma_r = 2;
Gamma_r = 6.5;
% Gamma_r = 7.0;
% Gamma_r = 4;
% Gamma_r = 6.0;
Theta_r = 0.5;

%%% Oil growth and one capital setting
sigma_tau = 0.01;
sigma_r = 0.03;
mu_k = mu_o;
sigma_k = sigma_o;

%% Damage

%%%% damage parameters differ:
%%%% depending on if damage comes from capital or consumption
gamma_C = 0.005;

gamma_K = 0.0003;

%%%% matthews: damage from consumption and capital
%%%% uncertainty no longer matters
xi_k = 1-alpha;
xi_tau_C = gamma_C;
xi_tau_K = xi_k.*gamma_K./delta;

%%%% castruccio: damage from consumption and capital
%%%% uncertainty matters for the parameter values
xi_d_C = -1;
xi_d_C_PU =(-delta+sqrt(delta^2-2*delta/kappa*sigma_d^2))/(sigma_d^2/kappa);

xi_d_K = -xi_k./delta;
xi_d_K_PU =(-delta+sqrt(delta^2-2*xi_k/kappa*sigma_d^2))/(sigma_d^2/kappa);




%% exogenous forcing
sigma_11 = 0.00031;
sigma_22 = -0.038;
mu_1 = -0.021;
mu_2 = -0.013;
%% parameters I cant remember when used before
lambda1 = 0;
lambda2 = 2.0;


% gamma_par = [0.001 0.005 0.01, 0.015, 0.02, 0.05, 0.1, 0.2, 0.5]';
gamma_par = [0.001 0.005 0.01, 0.015, 0.02, 0.05, 0.1]';

%% save

save(filename);

ik = (xi_k.*Gamma.*A_O-delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+xi_k.*Gamma);
(Alpha + Gamma.*log(1+ik./Theta))

% %%% This file combine all parameters so far for HJB, ODE, SPE and SCC
% clear all
% close all
% clc
% 
% filename = [pwd, '/climate_pars_high_uncert'];
% 
% %% Uncertainty Pars
% % theta = 0;
% % kappa = 5000;
% theta = 25;
% kappa = 0.1;
% if kappa>1000
%     RobDrift = 0;
% else
%     RobDrift = 1./(2.*kappa);
% end
% 
% %% Emission Parameters
% 
% %%% Castruccio:
% M = csvread('par_values.csv');
% par_rho = M(:,1);
% par_phi = M(:,2);
% par_beta0 = M(:,3);
% par_beta1 = M(:,4);
% par_sigma_z_2 = M(:,5);
% 
% phi = mean(par_phi);
% var_phi = var(par_phi);
% rho = mean(par_rho);
% var_rho = var(par_rho);
% beta0 = 0;
% par_beta1rho = par_beta1.*par_rho;
% beta1rho_ave = mean(par_beta1rho);
% beta1rho_min = min(par_beta1rho);
% beta1rho_max = max(par_beta1rho);
% var_beta1rho = var(par_beta1rho);
% cov_beta1rho_rho = cov(par_beta1rho,par_rho);
% cov_beta1rho_rho = cov_beta1rho_rho(1,2);
% sigma_z = mean(sqrt(par_sigma_z_2));
% 
% psi_PS = 0.07;
% var_psi_PS = 0.005;
% par_psi_PS = psi_PS.*ones(size(par_rho));
% 
% psi_hotel = 0.0064;
% var_psi_hotel = 0.0005;
% par_psi_hotel = psi_hotel.*ones(size(par_rho));
% 
% %%% Matthews
% %%% Use McD file for updated lambda parameters
% McD = csvread('TCRE_MacDougallEtAl2017_update.csv');
% McD = McD./1000;
% lambda_McD = mean(McD(:,1));
% var_lambda_McD = var(McD(:,1));
% lambda_max_McD = max(McD(:,1));
% skew_lambda_McD = skewness(McD(:,1)).*std(McD(:,1)).^3;
% kur_lambda_McD = kurtosis(McD(:,1)).*std(McD(:,1)).^4;
% par_lambda_McD = McD(:,1);
% 
% %%% Old version of lambda parameters
% N = csvread('matthews_param_estimates_update.csv');
% lambda_CCR = mean(N(:,3));
% var_lambda_CCR = var(N(:,3));
% lambda_max_CCR = max(N(:,3));
% skew_lambda_CCR = skewness(N(:,3)).*std(N(:,3)).^3;
% kur_lambda_CCR = kurtosis(N(:,3)).*std(N(:,3)).^4;
% par_lambda_CCR = N(:,3);
% 
% lambda_CCR_ppm = mean(N(:,4));
% var_lambda_CCR_ppm = var(N(:,4));
% lambda_max_CCR_ppm = max(N(:,4));
% skew_lambda_CCR_ppm = skewness(N(:,4)).*std(N(:,4)).^3;
% kur_lambda_CCR_ppm = kurtosis(N(:,4)).*std(N(:,4)).^4;
% par_lambda_CCR_ppm = N(:,4);
% 
% lambda_TCR = mean(N(:,7));
% var_lambda_TCR = var(N(:,7));
% lambda_max_TCR = max(N(:,7));
% skew_lambda_TCR = skewness(N(:,7)).*std(N(:,7)).^3;
% kur_lambda_TCR = kurtosis(N(:,7)).*std(N(:,7)).^4;
% par_lambda_TCR = N(:,7);
% 
% lambda_TCR_ppm = mean(N(:,8));
% var_lambda_TCR_ppm = var(N(:,8));
% lambda_max_TCR_ppm = max(N(:,8));
% skew_lambda_TCR_ppm = skewness(N(:,8)).*std(N(:,8)).^3;
% kur_lambda_TCR_ppm = kurtosis(N(:,8)).*std(N(:,8)).^4;
% par_lambda_TCR_ppm = N(:,8);
% 
% lambda_TCR_hotel = mean(N(:,5));
% % lambda_TCR_hotel = min(N(:,5));
% % lambda_TCR_hotel = max(N(:,5));
% var_lambda_TCR_hotel = var(N(:,5));
% lambda_max_TCR_hotel = max(N(:,5));
% skew_lambda_TCR_hotel = skewness(N(:,5)).*std(N(:,5)).^3;
% kur_lambda_TCR_hotel = kurtosis(N(:,5)).*std(N(:,5)).^4;
% par_lambda_TCR_hotel = N(:,5);
% 
% % %%% default using TCR values for now
% % lambda_PS = lambda_TCR;
% % var_lambda_PS = var_lambda_TCR;
% % par_lambda_PS = par_lambda_TCR;
% % 
% % lambda_hotel = lambda_TCR_hotel;
% % var_lambda_hotel = var_lambda_TCR_hotel;
% % par_lambda_hotel = par_lambda_TCR_hotel;
% 
% %% Economic/capital
% % delta = 0.03;
% delta = 0.02;
% % alphaO = 0.5;
% % alphaO = 0.1;
% alphaO= 0.05;
% xi_m = alphaO;
% xi_m_3state = 1;
% xi_o = alphaO;
% xi_g = 1-alphaO;
% sigma_g = 0.02;
% sigma_o = 0.02;
% % sigma_d = 0.005;
% % sigma_d = 0.01;
% sigma_d = 0.005;
% sigma_m = 0.01;
% sigma_n = 0.01;
% A_O = 0.1;
% A_G = 0.1;
% 
% alpha = alphaO;
% 
% %%% capital
% %%% quadratic adj cost
% delt = 0.02;
% nu2 = 10;
% % EW
% nu1Oew = 1-nu2.*A_O;
% nu1Gew = 1-nu2.*A_G;
% % BWY
% nu1Obwy = 1+nu2.*(delt-A_O);
% nu1Gbwy = 1+nu2.*(delt-A_G);
% % EW
% mu_oew = A_O.*nu1Oew + nu2.*0.5.*(A_O.^2)-delt;
% mu_gew = A_G.*nu1Gew + nu2.*0.5.*(A_G.^2)-delt;
% % BWY
% mu_obwy = A_O.*nu1Obwy + nu2.*0.5.*(A_O.^2-delt.^2)-delt;
% mu_gbwy = A_G.*nu1Gbwy + nu2.*0.5.*(A_G.^2-delt.^2)-delt;
% 
% nu1 = nu1Obwy;
% mu_o = mu_obwy;
% mu_g = mu_gbwy;
% 
% %%% Log adj cost
% Gamma = 0.2;
% Theta = 0.2;
% % Alpha = -0.03;
% % % Alpha = -0.045;
% Alpha = -0.05;
% % Alpha = -0.07133;
% 
% %%% for Expl model
% % Gamma_r = 2;
% Gamma_r = 6.5;
% % Gamma_r = 6.0;
% % Gamma_r = 5;
% % Gamma_r = 4;
% % Gamma_r = 6.0;
% Theta_r = 0.5;
% 
% %%% Oil growth and one capital setting
% sigma_tau = 0.01;
% sigma_r = 0.03;
% mu_k = mu_o;
% sigma_k = sigma_o;
% 
% %% Damage
% 
% %%%% damage parameters differ:
% %%%% depending on if damage comes from capital or consumption
% gamma_C = 0.005;
% 
% gamma_K = 0.0003;
% 
% %%%% matthews: damage from consumption and capital
% %%%% uncertainty no longer matters
% xi_k = 1-alpha;
% xi_tau_C = gamma_C;
% xi_tau_K = xi_k.*gamma_K./delta;
% 
% %%%% castruccio: damage from consumption and capital
% %%%% uncertainty matters for the parameter values
% xi_d_C = -1;
% xi_d_C_PU =(-delta+sqrt(delta^2-2*delta/kappa*sigma_d^2))/(sigma_d^2/kappa);
% 
% xi_d_K = -xi_k./delta;
% xi_d_K_PU =(-delta+sqrt(delta^2-2*xi_k/kappa*sigma_d^2))/(sigma_d^2/kappa);
% 
% 
% 
% 
% %% exogenous forcing
% sigma_11 = 0.00031;
% sigma_22 = -0.038;
% mu_1 = -0.021;
% mu_2 = -0.013;
% %% parameters I cant remember when used before
% lambda1 = 0;
% lambda2 = 2.0;
% 
% 
% % gamma_par = [0.001 0.005 0.01, 0.015, 0.02, 0.05, 0.1, 0.2, 0.5]';
% gamma_par = [0.001 0.005 0.01, 0.015, 0.02, 0.05]';
% 
% %% save
% 
% save(filename);
% 
% 
% ik = (xi_k.*Gamma.*A_O-delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+xi_k.*Gamma);
% (Alpha + Gamma.*log(1+ik./Theta))