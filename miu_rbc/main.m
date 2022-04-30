clear all
clc
%% money in utility in Dynare

%% set parameters


nss= 1/3;
alpha = 0.33;
beta = 0.99;
delta = 0.025;
rho_t = 0.95;
sigma_theta = 0.09;
rho_m = 0.5;
sigma_m = 0.01;
pi_ss = 0.02;
psi = 1;
zeta = 1;


%% find steady state values and calibration
options = optimset('Display','iter','MaxFunEvals',1000000,'TolFun',1e-8,'MaxIter',10000);
% steady state k 
fun_k = @(k) beta * (1-delta + alpha*k^(alpha-1) * nss^(1-alpha)) -1;
kss = fsolve(fun_k,0.1,options);
% steady state c

css = (1-delta + alpha*kss^(alpha-1)*nss^(1-alpha))*kss + (1-alpha)*kss^(alpha)*nss^(1-alpha) - kss;

% Fisher equation 
iss= (1+pi_ss)/beta - 1;

% steady state m
mss=((iss/(1+iss))*(css^(-1)/psi))^(-1/zeta);

% find A
A = (1-alpha)*kss^(alpha)*nss^(-alpha) * (1-nss) / css;



param = struct("alpha",alpha,"beta",beta,"delta",delta,"rho_t",rho_t,"sigma_theta",sigma_theta...
    ,"rho_m",rho_m,"sigma_m",sigma_m,"pi_ss",pi_ss,"psi",psi,"zeta",zeta,"A",A);
steady_val = struct("css",css,"kss",kss,"nss",nss,"mss",mss, "iss", iss);



steady_state_names=["c";"k";"n";"m";"i"];
steady_state_vals=[steady_val.css; steady_val.kss; steady_val.nss; steady_val.mss; steady_val.iss];
writematrix(strcat(strcat(steady_state_names,"="),...
    strcat(string(steady_state_vals),";")),"steady_state.txt");

params_names=["beta";"alpha";"delta";"rho_t";"sigma_theta";"rho_m";"sigma_m";"pi_ss";"psi";"zeta";"A"];
params_vals=[param.beta;param.alpha;param.delta;param.rho_t;param.sigma_theta;param.rho_m;param.sigma_m;param.pi_ss;param.psi;param.zeta;param.A];
writematrix(strcat(strcat(params_names,"="),...
    strcat(string(params_vals),";")),"params.txt");


dynare miu_RBC.mod




