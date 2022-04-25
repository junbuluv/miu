var c n k m i theta pi;
varexo epsit epsim;

parameters beta alpha delta rho_t sigma_theta rho_m sigma_m pi_ss psi zeta A;
@#include "params.txt"


model;
    [name='Intertemporal Euler Equation']
    c^(-1) = beta * c(+1)^(-1) * (1- delta + alpha * exp(theta)* k(-1)^(alpha-1) * n^(1-alpha));

    [name='Intratemporal Euler Equation']
    (1-alpha)*exp(theta)*k(-1)^(alpha)*n^(-alpha) = (A*c) / (1-n);
    
    [name='Money Demand']
    m = ((i/(1+i)) * (c^(-1)/psi))^(-1/zeta)  ;  
    
    [name='Euler Equation - Bond']
    
    c^(-1)/(1+i) = beta*(c(+1)^(-1) * exp(-pi(+1)));

    [name='Budget Constraint']
    c + k = (1-delta)*k(-1) + exp(theta)*k(-1)^(alpha)*n^(1-alpha); 

    [name='Money Supply']
    log(m) - log(m(-1)) = (1-rho_m)*pi_ss - pi + rho_m * log(m(-1)) + rho_m*pi(-1) + epsim ;
    
    [name='Production shock']
    theta = rho_t*theta(-1)+epsit;
    
  
end;

initval;
@#include "steady_state.txt"
end;

steady;
check;

shocks;
var epsit; stderr sigma_theta;
var epsim; stderr sigma_m;
end;

stoch_simul(order=1, periods = 500, irf = 100);
