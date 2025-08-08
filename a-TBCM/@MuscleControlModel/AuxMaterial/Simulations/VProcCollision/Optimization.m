close all; clear all; clc;

Variable.xCAJ = 10e-3;
Variable.yCAJ = -10e-3;
Variable.x__02 = 4e-3;
Variable.xi_a = 0e-3;
  
% Optimization parameters
mu = 100;
eps = 1e-6; % Optimization threshold
theta_a_opt = 0.8;
grad = 0;
iter=0;

while true
    f_val = func_aux(theta_a_opt, Variable);
    fprintf('f_val=%+2.4f mm, grad=%+1.6f, theta_a=%+2.4f rad, iter %u.\n',[f_val*1e3,grad,theta_a_opt,iter])
    
    if abs(f_val)<eps
        break
    end
    grad = 2*f_val*((Variable.xCAJ-Variable.x__02)*sin(theta_a_opt)+Variable.yCAJ*cos(theta_a_opt));
    theta_a_opt = theta_a_opt - mu*grad;
    iter = iter + 1;
end
