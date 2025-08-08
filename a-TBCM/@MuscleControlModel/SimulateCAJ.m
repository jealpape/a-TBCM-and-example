%%
% SimulateCAJ: Function for modeling the dynamics of CA joint as a function
% of the coupled activation of all the five intrinsic laryngeal (PCA, IA,
% LCA, TA and CT) muscles and the passive (vocal ligament and VF mucosa).
% fibrous tissue.
%
% Structure: SimulateCAJ(MCObj)
% where
%
% BCMObj: is an object from MuscleControlModel (handle) class.
%
% References:
% [1] I. R. Titze and E. J. Hunter, “A two-dimensional biomechanical model
%     of vocal fold posturing,” J. Acoust. Soc. Am., vol. 121, no. 4, pp.
%     2254–2260, Apr. 2007. 
%
% [2] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation, 1st
%     editio. National Center for Voice and Speech, 2006.
%
% Coded by Gabriel Alzamendi, February 2020.
function SimulateCAJ(MCObj)
  % CAJ variables
  xi_a = MCObj.xData_CAJ(1);
  dxi_a = MCObj.xData_CAJ(2);
  psi_a =  MCObj.xData_CAJ(4);
  dpsi_a =  MCObj.xData_CAJ(5);
  theta_a = MCObj.xData_CAJ(7);
  dtheta_a = MCObj.xData_CAJ(8);
  
  % Stiffness constants
  [k_ax, k_ay] = MCObj.CalcCAJoinTranslationalStiffness(xi_a, psi_a);
  kappa_a = MCObj.CalcCAJoinRotationalStiffness(theta_a);
  
  % Reaction moment due to collision between left/right arytenoids
%   tau_col = MCObj.TorqueCol;
  [Fx_col,Fy_col,tau_col] = MCObj.ColReaction;
  
  % *** x-acceleration of CAJ center Eq. (3.22) in [2] ***
  ddxi_a = 0;
  for cont_m = MCObj.LarMusIndex%(MCObj.LarMusIndex~=4)
    ddxi_a = ddxi_a + MCObj.LarMuscObj(cont_m).alpha*MCObj.LarMuscObj(cont_m).getMuscForce;  
  end
  ddxi_a = ddxi_a - k_ax*(xi_a+80/k_ax*MCObj.t_ax*dxi_a) + Fx_col;
  ddxi_a = ddxi_a/MCObj.Ma;
%   ddxi_a = (MCObj.xData_CAJ(3) + ddxi_a)/2;
  % *** End Eq. (3.22) ***
  
  % *** y-acceleration of CAJ center Eq. (3.23) in [2] ***
  ddpsi_a = 0;
  for cont_m = MCObj.LarMusIndex%(MCObj.LarMusIndex~=4)
    ddpsi_a = ddpsi_a + MCObj.LarMuscObj(cont_m).beta*MCObj.LarMuscObj(cont_m).getMuscForce;  
  end
  ddpsi_a = ddpsi_a - k_ay*(psi_a+250/k_ay*MCObj.t_ay*dpsi_a) + Fy_col;
  ddpsi_a = ddpsi_a/MCObj.Ma;
%   ddpsi_a = (MCObj.xData_CAJ(6) + ddpsi_a)/2;  % To smooth dynamic response
  % *** End Eq. (3.23) ***
  
  % *** Rotational acceleration of CAJ center Eq. (3.24) in [2] ***
  ddtheta_a = 0;
  for cont_m = MCObj.LarMusIndex(MCObj.LarMusIndex~=4)
    ddtheta_a = ddtheta_a + MCObj.LarMuscObj(cont_m).gamma*MCObj.LarMuscObj(cont_m).getMuscForce;  
  end
  ddtheta_a = ddtheta_a - kappa_a*(theta_a+0.01/kappa_a*MCObj.t_ar*dtheta_a) + tau_col;
  ddtheta_a = ddtheta_a/MCObj.Ia;
%   ddtheta_a = (MCObj.xData_CAJ(9) + ddtheta_a)/2;
  % *** End Eq. (3.24) ***
  
  %% ODE Solve: x-component of CAJ center
  dT = MCObj.Ts;
  MCObj.xData_CAJ(1) = MCObj.xData_CAJ(1) + dT*MCObj.xData_CAJ(2) + ...
                        0.5*dT^2*ddxi_a;
  MCObj.xData_CAJ(2) = MCObj.xData_CAJ(2) + dT*ddxi_a;
  MCObj.xData_CAJ(3) = ddxi_a;
  %% ODE Solve: y-component of CAJ center
  MCObj.xData_CAJ(4) = MCObj.xData_CAJ(4) + dT*MCObj.xData_CAJ(5) + ...
                        0.5*dT^2*ddpsi_a;
  MCObj.xData_CAJ(5) = MCObj.xData_CAJ(5) + dT*ddpsi_a;
  MCObj.xData_CAJ(6) = ddpsi_a;
  %% ODE Solve: Rotational component of CAJ center
  MCObj.xData_CAJ(7) = MCObj.xData_CAJ(7) + dT*MCObj.xData_CAJ(8) + ...
                        0.5*dT^2*ddtheta_a;
  MCObj.xData_CAJ(8) = MCObj.xData_CAJ(8) + dT*ddtheta_a;
  MCObj.xData_CAJ(9) = ddtheta_a;
end