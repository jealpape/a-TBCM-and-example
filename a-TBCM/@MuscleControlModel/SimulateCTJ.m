%%
% simulateCTJ: Function modeling the dynamics of CT joint as a function of
% the coupled activation of TA and CT muscles, and the passive (vocal
% ligament and VF mucosa). 
%
% Structure: simulateCTJ(MCObj)
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
function SimulateCTJ(MCObj)
    
%   delta_r = -MCObj.w_CT/MCObj.h_TA*MCObj.Lg0*MCObj.getRotStrain; % Eq. (3.61)
%   delta_t = -MCObj.Lg0*MCObj.getTranStrain/MCObj.cos_phi; % Eq. (3.62)
  delta_r = MCObj.Lg0*MCObj.getRotStrain; % Eq. (3.50)
  delta_t = MCObj.Lg0*MCObj.getTranStrain; % Eq. (3.51)
  k_t = MCObj.CalcCTJoinTranslationalStiffness(delta_t, delta_r);
  k_r = MCObj.CalcCTJoinRotationalStiffness(delta_t, delta_r);
  
  % Intrinsic muscle and passive tissue forces
  F_CT = MCObj.LarMuscObj(4).getMuscForce;
  F_TA = MCObj.LarMuscObj(5).getMuscForce;
  F_Lig = MCObj.LarMuscObj(6).getMuscForce;
  F_Muc = MCObj.LarMuscObj(7).getMuscForce;
  % CTJ variables
  eps_t = MCObj.xData_CTJ(1);
  deps_t =  MCObj.xData_CTJ(2);
  eps_r = MCObj.xData_CTJ(4);
  deps_r =  MCObj.xData_CTJ(5);
  
  % *** Rotational acceleration Eq. (3.56) in [2] ***
  ddeps_r = MCObj.w_CT*F_CT ...
            - MCObj.h_TA*(F_TA+F_Lig+F_Muc ) ...
            - k_r*MCObj.Lg0/MCObj.h_TA*(eps_r+0.1/k_r*MCObj.t_r*deps_r) ;
  ddeps_r = MCObj.h_TA*ddeps_r/(MCObj.Ir*MCObj.Lg0);
%   ddeps_r = (MCObj.xData_CTJ(6) + ddeps_r)/2;
  % *** End Eq. (3.56) ***
  
  % *** Translational acceleration Eq. (3.57) in [2] ***
  ddeps_t = F_CT*MCObj.cos_phi ...
            - (F_TA+F_Lig+F_Muc) ...
            - k_t*MCObj.Lg0*(eps_t+1000/k_t*MCObj.t_t*deps_t);
  ddeps_t = ddeps_t/(MCObj.Mt*MCObj.Lg0);
%   ddeps_t = (MCObj.xData_CTJ(3) + ddeps_t)/2;
  % *** End Eq. (3.57) ***
  
  %% ODE Solve: Translational VF strain actualization
  dT = MCObj.Ts;
  MCObj.xData_CTJ(1) = MCObj.xData_CTJ(1) + dT*MCObj.xData_CTJ(2) + ...
                        0.5*dT^2*ddeps_t;
  MCObj.xData_CTJ(2) = MCObj.xData_CTJ(2) + dT*ddeps_t;
  MCObj.xData_CTJ(3) = ddeps_t;
  %% ODE Solve: Rotational VF strain actualization
  MCObj.xData_CTJ(4) = MCObj.xData_CTJ(4) + dT*MCObj.xData_CTJ(5) + ...
                        0.5*dT^2*ddeps_r;
  MCObj.xData_CTJ(5) = MCObj.xData_CTJ(5) + dT*ddeps_r;
  MCObj.xData_CTJ(6) = ddeps_r;
end