%%
% TorqueCol: Function for computing the raction moment due to the
% collission between left and right arytenoid cartilages.
%
% Structure: tau_col = TorqueCol(MCObj),
% where
%
% BCMObj: is an object from MuscleControlModel (handle) class,
% tau_col: is the value of the reaction moment.
%
% References:
%
% Coded by Gabriel Alzamendi, February 2020.
function tau_col = TorqueCol(MCObj)
  if MCObj.LimitCollision
    xi_a = MCObj.xData_CAJ(1);
    theta_a = MCObj.xData_CAJ(7);
    
    cos_aux = 0.95;
    theta_col = (0e-5 - (1-cos_aux)*MCObj.x_CAJ - cos_aux*MCObj.xbar_02 - xi_a)/MCObj.y_CAJ;
    kappa_col = 10*MCObj.CalcCAJoinRotationalStiffness(theta_a);
    
    tau_col = -kappa_col*((theta_a-theta_col)+500*(theta_a-theta_col)^3)*(theta_a>theta_col);
  else
    tau_col=0;
  end
  
end