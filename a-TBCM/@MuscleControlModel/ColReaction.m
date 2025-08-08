%%
% ColReaction: Function for computing the raction forces (Fx_col, Fy_col)
% and the reaction moment Tau_col due to the collission between vocal
% processes from the left and right arytenoid cartilages. 
%
% Structure: [Fx_col,Fy_col,Tau_col] = ColReaction(MCObj),
% where
%
% BCMObj: is an object from MuscleControlModel (handle) class,
% Fx_col: is the horizontal reaction force component,
% Fy_col: is the vertcial reaction force component,
% tau_col: is the value of the reaction moment.
%
% References:
%
% Coded by Gabriel Alzamendi, May 2020.
function [Fx_col,Fy_col,Tau_col] = ColReaction(MCObj)
  [xi_02,~] = MCObj.getVocalProcessCoord;
  if (MCObj.LimitCollision)&&(xi_02<=0)
    theta_a = MCObj.xData_CAJ(7);
    theta_col = MCObj.ArithAngleCollision;
    kappa_col = 20*MCObj.CalcCAJoinRotationalStiffness(theta_a)/MCObj.R_CAJ;
    
    Tau_col = -kappa_col*((theta_a-theta_col)+500*(theta_a-theta_col)^3); % *(theta_a>=theta_col);
    
    Fx_col = Tau_col*(-(MCObj.x_CAJ-MCObj.xbar_02)*cos(theta_a) + MCObj.y_CAJ*sin(theta_a));
    Fy_col = Tau_col*(-MCObj.y_CAJ*cos(theta_a) - (MCObj.x_CAJ-MCObj.xbar_02)*sin(theta_a));
    Tau_col = Tau_col*MCObj.R_CAJ;
  else
    Fx_col = 0;
    Fy_col = 0;
    Tau_col = 0;
  end
  
end