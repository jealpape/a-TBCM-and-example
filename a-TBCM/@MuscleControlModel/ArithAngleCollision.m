%%
% ArithAngleCollision: Function for computing the arithenoid angle for the
% case of vocal process medial collision (xi_02=0), given the current
% laryngeal state.
%
% Structure: theta_a_opt = ArithAngleCollision(MCObj),
% where
%
% MCObj: is an object from MuscleControlModel (handle) class,
% theta_a_opt: is the value of the arithenoid angle for vocal process collision.
%
% References:
%
% Coded by Gabriel Alzamendi, May 2020.
function theta_a_opt = ArithAngleCollision(MCObj)

    % Optimization hyperparameters
    mu = 100;
    eps = 1e-6; % Optimization threshold
    % Initialization
    theta_a_opt = MCObj.xData_CAJ(7);
    grad = 0;
    iter=0;

    while true
        [x02_val,~] = MCObj.getVocalProcessCoord('theta_a',theta_a_opt);
        
        if ((x02_val<eps)&&(x02_val>0))||(iter>100)
            break
        end
        grad = 2*(x02_val-eps)*((MCObj.x_CAJ-MCObj.xbar_02)*sin(theta_a_opt)+MCObj.y_CAJ*cos(theta_a_opt));
    
        theta_a_opt = theta_a_opt - mu*grad;
        iter = iter + 1;
    end
    
end