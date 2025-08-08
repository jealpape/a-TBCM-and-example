classdef TwoFolds < handle
% Handle class for modeling vocal fold oscillations according to the
% triangular body cover model (TBCM) of the vocal folds [2]. This model is 
% an extension of  TBCM for an assymetric configuration.
properties (SetAccess = protected)
    yc_l = 0;
    yc_u = 0;
    
end


methods
     % Class constructor
    function TFObj = TwoFolds
        
    end
    
    function InitModel(TFObj)
    % Function for initializing the dynamic state
      TFObj.yc_u = 0; % Initializing CAJ data
      TFObj.yc_l = 0; % Initializing CTJ data
    end
    
    CalcuParam(TFObj,VFObj,VFRObj)

    [Fu_col,Fl_col] = CollisionForce(TFObj,VFLObj,VFRObj)

    Simulate(TFObj,VFLObj,VFRObj, Ps, Pe, varargin)    
    
end
     
% - END CLASS DEFINITION -
end