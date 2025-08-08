function Simulate(TFObj,VFLObj,VFRObj, Ps, Pe, varargin)

%% Left VF

% Sampling time
  Ts = VFLObj.Ts;
  
  % Past dynamic state variables
  xu = VFLObj.xData(1); % [m] upper mass displacement
  xl = VFLObj.xData(2); % [m] lower mass displacement
  xb = VFLObj.xData(3); % [m] body mass displacement
  vu = VFLObj.xData(4); % [m/s] upper mass velocity
  vl = VFLObj.xData(5); % [m/s] lower mass velocity
  vb = VFLObj.xData(6); % [m/s] body mass velocity
  
  % Computation of internal and external forces
  
  [Fku,Fkl,Fkb,Fkc] = VFLObj.ElasticForces; % Elastic forces  
  [Fdu,Fdl,Fdb] = VFLObj.DampingForces; % Damping forces 
  if (nargin == 5)
    [Feu, Fel] = VFLObj.AeroPressure2DrivingForces(Ps, Pe); % Externar driving forces
  elseif (nargin==6)
    if isnumeric(varargin{1})&&(varargin{1} > 0)
      Ae = varargin{1}; % Epilarynx tube area
      [Feu, Fel] = VFLObj.AeroPressure2DrivingForces(Ps, Pe, Ae); % Externar driving forces
    else
      error('Epilarynx tube area must be a positive number!')
    end
  else
    error('Incorrect number of function input variables! See the correct function estructures in the m file.')  
  end
  
  % Collision force
  [Fu_col, Fl_col] = TFObj.CollisionForce(VFLObj,VFRObj);
  %[Fu_col,Fl_col] = VFLObj.CollisionForces; % Collision forces
  
  
  % Total Force Equations
  Fu = Fku + Fdu - Fkc + Feu + Fu_col;         % Total Force in the Upper Mass [N] 
  Fl = Fkl + Fdl + Fkc + Fel + Fl_col;         % Total Force in the Lower Mass [N]
  Fb = Fkb + Fdb - (Fku + Fdu + Fkl + Fdl);  % Total Force in the Body Mass [N] (Typo in the original paper)
  
  % Acceleration
  au = Fu/VFLObj.mu;
  al = Fl/VFLObj.ml;
  ab = Fb/VFLObj.mb;

  % computation of present dynamic state variables
  % Truncated Taylor series method is applied for discretizing and solving
  % the body cover model of the vocal folds.
  xData = zeros(9,1);
  xData(1) = xu + Ts*vu + 0.5*Ts^2*au;      % Position of upper mass [m]
  xData(2) = xl + Ts*vl + 0.5*Ts^2*al;      % Position of lower mass [m]
  xData(3) = xb + Ts*vb + 0.5*Ts^2*ab;      % Position of body mass [m]
  xData(4) = vu + Ts*au;      % Velocity of upper mass [m s^{-1}]
  xData(5) = vl + Ts*al;      % Velocity of lower mass [m s^{-1}]
  xData(6) = vb + Ts*ab;      % Velocity of body mass [m s^{-1}]
  xData(7) = au; 
  xData(8) = al; 
  xData(9) = ab; 
  
  % Actuallize the dynamic state for TBCMobj
  VFLObj.xData = xData;
  
  %% Right VF
  
    % Past dynamic state variables
  xu = VFRObj.xData(1); % [m] upper mass displacement
  xl = VFRObj.xData(2); % [m] lower mass displacement
  xb = VFRObj.xData(3); % [m] body mass displacement
  vu = VFRObj.xData(4); % [m/s] upper mass velocity
  vl = VFRObj.xData(5); % [m/s] lower mass velocity
  vb = VFRObj.xData(6); % [m/s] body mass velocity
  
  % Computation of internal and external forces
  
  [Fku,Fkl,Fkb,Fkc] = VFRObj.ElasticForces; % Elastic forces  
  [Fdu,Fdl,Fdb] = VFRObj.DampingForces; % Damping forces 
  if (nargin == 5)
    [Feu, Fel] = VFRObj.AeroPressure2DrivingForces(Ps, Pe); % Externar driving forces
  elseif (nargin==6)
    if isnumeric(varargin{1})&&(varargin{1} > 0)
      Ae = varargin{1}; % Epilarynx tube area
      [Feu, Fel] = VFRObj.AeroPressure2DrivingForces(Ps, Pe, Ae); % Externar driving forces
    else
      error('Epilarynx tube area must be a positive number!')
    end
  else
    error('Incorrect number of function input variables! See the correct function estructures in the m file.')  
  end
  
  %[Fu_col,Fl_col] = VFRObj.CollisionForces; % Collision forces
  
  
  % Total Force Equations
  Fu = Fku + Fdu - Fkc + Feu + Fu_col;         % Total Force in the Upper Mass [N] 
  Fl = Fkl + Fdl + Fkc + Fel + Fl_col;         % Total Force in the Lower Mass [N]
  Fb = Fkb + Fdb - (Fku + Fdu + Fkl + Fdl);  % Total Force in the Body Mass [N] (Typo in the original paper)
  
  % Acceleration
  au = Fu/VFRObj.mu;
  al = Fl/VFRObj.ml;
  ab = Fb/VFRObj.mb;

  % computation of present dynamic state variables
  % Truncated Taylor series method is applied for discretizing and solving
  % the body cover model of the vocal folds.
  xData = zeros(9,1);
  xData(1) = xu + Ts*vu + 0.5*Ts^2*au;      % Position of upper mass [m]
  xData(2) = xl + Ts*vl + 0.5*Ts^2*al;      % Position of lower mass [m]
  xData(3) = xb + Ts*vb + 0.5*Ts^2*ab;      % Position of body mass [m]
  xData(4) = vu + Ts*au;      % Velocity of upper mass [m s^{-1}]
  xData(5) = vl + Ts*al;      % Velocity of lower mass [m s^{-1}]
  xData(6) = vb + Ts*ab;      % Velocity of body mass [m s^{-1}]
  xData(7) = au; 
  xData(8) = al; 
  xData(9) = ab; 
  
  % Actuallize the dynamic state for TBCMobj
  VFRObj.xData = xData;
  
  VFRObj.n_IterCont = VFRObj.n_IterCont + 1;
  VFLObj.n_IterCont = VFLObj.n_IterCont + 1;

end