%%
% Simulate: Function for produce the one-step simulation of the body 
% cover model of the vocal folds.
%
% Structure: BCMSimulate(TBCMobj, Ps, Pe)
%            BCMSimulate(TBCMobj, Ps, Pe, Ae)
% where
%
% TBCMobj: is an object from BodyCoverModel (handle) class,
% Ps: is the subglottal pressure in the trachea in Pascals,
% Pe: is the supraglottal pressure in the epilarynx in Pascals,
% Ae: Epilarynx tube area in m^2.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function Simulate(TBCMobj, Ps, Pe, varargin)
  if ~TBCMobj.SimParamOK % Simulation parameters are missing
    error('The simulation parameter are missing! See ''setSimulationParameter'' function.')
  end
  
  % Sampling time
  Ts = TBCMobj.Ts;
  
  % Past dynamic state variables
  xu = TBCMobj.xData(1); % [m] upper mass displacement
  xl = TBCMobj.xData(2); % [m] lower mass displacement
  xb = TBCMobj.xData(3); % [m] body mass displacement
  vu = TBCMobj.xData(4); % [m/s] upper mass velocity
  vl = TBCMobj.xData(5); % [m/s] lower mass velocity
  vb = TBCMobj.xData(6); % [m/s] body mass velocity
  
  % Computation of internal and external forces
  [Fku,Fkl,Fkb,Fkc] = TBCMobj.ElasticForces; % Elastic forces  
  [Fdu,Fdl,Fdb] = TBCMobj.DampingForces; % Damping forces 
  [Fu_col,Fl_col] = TBCMobj.CollisionForces; % Collision forces
  if (nargin == 3)
    [Feu, Fel] = TBCMobj.AeroPressure2DrivingForces(Ps, Pe); % Externar driving forces
  elseif (nargin==4)
    if isnumeric(varargin{1})&&(varargin{1} > 0)
      Ae = varargin{1}; % Epilarynx tube area
      [Feu, Fel] = TBCMobj.AeroPressure2DrivingForces(Ps, Pe, Ae); % Externar driving forces
    else
      error('Epilarynx tube area must be a positive number!')
    end
  else
    error('Incorrect number of function input variables! See the correct function estructures in the m file.')  
  end
  
  % Total Force Equations
  Fu = Fku + Fdu - Fkc + Feu + Fu_col;         % Total Force in the Upper Mass [N] 
  Fl = Fkl + Fdl + Fkc + Fel + Fl_col;         % Total Force in the Lower Mass [N]
  Fb = Fkb + Fdb - (Fku + Fdu + Fkl + Fdl);  % Total Force in the Body Mass [N] (Typo in the original paper)
  
  % Acceleration
  au = Fu/TBCMobj.mu;
  al = Fl/TBCMobj.ml;
  ab = Fb/TBCMobj.mb;

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
  TBCMobj.xData = xData;
  TBCMobj.n_IterCont = TBCMobj.n_IterCont + 1;
end