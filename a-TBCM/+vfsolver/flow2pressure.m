function [Pe_plus,Ps_minus] = flow2pressure(Qg,Ps_plus,Pe_minus,Am,constFlow)
%   config = inputParser;
%   config.KeepUnmatched = true;
%  
%   config.addRequired('Qg',@isnumeric);  % Supra-glottal area [m^2]
%   config.addRequired('Ps_plus',@isnumeric); % Incident sub-glottal pressure [Pa]
%   config.addRequired('Pe_minus',@isnumeric);  % Incident supra-glottal pressure [Pa]
%   config.addRequired('Am',@isnumeric); % Membranous glottal area [m^2]
%   config.addOptional('Ae',0,@isnumeric);  % Supra-glottal area [m^2]
%   config.addOptional('As',0,@isnumeric);  % Sub-glottal area [m^2]
%   config.addOptional('rho',1.14,@isnumeric);% air density [kg m^-3]
%   config.addOptional('c',350,@isnumeric); % sound velocity [m s^-1]
%   config.addOptional('PGO',0,@isnumeric); % Posterior glottal area [m^2]
%   config.addOptional('solver','galindo2016',@ischar); % solver to be used)
%  
%   config.parse(varargin{:});
%   
%   
%   Qg = config.Results.Qg;
%   Am = config.Results.Am;
%   Ap = config.Results.PGO;
%   Ae = config.Results.Ae;
%   As = config.Results.As;  
%   rho = config.Results.rho;
%   c = config.Results.c;  
%   Ps_plus = config.Results.Ps_plus;
%   Pe_minus = config.Results.Pe_minus;
  
  Ap  = constFlow.PGO;
  Ae  = constFlow.Ae;
  As  = constFlow.As;  
  rho = constFlow.rho;
  c   = constFlow.c;  
  
  Ag = Am + Ap;
  
  r_s = (As - Ag)/(As + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing
  r_e = (Ae - Ag)/(Ae + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing

  % Supraglottal departing wave (Lucero JASA 2015 - smoothness...
    Pe_plus = r_e*Pe_minus + Qg*(rho*c/Ae);
  % Subglottal departing wave
    Ps_minus = r_s*Ps_plus - Qg*(rho*c/As);  
end