function Qn = solveNoise(Qg,Am,Pm,constNoise)
%   config = inputParser;
% 
%   config.addRequired('Qg',@isnumeric);  % Glottal flow [m^3/s]
%   config.addRequired('Am',@isnumeric);  %  Membranous glottal area [m^2]
%   config.addRequired('Pm',@isnumeric); % Membranous glottal perimeter [m]
% %   config.addOptional('PGO',0,@isnumeric);  % Sub-glottal area [m^2]
% %   config.addOptional('rho',1.14,@isnumeric);% air density [kg m^-3]
% %   config.addOptional('mu',18.36922e-6,@isnumeric); % air viscosity [Pa s]
% %   config.addOptional('Re',1600,@isnumeric); % Reinolds number treshold [-]
%   
%   config.parse(varargin{:});
% 
%   Qg = config.Results.Qg;
%   Am = config.Results.Am;
%   Pm = config.Results.Pm;
%   Ap = config.Results.PGO;
%   rho = config.Results.rho;
%   mu = config.Results.mu;
%   Re_c = config.Results.Re;

  Ap   = constNoise.PGO;
  rho  = constNoise.rho;
  mu   = constNoise.mu;
  Re_c = constNoise.Re;
 
  % Hidraulic diameter (4A/P)
  Ag = Am + Ap;
  Pg = Pm + 2*sqrt(pi*Ap);
  if (Pg > 0)
    D_h = 4*Ag/Pg;
  else
    D_h = 0;
  end

  Re = Qg*(D_h/Ag)*(rho/mu);
  
  Qn = (Re > Re_c)*normrnd(0,max(0,(Re^2 - Re_c^2))*1e-14); % Titze 2006 - Myoelastic theory of phonation P.263 [m^3/s]
end
