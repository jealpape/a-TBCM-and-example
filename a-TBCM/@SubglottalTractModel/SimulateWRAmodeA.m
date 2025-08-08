%%
% SimulateWRAmodeA: Function for simulating acoustic wave propagation throughout
% the subglottal tract representation accordign to a particular
% area function. The function implements the half-time Wave Reflection
% Analogue (WRA) solver taking into account the attenuation effects due to
% wave propagation in the tubulets sections, according to Sec. 6.4 in [1].
%
% Structure: SimulateWRAmodeA(SGTObj,Ug_n)
%            SimulateWRAmodeA(SGTObj,Ug_n,rs_n)
%            PressureWaves = SimulateWRAmodeA(...)
%
% where
%
% SGTObj: is an object from VocalTractModel (handle) class,
% Ug_n: is the air flow (volume velocity) value for the n-th instant,
% rs_n: is the subglottal reflection coefficient (=1 by default),
% PressureWaves: is vector gathering the backward (Bn) and forward (Fn)
%                acoustic wave components and the radiated acoustic
%                pressure Pout
%                PressureWaves = [B1 B2 ... BL F1 F2 ... FL FL_p Pout],
%                with FL_p an auxiliari variable.
%
% References:
% [1] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation, 1st
%     editio. National Center for Voice and Speech, 2006. 
%
% Coded by Gabriel Alzamendi, January 2020.
function varargout = SimulateWRAmodeA(SGTObj,Ug_n,varargin)
    rs_n = 1;
    if (nargin == 3)
      if isnumeric(varargin{1})&&(abs(varargin{1})<=1)
        rs_n = varargin{1};
      else
        error('Incorrect ''rs_n'' subglottal reflection coefficient. Correct value mast fullf rs_n need to be |rs_n|<=1.')
      end
    end

    % Definition of simulation parameters
    rho = SGTObj.RHO_AIR; % [kg m^-3] Density of the air
    c = SGTObj.C_AIR; % [m/s] speed of sound
    fs = SGTObj.fs;
    Delta_z = SGTObj.Delta_z; % c/(2*fs); % lenght of each tube section [m]
    N_AreaSection=SGTObj.N_AreaSection;
    
    % Attenuations factors
    A_att  = 1 - (11.2e-3./sqrt(SGTObj.AreaFunction))*Delta_z; 
    
    % Reflections coefficients
    r_end = SGTObj.r_end;
    r_coef = (SGTObj.AreaFunction(1:N_AreaSection-1) - SGTObj.AreaFunction(2:N_AreaSection)) ./ ...
             (SGTObj.AreaFunction(1:N_AreaSection-1) + SGTObj.AreaFunction(2:N_AreaSection));
    
    % Input acoustic impedance
    Z_Ug = -rho*c/SGTObj.AreaFunction(1);
    
    % Initialization
    B_p = SGTObj.xData(1:N_AreaSection);
    F_p = SGTObj.xData(N_AreaSection+1:2*N_AreaSection);
        
    %% Sound wave propagation
    % Forward pressure wave at odd junctions and half sample computation
    F_p(1) = Z_Ug*Ug_n + rs_n*A_att(1)*B_p(1);
    F_p(3:2:end-1) = (1+r_coef(2:2:end-1)).*A_att(2:2:end-1).*F_p(2:2:end-1) ...
                            - r_coef(2:2:end-1).*A_att(3:2:end).*B_p(3:2:end);

    % Backward pressure wave at even junctions and half sample computation
    B_p(2:2:end-1) = r_coef(2:2:end-1).*A_att(2:2:end-1).*F_p(2:2:end-1) ...
                        + (1-r_coef(2:2:end-1)).*A_att(3:2:end-1).*B_p(3:2:end-1);
    B_p(end) = r_end*A_att(end)*F_p(end);

    % Forward pressure wave at even junctions and integer sample computation
    F_p(2:2:end) = (1+r_coef(1:2:end)).*A_att(1:2:end-1).*F_p(1:2:end-1) ...
                            - r_coef(1:2:end).*A_att(2:2:end).*B_p(2:2:end);
    % Backward pressure wave at odd junctions and integer sample computation
    B_p(1:2:end-1) = r_coef(1:2:end).*A_att(1:2:end-1).*F_p(1:2:end-1) ...
                        + (1-r_coef(1:2:end)).*A_att(2:2:end).*B_p(2:2:end);

    %% Storing resulting acoustic waves
    PressureWaves = [B_p; F_p]; 
    SGTObj.xData = PressureWaves;
    
    if nargout == 1
      varargout{1} =  PressureWaves;
    elseif (nargout>1)
      error('It is requested more output varaibles than allowed!')  
    end
end