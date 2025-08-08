%%
% SimulateWRA: Function for simulating acoustic wave propagation throughout
% the (supraglottal) vocal tract representation accordign to a particular
% area function. The function implements the half-time Wave Reflection
% Analogue (WRA) solver taking into account the attenuation effects due to
% wave propagation in the tubulets sections, according to Sec. 6.4 in [1].
%
% Structure: SimulateWRA(VTobj,Ug_n)
%            SimulateWRA(VTobj,Ug_n,re_n)
%            PressureWaves = SimulateWRA(...)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% Ug_n: is the air flow (volume velocity) value for the n-th instant,
% re_n: is the supraglottal reflection coefficient (=1 by default),
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
function varargout = SimulateWRA(VTobj,Ug_n,varargin)
    re_n = 1;
    if (nargin == 3)
      re_n = varargin{1};
    end

    % Definition of simulation parameters
    rho = VTobj.RHO_AIR; % [kg m^-3] Density of the air
    c = VTobj.C_AIR; % [m/s] speed of sound
    fs = VTobj.fs;
    Delta_z = VTobj.Delta_z; % c/(2*fs); % lenght of each tube section [m]
    N_AreaSection=VTobj.N_AreaSection;
    
    % Attenuations factors
    A_att  = 1 - 2*(3.8e-3./sqrt(VTobj.AreaFunction))*Delta_z; 
    
    % Reflections coefficients
    r_coef = (VTobj.AreaFunction(1:N_AreaSection-1) - VTobj.AreaFunction(2:N_AreaSection)) ./ ...
             (VTobj.AreaFunction(1:N_AreaSection-1) + VTobj.AreaFunction(2:N_AreaSection));
         
    % Radiation parameters
    I = (2*fs/c)*(8/3)*sqrt(VTobj.AreaFunction(end)/pi^3);
    R = 128/(9*pi^2);
    b1 = + I - R + R*I;
    b2 = + I + R + R*I;
    c1 = + I - R - R*I;
    c2 = - I - R + R*I;
    lambda1 = (b2+c2)*A_att(end)/b2;
    lambda2 = (-b1+c1)*A_att(end)/b2;
    lambda3 = b1/b2;
    gamma1 = c2/b2*A_att(end);
    gamma2 = c1/b2*A_att(end);
    gamma3 = b1/b2;
    
    % Input acoustic impedance
    Z_Ug = rho*c/VTobj.AreaFunction(1);
    
    % Initialization
    B_p = VTobj.xData(1:N_AreaSection);
    F_p = VTobj.xData(N_AreaSection+1:2*N_AreaSection);
    F_end_p = VTobj.xData(2*N_AreaSection);
    F_end_pp = VTobj.xData(2*N_AreaSection+1);
    Pout_p = VTobj.xData(2*N_AreaSection+2);
        
    %% Sound wave propagation
    % Forward pressure wave at odd junctions and half sample computation
    F_p(1) = Z_Ug*Ug_n + re_n*A_att(1)*B_p(1);
    F_p(3:2:end-1) = (1+r_coef(2:2:end-1)).*A_att(2:2:end-1).*F_p(2:2:end-1) ...
                            - r_coef(2:2:end-1).*A_att(3:2:end).*B_p(3:2:end);

    % Backward pressure wave at even junctions and half sample computation
    B_p(2:2:end-1) = r_coef(2:2:end-1).*A_att(2:2:end-1).*F_p(2:2:end-1) ...
                        + (1-r_coef(2:2:end-1)).*A_att(3:2:end-1).*B_p(3:2:end-1);
    B_p(end) = gamma2*F_end_pp+gamma1*F_end_p + gamma3*B_p(end);

    % Forward pressure wave at even junctions and integer sample computation
    F_p(2:2:end) = (1+r_coef(1:2:end)).*A_att(1:2:end-1).*F_p(1:2:end-1) ...
                            - r_coef(1:2:end).*A_att(2:2:end).*B_p(2:2:end);
    % Backward pressure wave at odd junctions and integer sample computation
    B_p(1:2:end-1) = r_coef(1:2:end).*A_att(1:2:end-1).*F_p(1:2:end-1) ...
                        + (1-r_coef(1:2:end)).*A_att(2:2:end).*B_p(2:2:end);

    %% Radiated pressure
    Pout = lambda2*F_end_pp + lambda1*F_end_p + lambda3*Pout_p;

    %% Storing resulting acoustic waves
    PressureWaves = [B_p; F_p; F_end_p; Pout]; 
    VTobj.xData = PressureWaves;
    
    if nargout == 1
      varargout{1} =  PressureWaves;
    elseif (nargout>1)
      error('It is requested more output varaibles than allowed!')  
    end
end