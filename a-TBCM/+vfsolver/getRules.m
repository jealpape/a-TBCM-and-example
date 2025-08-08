% const_rulesBCM: Function of constants for the body cover model from Titze
% paper "Rules for controlling low-dimensional vocal fold models" JASA,
% 2002
%
% Output units are SI
% 
% OUT = const_rulesBCM(Act,Ata,Alc,varargin)
% 
%   Act [-] : Activation of the cricothyroid Muscle [0 to 1]
%   Ata [-] : Activation of the thyroarytenoid Muscle [0 to 1]
%   Alc [-] : Activation of the cricoarytenoid Muscle [-1 to 1]
% 
%   VARARGIN:
%     'male'(default) or 'female' : Identify the sex of the parameters
% 
%   OUT:
%     OUT.L    [m]  : Vocal folds length
%     OUT.T1   [m]  : Thickness lower mass
%     OUT.T2   [m]  : Thickness upper mass
%     OUT.k1   [N/m]: Spring constant of lower mass
%     OUT.k2   [N/m]: Spring constant of upper mass
%     OUT.kc   [N/m]: Coupling spring constant
%     OUT.kb   [N/m]: Spring constant of body mass
%     OUT.m1   [kg] : Weight of lower mass
%     OUT.m2   [kg] : Weight of upper mass
%     OUT.mb   [kg] : Weight of body mass
%     OUT.x_01 [m]  : rest position of lower mass
%     OUT.x_02 [m]  : rest position of upper mass 
% 
% Observations: This code has enmendations to the actual paper, based on
% direct comunication with Brad Story.
%
% GGF 02.2016

% function out = const_rulesBCM(Act,Ata,Alc,varargin)
function out = getRules(varargin)
  % This is a modification to make it compatible with the old rules file
    % Configuration
      if nargin == 1
        if ~isstruct(varargin{1})
          error('Invalid configuration variable, See vfconfig.m for more info');
        else
          config = varargin{1};
        end
      elseif nargin >= 3
        % Unisex Constants
          config.Act   = varargin{1};  % Activation of the Cricothyroid [-]
          config.Ata   = varargin{2};  % Activation of the Thyroarytenoid [-]
          config.Alc   = varargin{3};  % Activation of the Lateral Cricoarytenoid [-]
          
          config.G     = 0.2;  % Gain of elongation [-]
          config.R     = 3.0;  % Torque ratio [-]
          config.H     = 0.3;  % Adductory strain factor [-]
          config.rho   = 1.040e-3;  % Tissue density (0.001040 [kg/cm^3])
          config.mu_c  = 0.500;  % Shear moduli of cover layers (0.500 [kPa])
          config.mu_b  = 0.100;  % Shear moduli of body layer (1 [kPa])

        % Sex dependent Constants
          if nargin >= 4
            config.sex = varargin{4};
          else
            config.sex = 'male';
          end
          
          if strcmpi(config.sex,'female')
            config.L0   = 1.0; % Lenght of the vocal fold (1.6 [cm] for males and 1.0 [cm] for females)
            config.T0   = 0.2; % Thickness of the vocal folds (0.3 [cm] for males and 0.2 [cm] for females)
            config.Dmuc = 0.15; % Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
            config.Dlig = 0.15; % Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
            config.Dmus = 0.3; % Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
          elseif strcmpi(config.sex,'male')
            config.L0   = 1.6; % Lenght of the vocal fold (1.6 [cm] for males and 1.0 [cm] for females)
            config.T0   = 0.3; % Thickness of the vocal folds (0.3 [cm] for males and 0.2 [cm] for females)
            config.Dmuc = 0.2; % Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
            config.Dlig = 0.2; % Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
            config.Dmus = 0.4; % Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
          else
            error('Unrecognized sex specification');
          end  
          
          if nargin == 5 % adjustement of the proportion parameter
            config.L0   = config.L0*varargin{5};
            config.T0   = config.T0*varargin{5};
            config.Dmuc = config.Dmuc*varargin{5};
            config.Dlig = config.Dlig*varargin{5};
            config.Dmus = config.Dmus*varargin{5};
          end
      else
        error('Invalid number of parameters for configuration');
      end
  
  % Paper Code
  Act   = config.Act;
  Ata   = config.Ata;
  Alc   = config.Alc;
  
  G     = config.G; % Gain of elongation [-]
  R     = config.R; % Torque ratio [-]
  H     = config.H; % Adductory strain factor [-]
  rho   = config.rho; % Tissue density (0.001040 [kg/cm^3])
  mu_c  = config.mu_c; % Shear moduli of cover layers (0.500 [kPa])
  mu_b  = config.mu_b; % Shear moduli of body layer (1 [kPa])

  L0    = config.L0; % Lenght of the vocal fold (1.6 [cm] for males and 1.0 [cm] for females)
  T0    = config.T0; % Thickness of the vocal folds (0.3 [cm] for males and 0.2 [cm] for females)
  Dmuc  = config.Dmuc; % Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
  Dlig  = config.Dlig; % Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
  Dmus  = config.Dmus; % Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
    
  % Elongation rule [-]
    epsilon = G*(R*Act - Ata) - H*Alc;
  % Thickness rule [cm]
    T = T0/(1+0.8*epsilon);
  % Nodal Point rule [cm]
    Zn = (1+Ata)*T/3;
  % Depth rules [cm]
    Db = (Ata*Dmus+0.5*Dlig)/(1+0.2*epsilon);
    Dc = (Dmuc+0.5*Dlig)/(1+0.2*epsilon);
  % Adduction rule [cm]
    xi_02 = 0.25*L0*(1-2.0*Alc);
  % Convergence rule [cm]
%     xi_c = T*(0.05-0.15*Ata);
%     % Nearly rectangular approach
      xi_c = T*tan(0.0001);

    xi_01 = xi_c + xi_02; % <- errata of Version 26.01.B (old = xi_c - xi_02);

  % Elongation
    L = L0/(1+epsilon);
%     L = L0*(1+epsilon); % [cm] (EN EL PAPER DICE L = L0*(1+epsilon) pero en la enmienda aparece L = L0/(1+epsilon);!!!!!) <------!!!!
    
  % Stress [kPa]
    epsilon_1_muc = -0.5; % [-]
    epsilon_2_muc = -0.35; % [-] (EN EL PAPER DICE +0,35!!!!!) <------!!!!
    sigma_0_muc   = 0.5; %  0.5 [kPa]
    sigma_2_muc   = 20.0; %  30.0 [kPa]
    C_muc         = 4.4; % [-]

    sigma_muc = sigma_p(epsilon,epsilon_1_muc,epsilon_2_muc,sigma_0_muc,sigma_2_muc,C_muc); % [kPa]

    epsilon_1_lig = -0.5; % [-]
    epsilon_2_lig = -0.00; % [-]
    sigma_0_lig   = 0.4; %  <-- [dyne/cm2] original 0.4 [kPa]
    sigma_2_lig   = 1.393; %  <-- [dyne/cm2] original 1.393 [kPa]
    C_lig         = 17.0; % [-]

    sigma_lig = sigma_p(epsilon,epsilon_1_lig,epsilon_2_lig,sigma_0_lig,sigma_2_lig,C_lig); % [kPa]

    epsilon_1_mus = -0.5; % [-]
    epsilon_2_mus = -0.05; % [-]
    sigma_0_mus   = 1.0; %  1.00 [kPa]
    sigma_2_mus   = 1.50; % 1.50 [kPa]
    C_mus         = 6.5; % [-]
    sigma_m       = 105; % maximum active stress in the TA muscle <-- [dyne/cm2] original 105 [kPa]
    epsilon_m     = 0.4; % [-]
    b             = 1.07; % [-]
  
    sigma_mus = Ata*sigma_m*max(0,1-b*(epsilon-epsilon_m)^2)+sigma_p(epsilon,epsilon_1_mus,epsilon_2_mus,sigma_0_mus,sigma_2_mus,C_mus); % [kPa]
    
    sigma_b = (0.5*sigma_lig*Dlig+sigma_mus*Dmus)/Db; % [kPa]
    sigma_c = (0.5*sigma_lig*Dlig+sigma_muc*Dmuc)/Dc; % [kPa]

  % Parameters
    T1 = Zn; % [cm]
    T2 = T-Zn; % [cm]
    % Cover 
%       Ic = rho*L*T*Dc*(T^2)*((1/3) - (Zn/T)*(1-(Zn/T))); % Moment of inertia
      k1 = 2*mu_c*(L*T/Dc)*(Zn/T) + (pi^2)*sigma_c*(Dc/L)*Zn; % [kPa*cm]
      k2 = 2*mu_c*(L*T/Dc)*(1-(Zn/T)) + (pi^2)*sigma_c*(Dc/L)*T*(1-(Zn/T)); % [kPa*cm]
      kc = ((1/2)*mu_c*(L*Dc/T)*((1/3)-(Zn/T)*(1-(Zn/T)))^(-1) - 2*mu_c*(L*T/Dc))*(Zn/T)*(1-(Zn/T)); % [kPa*cm]
      m1 = rho*L*T*Dc*(Zn/T); % [kg]
      m2 = rho*L*T*Dc*(1 - (Zn/T)); % [kg]
    % Body
      K = 2*mu_b*(L*T/Db) + (pi^2)*sigma_b*(Db/L)*T; % [kPa*cm]
      M = rho*L*T*Db; % [kg]

  % Output [SI]
    out.L = L*1e-2; % [cm] -> [m]
    out.T1 = T1*1e-2; % [cm] -> [m]
    out.T2 = T2*1e-2; % [cm] -> [m]
    out.k1 = k1*1e1; % [kPa cm] -> [N/m]
    out.k2 = k2*1e1; % [kPa cm] -> [N/m]
    out.kc = kc*1e1; % [kPa cm] -> [N/m]
    out.kb = K*1e1; % [kPa cm] -> [N/m]
    out.m1 = m1*1e0; % [kg] -> [kg]
    out.m2 = m2*1e0; % [kg] -> [kg]
    out.mb = M*1e0; % [kg] -> [kg]
    out.x_01 = xi_01*1e-2; % [cm] -> [m]
    out.x_02 = xi_02*1e-2; % [cm] -> [m]
end

% Passive stress formula
function out = sigma_p (epsilon,epsilon_1,epsilon_2,sigma_0,sigma_2,C)
  if epsilon < epsilon_1
    out = 0;
  elseif epsilon <= epsilon_2
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1);
  else
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1) + sigma_2*(exp(C*(epsilon-epsilon_2))-C*(epsilon-epsilon_2)-1);
  end
end
