%%
% CalcBodyCoverParameters: Function for implementing the computation of the 
% biomechanical parameters of the Body Cover Model resulting from the
% laryngeal posture and the activation of the instrinsic muscles. This
% method  allows for the dynamical simulation of muscle control of vocal
% folds oscillations.
%
% Structure: BCMParam = CalcBodyCoverParameters(BCMObj), 
% where
%
% BCMObj: is an object from BodyCoverModel o TriangularBodyCoverModel
%           (handle) classes,  
% BCMParam: Struct gathering the computed model parameters.
%
% References:
% [1] I. R. Titze and B. H. Story, “Rules for controlling low-dimensional 
%     vocal fold models with muscle activation,” J. Acoust. Soc. Am., 
%     vol. 112, p. 1064, 2002.
%
% Coded by Gabriel Alzamendi, February 2020.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function BCMParam = CalcBodyCoverParameters(MCObj,BCMObj,varargin)
  % Activation levels
  try
    a_CT = BCMObj.a_CT;
    a_TA = BCMObj.a_TA;
    a_LC = BCMObj.a_LC;
  catch
    error('An object of BodyCoverModel or TriangularBodyCoverModel classes is required! See function definition.')
  end
 
  % Check the input: Normalized activation levels
  if (nargin>2)
    if (nargin==3)&&isnumeric(varargin{1})&&(varargin{1}>=0)&&(varargin{1}<=1)
      a_TA = varargin{1};
    else
      error('Activation level for TA muscle must be a real number in the range [0, 1]!')
    end
  end
  propLig_body=0.5;
  propLig_cover=1-propLig_body;
  
%   if (nargin>=2)
%     if isnumeric(varargin{1})&&(varargin{1}>=0)&&(varargin{1}<=1)
%       a_CT = varargin{1};
%     else
%       error('Activation level for CT muscle must be a real number in the range [0, 1]!')
%     end
%     
%     if (nargin>=3)
%       if isnumeric(varargin{2})&&(varargin{2}>=0)&&(varargin{2}<=1)
%         a_TA = varargin{2};
%       else
%         error('Activation level for TA muscle must be a real number in the range [0, 1]!')
%       end
%       
%       if (nargin==4)
%         if isnumeric(varargin{3})&&(varargin{3}>=0)&&(varargin{3}<=1)
%           a_LC = varargin{3};
%         else
%           error('Activation level for LC muscle must be a real number in the range [0, 1]!')
%         end
%       elseif  (nargin>4)
%         error('More input variables than allowed. Read function description!')  
%       end
%     end
%   end
      
  % Elongation rule [-]
%   epsilon = MuscleActivation.G*(MuscleActivation.R*a_CT - a_TA) - MuscleActivation.H*a_LC;
  epsilon = MCObj.getStrainVF; % + 0.12;
  
  % Fold elongation 
  Lg = BCMObj.Lg0*(1+epsilon); % (EN EL PAPER DICE L = L0*(1+epsilon) pero en la enmienda aparece L = L0/(1+epsilon);!!!!!) <------!!!!
%   Lg = BCMObj.Lg0/(1+epsilon); % (EN EL PAPER DICE L = L0*(1+epsilon) pero en la enmienda aparece L = L0/(1+epsilon);!!!!!) <------!!!! THIS IS WRONG!!
  
  % Thickness rule
  Tg = BCMObj.Tg0/(1+0.8*epsilon);
  
  % - Point rule
  Znodal = (1+a_TA)*Tg/3;
  
  % Depth rules
%   if strcmpi(BCMObj.sex, 'male')
%     Dmuc = MCObj.DEPTH_MUC_MALE;
%     Dlig = MCObj.DEPTH_LIG_MALE;
%     Dmus = MCObj.DEPTH_MUS_MALE;
%   elseif strcmpi(BCMObj.sex, 'female')
%     Dmuc = MCObj.DEPTH_MUC_FEMALE;
%     Dlig = MCObj.DEPTH_LIG_FEMALE;
%     Dmus = MCObj.DEPTH_MUS_FEMALE;
%   else
%     error('Error in sex especification in the Body Cover Model')  
%   end
  
  Dmuc = BCMObj.DMuc0;
  Dlig = BCMObj.DLig0;
  Dmus = BCMObj.DMus0;
  D_body = (a_TA*Dmus+propLig_body*Dlig)/(1+0.2*epsilon);
  D_cover = (Dmuc+propLig_cover*Dlig)/(1+0.2*epsilon);
  % Adduction rule
%   xi_02 = 0.25*BCMObj.Lg0*(1-2.0*a_LC); % [m] Upper mass
  [xi_02, ~] = MCObj.getVocalProcessCoord; % [m] Upper mass
  % Convergence rule [cm]
%     xi_c = T*(0.05-0.15*Ata);
%     % Nearly rectangular approach
  xi_c = Tg*tan(0.0001);
%   xi_01 = max([xi_c, xi_c + xi_02]); % [m] Lower mass <- errata of Version 26.01.B (old = xi_c - xi_02);
  xi_01 = xi_c + xi_02;
  
  % Cover and Body Stress
%   sigma_muc = PassiveStress(epsilon,MCObj.EPSILON_1_MUC,MCObj.EPSILON_2_MUC,...
%                     MCObj.SIGMA_0_MUC,MCObj.SIGMA_2_MUC,MCObj.C_MUC); % [Pa]
%                 
%   sigma_lig = PassiveStress(epsilon,MCObj.EPSILON_1_LIG,MCObj.EPSILON_2_LIG,...
%                     MCObj.SIGMA_0_LIG,MCObj.SIGMA_2_LIG,MCObj.C_LIG); % [Pa]
%                 
%   sigma_mus = a_TA*MCObj.SIGMA_M_TA*max(0,1-MCObj.B_TA*(epsilon-MCObj.EPSILON_M_TA)^2) + ...
%               PassiveStress(epsilon,MCObj.EPSILON_1_TA,MCObj.EPSILON_2_TA,...
%                     MCObj.SIGMA_0_TA,MCObj.SIGMA_2_TA,MCObj.C_TA); % [Pa]
  K_mult = 1.0;
  sigma_muc = K_mult*MCObj.LarMuscObj(7).getMuscStress; % [Pa]
                
  sigma_lig = K_mult*MCObj.LarMuscObj(6).getMuscStress; % [Pa]
                
  sigma_mus = K_mult*MCObj.LarMuscObj(5).getMuscStress; % [Pa]
  
  sigma_body = (propLig_body*sigma_lig*Dlig+sigma_mus*Dmus)/D_body; % [Pa]
  sigma_cover = (propLig_cover*sigma_lig*Dlig+sigma_muc*Dmuc)/D_cover; % [Pa]
  
  % Computation of Body cover model parameters
  Tl = Znodal; % [m]
  Tu = Tg-Znodal; % [m]
  % Geometric and elastic parameters for the Cover (l: lower mass, u: upper mass, c: coupling)
  kl = 2*MCObj.SHEARMODULUS_COVER*(Lg*Tg/D_cover)*(Znodal/Tg) ...
        + (pi^2)*sigma_cover*(D_cover/Lg)*Znodal; % [Pa*m] or [N/m]
  ku = 2*MCObj.SHEARMODULUS_COVER*(Lg*Tg/D_cover)*(1-(Znodal/Tg)) ...
        + (pi^2)*sigma_cover*(D_cover/Lg)*Tg*(1-(Znodal/Tg)); % [Pa*m] or [N/m]
  kc = ((1/2)*MCObj.SHEARMODULUS_COVER*(Lg*D_cover/Tg)*((1/3)-(Znodal/Tg)*(1-(Znodal/Tg)))^(-1) ...
        - 2*MCObj.SHEARMODULUS_COVER*(Lg*Tg/D_cover))*(Znodal/Tg)*(1-(Znodal/Tg)); % [Pa*m] or [N/m]
  ml = MCObj.TISSUE_DENS*Lg*Tg*D_cover*(Znodal/Tg); % [kg]
  mu = MCObj.TISSUE_DENS*Lg*Tg*D_cover*(1-(Znodal/Tg)); % [kg]
  % Geometric and elastic parameters for the Body
  kb = 2*MCObj.SHEARMODULUS_BODY*(Lg*Tg/D_body) ...
        + (pi^2)*sigma_body*(D_body/Lg)*Tg; % [Pa*m] or [N/m]
  mb = MCObj.TISSUE_DENS*Lg*Tg*D_body; % [kg]
      
  % Updating the parameters of the Body Cover object
%   BCMObj.a_CT = a_CT;
%   BCMObj.a_TA = a_TA;
%   BCMObj.a_LC = a_LC;
  BCMObj.epsilon = epsilon;
  BCMObj.Lg = Lg;
  BCMObj.Tg = Tg;
  BCMObj.Znodal = Znodal;
  BCMObj.Tl = Tl;
  BCMObj.Tu = Tu;
  BCMObj.kl = kl;
  BCMObj.ku = ku;
  BCMObj.kc = kc;
  BCMObj.ml = ml;
  BCMObj.mu = mu;
  BCMObj.kb = kb;
  BCMObj.mb = mb;
  BCMObj.xu0 = xi_02/2; % max([xi_02 xi_02/2])/2;
  BCMObj.xl0 = xi_01/2; % max([xi_01 xi_01/2])/2;
  BCMObj.xi_01 = xi_01;
  BCMObj.xi_02 = xi_02;
  BCMObj.ParamSet_MuscleRules = true;
  
  % Gathering the resulting parameters in the output struct
  BCMParam.epsilon = epsilon;
  BCMParam.Lg = Lg;
  BCMParam.Tg = Tg;
  BCMParam.Znodal = Znodal;
  BCMParam.Tl = Tl;
  BCMParam.Tu = Tu;
  BCMParam.kl = kl;
  BCMParam.ku = ku;
  BCMParam.kc = kc;
  BCMParam.ml = ml;
  BCMParam.mu = mu;
  BCMParam.kb = kb;
  BCMParam.mb = mb;
  BCMParam.xi_01 = xi_01;
  BCMParam.xi_02 = xi_02;
end

% Passive stress formula
function out = PassiveStress(epsilon,epsilon_1,epsilon_2,sigma_0,sigma_2,C)
  if epsilon < epsilon_1
    out = 0;
  elseif epsilon <= epsilon_2
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1);
  else
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1) + sigma_2*(exp(C*(epsilon-epsilon_2))-C*(epsilon-epsilon_2)-1);
  end
end