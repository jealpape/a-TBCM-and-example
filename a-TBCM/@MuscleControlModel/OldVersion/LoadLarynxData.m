% G. A. Alzamendi, October 2018
% Code for setting the laringeal biomedical parameters for intinsic muscles
% (TA, CT, LCA, IA, PCA), mucosa and, vocal ligament. It is also setting 
% the geometrical variables for cricoaryrtenoid joint (CAJ) and 
% cricothyroid joint (CTJ).
%
% Function structure
%
% LarynxData = LoadLarynxData(L0_VF), where L0_VF is the vocal fold length
%
% LarynxData = LoadLarynxData(L0_VF,'HumanVar',boolval), where 
% L0_VF is the vocal fold length, and boolval is 1 (0) for choosing male 
% (canine) biomedical data.
function LarynxData = LoadLarynxData(L0_VF,varargin)
if nargin==1
    HumanVar=1;
elseif (nargin==2)&&strcmp(varargin{1},'HumanVar')
    if varargin{2}
        HumanVar=1;
    else
        HumanVar=0;
    end
end

    % Initialization
    LarynxData=struct;
    
    % Vocal fold length
    LarynxData.L0=L0_VF;
    
    %% Thyroarytenoid muscle biomedical Data
    %  Extracted from (Titze, 2006), Table 2.3.
%     if HumanVar
%         % Parameters  for human males
%         LarynxData.TA.Length = 18.3*1e-3; % [m] Length at rest
%         LarynxData.TA.CrossArea = 40.9*1e-6; % [m^2] Cross sectional areal
%     else
        % Parameters  for canine
        LarynxData.TA.Length = 21.9*1e-3; % [m] Cricothyroid length at rest
        LarynxData.TA.CrossArea = 63.8*1e-6; % [m^2] Cross sectional areal
%     end
    % Parameters  for canine
    LarynxData.TA.sigma0 = 1.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.TA.sigma2 = 1.5*1e3; % [Pa] Scaling of exponential stress
    LarynxData.TA.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.TA.epsilon2 = -0.05; % Strain at exponential stress
    LarynxData.TA.B = 6.5; % Exponential strain constant
    LarynxData.TA.sigmam = 105.0*1e3; % [Pa] Maximum active stress
    LarynxData.TA.epsilonm = 0.2; % Strain at maximum active stress
    LarynxData.TA.b = 1.07; % Coeff. for active stress-strain
    LarynxData.TA.alpha = 0.015; % Direction coseno for horizontal component
    LarynxData.TA.beta = 0.990; % Direction coseno for vertical component
    LarynxData.TA.gamma = 0.8*1e-3; % [m] Directional moment arm
    
    %% Cricothyroid muscle biomedical Data
    %  Extracted from (Titze, 2006), Table 2.4.
%     if HumanVar
%         % Parameters  for human males
%         LarynxData.CT.Length = 13.8*1e-3; % [m] Length at rest
%         LarynxData.CT.CrossArea = 73.8*1e-6; % [m^2] Cross sectional areal (including 
%                                   %       (including CT fibers that incert 
%                                   %        on the thyroid cartilage)
%     else
        % Parameters  for canine
        LarynxData.CT.Length = 15.2*1e-3; % [m] Cricothyroid length at rest
        LarynxData.CT.CrossArea = 105.3*1e-6; % [m^2] Cross sectional areal
%     end
    % Parameters  for canine
    LarynxData.CT.sigma0 = 2.2*1e3; % [Pa] Passive stress at rest position
    LarynxData.CT.sigma2 = 5.0*1e3; % [Pa] Scaling of exponential stress
    LarynxData.CT.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.CT.epsilon2 = -0.06; % Strain at exponential stress
    LarynxData.CT.B = 7.0; % Exponential strain constant
    LarynxData.CT.sigmam = 3*87.0*1e3; % [Pa] Maximum active stress
    LarynxData.CT.epsilonm = 0.2; % Strain at maximum active stress
    LarynxData.CT.b = 2.4; % Coeff. for active stress-strain
    LarynxData.CT.alpha = 0.0; % Direction coseno for horizontal component
    LarynxData.CT.beta = 0.0; %-1.0; % Direction coseno for vertical component
    LarynxData.CT.gamma = 0.0*1e-3; % [m] Directional moment arm
    
    %% Lateral Cricoarytenoid muscle biomedical Data
    %  Extracted from (Titze, 2006), Table B.1.
    % Parameters  for canine
    LarynxData.LCA.Length = 14.4*1e-3; % [m] Cricothyroid length at rest
%     if HumanVar
%         % Parameters  for human males
%         LarynxData.LCA.CrossArea = 11.9*1e-6; % [m^2] Cross sectional areal
%     else
        % Parameters  for canine
        LarynxData.LCA.CrossArea = 21.2*1e-6; % [m^2] Cross sectional areal
%     end
    % Parameters  for canine
    LarynxData.LCA.sigma0 = 3.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.LCA.sigma2 = 59.0*1e3; % [Pa] Scaling of exponential stress
    LarynxData.LCA.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.LCA.epsilon2 = 0.05; % Strain at exponential stress
    LarynxData.LCA.B = 4.0; % Exponential strain constant
    LarynxData.LCA.sigmam = 1.5*96.0*1e3; % [Pa] Maximum active stress
    LarynxData.LCA.epsilonm = 0.4; % Strain at maximum active stress
    LarynxData.LCA.b = 2.37; % Coeff. for active stress-strain
    LarynxData.LCA.alpha = -0.198; % Direction coseno for horizontal component
    LarynxData.LCA.beta = 0.886; % Direction coseno for vertical component
    LarynxData.LCA.gamma = 3.915*1e-3; % [m] Directional moment arm
    
    %% Interarytenoid muscle biomedical Data
    %  Extracted from (Titze, 2006), Table B.2.
    % Parameters  for canine
    LarynxData.IA.Length = 9.3*1e-3; % [m] Cricothyroid length at rest
%     if HumanVar
%         % Parameters  for human males
%         LarynxData.IA.CrossArea = 24.5*1e-6; % [m^2] Cross sectional areal
%     else
        % Parameters  for canine
        LarynxData.IA.CrossArea = 12.2*1e-6; % [m^2] Cross sectional areal
%     end
    % Parameters  for canine
    LarynxData.IA.sigma0 = 2.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.IA.sigma2 = 30.0*1e3; % [Pa] Scaling of exponential stress
    LarynxData.IA.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.IA.epsilon2 = 0.00; % Strain at exponential stress
    LarynxData.IA.B = 3.5; % Exponential strain constant
    LarynxData.IA.sigmam = 0.6*96.0*1e3; % [Pa] Maximum active stress
    LarynxData.IA.epsilonm = 0.4; % Strain at maximum active stress
    LarynxData.IA.b = 1.25; % Coeff. for active stress-strain
    LarynxData.IA.alpha = -0.697; % Direction coseno for horizontal component
    LarynxData.IA.beta = -0.644; % Direction coseno for vertical component
    LarynxData.IA.gamma = -3.30*1e-3; % [m] Directional moment arm
    
    %% Posterior Cricoarytenoid muscle biomedical Data
    %  Extracted from (Titze, 2006), Table B.3.
    % Parameters  for canine
    LarynxData.PCA.Length = 15.0*1e-3; % [m] Cricothyroid length at rest
%     if HumanVar
%         % Parameters  for human males
%         LarynxData.PCA.CrossArea = 48.1*1e-6; % [m^2] Cross sectional areal
%     else
        % Parameters  for canine
        LarynxData.PCA.CrossArea = 34.4*1e-6; % [m^2] Cross sectional areal
%     end
    % Parameters  for canine
    LarynxData.PCA.sigma0 = 5.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.PCA.sigma2 = 55.0*1e3; % [Pa] Scaling of exponential stress
    LarynxData.PCA.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.PCA.epsilon2 = 0.1; % Strain at exponential stress
    LarynxData.PCA.B = 5.3; % Exponential strain constant
    LarynxData.PCA.sigmam = 96.0*1e3; % [Pa] Maximum active stress
    LarynxData.PCA.epsilonm = 0.4; % Strain at maximum active stress ** This value must be checked **
    LarynxData.PCA.b = 1.86; % Coeff. for active stress-strain
    LarynxData.PCA.alpha = -0.1; % Direction coseno for horizontal component
    LarynxData.PCA.beta = -0.8; % Direction coseno for vertical component
    LarynxData.PCA.gamma = -5.49*1e-3; % [m] Directional moment arm
    
    %% Vocal Fold Mucosal tissue biomedical Data
    %  Extracted from (Titze, 2006), Table 2.5.
    if HumanVar
        % Parameters  for human males
        LarynxData.Muc.Length = 16*1e-3; % [m] Length at rest
    else
        % Parameters  for canine
        LarynxData.Muc.Length = 16*1e-3; % [m] Cricothyroid length at rest
    end
    % Parameters  for human males
    LarynxData.Muc.CrossArea = 5*1e-6; % [m^2] Cross sectional areal
    % Parameters  for canine
    LarynxData.Muc.sigma0 = 1.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.Muc.sigma2 = 9.0*1e3; % [Pa] Scaling of exponential stress
    LarynxData.Muc.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.Muc.epsilon2 = -0.35; % Strain at exponential stress
    LarynxData.Muc.B = 4.4; % Exponential strain constant
    LarynxData.Muc.alpha = 0.0; % Direction coseno for horizontal component
    LarynxData.Muc.beta = 0.0; % Direction coseno for vertical component
    LarynxData.Muc.gamma = 0.0*1e-3; % [m] Directional moment arm
    
    %% Vocal Fold Ligament tissue biomedical Data
    %  Extracted from (Titze, 2006), Table 2.6.
    if HumanVar
        % Parameters  for human males
        LarynxData.Lig.Length = 16*1e-3; % [m] Length at rest
    else
        % Parameters  for canine
        LarynxData.Lig.Length = 16*1e-3; % [m] Cricothyroid length at rest
    end
    % Parameters  for human males
    LarynxData.Lig.CrossArea = 6.1*1e-6; % [m^2] Cross sectional areal
    LarynxData.Lig.sigma0 = 1.0*1e3; % [Pa] Passive stress at rest position
    LarynxData.Lig.sigma2 = 1.4*1e3; % [Pa] Scaling of exponential stress
    LarynxData.Lig.epsilon1 = -0.5; % Strain at zero stress
    LarynxData.Lig.epsilon2 = -0.0; % Strain at exponential stress
    LarynxData.Lig.B = 17.; % Exponential strain constant
    LarynxData.Lig.alpha = 0.0; % Direction coseno for horizontal component
    LarynxData.Lig.beta = 0.0; % Direction coseno for vertical component
    LarynxData.Lig.gamma = 0.0*1e-3; % [m] Directional moment arm
    

    %% Cricothyroid Joint (CTJ) biomechanical constants
    % Extracted from (Titze, 2006), Table 3.4
    if HumanVar
        % Parameters  for human males
        LarynxData.CTJ.cos_phi=0.76; % cosine of CT angle
        LarynxData.CTJ.h_TA=16.1*1e-3; % [m] TA moment arm
        LarynxData.CTJ.w_CT=11.1*1e-3; % [m] CT moment arm
    else
        % Parameters  for canine
        LarynxData.CTJ.cos_phi=0.281; % cosine of CT angle
        LarynxData.CTJ.h_TA=16.25*1e-3; % [m] TA moment arm
        LarynxData.CTJ.w_CT=16.29*1e-3; % [m] CT moment arm
    end
    LarynxData.CTJ.k_t=1500; % [N/m] Translational stiffness (linearized). There are a more accurate, non-linear rule on (Titze, 2006) Pag. 127.
    LarynxData.CTJ.k_r=0.082; % [N m/rad] Rotational stiffness (linearized). There are a more accurate, non-linear rule on (Titze, 2006) Pag. 125.
    
    % Parameters for a Arytenoid cartilage cadaveric position based on
    % canine probes
    LarynxData.x_CAJ=7.1*1e-3; % [m] Center of CAJ x coordinate
    LarynxData.y_CAJ=-7.1*1e-3; % [m] Center of CAJ y coordinate
    LarynxData.xbar_02=4e-3; % [m] cadaveric x position of the vocal process
%     R_CA=12.8*1e-3; % [m] Distance from CAJ to vocal process tip
    LarynxData.R_CA=sqrt(LarynxData.y_CAJ^2+(LarynxData.x_CAJ-LarynxData.xbar_02)^2); % [m] Distance from CAJ to vocal process tip

end