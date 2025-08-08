classdef MuscleControlModel < handle
% Implements muscular control of prephonatory laryngeal posture and vocal fold strain-stress state.
% Based on [1], [2] Chapter 3, and [3] for biomechanical parameters of the body cover model.
% Simulates dynamic posture due to activation of intrinsic muscles (CT, TA, LCA, PCA, IA) and tissue layers.
% Includes static methods for computing geometry and mechanics as a function of CT, TA, LC activations.

% References:
% [1] Titze & Hunter (2007), JASA, 121(4), 2254–2260.  
% [2] Titze (2006), The Myoelastic Aerodynamic Theory of Phonation.  
% [3] Titze & Story (2002), JASA, 112, 1064.
%
% Coded by Gabriel Alzamendi, Januaru 2020.
%
  properties (Constant, Hidden)
    % Unisex Constants
    G = 0.2; % Gain of elongation [-]
    R = 3.0; % Torque ratio [-]
    H = 0.2; % Adductory strain factor [-]
    TISSUE_DENS = 1040; % [kg/m^3] Tissue density
    SHEARMODULUS_COVER = 300; % [Pa] Shear modulus of cover layers (0.500 [kPa])
    SHEARMODULUS_BODY = 600; % [Pa] Shear modulus of body layer (1 [kPa])

    % Stress-strain Constants
    % Mucosa
    EPSILON_1_MUC = -0.5; % [-]
    EPSILON_2_MUC = -0.35; % [-] (EN EL PAPER DICE +0,35!!!!!) <------!!!!
    SIGMA_0_MUC   = 500; % [Pa]       0.5 [kPa]
    SIGMA_2_MUC   = 30000; % [Pa]     30 [kPa]
    C_MUC         = 4.4; % [-]
    
    % Ligament
    EPSILON_1_LIG = -0.5; % [-]
    EPSILON_2_LIG = -0.00; % [-]
    SIGMA_0_LIG   = 400; % [Pa]      0.4 [kPa]
    SIGMA_2_LIG   = 1393; % [Pa]     1.393 [kPa] 
    C_LIG         = 17.0; % [-]
    
    % TA Muscle
    EPSILON_1_TA = -0.5; % [-]
    EPSILON_2_TA = -0.05; % [-]
    SIGMA_0_TA   = 1000; % [Pa]     1.0 [kPa]
    SIGMA_2_TA   = 1500; % [Pa]     1.5 [kPa]
    C_TA         = 6.5; % [-]
    SIGMA_M_TA   = 105000; % [Pa] maximum active stress in the TA muscle 105 [kPa]
    EPSILON_M_TA = 0.4; % [-]
    B_TA         = 1.07; % [-]
    
    % CAJ Joint
    t_ar0 = 0.02; % [s] Rotational time constant for CAJ movement
    t_ax = 0.02; % [s] x-translational time constant for CAJ movement
    t_ay = 0.02; % [s] y-translational time constant for CAJ movement
    Ia = 1.6e-6; % [kg m^2] Arytenoid moment of inertia
    Ma = 1.37e-3; % [kg] Arytenoid mass
    y_CAJ = -10.1e-3; % [m] y coordinate of the CAJ center
    x_CAJ = 10.1e-3; % [m] x coordinate of the CAJ cente
    xbar_02 = 4.0e-3; % [m] Cadaveric x coordinate of the vocal process
    R_CAJ = sqrt(MuscleControlModel.y_CAJ^2+(MuscleControlModel.x_CAJ-MuscleControlModel.xbar_02)^2);
    % CT Joint
    cos_phi=0.76; % cosine of CT angle
    h_TA=16.1*1e-3; % [m] TA moment arm
    w_CT=11.1*1e-3; % [m] CT moment arm
    Ir = 10e-6; % [kg m^2] moment of inertia associated with CTJ rotation
    Mt = 10e-3; % [kg] mass associated with CTJ translation
    t_r = 0.04; % [s] Rotational time constant for CTJ movement
    t_t = 0.04; % [s] Translational time constant for CTJ movement
    
  end
  
  properties (SetAccess = protected)
    % Dynamic CAJ and CTJ state variables
    Model = 'MuscleControl';
    xData_CAJ = zeros(9,1); % State variable descring the cricoarytenoid joint (CAJ):
                    % xData_CAJ = [Xia dXia ddXia Psia dPsia ddPsia Thetaa dThetaa ddThetaa], where
                    %  Xi_a: x-displacement of CAJ center in [m]
                    %  dXi_a: x-velocity of CAJ center in [m/s]
                    %  ddXi_a: x-acceleration of CAJ center in [m/s^2]
                    %  Psi_a: y-displacement of CAJ center in [m]
                    %  dPsi_a: y-velocity of CAJ center in [m/s]
                    %  ddPsi_a: y-acceleration of CAJ center in [m/s^2]
                    %  Theta_a: Rotationnal angle of CAJ center in [rad]
                    %  dTheta_a: Rotationnal velocity of CAJ center in [rad/s]
                    %  ddTheta_a: Rotationnal acceleration of CAJ center in [rad/s^2]

    xData_CTJ = zeros(6,1); % State variable descring the cricothyroid joint (CTJ):
                    % xData_CTJ = [Epst dEpst ddEpst Epsr dEpsr ddEpsr], where
                    %  Eps_t: translational VF strain center in [-]
                    %  dEps_t: translational VF strain velocity in [-/s]
                    %  ddEps_t: translational VF strain acceleration in [-/s^2]
                    %  Eps_r: rotational VF strain center in [-]
                    %  dEps_r: rotational VF strain velocity in [-/s]
                    %  ddEps_r: rotational VF strain acceleration in [-/s^2]
                    
    % General object description
    Lg0 = 16e-3; % [m] Vocal fold length at rest
    sex = '';
    
    % Intrinsic muscles and conective tissue
    LarMuscObj = []; % Array of Laryngeal muscles/tissues
    LarMusIndex = 1:7; % Index for accesing Laryngeal muscle data, where:
                % =1 for LCA muscle, =2 for IA muscle, =3 for PCA muscle,
                % =4 for the CT muscle, =5 for the TA muscle;
                % =6 for the vocal ligament, and =7 for the mucosa.
	deps_VF; % Estimate of VF strain rate
    deps_CT; % Estimate of CT strain rate
    
    % Parameter describing the posterior glottal opening
    Xi_PGO = 3.0e-3; % Poterior wall (half-width) x-position
    Psi_PGO = -2.0e-3; % Posterior wall y-position
     
    % Simulation parameters
    n_IterCont = 0; % Simulation time index
    fs = 0; % [Hz] Sampling frequency for the numerical simulation
    Ts = 0; % [s] Sampling period for the numerical simulation
    SimParamOK = false;  % If 'true' simulation parameters are set,
                         % otherwise simulation parameters are missing
    LimitCollision = 1;
  end
  
  properties (Dependent)
    t_ar;
    strain_CT; % [-] Total (translational + rotational) CT muscle strain
    aPGO; % [m2] area of the Posterior Glottal Opening (PGO)
    angleGlottic; % [deg] glottic angle 
  end
  
  methods
    % Class constructor
    function MCObj = MuscleControlModel
      MCObj.LarMuscObj = [IntrinsicMuscles.LateralCricoarytenoidModel; ... % LCA muscle
                          IntrinsicMuscles.InterarytenoidModel; ...        % IA muscle 
                          IntrinsicMuscles.PosteriorCricoarytenoidModel; ... % PCA muscle
                          IntrinsicMuscles.CricothyroidModel; ...          % CT muscle 
                          IntrinsicMuscles.ThyroarytenoidModel; ...        % TA muscle 
                          IntrinsicMuscles.LigamentModel; ...              % Vocal Ligament 
                          IntrinsicMuscles.MucosaModel ];                  % VF Mucosa 
      MCObj.sex = 'male-canine';
      MCObj.Lg0 = MCObj.LarMuscObj(5).L0;
    end
    
    function InitModel(MCObj)
    % Function for initializing the dynamic state and simulation time index
    % prior to run the simulation of the Control Larynx pre-phonatory posture
      MCObj.xData_CAJ = zeros(9,1); % Initializing CAJ data
      MCObj.xData_CTJ = zeros(6,1); % Initializing CTJ data
      
      for cont_m = MCObj.LarMusIndex  % Initializing muscle interval variables
        MCObj.LarMuscObj(cont_m).InitModel;
      end
      MCObj.deps_VF = 0;
      MCObj.deps_CT = 0;
      
      MCObj.n_IterCont = 0; % Simulation time index
    end
    
    function setSimulationParameter(MCObj,Param)
    % setSimulationParameter: Function for setting the simulation parameters.
    %
    % Structure 1: setSimulationParameter(fs)
    % where fs is the sampling frequency in Hertz (real positive numbre, >=1).
    %
    % Structure 2: setSimulationParameter(Ts)
    % where Ts is the sampling period in seconds (real positive numbre, <1).
    %
      if isnumeric(Param)&&(Param>0)
        if (Param>=1) % Param is sampling frequency in Hertz
          MCObj.fs = Param;
          MCObj.Ts = 1/Param;
        else % Param is sampling period in seconds
          MCObj.Ts = Param;
          MCObj.fs = 1/Param;
        end
        
        % Seting simulation parameters for the intrinsic muscle models
        for cont_m = MCObj.LarMusIndex 
          MCObj.LarMuscObj(cont_m).setSimulationParameter(Param);
        end 
      
        MCObj.SimParamOK = true;
      else
        error('Incorrect input arguments! See function structures.')
      end
    end
    
    function SelectModelParameters(MCObj,varargin)
    % Function for setting the solver for the set of muscle model
    % pareameters to use in the simulations.
      for cont_m = MCObj.LarMusIndex  % Set muscle model parameters
        if nargin == 1
          MCObj.LarMuscObj(cont_m).SelectModelParameters;
        elseif (nargin == 2)&&(ischar(varargin{1}))
          MCObj.LarMuscObj(cont_m).SelectModelParameters(varargin{1});
        end
      end
    end
    
    function setArytenoidCollisionMode(MCObj,ModeState)
    % Function for setting on (ModeState = 'on') or off (ModeState = 'off') 
    % the physical constrain due to the collision between the left and
    % right arytenoid cartilages.
      if strcmp(ModeState,'on')
          MCObj.LimitCollision = 1;
      elseif strcmp(ModeState,'off')
          MCObj.LimitCollision = 0;
      else
          error('Available options: ''on'' for modeling arytenoid collision, ''off'' for disregarding the collision constraint.')
      end
    end
    
    function strain_CT = get.strain_CT(MCObj)
    % Function for computing strain_CT the total (translational +
    % rotational) CT muscle strain, following Eq. (3.63) in [2].
      L0_CT = MCObj.LarMuscObj(4).L0;
      strain_CT = -MCObj.Lg0/L0_CT*(MCObj.getTranStrain/MCObj.cos_phi + ...
                                 MCObj.w_CT/MCObj.h_TA*MCObj.getRotStrain);
    end
    
    function eps = getStrainVF(MCObj)
    % Function for computing the total (translational + rotational +
    % adductory) VF strain, following Eqs. (3.60) and (3.73) in [2]. 
      eps_r = MCObj.getRotStrain; % Rotational strain
      eps_t = MCObj.getTranStrain; % Translational strain
      eps_a = MCObj.getAddStrain; % Adductory strain
      
      eps = eps_r + eps_t + eps_a;
    end
    
    function [xi_02, psi_02] = getVocalProcessCoord(MCObj,varargin)
    % Function for computing the vocal process (x_02,y_02) coordinates,
    % following Eqs. (3.71) and (3.72) in [2]. 
      xi_a = MCObj.xData_CAJ(1);
      psi_a = MCObj.xData_CAJ(4);
      theta_a = MCObj.xData_CAJ(7);
      if (nargin>1)
        if strcmpi('theta_a',varargin{1})
            theta_a = varargin{2};
        else
            error('The input option is not correct!')
        end
      end
      
      xi_02 = MCObj.x_CAJ - (MCObj.x_CAJ-MCObj.xbar_02)*cos(theta_a) ...  % Eq. (3.71)
                + MCObj.y_CAJ*sin(theta_a) + xi_a;
      psi_02 = MCObj.y_CAJ*(1-cos(theta_a)) ...  % Eq. (3.72)
                - (MCObj.x_CAJ-MCObj.xbar_02)*sin(theta_a) + psi_a ...
                - MCObj.Lg0*(MCObj.getRotStrain+MCObj.getTranStrain);
    end
    
    function aPGO = get.aPGO(MCObj)
    % Function for computing area of the Posterior Glottal Opening (PGO), 
    % following [2], pp. 210-211.
    % PGO has a trapeizodal shape with vertices:
    %           (+/-Xi_02, Psi_02) and (+/-Xi_PGO, Psi_PGO)
    
      Tg = 5e-3; % [m] glottal thickness at the posterior cartilaginous portion
      xi_c = Tg*tan(0.0001);
      [Xi_02, Psi_02] = MCObj.getVocalProcessCoord;
      
      xi_a = MCObj.xData_CAJ(1);
      xi_pnt = MCObj.Xi_PGO+xi_a;
      aPGO = max([0, min([Xi_02+xi_pnt, Xi_02+xi_pnt+xi_c]) * (Psi_02-MCObj.Psi_PGO) ]);
    end
    
    function angGlot = get.angleGlottic(MCObj)
    % Function for computing the glottic angle in degrees
      [xi_02, ~] = MCObj.getVocalProcessCoord;
      Lg = MCObj.Lg0*(1+MCObj.getStrainVF);
      angGlot = 2*atand(xi_02/Lg);
    end
    
    function t_ar = get.t_ar(MCObj)
        [xi_02,~] = MCObj.getVocalProcessCoord;
        if xi_02>0
            t_ar = MCObj.t_ar0;
        else
            t_ar = 10*MCObj.t_ar0;
        end
    end
    
    SimulatePosture(MCObj,a_Act)
    
    BCMParam = CalcBodyCoverParameters(MCObj,BCMObj,varargin)
  end
  
  methods (Hidden)
    
    function eps_t = getTranStrain(MCObj)
    % Function for getting the VF translational strain
      eps_t = MCObj.xData_CTJ(1);
    end
    
    function eps_r = getRotStrain(MCObj)
    % Function for getting the VF rotational strain
      eps_r = MCObj.xData_CTJ(4);
    end
    
    function eps_a = getAddStrain(MCObj)
    % Function for computing the VF adductory strain, following Eq. (3.74)
    % in [2].
      psi_a = MCObj.xData_CAJ(4);
      theta_a = MCObj.xData_CAJ(7);
      
      eps_a = MCObj.y_CAJ*(1-cos(theta_a)) - ...
              (MCObj.x_CAJ-MCObj.xbar_02)*sin(theta_a) + psi_a;
      eps_a = -eps_a/MCObj.Lg0;
    end
    
    SimulateCTJ(MCObj)
    
    SimulateCAJ(MCObj)
    
%     tau_col = TorqueCol(MCObj)
    
    [Fx_col,Fy_col,Tau_col] = ColReaction(MCObj)
    
    theta_a_opt = ArithAngleCollision(MCObj)

  end
  
  methods (Static)
    function Param = getBodyCoverModelParameters(BCMObj)
    % Function for showing the Body Cover parameters handled by the
    % MuscleActivation methods.
      Param=struct;
      try
        Param.Lg = BCMObj.Lg; % [m] Dynamic vocal fold length
        Param.Tg = BCMObj.Tg; % [m] Dynamic vocal fold thickness
        Param.Tu = BCMObj.Tu; % [m] Dynamic thickness for the upper mass
        Param.Tl = BCMObj.Tl; % [m] Dynamic thickness for the lower mass
        Param.epsilon = BCMObj.epsilon; % [-] Vocal fold elongation
        Param.Znodal = BCMObj.Znodal;  % [m] Nodal point on the medial surface
        Param.mu = BCMObj.mu; % [kg] Mass of the upper cover mass
        Param.ml = BCMObj.ml; % [kg] Mass of the lower cover mass
        Param.mb = BCMObj.mb; % [kg] Mass of the body mass
        Param.ku = BCMObj.ku; % [N/m] Linear spring constant for the upper mass
        Param.kl = BCMObj.kl; % [N/m] Linear spring constant for the lower mass
        Param.kc = BCMObj.kc; % [N/m] Coupling spring constant for the cover layer
        Param.kb = BCMObj.kb; % [N/m] Linear spring constant for the body mass
        Param.a_CT = BCMObj.a_CT; % [-] Activity level for cricothyroid (CT) muscle
        Param.a_TA = BCMObj.a_TA; % [-] Activity level for thyroarytenoid (TA) muscle
        Param.a_LC = BCMObj.a_LC; % [-] Activity level for lateral cricoarytenoid (LC) muscle
      catch
        error('The input is not a BodyCoverModel object!')
      end
    end
    
  end
  
  methods (Static, Hidden)
    function [k_x, k_y] = CalcCAJoinTranslationalStiffness(xi_a, psi_a)
    % Function for implementing the calculation of the translational 
    % stiffnesses in the CA joint according to [2] Sec. 3.1.2
    
      % Horizontal translational stiffness k_x in [N/m]
      CAJdisp_max = 3e-3;
      k_x = 60*( 1 + 2e5*xi_a^2 + ...                        % eq. (3.14)
                 (abs(xi_a)>CAJdisp_max)*2e13*(abs(xi_a)-CAJdisp_max)^4 ); % eq. (3.15) if |xi_a|> 3mm
      % Vertical translational stiffness k_y in [N/m]  
      k_y = 200*( 1 + 2e5*psi_a^2 + ...                         % eq. (3.16)
                  (abs(psi_a)>CAJdisp_max)*2e13*(abs(psi_a)-CAJdisp_max)^4 ); % eq. (3.17) if |psi_a|> 3mm
    end
    
    function kappa = CalcCAJoinRotationalStiffness(theta_a)
    % Function for implementing the calculation of the rotational 
    % stiffness in the CA joint according to [2] Sec. 3.1.2
    
      % Rotational stiffness kappa in [Nm/rad]  
      kappa = 0.005*( 1 + 10*theta_a^2 + ...                         % eq. (3.18)
                      (abs(theta_a)>0.3)*1e4*(abs(theta_a)-0.3)^4 ); % eq. (3.19) if |theta_a|> 0.3 rad
    end
    
    function k_t = CalcCTJoinTranslationalStiffness(delta_t, delta_r) % delta_t = L0*epsilon_t, delta_r = L0*epsilon_r
    % Function for implementing the calculation of the translational 
    % stiffnesses in the CT joint according to [2] Sec. 3.2.2
    
      % Translational stiffness k_t in [N/m]
      h = MuscleControlModel.h_TA;
      CTJdisp_max = 1.5e-3;
      k_t = 500 * ( 1 + 30/(pi*h)*abs(delta_r) ) * ...
            ( 1 + 1e6*delta_t^2 + ...                            % eq. (3.52)
              (abs(delta_t)>CTJdisp_max)*2e14*(abs(delta_t)-CTJdisp_max)^4 ); % eq. (3.53) if |L_0 epsilon_t|> 1.5mm
    end
    
    function k_r = CalcCTJoinRotationalStiffness(delta_t, delta_r) % delta_t = L0*epsilon_t, delta_r = L0*epsilon_r = h_TA theta
    % Function for implementing the calculation of the rotational 
    % stiffnesses in the CT joint according to [2] Sec. 3.2.2
    
      % Rotational stiffness k_r in [N/m]
      h = MuscleControlModel.h_TA;
      theta = delta_r/h;
      CTJang_max = 0.2;
      k_r = 0.05 * ( 1 + 1500*abs(delta_t) ) * ...
            ( 1 + 40*theta^2 + ...                         % eq. (3.47)
              (abs(theta)>CTJang_max)*2e5*(abs(theta)-CTJang_max)^4 ); % eq. (3.48) if |delta_r/h_TA| = |theta|> 0.2 rad
    end
    
    plotCAJointSpringConstants
    plotCTJointSpringConstants
    
  end
end

