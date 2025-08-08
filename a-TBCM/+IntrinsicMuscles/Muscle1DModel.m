classdef Muscle1DModel < handle & matlab.mixin.Heterogeneous
% Handle generic class for the implementation of Kelvin model for 
% simulating 1D (longitudinal) muscle contraction. This class is pretended
% to be the foundation for the especific classes for the computational
% models for the five intrinsic laryngeal muscles.
% This implementation is based on theoretical material in [1], Chapter 2.
%
% Reference:
% [1] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation,  
%     1st edition.  National Center for Voice and Speech, 2006.
%
% Coded by Gabriel Alzamendi, February 2020
  
  properties (Constant, Hidden)
   PARAMASET_LIST = {'Palaparthi2019','Titze2006','Alzamendi2020'}; % List of available model parameter sets
  end
  
  properties  (Constant, Hidden, Abstract)
    % Definition according to biomechanical constant tables in [1]. The
    % values below have none physiological/physical meaning for this
    % generic class.
    L0; % [m] Resting length
    T0; % [m] Resting thickness
    D0; % [m] Resting depth
    Ac; % [m] Cross section area
    M; % [kg] Mass
    rho; % [kg/m^3] Density
%     Sigma_0; % [Pa] Passive stress at resting length
%     Sigma_2; % [Pa] Scalling of exponential stress
%     Epsilon_1; % [-] Strain at zero stress
%     Epsilon_2; % [-] Strain at exponential stress 
%     B; % [-] Exponential strain constant
%     Sigma_m; % [Pa] Maximum active stress
    Epsilon_m; % [-] Strain at maximum active stress
    b; % [-] Coefficient for active stress-strain
    d_Sigma_m; % [s^-1] Maximum strain rate
    t_i; % [s] Activation time
    t_p; % [s] Parallel contraction time
    t_s; % [s] Series contraction time
    % Direction cosines and directional moment arms
    alpha; % [-] Direction cosine of the horizontal force
    beta; % [-] Direction cosine of the vertical force
    gamma; % [m] Directional moment arms
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006; % [-] Strain at zero stress
    Epsilon_2_Tit2006; % [-] Strain at exponential stress
    B_Tit2006; % [-] Exponential strain constant
    Sigma_m_Tit2006; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019; % [-] Strain at zero stress
    Epsilon_2_Pal2019; % [-] Strain at exponential stress
    B_Pal2019; % [-] Exponential strain constant
    Sigma_m_Pal2019; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020; % [-] Strain at zero stress
    Epsilon_2_Alz2020; % [-] Strain at exponential stress
    B_Alz2020; % [-] Exponential strain constant
    Sigma_m_Alz2020; % [Pa] Maximum active stress
  end
  
  properties  (SetAccess = protected, Hidden)
    % Definition according to biomechanical constant tables in [1]. The
    % values below have none physiological/physical meaning for this
    % generic class.
    Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length
    Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1 = 0.0; % [-] Strain at zero stress
    Epsilon_2 = 0.0; % [-] Strain at exponential stress 
    B = 0.0; % [-] Exponential strain constant
    Sigma_m = 0.0e3; % [Pa] Maximum active stress
  end
  
  properties (SetAccess = protected)
    % Dynamic state variable
    Model = 'Muscle1DModel';
    xData = zeros(4,1); % State variable describinng 1D longitudinal muscle 
                   % contraction according to Kelvin model:
                   % xData = [Sigma_y, d_Sigma_y, Sigma_i, d_Sigma_i]^T,
                   % where:
                   % Sigma_y [Pa] is the resulting muscle stress,
                   % d_Sigma_y [Pa s^-1] is the muscle stress rate,
                   % Sigma_i [Pa] is the internal active muscle stress,
                   % d_Sigma_y [Pa s^-1] is the active muscle stress rate.
    ModelParamSet = 'Palaparthi2019';
    n_IterCont = 0; % Simulation time index
    % Simulation parameters
    fs = 0; % [Hz] Sampling frequency for the numerical simulation
    Ts = 0; % [s] Sampling period for the numerical simulation
    SimParamOK = false;  % If 'true' simulation parameters are set,
                         % otherwise simulation parameters are missing
  end
  
  methods (Sealed)
    % Class constructor
    function MuscObj = Muscle1DModel
      MuscObj.SelectModelParameters;
      MuscObj.InitModel;
    end
      
    function InitModel(MuscObj)
    % Function for initializing the internal muscle state and time index
    % prior to run the simulation of the longitudinal muscle model
      MuscObj.xData = zeros(4,1);
      
      MuscObj.n_IterCont = 0; % Simulation time index
    end
    
    function setSimulationParameter(MuscObj,Param)    
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
          MuscObj.fs = Param;
          MuscObj.Ts = 1/Param;
        else % Param is sampling period in seconds
          MuscObj.Ts = Param;
          MuscObj.fs = 1/Param;
        end
        MuscObj.SimParamOK = true;
      else
        error('Incorrect input arguments! See function structures.')
      end
    end
    
    function SelectModelParameters(MuscObj,varargin)
    % Function for setting the solver for the set of muscle model
    % pareameters to use in the simulations.
      if nargin == 1
        MuscObj.ModelParamSet = MuscObj.PARAMASET_LIST{1}; % Default setting
      elseif (nargin == 2)&&(ischar(varargin{1}))&&any(strcmp(varargin{1},MuscObj.PARAMASET_LIST))
        MuscObj.ModelParamSet = varargin{1}; % solver setting
      else
        ModelSetOpts = '';
        for OptsAux = MuscObj.PARAMASET_LIST
          ModelSetOpts = strcat(ModelSetOpts, ' ''', char(OptsAux), ''', ');
        end
        ModelSetOpts = ModelSetOpts(1:end-1);
        error('The acceptable set of model parameters are:%s. \n',ModelSetOpts);
      end
      MuscObj.UpdateModelParameters;
    end
      
    function Sigma_p = PassiveStress(MuscObj, Strain)
    % Function for computing the passive stress in fibrous tissues
      Sigma_p = -(MuscObj.Sigma_0/MuscObj.Epsilon_1)*(Strain-MuscObj.Epsilon_1) + ...
                MuscObj.Sigma_2*( exp(MuscObj.B*(Strain-MuscObj.Epsilon_2)) - 1 ...
                - MuscObj.B*(Strain-MuscObj.Epsilon_2) ) * (Strain > MuscObj.Epsilon_2);
      Sigma_p = Sigma_p*(Strain>=MuscObj.Epsilon_1);
    end
      
    function E_young = TangentYoungModulus(MuscObj, Strain)
    % Function for computing the tangent Young's modulus E = \frac{d Sigma_p}{d Epsilon_y}
      E_young = -(MuscObj.Sigma_0/MuscObj.Epsilon_1) + ...
                MuscObj.B*MuscObj.Sigma_2*( exp(MuscObj.B*(Strain-MuscObj.Epsilon_2)) - 1 ) * (Strain > MuscObj.Epsilon_2);
      E_young = E_young*(Strain>=MuscObj.Epsilon_1);      
    end
    
    function Sigma_a = ActiveStress(MuscObj, ActLevel, Strain, d_Strain)
    % Function for computing the internal active muscle stress.
      f_eps_y = StressStrainFactor(MuscObj,Strain);
      g_eps_y = StressVelocityFactor(MuscObj,d_Strain);
      Sigma_a = ActLevel*MuscObj.Sigma_m*f_eps_y*g_eps_y;
    end
    
    function SimulateStress(MuscObj, ActLevel, Strain, d_Strain)
    % Function for simulating muscle stress development according to Kelvin
    % model considering active contraction and passive tissue effects 
    % (see [1], Section 2.2).
      Sigma_p = PassiveStress(MuscObj, Strain);
      Sigma_a = ActiveStress(MuscObj, ActLevel, Strain, d_Strain);
      E_young = TangentYoungModulus(MuscObj, Strain);
      
      Sigma_y_prev = MuscObj.xData(1);
%       d_Sigma_y_prev = MuscObj.xData(2);
      Sigma_i_prev = MuscObj.xData(3);
%       d_Sigma_i_prev = MuscObj.xData(4);
      
      diff_Sigma_y_ = 1/MuscObj.t_s*(-Sigma_y_prev+Sigma_i_prev+ ...   % Eq (2.64)
                                    Sigma_p+MuscObj.t_p*E_young*d_Strain);
      diff_Sigma_i_ = 1/MuscObj.t_i*(-Sigma_i_prev+Sigma_a); % Eq (2.65)
      
      % Solving the differential equation iteratively by using truncated
      % Taylor series.
      MuscObj.xData(1) = MuscObj.xData(1) + MuscObj.Ts*diff_Sigma_y_;
      MuscObj.xData(2) = diff_Sigma_y_;
      MuscObj.xData(3) = MuscObj.xData(3) + MuscObj.Ts*diff_Sigma_i_;
      MuscObj.xData(4) = diff_Sigma_i_;
      
      % Incrementing time simulation index
      MuscObj.n_IterCont = MuscObj.n_IterCont+1;
    end
    
    function Stress_y = getMuscStress(MuscObj)
    % Function for gettin the muscle stress.    
      Stress_y = MuscObj.xData(1);
    end
    
    function Stress_y = getMuscForce(MuscObj)
    % Function for gettin the muscle stress.    
      Stress_y = MuscObj.Ac*MuscObj.xData(1);
    end
    
    function Stress_i = getMuscActiveStress(MuscObj)
    % Function for gettin the internal active muscle stress.
      Stress_i = MuscObj.xData(3);
    end
    
    function [alp, bet, gam] = getDirCosDirMom(MuscObj)
    % Function for gettin the direction cosines and directional moment arm
    % corresponding to a given intrinsic muscle.
      alp = MuscObj.alpha;
      bet = MuscObj.beta;
      gam = MuscObj.gamma;
    end
  end
  
  methods (Sealed,Hidden)
      
    function  g_eps_y = StressVelocityFactor(MuscObj,d_eps_y)
    % Function to compute the (piecewise continuous) active stress-velocity
    % factor.
      Val_aux = d_eps_y/MuscObj.d_Sigma_m;
      g_eps_y = max([0, (Val_aux+1)/(1-3*Val_aux)])*(d_eps_y<=0) + ...
                (Val_aux+1/9)/(5/9*Val_aux+1/9)*(d_eps_y>0);
    end
    
    function  f_eps_y = StressStrainFactor(MuscObj,eps_y)
    % Function to compute the (piecewise continuous) active stress-velocity
    % factor.
      f_eps_y = max([0, 1 - MuscObj.b*(eps_y-MuscObj.Epsilon_m)^2]);
    end
    
    function  UpdateModelParameters(MuscObj)
    % Function to update the set of muscle model parameters.
      switch MuscObj.ModelParamSet
        case MuscObj.PARAMASET_LIST{1} % Set: 'Palaparthi2019'
          MuscObj.Sigma_0 = MuscObj.Sigma_0_Pal2019; 
          MuscObj.Sigma_2 = MuscObj.Sigma_2_Pal2019; 
          MuscObj.Epsilon_1 = MuscObj.Epsilon_1_Pal2019; 
          MuscObj.Epsilon_2 = MuscObj.Epsilon_2_Pal2019;
          MuscObj.B = MuscObj.B_Pal2019;
          MuscObj.Sigma_m = MuscObj.Sigma_m_Pal2019;
        case MuscObj.PARAMASET_LIST{2} % Set: 'Titze2006'
          MuscObj.Sigma_0 = MuscObj.Sigma_0_Tit2006; 
          MuscObj.Sigma_2 = MuscObj.Sigma_2_Tit2006; 
          MuscObj.Epsilon_1 = MuscObj.Epsilon_1_Tit2006; 
          MuscObj.Epsilon_2 = MuscObj.Epsilon_2_Tit2006;
          MuscObj.B = MuscObj.B_Tit2006;
          MuscObj.Sigma_m = MuscObj.Sigma_m_Tit2006;
        case MuscObj.PARAMASET_LIST{3} % Set: 'Alzamendi2020'
          MuscObj.Sigma_0 = MuscObj.Sigma_0_Alz2020; 
          MuscObj.Sigma_2 = MuscObj.Sigma_2_Alz2020; 
          MuscObj.Epsilon_1 = MuscObj.Epsilon_1_Alz2020; 
          MuscObj.Epsilon_2 = MuscObj.Epsilon_2_Alz2020;
          MuscObj.B = MuscObj.B_Alz2020;
          MuscObj.Sigma_m = MuscObj.Sigma_m_Alz2020;
        otherwise
          error('Invalid set of model parameters!')
      end
    end
    
  end
end
