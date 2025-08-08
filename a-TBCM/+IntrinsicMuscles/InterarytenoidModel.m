classdef InterarytenoidModel < IntrinsicMuscles.Muscle1DModel
% Handle class for modeling of InterarytenoidModel muscle. The model
% simulates the longitudinal contraction of this laryngeal model according
% to the Kelvin model of muscle tissue [1].
% This implementation is based on theoretical material in [1], Chapter 2.
%
% Reference:
% [1] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation,  
%     1st edition.  National Center for Voice and Speech, 2006.
%
% Coded by Gabriel Alzamendi, February 2020
  properties  (Constant, Hidden)
    % Definition according to biomechanical constant Table 2.3 (page 86  in [1]).
    L0 = 9.3e-3; % [m] Resting length (canine)
    T0 = 0.0e-3; % [m] Resting thickness (not available)
    D0 = 0.0e-3; % [m] Resting depth (not available)
    Ac = 12.2e-6; % [m^2] Cross section area (canine)
    M = 0.121e-3; % [kg] Mass (canine)
    rho = 1.04e3; % [kg/m^3] Density (male)
%     Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length (canine)
%     Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress (canine)
%     Epsilon_1 = 0.0; % [-] Strain at zero stress (canine)
%     Epsilon_2 = 0.0; % [-] Strain at exponential stress (canine)
%     B = 0.0; % [-] Exponential strain constant (canine)
%     Sigma_m = 0.0e3; % [Pa] Maximum active stress (canine)
    Epsilon_m = 0.4; % [-] Strain at maximum active stress (canine)
    b = 1.25; % [-] Coefficient for active stress-strain (canine)
    d_Sigma_m = 4.0; % [s^-1] Maximum strain rate (canine)
    t_i = 0.010; % [s] Activation time (canine)
    t_p = 0.100; % [s] Parallel contraction time (canine)
    t_s = 0.060; % [s] Series contraction time (canine)
    % Direction cosines and directional moment arms
    alpha = -0.697; % [-] Direction cosine of the horizontal force
    beta = -0.644; % [-] Direction cosine of the vertical force
    gamma = -3.30e-3; % [m] Directional moment arms (-2.27e-3)
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006 = 2.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006 = 30.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Tit2006 = 0.0; % [-] Strain at exponential stress
    B_Tit2006 = 3.5; % [-] Exponential strain constant
    Sigma_m_Tit2006 = 100e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019 = 2.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019 = 30.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019 = -0.9; % [-] Strain at zero stress
    Epsilon_2_Pal2019 = -0.0; % [-] Strain at exponential stress
    B_Pal2019 = 3.5; % [-] Exponential strain constant
    Sigma_m_Pal2019 = 140e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020 = 2.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020 = 30.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Alz2020 = -0.06; % [-] Strain at exponential stress
    B_Alz2020 = 3.5; % [-] Exponential strain constant
    Sigma_m_Alz2020 = 100e3; % [Pa] Maximum active stress
  end
  
  methods
    % Class constructor
    function MuscObj = InterarytenoidModel
      MuscObj.Model = "InterarytenoidModel";
      MuscObj.InitModel;
    end
  end
end
