classdef CricothyroidModel < IntrinsicMuscles.Muscle1DModel
% Handle class for modeling of Cricothyroid muscle. The model simulates
% the longitudinal contraction of this laryngeal model according to the
% Kelvin model of muscle tissue [1].
% This implementation is based on theoretical material in [1], Chapter 2.
%
% Reference:
% [1] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation,  
%     1st edition.  National Center for Voice and Speech, 2006.
%
% Coded by Gabriel Alzamendi, February 2020
  properties  (Constant, Hidden)
    % Definition according to biomechanical constant Table 2.3 (page 86  in [1]).
    L0 = 13.8e-3; % [m] Resting length (male)
    T0 = 0.0e-3; % [m] Resting thickness (not available)
    D0 = 0.0e-3; % [m] Resting depth (not available)
    Ac = 73.8e-6; % [m^2] Cross section area (male)
    M = 0.9423e-3; % [kg] Mass (male)
    rho = 1.04e3; % [kg/m^3] Density (male)
%     Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length (canine)
%     Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress (canine)
%     Epsilon_1 = 0.0; % [-] Strain at zero stress (canine)
%     Epsilon_2 = 0.0; % [-] Strain at exponential stress (canine)
%     B = 0.0; % [-] Exponential strain constant (canine)
%     Sigma_m = 0.0e3; % [Pa] Maximum active stress (canine)
    Epsilon_m = -0.0; % [-] Strain at maximum active stress (canine)
    b = 2.4; % [-] Coefficient for active stress-strain (canine)
    d_Sigma_m = 2.2; % [s^-1] Maximum strain rate (canine)
    t_i = 0.010; % [s] Activation time (canine)
    t_p = 0.100; % [s] Parallel contraction time (canine)
    t_s = 0.080; % [s] Series contraction time (canine)
    % Direction cosines and directional moment arms
    alpha = 0; % [-] Direction cosine of the horizontal force
    beta = 0; % [-] Direction cosine of the vertical force
    gamma = 0; % [m] Directional moment arms
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006 = 2.2e3; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006 = 5.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Tit2006 = -0.06; % [-] Strain at exponential stress
    B_Tit2006 = 7.0; % [-] Exponential strain constant
    Sigma_m_Tit2006 = 90e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019 = 2.2e3; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019 = 5.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019 = -0.9; % [-] Strain at zero stress
    Epsilon_2_Pal2019 = -0.06; % [-] Strain at exponential stress
    B_Pal2019 = 7.0; % [-] Exponential strain constant
    Sigma_m_Pal2019 = 400e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020 = 2.2e3; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020 = 5.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Alz2020 = -0.0; % [-] Strain at exponential stress
    B_Alz2020 = 7.0; % [-] Exponential strain constant
    Sigma_m_Alz2020 = 300e3; % [Pa] Maximum active stress
  end
  
  methods
    % Class constructor
    function MuscObj = CricothyroidModel
      MuscObj.Model = "CricothyroidModel";
      MuscObj.InitModel;
    end
  end
end
