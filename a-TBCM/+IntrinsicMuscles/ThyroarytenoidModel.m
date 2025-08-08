classdef ThyroarytenoidModel < IntrinsicMuscles.Muscle1DModel
% Handle class for modeling of Thyroarytenoid muscle. The model simulates
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
    L0 = 17.3e-3; % [m] Resting length (male)
    T0 = 7.0e-3; % [m] Resting thickness (canine)
    D0 = 4.0e-3; % [m] Resting depth (canine)
    Ac = 40.9e-6; % [m^2] Cross section area (male)
    M = 0.8232e-3; % [kg] Mass (male)
    rho = 1.04e3; % [kg/m^3] Density (male)
%     Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length (canine)
%     Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress (canine)
%     Epsilon_1 = 0.0; % [-] Strain at zero stress (canine)
%     Epsilon_2 = 0.0; % [-] Strain at exponential stress (canine)
%     B = 0.0; % [-] Exponential strain constant (canine)
%     Sigma_m = 0.0e3; % [Pa] Maximum active stress (canine)
    Epsilon_m = 0.2; % [-] Strain at maximum active stress (canine)
    b = 1.07; % [-] Coefficient for active stress-strain (canine)
    d_Sigma_m = 6.0; % [s^-1] Maximum strain rate (canine)
    t_i = 0.010; % [s] Activation time (canine)
    t_p = 0.100; % [s] Parallel contraction time (canine)
    t_s = 0.060; % [s] Series contraction time (canine)
    % Direction cosines and directional moment arms
    alpha = 0.015; % [-] Direction cosine of the horizontal force
    beta = 0.990; % [-] Direction cosine of the vertical force
    gamma = 0.8e-3; % [m] Directional moment arms (-1.23e-3)
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006 = 1.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006 = 1.5e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Tit2006 = -0.05; % [-] Strain at exponential stress
    B_Tit2006 = 6.5; % [-] Exponential strain constant
    Sigma_m_Tit2006 = 105e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019 = 2.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019 = 1.5e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019 = -0.9; % [-] Strain at zero stress
    Epsilon_2_Pal2019 = -0.5; % [-] Strain at exponential stress
    B_Pal2019 = 6.5; % [-] Exponential strain constant
    Sigma_m_Pal2019 = 180e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020 = 1.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020 = 1.5e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Alz2020 = -0.05; % [-] Strain at exponential stress
    B_Alz2020 = 6.5; % [-] Exponential strain constant
    Sigma_m_Alz2020 = 150e3; % [Pa] Maximum active stress
  end
  
  methods
    % Class constructor
    function MuscObj = ThyroarytenoidModel
      MuscObj.Model = "ThyroarytenoidModel";
      MuscObj.InitModel;
    end
  end
end
