classdef PosteriorCricoarytenoidModel < IntrinsicMuscles.Muscle1DModel
% Handle class for modeling of Posterior Cricoarytenoid muscle. The model
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
    L0 = 15.0e-3; % [m] Resting length (canine)
    T0 = 0.0e-3; % [m] Resting thickness (not available)
    D0 = 0.0e-3; % [m] Resting depth (not available)
    Ac = 34.4e-6; % [m^2] Cross section area (canine)
    M = 0.544e-3; % [kg] Mass (male)
    rho = 1.04e3; % [kg/m^3] Density (male)
%     Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length (canine)
%     Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress (canine)
%     Epsilon_1 = 0.0; % [-] Strain at zero stress (canine)
%     Epsilon_2 = 0.0; % [-] Strain at exponential stress (canine)
%     B = 0.0; % [-] Exponential strain constant (canine)
%     Sigma_m = 0.0e3; % [Pa] Maximum active stress (canine)
    Epsilon_m = 0.4; % [-] Strain at maximum active stress (canine)
    b = 1.86; % [-] Coefficient for active stress-strain (canine)
    d_Sigma_m = 4.0; % [s^-1] Maximum strain rate (canine)
    t_i = 0.010; % [s] Activation time (canine)
    t_p = 0.100; % [s] Parallel contraction time (canine)
    t_s = 0.060; % [s] Series contraction time (canine)
    % Direction cosines and directional moment arms
    alpha = -0.1; % [-] Direction cosine of the horizontal force (-0.639)
    beta = -0.95; % [-] Direction cosine of the vertical force (-0.253)
    gamma = -5.49e-3; % [m] Directional moment arms
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006 = 5.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006 = 55.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Tit2006 = 0.1; % [-] Strain at exponential stress
    B_Tit2006 = 5.3; % [-] Exponential strain constant
    Sigma_m_Tit2006 = 100e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019 = 5.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019 = 55.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019 = -0.9; % [-] Strain at zero stress
    Epsilon_2_Pal2019 = 0.1; % [-] Strain at exponential stress
    B_Pal2019 = 5.3; % [-] Exponential strain constant
    Sigma_m_Pal2019 = 100e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020 = 5.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020 = 55.0e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Alz2020 = -0.05; % 0.0; % [-] Strain at exponential stress
    B_Alz2020 = 5.3; % [-] Exponential strain constant
    Sigma_m_Alz2020 = 100e3; % [Pa] Maximum active stress
  end
  
  methods
    % Class constructor
    function MuscObj = PosteriorCricoarytenoidModel
      MuscObj.Model = "PosteriorCricoarytenoidModel";
      MuscObj.InitModel;
    end
  end
end
