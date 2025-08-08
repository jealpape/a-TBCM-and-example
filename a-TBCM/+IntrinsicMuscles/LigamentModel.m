classdef LigamentModel < IntrinsicMuscles.Muscle1DModel
% Handle class for modeling the Vocal Ligament tissue, which is made of
% passive (non-contractile) fibers and, therefore, no longitudinal
% contraction is developed. Nevertheless, Kelvin model can be applied by
% simulating only the passive stress-strain dynamics [1].
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
    T0 = 0.0e-3; % [m] Resting thickness (not available)
    D0 = 2.0e-3; % [m] Resting depth (not available)
    Ac = 6.1e-6; % [m^2] Cross section area (not available)
    M = 0.0e-3; % [kg] Mass (not available)
    rho = 0.0e3; % [kg/m^3] Density (not available)
%     Sigma_0 = 0.0e3; % [Pa] Passive stress at resting length (canine)
%     Sigma_2 = 0.0e3; % [Pa] Scalling of exponential stress (canine)
%     Epsilon_1 = 0.0; % [-] Strain at zero stress (canine)
%     Epsilon_2 = 0.0; % [-] Strain at exponential stress (canine)
%     B = 0.0; % [-] Exponential strain constant (canine)
%     Sigma_m = 0.0e3; % [Pa] Maximum active stress (canine)
    Epsilon_m = 1.0; % [-] Strain at maximum active stress (passive only)
    b = 1.0; % [-] Coefficient for active stress-strain (passive only)
    d_Sigma_m = 1.0; % [s^-1] Maximum strain rate (passive only)
    t_i = 0.010; % [s] Activation time (passive only)
    t_p = 0.100; % [s] Parallel contraction time (canine)
    t_s = 0.080; % [s] Series contraction time (canine)
    % Direction cosines and directional moment arms
    alpha = 0.0; % [-] Direction cosine of the horizontal force
    beta = 0.0; % [-] Direction cosine of the vertical force
    gamma = 0.0e-3; % [m] Directional moment arms (-1.23e-3)
    
    % Muscular parameters from Titze (2006)
    Sigma_0_Tit2006 = 1.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Tit2006 = 1.4e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Tit2006 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Tit2006 = -0.0; % [-] Strain at exponential stress
    B_Tit2006 = 17.0; % [-] Exponential strain constant
    Sigma_m_Tit2006 = 0.0e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Palaparthi et al. (2019)
    Sigma_0_Pal2019 = 2.0e3; % [Pa] Passive stress at resting length
    Sigma_2_Pal2019 = 0.15e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Pal2019 = -0.9; % [-] Strain at zero stress
    Epsilon_2_Pal2019 = -0.5; % [-] Strain at exponential stress
    B_Pal2019 = 12.8; % [-] Exponential strain constant
    Sigma_m_Pal2019 = 0.0e3; % [Pa] Maximum active stress
    
    % Muscular parameters from Alzamendi (2020)
    Sigma_0_Alz2020 = 0.5e3; % [Pa] Passive stress at resting length
    Sigma_2_Alz2020 = 1.4e3; % [Pa] Scalling of exponential stress
    Epsilon_1_Alz2020 = -0.5; % [-] Strain at zero stress
    Epsilon_2_Alz2020 = -0.3; % [-] Strain at exponential stress 0.3
    B_Alz2020 = 13.0; % [-] Exponential strain constant
    Sigma_m_Alz2020 = 0.0e3; % [Pa] Maximum active stress
  end
  
  methods
    % Class constructor
    function MuscObj = LigamentModel
      MuscObj.Model = "LigamentModel";
      MuscObj.InitModel;
    end
  end
end
