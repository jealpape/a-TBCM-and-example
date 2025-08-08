function [dXR,dXL,Lg,X_L,X_R,a_g,Ut,Pgo,Pout,Psub,Pcol,tm] = Simu_aTBCM(a_R,a_L,Pl,vow,gender)
%% Fuction for aTBCM simulation:

% Input:
% a_R = [lca;ia;pca;ct;ta];  Right muscle activation
% a_L = [lca;ia;pca;ct;ta];  Left muscle activation
% Pl = number; Lung Pressure
% vow = /a/; vowel tract
% gender; male or female
%
% Output:
% tm = time;
% dXR = Right VF Posterior desplacement mm; 
% dXL = Right VF Posterior desplacement mm;
% Lg =  VFs length mm; 
% 
% X_L = mm;
% X_R = mm;
% 
% a_g = mm 
% Ut =  mL
% Pgo = mm^2
% 
% Pout = Pa
% Psub = Pa
% Pcol = Pa

%% Simulation
Vowel = vow;
PL = Pl;

% Simulation information
fs = 44100; % [Hz] sampling frequency
Ns_Larynx = 5;
fs_Larynx = fs/Ns_Larynx;
T_pre = 0.2; N_pre = ceil(T_pre*fs);
T_sim = 1.0; N_Sim = ceil(T_sim*fs); N_tot = N_pre + N_Sim; t = (1:N_tot)/fs;
N_tran = ceil(0.5*T_pre*fs);

% Filter design
[b1,a1] = butter(4,50/(44.1e3/2),'high');
[b2,a2] = butter(6,5.5e3/(44.1e3/2));
bf = conv(b1,b2); af = conv(a1,a2);

%% Muscle activation vector
% a_R, a_L
act_cerrada = [0.7;0.7;0;0.5;0.5];

%% Larynx definition
% LARYNX DEFINITION (LEFT)
LCLObj = MuscleControlModel;
LCLObj.setSimulationParameter(fs_Larynx);
LCLObj.SelectModelParameters('Alzamendi2020'); % ('Palaparthi2019')  ('Titze2006')  ('Alzamendi2020')

% LARYNX DEFINITION (RIGHT)
LCRObj = MuscleControlModel;
LCRObj.setSimulationParameter(fs_Larynx);
LCRObj.SelectModelParameters('Alzamendi2020'); % ('Palaparthi2019')  ('Titze2006')  ('Alzamendi2020')

%% Vocal Fold definition
% VOCAL FOLD DEFINITION (LEFT)
VFLObj = TriangularBodyCoverModel(gender);
VFLObj.setSimulationParameter(fs)
VFLObj.setDrivingForceSolver('new'); % ('original')  ('new')
VFLObj.setNonLinMode('on')
VFLObj.setMuscleActivity(a_L(5),a_L(4),a_L(1)-a_L(3));

% VOCAL FOLD DEFINITION (RIGHT)
VFRObj = TriangularBodyCoverModel(gender);
VFRObj.setSimulationParameter(fs)
VFRObj.setDrivingForceSolver('new'); % ('original')  ('new')
VFRObj.setNonLinMode('on')
VFRObj.setMuscleActivity(a_R(5),a_R(4),a_R(1)-a_R(3));

LCLObj.CalcBodyCoverParameters(VFLObj);
LCRObj.CalcBodyCoverParameters(VFRObj);

%% TWO FOLDS DEFINITION= Definición de Colisión + �?rea Glotal+ Alpha(Penetración)
TFObj = TwoFolds;

%% VOCAL TRACT DEFINITION

VTObj = VocalTractModel(gender);
VTObj.SetSolver('WRA')
if gender(1:4)=='male'
    VTObj.getMaleVocalTract_Story2008(['/' Vowel '/']); % for male
else
    VTObj.getFemaleVocalTract_Story1998(['/' Vowel '/']) % for female
end
VTObj.setSimulationParameter(fs)

%% SUBGLOTTAL TRACT DEFINITION
SGTObj = SubglottalTractModel(gender);
SGTObj.SetSolver('WRA')
% SGTObj.getSubglottalTract('no-interaction')
SGTObj.getSubglottalTract('Story_smooth')
VTObj.setSimulationParameter(fs)

%% VT Impulse response computation = Equivalente a segundos en la frecuencia de muestreo
a_R = a_R*heaviside(t-T_pre)+act_cerrada*(heaviside(t-0.00)-heaviside(t-T_pre));
a_L = a_L*heaviside(t-T_pre)+act_cerrada*(heaviside(t-0.00)-heaviside(t-T_pre));

%Variables a guardar para el calculo de las salidas del modelo
a_g = zeros(N_tot,1); X_L = zeros(N_tot,9); X_R = zeros(N_tot,9);

Psub = zeros(N_tot,1); Psup = zeros(N_tot,1); Pout = zeros(N_tot,1);
Pcol = zeros(N_tot,1); Ang = zeros(N_tot,1);
Ug = zeros(N_tot,1); Un = zeros(N_tot,1); Ut = zeros(N_tot,1);
Fu_col = zeros(N_tot,1); Fl_col = zeros(N_tot,1);

alpha_u = zeros(N_tot,1); alpha_l = zeros(N_tot,1); Pgo = zeros(N_tot,1);
dXR = zeros(N_tot,2); dXL = zeros(N_tot,2);  Lg = zeros(N_tot,2); % left, rigth

%% INIT ELEMENTS
VFLObj.InitModel
VFRObj.InitModel

LCLObj.InitModel
LCRObj.InitModel

TFObj.InitModel

VTObj.InitModel
SGTObj.InitModel

% FLOW SOLVER DEFINITION
PGO = 0.0e-6;
constFlow.PGO = PGO;
constFlow.Ae = VTObj.AreaFunction(1);
constFlow.As = SGTObj.AreaFunction(1);
constFlow.mu = 18.36922e-6;  % Air Viscosity [Pa s]
constFlow.rho = 1.146; % [kg m^-3]  air density
constFlow.c = 350; % speed of sound
constFlow.solver = 'LUCERO2015';
constFlow.L = (VFLObj.Lg + VFRObj.Lg)/2;
constFlow.T = (VFLObj.Tg + VFRObj.Tg)/2;

Re_c = 1200;
[a_bp,b_bp] = butter(4,[300 3000]/(fs/2));
Nf = rand(N_tot,1)-0.5;
Nf = filter(a_bp,b_bp,Nf);
r_e = 1.0;
r_s = 1.0;

for cont_sim = 1:N_tot
    %% Subglottal pressure transition
    PL_sim = PL * (sin(pi/2*cont_sim/N_tran)*(heaviside(cont_sim)-heaviside(cont_sim-N_tran)) + ...
                 heaviside(cont_sim-N_tran) );

    %% LARYNX UPDATE
    if rem(cont_sim,Ns_Larynx)==0
        LCLObj.SimulatePosture(a_L(:,cont_sim));
        VFLObj.setMuscleActivity(a_L(5,cont_sim),a_L(4,cont_sim),a_L(1,cont_sim)-a_L(3,cont_sim));
        LCLObj.CalcBodyCoverParameters(VFLObj);

        LCRObj.SimulatePosture(a_R(:,cont_sim));
        VFRObj.setMuscleActivity(a_R(5,cont_sim),a_R(4,cont_sim),a_R(1,cont_sim)-a_R(3,cont_sim));
        LCRObj.CalcBodyCoverParameters(VFRObj);
        %       VFObj.setPGO(LarynxControl.aPGO);
    end

    %% VOCAL FOLD STATE UPDATE
    Ps_plus = SGTObj.xData(1);
    Pe_minus = VTObj.xData(1);
    TFObj.CalcuParam(VFLObj,VFRObj)
    TFObj.Simulate(VFLObj,VFRObj, SGTObj.Pressure(1),VTObj.Pressure(1),VTObj.AreaFunction(1))
    TFObj.CalcuParam(VFLObj,VFRObj)

    X_L(cont_sim,:) = VFLObj.xData';
    X_R(cont_sim,:) = VFRObj.xData';

    a_g(cont_sim) = VFLObj.ag;

    constFlow.PGO = (LCLObj.aPGO + (LCRObj.aPGO))/2;
    
%     constFlow.T = 0.5*(VFLObj.Tg+VFRObj.Tg);
%     constFlow.L = (LCLObj.Lg0 + (LCRObj.Lg0))/2 + 0.5*(abs(LCRObj.Psi_PGO)+abs(LCLObj.Psi_PGO));
%     L_tot = constFlow.L;
    L_tot = 21e-3;
    
    Ug(cont_sim) = vfsolver.solveFlow(Ps_plus,Pe_minus,a_g(cont_sim),constFlow);
    a_g(cont_sim) = VFLObj.ag + (LCLObj.aPGO + LCRObj.aPGO)/2;
    Pgo(cont_sim) = (LCLObj.aPGO + LCRObj.aPGO)/2;

    Re_n = constFlow.rho*Ug(cont_sim)/(L_tot*constFlow.mu);
    Un(cont_sim) = Nf(cont_sim)*(Re_n^2 - Re_c^2)*1e-12*(Re_n>Re_c);
    Ut(cont_sim) = Ug(cont_sim) + Un(cont_sim);
    SGTObj.Simulate(Ut(cont_sim),'PL',PL_sim,r_s)
    VTObj.Simulate(Ut(cont_sim),r_e)
    Pout(cont_sim) = VTObj.xData(end);

    Psubtract = SGTObj.Pressure;
    Psub(cont_sim) = Psubtract(1);
    Psuptract = VTObj.Pressure;
    Psup(cont_sim) = Psuptract(1);
    %cols(cont_sim,:) = TotalCollisionForces(VFObj);
    % get the forces.
    

    [Fu_col(cont_sim), Fl_col(cont_sim)] = TFObj.CollisionForce(VFLObj,VFRObj);
    alpha_u(cont_sim) = VFLObj.alpha_u;
    alpha_l(cont_sim) = VFLObj.alpha_l;
    dXR(cont_sim,:) = [max([0, VFRObj.xi_02]),max([0, VFRObj.xi_01])]; % dxu,dxl
    dXL(cont_sim,:) = [max([0, VFLObj.xi_02]),max([0, VFLObj.xi_02])];
    Lg(cont_sim,:) = [VFLObj.Lg, VFRObj.Lg]; % Left, rigth
    

    if (VFLObj.acont~=0)
    Pcol(cont_sim) = (Fu_col(cont_sim)+Fl_col(cont_sim))/VFLObj.acont;
    end

    Ang(cont_sim) = max(0,LCLObj.angleGlottic/2) + max(0,LCRObj.angleGlottic/2);

end

% Output vectors
tm = round(fs*0.1);
dXR = 1e3*dXR(end-tm:end,:); % dxu,dxl mm
dXL = 1e3*dXL(end-tm:end,:); % mm
Lg = 1e3*Lg(end-tm:end,:); % Left, rigth mm

X_L = 1e3*X_L(end-tm:end,1:2); % mm
X_R = 1e3*X_R(end-tm:end,1:2); % mm

a_g = 1e6*a_g(end-tm:end); %mm 
Ut = 1e6*Ut(end-tm:end); % mL
Pgo = 1e6*Pgo(end-tm:end); % mm²

Pout = Pout(end-tm:end);
Psub = Psub(end-tm:end);
Pcol = Pcol(end-tm:end);

end