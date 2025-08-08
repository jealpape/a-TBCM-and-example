classdef TriangularBodyCoverModel < handle
% Handle class for modeling vocal fold oscillations according to the
% triangular body cover model (TBCM) of the vocal folds [2]. This model is 
% an extension of  the three-mass lumped-element body cover model [1] allowing
% a triangle-shape glottal posture due to the prephonatory position of the
% arytenoids cartilages. It is built around on the modified vocal fold
% introduced in [3,4] but taking into account the body-cover layered
% structure of the vocal folds.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body�?cover 
%     model of the vocal folds,�? J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
% [2] G. E. Galindo, S. D. Peterson, B. D. Erath, C. Castro, R. E. Hillman,
%     and M. Zañartu, “Modeling the Pathophysiology of Phonotraumatic Vocal
%     Hyperfunction With a Triangular Glottal Model of the Vocal Folds,�? 
%     J. Speech Lang. Hear. Res., vol. 60, no. 9, p. 2452, Sep. 2017.
% [3] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Articulatory 
%     synthesis ofwords in six voice qualities using a modified two-mass 
%     model ofthe vocal folds.�? Paper presented at the First International
%     Workshop on Performative Speech and Singing Synthesis, Vancouver,
%     British Columbia, Canada, 2011.
% [4] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Synthesis of 
%     breathy, normal, and pressed phonation using a two-mass model with a 
%     triangular glottis.�? In P. Cosi, R. De Mori, G. Di Fabbrizio, & R.
%     Pieraccini (Eds.), Interspeech 2011: 12th Annual Conference of the
%     International Speech Communication Association (pp. 2681–2684).
%     Baixas, France: Interna- tional Speech Communication Association. 2011    
%
% Coded by Gabriel Alzamendi, January 2020.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
  
  properties (Constant, Hidden)
    % Constant parameters describign the anatomical average dimension for
    % normal male and female vocal folds.
    LREST_MALE = 1.6e-2; % [m] Vocal fold length at rest for male subject
    TREST_MALE = 0.3e-2; % [m] Vocal fold thickness at rest for male subject
    DEPTH_MUC_MALE = 0.2e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_LIG_MALE = 0.2e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_MUS_MALE = 0.4e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
      
    LREST_FEMALE = 1.0e-2; % [m] Vocal fold length at rest for female subject
    TREST_FEMALE = 0.2e-2; % [m] Vocal fold thickness at rest for female subject
    DEPTH_MUC_FEMALE = 0.15e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_LIG_FEMALE = 0.15e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_MUS_FEMALE = 0.3e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
  end
  
  properties (SetAccess = {?TwoFolds})
    % Dynamic state variable
    Model = 'TBCM';
    xData = zeros(9,1); % State variable descring vocal fold oscillations:
                   % xData = [xu xl xb vu vl vb au al ab], where
                   %   x_: mass displacement in [m]
                   %   v_: mass velocity in [m/s]
                   %   a_: mass acceleration in [m/s^2]
    n_IterCont = 0; % Simulation time index
    % General object description
    Lg0 = 0; % [m] Vocal fold length at rest
    Tg0 = 0;% [m] Vocal fold thickness at rest
    sex = '';
    scale = 1.0; % Scale factor for converting the volume of the vocal fold (VFvol_out = scale * VFvol_rest)
    A_PGO = 0; % [m^2]
    
    % Simulation parameters
    fs = 0; % [Hz] Sampling frequency for the numerical simulation
    Ts = 0; % [s] Sampling period for the numerical simulation
    SimParamOK = false;  % If 'true' simulation parameters are set,
                         % otherwise simulation parameters are missing
    NonLinMode = 1; % Binary factor setting on (=1) or off (=0) the 
                    % non-linear terms in the mechanical rules.
    NonLinDamping = 1; % Binary factor setting on (=1) or off (=0) the 
                    % non-linear terms in the damping rules.
  end
    
  properties (GetAccess = public)
    % Structural parameters
    Lg_init = 0; % [m] Vocal fold length at rest
    Tg_init = 0;% [m] Vocal fold thickness at rest
    DMuc_init = 0;% [m] Mucosa depth at rest
    DLig_init = 0;% [m] Ligament depth at rest
    DMus_init = 0;% [m] TA muscle depth at rest
    
    % Biomechanical parameters
    etau = 100*1e4; % [1/m^2] Non-linear spring coefficient for the upper mass
    etal = 100*1e4; % [1/m^2] Non-linear spring coefficient for the lower mass
    etab = 100*1e4; % [1/m^2] Non-linear spring coefficient for the body mass
    zetau0 = 0.6; % 0.4; % [-] Basic damping factor for the upper mass
    zetal0 = 0.1; % 0.4; % [-] Basic damping factor for the lower mass
    zetab = 0.15; % 0.2; % [-] Damping factor for the body mass
    zetauCol = 1.0; % 0.4; % [-] Collision damping factor for the upper mass
    zetalCol = 1.0; % 0.4; % [-] Collision damping factor for the lower mass
    
    % Collision biomechanical parameters
    xu_col = 0; % [m] Displacement where collision occurs for upper mass
    xl_col = 0; % [m] Displacement where collision occurs for lower mass
    etau_col = 500*1e4; % [1/m^2] Non-linear spring coefficient during collision for the upper mass
    etal_col = 500*1e4; % [1/m^2] Non-linear spring coefficient during collision for the lower mass
    
    % Boolean variables controlling TriangularBodyCoverModel object behavior
    symmetric = true; % If 'true' symmetric lateral displacements are considered, 
                      %  otherwise asymmetric displacement are assumed (INCOMPLETE)
    UseUpdated_ContactRules = false; % If 'true' new contact rules are used,
                                     %  otherwise original expressions are applied
    UseUpdated_AeroDrivingForces = false; % If 'true' new rules for computing the aerodynamic driving pressure are used, 
                                          % otherwise original expression are applied
  end
  
  properties (GetAccess = public, SetAccess = {?MuscleActivation,?MuscleControlModel,?TwoFolds})
    % Structural parameters  
    Lg = 0; % [m] Dynamic vocal fold length
    Tg = 0; % [m] Dynamic vocal fold thickness
    Tu = 0; % [m] Dynamic thickness for the upper mass
    Tl = 0; % [m] Dynamic thickness for the lower mass
    DMuc0 = 0;% [m] Mucosa depth at rest
    DLig0 = 0;% [m] Ligament depth at rest
    DMus0 = 0;% [m] TA muscle depth at rest
    epsilon = 0; % [-] Vocal fold elongation
    Znodal = 0; % [m] Nodal point on the medial surface
        
    % Normalized activity levels for the laryngeal muscles controlling the 
    % oscillations of the body Cover model
    a_CT = 0; % [-] Activity level for cricothyroid (CT) muscle
    a_TA = 0; % [-] Activity level for thyroarytenoid (TA) muscle
    a_LC = 0.5; % [-] Activity level for lateral cricoarytenoid (LC) muscle
    
    % Biomechanical parameters
    mu = 0; % [kg] Mass of the upper cover mass
    ml = 0; % [kg] Mass of the lower cover mass
    mb = 0; % [kg] Mass of the body mass
    xu0 = 0; % [m] Initial position of upper mass
    xl0 = 0; % [m] Initial position of lower mass
    xb0 = 3e-3; % [m] Initial position of body mass
    ku = 0; % [N/m] Linear spring constant for the upper mass
    kl = 0; % [N/m] Linear spring constant for the lower mass
    kc = 0; % [N/m] Coupling spring constant for the cover layer
    kb = 0; % [N/m] Linear spring constant for the body mass
    xi_01 = 0; % [m] Vocal process horizontal displacement for the lower mass 
    xi_02 = 0; % [m] Vocal process horizontal displacement for the upper mass
    % Boolean variables controlling TriangularBodyCoverModel object behavior
    
    alpha_u = 0;
    alpha_l = 0;
    au = 0;
    al = 0;
    
    ParamSet_MuscleRules = false;  % If 'true' model parameters are set trough muscle rules, 
                                   % otherwise model parameters need to be computed
  end
  
  properties (Dependent)%, Access = private)
    % Biomechanical parameters
    du % [Ns/m] Damping coefficient for the upper mass
    dl % [Ns/m] Damping coefficient for the lower mass
    db % [Ns/m] Damping coefficient for the body mass
    zetau % [-] Damping factor for the upper mass
    zetal % [-] Damping factor for the lower mass
    
    % Collision biomechanical parameters
    hu_col % [N/m] Linear spring constant during collision for the upper mass
    hl_col % [N/m] Linear spring constant during collision for the upper mass
    
    % Area variables
    ag % [m2] Proyected glottal area
    acont % [m2] Glottal contact area
  end
  
  methods
    % Class constructor
    function TBCMobj = TriangularBodyCoverModel(varargin)
      if nargin ==0
        TBCMobj.sex = 'male';
      elseif (nargin == 1)&&(ischar(varargin{1}))
        if (strcmpi(varargin{1},'male'))
          TBCMobj.sex = 'male';
          TBCMobj.Lg_init = TriangularBodyCoverModel.LREST_MALE;
          TBCMobj.Tg_init = TriangularBodyCoverModel.TREST_MALE;
          TBCMobj.DMuc_init = TriangularBodyCoverModel.DEPTH_MUC_MALE;
          TBCMobj.DLig_init = TriangularBodyCoverModel.DEPTH_LIG_MALE;
          TBCMobj.DMus_init = TriangularBodyCoverModel.DEPTH_MUS_MALE;
          InitModel(TBCMobj);
          scalingVocalFold(TBCMobj);
        elseif (strcmpi(varargin{1},'female'))
          TBCMobj.sex = 'female';
          TBCMobj.Lg_init = TriangularBodyCoverModel.LREST_FEMALE;
          TBCMobj.Tg_init = TriangularBodyCoverModel.TREST_FEMALE;
          TBCMobj.DMuc_init = TriangularBodyCoverModel.DEPTH_MUC_FEMALE;
          TBCMobj.DLig_init = TriangularBodyCoverModel.DEPTH_LIG_FEMALE;
          TBCMobj.DMus_init = TriangularBodyCoverModel.DEPTH_MUS_FEMALE;
          InitModel(TBCMobj);
          scalingVocalFold(TBCMobj);
        else
          error('Acceptable ''Sex'' especification are ''male'' or ''female'' ')
        end
      elseif (nargin == 1)&&(isa(varargin{1},'TriangularBodyCoverModel'))
        BCMAux=TriangularBodyCoverModel(varargin{1}.sex);
        TBCMobj.sex = varargin{1}.sex;
        TBCMobj.Lg_init = BCMAux.Lg_init;
        TBCMobj.Tg_init = BCMAux.Tg_init;
        TBCMobj.DMuc_init = BCMAux.DMuc_init;
        TBCMobj.DLig_init = BCMAux.DLig_init;
        TBCMobj.DMus_init = BCMAux.DMus_init;
        TBCMobj.xData = varargin{1}.xData;
        TBCMobj.n_IterCont = varargin{1}.n_IterCont;
        TBCMobj.Lg0 = varargin{1}.Lg0;
        TBCMobj.Tg0 = varargin{1}.Tg0;
        TBCMobj.scale = varargin{1}.scale;
        TBCMobj.fs = varargin{1}.fs;
        TBCMobj.Ts = varargin{1}.Ts;
        TBCMobj.SimParamOK = varargin{1}.SimParamOK;
        TBCMobj.Lg = varargin{1}.Lg;
        TBCMobj.Tg = varargin{1}.Tg;
        TBCMobj.Tu = varargin{1}.Tu;
        TBCMobj.Tl = varargin{1}.Tl;
        TBCMobj.DMuc0 = varargin{1}.DMuc0;
        TBCMobj.DLig0 = varargin{1}.DLig0;
        TBCMobj.DMus0 = varargin{1}.DMus0;
        TBCMobj.epsilon = varargin{1}.epsilon;
        TBCMobj.Znodal = varargin{1}.Znodal;
        TBCMobj.a_CT = varargin{1}.a_CT;
        TBCMobj.a_TA = varargin{1}.a_TA;
        TBCMobj.a_LC = varargin{1}.a_LC;
        TBCMobj.mu = varargin{1}.mu;
        TBCMobj.ml = varargin{1}.ml;
        TBCMobj.mb = varargin{1}.mb;
        TBCMobj.xu0 = varargin{1}.xu0;
        TBCMobj.xl0 = varargin{1}.xl0;
        TBCMobj.xb0 = varargin{1}.xb0;
        TBCMobj.ku = varargin{1}.ku;
        TBCMobj.kl = varargin{1}.kl;
        TBCMobj.kc = varargin{1}.kc;
        TBCMobj.kb = varargin{1}.kb;
        TBCMobj.xi_01 = varargin{1}.xi_01;
        TBCMobj.xi_02 = varargin{1}.xi_02;
        TBCMobj.alpha_u = varargin{1}.alpha_u;
        TBCMobj.alpha_l = varargin{1}.alpha_l;
        TBCMobj.au = varargin{1}.au;
        TBCMobj.al = varargin{1}.al;
        TBCMobj.ParamSet_MuscleRules = varargin{1}.ParamSet_MuscleRules;
      else
        error('Incorrect input arguments! Valid options: ''male'', ''female'' or no argument at all!')
      end
      
    end
        
    function zetau = get.zetau(TBCMobj) % CHECK!
    % Function for computing zetau the damping factor for the upper mass 
      zetau = TBCMobj.zetau0 + TBCMobj.zetauCol*TBCMobj.alpha_u;
%       zetau = TBCMobj.zetau0 * (1.0+1.0*TBCMobj.alpha_u);
    end
    
    function zetal = get.zetal(TBCMobj) % CHECK!
    % Function for computing zetau the damping factor for the lower mass 
      zetal = TBCMobj.zetal0 + TBCMobj.zetalCol*TBCMobj.alpha_l;
%       zetal = TBCMobj.zetal0 * (1.0+3.0*TBCMobj.alpha_l); 
    end
    
    function hu_col = get.hu_col(TBCMobj) % CHECK!
    % Function for computing hu_col the linear spring constant during collision for the upper mass
      hu_col = 3*TBCMobj.ku; % [N/m]
    end
    
    function hl_col = get.hl_col(TBCMobj) % CHECK!
    % Function for computing hl_col the linear spring constant during collision for the lower mass
      hl_col = 3*TBCMobj.kl; % [N/m]
    end
    
    function du = get.du(TBCMobj)
    % Function for computing du the damping coefficient for the upper mass
      du = 2*TBCMobj.zetau*sqrt(TBCMobj.mu*TBCMobj.ku); % [Ns/m]
    end
    
    function dl = get.dl(TBCMobj)
    % Function for computing dl the damping coefficient for the lower mass
      dl = 2*TBCMobj.zetal*sqrt(TBCMobj.ml*TBCMobj.kl); % [Ns/m]
    end
    
    function db = get.db(TBCMobj)
    % Function for computing du the damping coefficient for the upper mass
      db = 2*TBCMobj.zetab*sqrt(TBCMobj.mb*TBCMobj.kb); % [Ns/m]
    end
    
    function ag = get.ag(TBCMobj)
    % Function for computing ag the proyected glottal area
      ag = max([0, min([TBCMobj.au, TBCMobj.al])]); % [m^2]
    end
    
    function acont = get.acont(TBCMobj)
    % Function for computing acont the contact area
      acont = max([0, TBCMobj.Lg*(TBCMobj.Tu*TBCMobj.alpha_u+TBCMobj.Tl*TBCMobj.alpha_l)]); % [m^2]
    end
    
    function set.a_CT(TBCMobj,act_val)
      if isnumeric(act_val)&&(act_val>=0)&&(act_val<=1)
        TBCMobj.a_CT = act_val;
      else
        error('Activation level for CT muscle must be a real number in the range [0, 1]!')
      end
    end
    
    function set.a_TA(TBCMobj,act_val)
      if isnumeric(act_val)&&(act_val>=0)&&(act_val<=1)
        TBCMobj.a_TA = act_val;
      else
        error('Activation level for TA muscle must be a real number in the range [0, 1]!')
      end
    end
    
    function set.a_LC(TBCMobj,act_val)
      if isnumeric(act_val)&&(act_val>=-1)&&(act_val<=1)
        TBCMobj.a_LC = act_val;
      else
        error('Activation level for LC muscle must be a real number in the range [-1, 1]!')
      end
    end
    
    function InitModel(TBCMobj)
    % Function for initializing the dynamic state and simulation time index
    % prior to run the simulation of the Body Cover Model
      TBCMobj.xData = zeros(9,1);
      TBCMobj.xData(1) = TBCMobj.xu0; % 0.5*TBCMobj.xi_02; % 
      TBCMobj.xData(2) = TBCMobj.xl0; % 0.5*TBCMobj.xi_01; % 
      TBCMobj.xData(3) = TBCMobj.xb0;
      
      TBCMobj.n_IterCont = 0; % Simulation time index
    end
    
    function setMuscleActivity(TBCMobj,a_TA,a_CT,a_LC)
      if (nargin==4)
        TBCMobj.a_TA = a_TA;
        TBCMobj.a_CT = a_CT;
        TBCMobj.a_LC = a_LC;
      else
        error('Incorrect number of imput variables!')
      end
    end
    
    function setState(TBCMobj,Xstate)
    % Function for setting the dynamic state of the body cover model
      Xstate = Xstate(:);
      [row,col] = size(Xstate);
      if (row==9)&&(col==1)
        TBCMobj.xData = Xstate;
      elseif (row==6)&&(col==1)
        TBCMobj.xData = [Xstate; zeros(3,1)];
      elseif (row==3)&&(col==1)
        TBCMobj.xData = [Xstate; zeros(6,1)];
      else
        error('Incorrect state vector dimension!')
      end
    end
    
    function setDrivingForceSolver(TBCMobj,version)
    % Function for setting the rules version for computing the driving
    % forces due to aerodinamic pressure
      if ischar(version)&&strcmpi(version,'original')
        TBCMobj.UseUpdated_ContactRules = false; 
        TBCMobj.UseUpdated_AeroDrivingForces = false; 
      elseif ischar(version)&&strcmpi(version,'new')
        TBCMobj.UseUpdated_ContactRules = true; 
        TBCMobj.UseUpdated_AeroDrivingForces = true; 
      else
        error('Incorrect ''version'' string variable! Valid options are: ''original'', ''new''.')
      end
    end
    
    function setLgInit(TBCMobj,LgVal)
    % Function for setting the rest length of the vocal fold
      TBCMobj.Lg_init = LgVal;
      TBCMobj.scalingVocalFold;
    end
    
    function setNonLinMode(TBCMobj,ModeState)
    % Function for setting on (ModeState = 'on') or off (ModeState = 'off') 
    % the non-linear terms in the rules involved in the Body-Cover 
    % vocal fold model.
      if strcmp(ModeState,'on')
          TBCMobj.NonLinMode = 1;
      elseif strcmp(ModeState,'off')
          TBCMobj.NonLinMode = 0;
      else
          error('Available options: ''on'' for using non-linear terms, ''off'' for disregardin non linear terms.')
      end
    end
    
    function setNonLinDamping(TBCMobj,ModeState)
    % Function for setting on (ModeState = 'on') or off (ModeState = 'off') 
    % the non-linear terms in the damping rules involved in the Body-Cover 
    % vocal fold model.
      if strcmp(ModeState,'on')
          TBCMobj.NonLinDamping = 1;
      elseif strcmp(ModeState,'off')
          TBCMobj.NonLinDamping = 0;
      else
          error('Available options: ''on'' for using non-linear damping terms, ''off'' for disregardin non linear damping terms.')
      end
    end
    
    function agcalc = calcGlottalArea(TBCMobj)
    % Function for computing the glottal area
      agcalc = TBCMobj.ag;
    end
    
    function aucalc = calcAreagu(TBCMobj)
    % Function for computing the glottal area for the upper mass
      aucalc = TBCMobj.au;
    end
    
    function alcalc = calcAreagl(TBCMobj)
    % Function for computing the glottal area for the lower mass
      alcalc = TBCMobj.al;
    end
    
    function acontcalc = calcContactArea(TBCMobj)
    % Function for computing the contact area
      acontcalc = TBCMobj.acont;
    end
    
    function setPGO(TBCMobj,PGO_val)
    % Function for setting the PGO area
      TBCMobj.A_PGO = PGO_val;
    end
    
    % Functions defined on separate files
    [Fku,Fkl,Fkb,Fkc] = ElasticForces(TBCMobj)
    
    [Fdu,Fdl,Fdb] = DampingForces(TBCMobj)
    
    [Fu_col,Fl_col] = CollisionForces(TBCMobj)
    
    [Feu, Fel] = AeroPressure2DrivingForces(TBCMobj, Ps, Pe, varargin)
    
    P_col = CollisionPressure(TBCMobj)
    
    setSimulationParameter(TBCMobj,Param)
    
    Simulate(TBCMobj, Ps, Pe, varargin)
    
    scalingVocalFold(TBCMobj,varargin)
        
  end
% - END CLASS DEFINITION -
end
