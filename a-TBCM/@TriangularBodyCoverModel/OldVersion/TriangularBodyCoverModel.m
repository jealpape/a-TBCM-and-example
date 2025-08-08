classdef TriangularBodyCoverModel < BodyCoverModel
% Handle class for modeling vocal fold oscillations according to the
% triangular body cover model (TBCM) of the vocal folds [2]. This model is 
% an extension of  the three-mass lump-element body cover model [1] allowing
% a triangle-shape glottal posture due to the prephonatory position of the
% arytenoids cartilages. It is built around on the modified vocal fold
% introduced in [3,4] but taking into account the body-cover layered
% structure of the vocal folds. The TriangularBodyCoverModel class 
% inherits all the constants, atributes and methods of the BodyCoverModel
% class. The modifications ans especial functions required for building
% the TBCM are here implemented. 
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
% [2] G. E. Galindo, S. D. Peterson, B. D. Erath, C. Castro, R. E. Hillman,
%     and M. Zañartu, “Modeling the Pathophysiology of Phonotraumatic Vocal
%     Hyperfunction With a Triangular Glottal Model of the Vocal Folds,” 
%     J. Speech Lang. Hear. Res., vol. 60, no. 9, p. 2452, Sep. 2017.
% [3] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Articulatory 
%     synthesis ofwords in six voice qualities using a modified two-mass 
%     model ofthe vocal folds.” Paper presented at the First International
%     Workshop on Performative Speech and Singing Synthesis, Vancouver,
%     British Columbia, Canada, 2011.
% [4] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Synthesis of 
%     breathy, normal, and pressed phonation using a two-mass model with a 
%     triangular glottis.” In P. Cosi, R. De Mori, G. Di Fabbrizio, & R.
%     Pieraccini (Eds.), Interspeech 2011: 12th Annual Conference of the
%     International Speech Communication Association (pp. 2681–2684).
%     Baixas, France: Interna- tional Speech Communication Association. 2011    
%
% Coded by Gabriel Alzamendi, December 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
  
  properties (SetAccess = protected)
    % Dynamic state variable
  end
  
  properties (GetAccess = protected)
    % Collision biomechanical parameters  
    zetau_col = 0.4; % [-] Collision damping factor for the upper mass
    zetal_col = 0.4; % [-] Collision damping factor for the lower mass
  end
  
  properties (Dependent)
    % Biomechanical parameters
    du % [Ns/m] Damping coefficient for the upper mass
    dl % [Ns/m] Damping coefficient for the lower mass
    db % [Ns/m] Damping coefficient for the body mass
    zetau % [-] Damping factor for the upper mass
    zetal % [-] Damping factor for the lower mass
    
    % Collision biomechanical parameters
    hu_col % [N/m] Linear spring constant during collision for the upper mass
    hl_col % [N/m] Linear spring constant during collision for the upper mass
    zpos_u; % [m] Posterior limit for the collision region in the uppper mass (0<=zpos_u<=Lg).
    zpos_l; % [m] Posterior limit for the collision region in the lower mass (0<=zpos_l<=Lg).
    xrest_u; % [m] Rest position in the middle of the collision region in the upper mass.
    xrest_l; % [m] Rest position in the middle of the collision region in the lower mass.
    
    % Area variables
    alpha_u % Portion of the upper mass in contact (0<=alpha_u<=1).
    alpha_l % Portion of the lower mass in contact (0<=alpha_l<=1).
    au % [m2] Glottal area for the upper portion
    al % [m2] Glottal area for the lower portion
    ag % [m2] Proyected glottal area
    aucont % [m2] Contact area for the upper mass
    alcont % [m2] Contact area for the lower mass
    acont % [m2] Total contact area
  end
  
  methods
    % Class constructor
    function TBCMObj = TriangularBodyCoverModel(varargin)
      if nargin ==0
        InputData = 'male';
      elseif (nargin == 1)&&(ischar(varargin{1}))
        if (strcmpi(varargin{1},'male'))
          InputData = 'male';
        elseif (strcmpi(varargin{1},'female'))
          InputData = 'female';
        else
          error('Acceptable ''Sex'' especification are ''male'' or ''female'' ')
        end
      elseif (nargin == 1)&&(isa(varargin{1},'BodyCoverModel')||isa(varargin{1},'TriangularBodyCoverModel'))
        InputData = varargin{1};
      else
        error('Incorrect input arguments! Valid options: ''male'', ''female'' or no argument at all!')
      end
      
      TBCMObj@BodyCoverModel(InputData);
      TBCMObj.Model = 'TBCM';
    end
    
    function alpha_u = get.alpha_u(TBCMObj) % CHECK!
    % Function for computing alpha_u the portion of the upper mass in
    % contact (0<=alpha_u<=1). 
      xu_f = TBCMObj.xData(1);
      xu_col_f = TBCMObj.xu_col;
      if (xu_f - xu_col_f) > 0
        alpha_u = 0;  
      else
        if TBCMObj.xi_02<=0
          alpha_u = 1;  
        else
          alpha_u = min([1, max([0,-(xu_f-xu_col_f)]) / TBCMObj.xi_02]);
        end
      end
    end
    
    function alpha_l = get.alpha_l(TBCMObj) % CHECK!
    % Function for computing alpha_l the portion of the lower mass in
    % contact (0<=alpha_l<=1). 
      xl_f = TBCMObj.xData(2);
      xl_col_f = TBCMObj.xl_col;
      if (xl_f - xl_col_f) > 0
        alpha_l = 0;
      else
        if TBCMObj.xi_01<=0
          alpha_l = 1;  
        else
          alpha_l = min([1, max([0,-(xl_f-xl_col_f)]) / TBCMObj.xi_01]);
        end
      end
    end
    
    function zpos_u = get.zpos_u(TBCMObj) % CHECK!
    % Function for computing zpos_u the posterior limit for the collision
    % region in the uppper mass
      zpos_u = TBCMObj.Lg*TBCMObj.alpha_u;
    end
    
    function zpos_l = get.zpos_l(TBCMObj) % CHECK!
    % Function for computing zpos_l the posterior limit for the collision
    % region in the lower mass
      zpos_l = TBCMObj.Lg*TBCMObj.alpha_l;
    end
    
    function xrest_u = get.xrest_u(TBCMObj) % CHECK!
    % Function for computing xrest_u the rest position in the middle of the
    % collision region in the upper mass. 
      z_rest = TBCMObj.zpos_u/2;
      xrest_u = (TBCMObj.xi_02/TBCMObj.Lg)*z_rest;
    end
    
    function xrest_l = get.xrest_l(TBCMObj) % CHECK!
    % Function for computing xrest_u the rest position in the middle of the
    % collision region in the upper mass. 
      z_rest = TBCMObj.zpos_l/2;
      xrest_l = (TBCMObj.xi_02/TBCMObj.Lg)*z_rest;
    end
        
    function zetau = get.zetau(TBCMObj) % CHECK!
    % Function for computing zetau the damping factor for the upper mass 
      xu_f = TBCMObj.xData(1);
      xu_col_f = TBCMObj.xu_col;
      if (xu_f > xu_col_f)
        zetau = TBCMObj.zetau0;
      else
%         zetau = 2*TBCMobj.zetau0 + 1.0; % 0.4;  
        zetau = TBCMObj.zetau0 + 0.4*TBCMObj.alpha_u; % 0.4; 
      end
    end
    
    function zetal = get.zetal(TBCMObj) % CHECK!
    % Function for computing zetau the damping factor for the lower mass 
      xl_f = TBCMObj.xData(2);
      xl_col_f = TBCMObj.xl_col;
      if (xl_f > xl_col_f)
        zetal = TBCMObj.zetal0;
      else
%         zetal = 2*TBCMObj.zetal0 + 1.0; % 0.4; 
        zetal = TBCMObj.zetal0 + 0.4*TBCMObj.alpha_l; % 0.4; 
      end
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
    
    function au = get.au(TBCMobj)
    % Function for computing au the glottal area for the upper portion
      deltau = (TBCMobj.xData(1)-TBCMobj.xu_col);
      if deltau<0
        au = max([0,(1-TBCMobj.alpha_u)*TBCMobj.Lg*(deltau + TBCMobj.xi_02)]); % [m^2]
      else
        au = max([0,TBCMobj.Lg*(2*deltau + TBCMobj.xi_02)]); % [m^2]
      end
    end
    
    function al = get.al(TBCMobj)
    % Function for computing al the glottal area for the lower portion
      deltal = (TBCMobj.xData(2)-TBCMobj.xl_col);
      if deltal<0
        al = max([0,(1-TBCMobj.alpha_l)*TBCMobj.Lg*(deltal + TBCMobj.xi_01)]); % [m^2]
      else
        al = max([0,TBCMobj.Lg*(2*deltal + TBCMobj.xi_01)]); % [m^2]
      end
    end
    
    function ag = get.ag(TBCMobj)
    % Function for computing ag the proyected glottal area
      ag = max([0, min([TBCMobj.au, TBCMobj.al])]); % [m^2]
    end
    
    function aucont = get.aucont(TBCMobj)
    % Function for computing aucont the contact area for the upper mass
      aucont = max([0,TBCMobj.alpha_u*TBCMobj.Lg*TBCMobj.Tu]); % [m^2]
    end
    
    function alcont = get.alcont(TBCMobj)
    % Function for computing aucont the contact area for the upper mass
      alcont = max([0,TBCMobj.alpha_l*TBCMobj.Lg*TBCMobj.Tl]); % [m^2]
    end
    
    function acont = get.acont(TBCMobj)
    % Function for computing acont the proyected glottal area
      acont = max([0, (TBCMobj.alcont+TBCMobj.aucont)]); % [m^2]
    end
    
    function agcalc = calcGlottalArea(TBCMobj)
    % Function for computing the glottal area
      agcalc = TBCMobj.ag;
    end
    
    function acontcalc = calcContactArea(TBCMobj)
    % Function for computing the glottal area
      acontcalc = TBCMobj.acont;
    end
    
    % Functions defined on separate files
    [Fu_col,Fl_col] = CollisionForces(TBCMObj)
    
    [Fdu,Fdl,Fdb] = DampingForces(TBCMObj)
    
    [Feu, Fel] = AeroPressure2DrivingForces(TBCMobj, Ps, Pe, varargin)
    
%     setSimulationParameter(TBCMobj,Param)
%     
%     Simulate(TBCMobj, Ps, Pe)
%     
%     scalingVocalFold(TBCMobj,varargin)
        
  end
% - END CLASS DEFINITION -
end
