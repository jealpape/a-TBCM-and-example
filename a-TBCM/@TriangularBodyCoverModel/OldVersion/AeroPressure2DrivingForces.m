%%
% AeroPressure2DrivingForces: Function implementing the computation of the 
% driving forces impiging on the cover masses due to the aerodynamic
% pressure in the glottis. The method modifies the expression
% introduced in the original formulation of the body cover model [1] or
% the new rules introduced later in [2] for the geometry of the triangular
% body cover model [3]. The switching between the solvers in [1] or in [2]
% is controlled by means of the function 'setDrivingForceSolver'. 
%
% Structure: [Feu, Fel] = AeroPressure2DrivingForces(TBCMobj, Ps, Pe, Ae, Ph)
% where
%
% TBCMobj: is an object from TriangularBodyCoverModel (handle) class,
% Ps: is the subglottal pressure in the trachea in Pascals,
% Pe: is the supraglottal pressure in the epilarynx in Pascals,
% Ae (optional, default 5e-4 [m^2]): Epilarynx (first supraglottal) tube area in meters,
% Ph (optional, default 0.0 [Pa]): hydrostatic pressure involved in collitions of cover masses in Pascals,
% Feu: is the resulting elastic force in the upper mass,
% Fel: is the resulting elastic force in the lower mass.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
% [2] I. R. Titze, “Regulating glottal airflow in phonation: Application of
%     the maximum power transfer theorem to a low dimensional phonation 
%     model,” J. Acoust. Soc. Am., vol. 111, no. 1, pp. 367–376, Jan. 2002.
% [3] G. E. Galindo, S. D. Peterson, B. D. Erath, C. Castro, R. E. Hillman,
%     and M. Zañartu, “Modeling the Pathophysiology of Phonotraumatic Vocal
%     Hyperfunction With a Triangular Glottal Model of the Vocal Folds,” 
%     J. Speech Lang. Hear. Res., vol. 60, no. 9, p. 2452, Sep. 2017.
%
% Coded by Gabriel Alzamendi, December 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Feu, Fel] = AeroPressure2DrivingForces(TBCMobj, Ps, Pe, varargin)
  % Check input variables
  if (nargin>=3)&&(nargin<=5)&&isnumeric(Ps)&&isnumeric(Pe)
    Ae = 5e-4; % [m^2] Epilarynx tube area
    Ph = 0.0; % [Pa] Hydrostatic pressure    
    if (nargin>=4)
      if isnumeric(varargin{1})&&(varargin{1}>=0)
        Ae = varargin{1}; % [m^2] Epilarynx tube area
      else
        error('Epilarynx area Ae must be a real positive numer.')
      end
    end
    if (nargin==5)
      if isnumeric(varargin{2})&&(varargin{2}>=0)
        Ph = varargin{2}; % [Pa] Hydrostatic pressure
      else
        error('Hydrostatic pressure must be a real positive numer.')
      end
    end
  else
    error('Incorrect input arguments! See function structure.')
  end
  
  % Current areas for the upper and lower masses
  au = TBCMobj.au; % Upper glottal area [m^2]
  al = TBCMobj.al; % Lower glottal area [m^2]
  am = max([0, min([al au])]);
  alpha_uc = TBCMobj.alpha_u;
  alpha_lc = TBCMobj.alpha_l;
  
  % Driving forces produced during contact phase of the cover masses
  if (am <= 0)
    if ~TBCMobj.UseUpdated_ContactRules % Original formulation introduced in [1]
      if (au <= 0) && (al <= 0)     % Both masses colliding
        Pu = 0;                        % Upper mass pressure [N m^{-2}]
        Pl = 0;                        % Lower mass pressure [N m^{-2}]
      elseif (au <= 0) && (al > 0)  % Upper mass colliding
        Pu = 0;                        % Upper mass pressure [N m^{-2}]
        Pl = (1-alpha_lc)*Ps;                      % Lower mass pressure [N m^{-2}]
      elseif (al <= 0) && (au > 0)  % Lower mass colliding
        Pu = (1-alpha_uc)*Pe;                      % Upper mass pressure [N m^{-2}]
        Pl = 0;                        % Lower mass pressure [N m^{-2}]
      end
    else % New expressions introduced in [2]
      if (au <= 0) && (al <= 0)     % Both masses colliding
        Pu = Ph;                        % Upper mass pressure [N m^{-2}]
        Pl = Ph;                        % Lower mass pressure [N m^{-2}]
      else
%         zc = abs(al)/(abs(au)+abs(al))*(TBCMobj.Tl+TBCMobj.Tu); % Check this expression!
        zc = abs(TBCMobj.xData(2)/(TBCMobj.xData(2)-TBCMobj.xData(1)))*(TBCMobj.Tl+TBCMobj.Tu); % improved expression!
        if (au <= 0) && (al > 0)  % Upper mass colliding
          if zc<TBCMobj.Tl
            Pl=1/TBCMobj.Tl*(zc*Ps+(TBCMobj.Tl-zc)*Ph);
            Pu=Ph;
          else
            Pl=Ps;
            Pu=1/TBCMobj.Tu*((zc-TBCMobj.Tl)*Ps+(TBCMobj.Tl+TBCMobj.Tu-zc)*Ph);
          end
%           Pu = Ph;                        % Upper mass pressure [N m^{-2}]
%           Pl = Ps;                      % Lower mass pressure [N m^{-2}]
        elseif (al <= 0) && (au > 0)  % Lower mass colliding
          if zc<TBCMobj.Tl
            Pl=1/TBCMobj.Tl*(zc*Ph + (TBCMobj.Tl-zc)*Pe);
            Pu=Pe;
          else
            Pl=Ph;
            Pu=1/TBCMobj.Tu*((zc-TBCMobj.Tl)*Ph + (TBCMobj.Tl+TBCMobj.Tu-zc)*Pe);
          end
%           Pu = Pe;                      % Upper mass pressure [N m^{-2}]
%           Pl = Ph;                        % Lower mass pressure [N m^{-2}]
        end
      end
    end
  % Driving forces produced during open phase of the glottis
  else
    if ~TBCMobj.UseUpdated_AeroDrivingForces % Original formulation introduced in [1]
      Pu = (1-alpha_uc)*Pe;                     % Upper mass pressure [N m^{-2}]
      Pl = (1-alpha_lc)*(Ps-(Ps-Pe)*(am/al)^2);   % Lower mass pressure [N m^{-2}]
    else % New expressions introduced in [2]
      % New aerodynamic terms
      spr= 1.2; %separation point ratio
      ad = max(0,min(spr*al,au));
      an = al+TBCMobj.Tl*(au-al)/(TBCMobj.Tl+TBCMobj.Tu);
      zd = min(TBCMobj.Tl+TBCMobj.Tu,max(0,(ad-al)*(TBCMobj.Tl+TBCMobj.Tu)/(au-al)));
      ke = 2*ad/Ae*(1-ad/Ae);
%       ke = 0;
      Pkd=(Ps-Pe)/(1-ke);
      if al>=au     %convergent
        Pl = Ps - Pkd*(au^2/(al*an));
        Pu = Ps - Pkd*(au/an);
      else          %divergent
        if zd>TBCMobj.Tl  %separation point above nodal point
          Pl = Ps - Pkd*(ad^2/(al*an));
          Pu = Ps - Pkd/TBCMobj.Tu*(TBCMobj.Tl+TBCMobj.Tu-zd + (ad/an)*(zd-TBCMobj.Tl));
        else      %separation point below nodal point
          Pl = Ps - Pkd/TBCMobj.Tl*(TBCMobj.Tl-zd+(ad/al)*zd);
          Pu = Ps - Pkd;
        end
      end
    end
  end
  
  % Aerodinamic pressure to driving forces on the cover masses
  Feu = Pu*TBCMobj.Lg*TBCMobj.Tu;  % Upper pressure force [N]
  Fel = Pl*TBCMobj.Lg*TBCMobj.Tl;  % Lower pressure force [N]
end