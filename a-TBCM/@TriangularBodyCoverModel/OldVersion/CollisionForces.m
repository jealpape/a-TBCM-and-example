%%
% CollisionForces: Function implementing the computation of the 
% collision forces produced during closure phase in accordance to the  
% triangular body-cover model of the vocal folds.
%
% Structure: [Fu_col,Fl_col] = CollisionForces(TBCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fu_col: is the collision force in the upper mass,
% Fl_col: is the collision force in the lower mass.
%
% References:
% [1] G. E. Galindo, S. D. Peterson, B. D. Erath, C. Castro, R. E. Hillman,
%     and M. Zañartu, “Modeling the Pathophysiology of Phonotraumatic Vocal
%     Hyperfunction With a Triangular Glottal Model of the Vocal Folds,” 
%     J. Speech Lang. Hear. Res., vol. 60, no. 9, p. 2452, Sep. 2017.
% [2] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Articulatory 
%     synthesis ofwords in six voice qualities using a modified two-mass 
%     model ofthe vocal folds.” Paper presented at the First International
%     Workshop on Performative Speech and Singing Synthesis, Vancouver,
%     British Columbia, Canada, 2011.
% [3] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Synthesis of 
%     breathy, normal, and pressed phonation using a two-mass model with a 
%     triangular glottis.” In P. Cosi, R. De Mori, G. Di Fabbrizio, & R.
%     Pieraccini (Eds.), Interspeech 2011: 12th Annual Conference of the
%     International Speech Communication Association (pp. 2681–2684).
%     Baixas, France: Interna- tional Speech Communication Association. 2011   
%
% Coded by Gabriel Alzamendi, December 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.

% IMPORTAN!! 
% There are two typos in [1] in the Eqs. (A19) and (A20)
% (A19) ... + ñ_{uCol} (x_u^3 + 1.5 x_u^2 \deltax_u \alpha_u + x_u \deltax_u^2 [+] \alpha_u^2 + 0.25 \deltax_u^3 \alpha_u^3)    -->>(CHANGE FOR)-->>
%       ... + ñ_{uCol} (x_u^3 + 1.5 x_u^2 \deltax_u \alpha_u + x_u \deltax_u^2 \alpha_u^2 + 0.25 \deltax_u^3 \alpha_u^3)  
% (A20) ... + ñ_{lCol} (x_l^3 + 1.5 x_l^2 \deltax_l \alpha_l + x_l \deltax_l^2 [+] \alpha_l^2 + 0.25 \deltax_l^3 \alpha_l^3)    -->>(CHANGE FOR)-->>
%       ... + ñ_{lCol} (x_l^3 + 1.5 x_l^2 \deltax_l \alpha_l + x_l \deltax_l^2 \alpha_l^2 + 0.25 \deltax_l^3 \alpha_l^3)  


function [Fu_col,Fl_col] = CollisionForces(TBCMobj)
  % Mass displacements
  xu = TBCMobj.xData(1);
  xl = TBCMobj.xData(2);
  xu_col = TBCMobj.xu_col;
  xl_col = TBCMobj.xl_col;
  deltau = (xu-xu_col)*(xu<=xu_col);
  deltal = (xl-xl_col)*(xl<=xl_col);
  alphau_aux = TBCMobj.alpha_u;
  alphal_aux = TBCMobj.alpha_l;
  
  % Collision forces
  Fu_col = - TBCMobj.hu_col*TBCMobj.Lg*alphau_aux*(...                 % [N] Collision force in the upper mass
              deltau + 0.5*TBCMobj.xi_02*alphau_aux ...
              + TBCMobj.NonLinMode*TBCMobj.etau_col*( deltau^3 + ...
                1.5*deltau^2*TBCMobj.xi_02*alphau_aux + ...
                deltau*TBCMobj.xi_02^2*alphau_aux^2 + ...
                0.25*TBCMobj.xi_02^3*alphau_aux^3 )            );

  Fl_col = - TBCMobj.hl_col*TBCMobj.Lg*alphal_aux*(...                 % [N] Collision force in the lower mass
              deltal + 0.5*TBCMobj.xi_01*alphal_aux ...
              + TBCMobj.NonLinMode*TBCMobj.etal_col*( deltal^3 + ...
                1.5*deltal^2*TBCMobj.xi_01*alphal_aux + ...
                deltal*TBCMobj.xi_01^2*alphal_aux^2 + ...
                0.25*TBCMobj.xi_01^3*alphal_aux^3 )            );  
end