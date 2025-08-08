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
% Coded by Gabriel Alzamendi, January 2019.
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
  xu_vproc = TBCMobj.xi_02;
  xl_vproc = TBCMobj.xi_01;
  
  deltau0 = xu_vproc-xu_col; % Posterior position of the upper mass in the TBCM.    
  deltal0 = xl_vproc-xl_col; % Posterior position of the lower mass in the TBCM.
  alphau = TBCMobj.alpha_u;
  alphal = TBCMobj.alpha_l;
  
  xu_neg = xu*(xu<=0);
  triu_col = alphau*deltau0;
  
  xl_neg = xl*(xl<=0);
  tril_col = alphal*deltal0;
  
  % Collision forces
  Fu_col = - TBCMobj.hu_col*alphau*(...                 % [N] Collision force in the upper mass   *TBCMobj.Lg
              xu_neg + 0.5*triu_col ...
              + TBCMobj.NonLinMode*TBCMobj.etau_col*( xu_neg^3 + ...
                1.5*xu_neg^2*triu_col + xu_neg*triu_col^2 + 0.25*triu_col^3 )            );

  Fl_col = - TBCMobj.hl_col*alphal*(...                 % [N] Collision force in the lower mass    *TBCMobj.Lg
              xl_neg + 0.5*tril_col ...
              + TBCMobj.NonLinMode*TBCMobj.etal_col*( xl_neg^3 + ...
                1.5*xl_neg^2*tril_col + xl_neg*tril_col^2 + 0.25*tril_col^3 )            );  
end


%   Fu_col = - TBCMobj.hu_col*TBCMobj.Lg*alphau*(...                 % [N] Collision force in the upper mass
%               deltau + 0.5*TBCMobj.xi_02*alphau ...
%               + TBCMobj.NonLinMode*TBCMobj.etau_col*( deltau^3 + ...
%                 1.5*deltau^2*TBCMobj.xi_02*alphau + ...
%                 deltau*TBCMobj.xi_02^2*alphau^2 + ...
%                 0.25*TBCMobj.xi_02^3*alphau^3 )            );
% 
%   Fl_col = - TBCMobj.hl_col*TBCMobj.Lg*alphal*(...                 % [N] Collision force in the lower mass
%               deltal + 0.5*TBCMobj.xi_01*alphal ...
%               + TBCMobj.NonLinMode*TBCMobj.etal_col*( deltal^3 + ...
%                 1.5*deltal^2*TBCMobj.xi_01*alphal + ...
%                 deltal*TBCMobj.xi_01^2*alphal^2 + ...
%                 0.25*TBCMobj.xi_01^3*alphal^3 )            );  