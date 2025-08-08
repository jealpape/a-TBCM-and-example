%%
% ElasticForces: Function implementing the computation of the (linear and
% nonlinear) elastic forces involved in the body cover model of the vocal
% folds.
%
% Structure: [Fku,Fkl,Fkb,Fkc] = ElasticForces(TBCMobj)
% where
%
% TBCMObj: is an object from BodyCoverModel (handle) class,
% Fku: is the resulting elastic force in the upper mass,
% Fkl: is the resulting elastic force in the lower mass,
% Fkb: is the resulting elastic force in the body mass,
% Fkc: is the resulting coupling force between the upper and lower mases.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Fku,Fkl,Fkb,Fkc] = ElasticForces(TBCMobj)
  
  % Mass displacements
  xb0 = TBCMobj.xb0;
  xu0 = 0.5*max([-1e-4, TBCMobj.xi_02]) - TBCMobj.xu_col; % Rest position for the upper mass according to the TBCM
  xl0 = 0.5*max([TBCMobj.xi_01-TBCMobj.xi_02, TBCMobj.xi_01]) - TBCMobj.xl_col; % Rest position for the lower mass according to the TBCM
  
  xiu_f = TBCMobj.xData(1)-xu0-TBCMobj.xu_col;
  xil_f = TBCMobj.xData(2)-xl0-TBCMobj.xl_col;
  xib_f = TBCMobj.xData(3)-xb0;
  
  % Elastic forces
  Fku = - TBCMobj.ku*((xiu_f-xib_f) + TBCMobj.NonLinMode*TBCMobj.etau*(xiu_f-xib_f)^3); % [N] Elastic force in the upper mass
  Fkl = - TBCMobj.kl*((xil_f-xib_f) + TBCMobj.NonLinMode*TBCMobj.etal*(xil_f-xib_f)^3); % [N] Elastic force in the lower mass
  Fkb = - TBCMobj.kb*(xib_f + TBCMobj.NonLinMode*TBCMobj.etab*xib_f^3); % [N] Elastic force in the body mass
  Fkc = - TBCMobj.kc*(xil_f - xiu_f); % [N] Coupling force between the upper and lower mases

end