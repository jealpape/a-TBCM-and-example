%%
% DampingForces: Function implementing the computation of the (linear and
% nonlinear) elastic forces involved in the body cover model of the vocal
% folds.
%
% Structure: [Fdu,Fdl,Fdb] = DampingForces(TBCMobj)
% TBCMobj: is an object from BodyCoverModel (handle) class,
% Fdu: is the resulting elastic force in the upper mass,
% Fdl: is the resulting elastic force in the lower mass,
% Fdb: is the resulting elastic force in the body mass.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
% [2] J. C. Lucero and L. L. Koenig, “Simulations of temporal patterns of
%     oral airflow in men and women using a two-mass model of the vocal
%     folds under dynamic control,” J. Acoust. Soc. Am., vol. 117, no. 3,
%     pp. 1362–1372, Mar. 2005.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Fdu,Fdl,Fdb] = DampingForces(TBCMobj)
  % Mass velocities
  xu = TBCMobj.xData(1);
  xl = TBCMobj.xData(2);
  xb = TBCMobj.xData(3);
  vu = TBCMobj.xData(4);
  vl = TBCMobj.xData(5);
  vb = TBCMobj.xData(6);
  
  xb0 = TBCMobj.xb0;
%   nonLinearTerm = 0;
  % Damping Equations, with the modification introduced in [2]
  deltaxu = xu*(xu>=0);
  deltaxl = xl*(xl>=0);

  Fdu = -TBCMobj.du*(1+TBCMobj.NonLinDamping*0.0e3*deltaxu)*(vu-vb);  % Damping force of upper mass [N]
  Fdl = -TBCMobj.dl*(1+TBCMobj.NonLinDamping*0.0e3*deltaxl)*(vl-vb);  % Damping force of lower mass [N]
  Fdb = -TBCMobj.db*(1+TBCMobj.NonLinDamping*0.0e3*abs(max(deltaxu,deltaxl)))*(vb);     % Damping force of body mass [N]
end