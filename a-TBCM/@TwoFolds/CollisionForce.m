function [Fu, Fl] = CollisionForce(TFObj,VFLObj,VFRObj)

%% Upper Vocal Folds
% Right
dxR = max([0, VFRObj.xi_02]);
xR = VFRObj.xData(1);
LR = VFRObj.Lg;
xbR = xR - 0.5*dxR;
nR = dxR/LR;

% Left
dxL = max([0, VFLObj.xi_02]);
xL = -VFLObj.xData(1);
LL = VFLObj.Lg;
xbL = xL + 0.5*dxL;
nL = -dxL/LL;

% Variables
kL = VFLObj.hu_col/LL;
kR = VFRObj.hu_col/LR;

kT = kL*kR/(kL+kR);
xT = xbL-xbR;
mT = nL-nR;

%yc = TFObj.yc_u;
yc = VFRObj.alpha_u*LR; 

% --------- new version --------------
Fu = force(yc,xT,mT,kT,0*VFRObj.NonLinMode*VFRObj.etau_col);
  
%% Lower Vocal Folds
% Right
dxR = max([0, VFRObj.xi_01]);
xR = VFRObj.xData(2);
LR = VFRObj.Lg;
xbR = xR - 0.5*dxR;
nR = dxR/LR;

% Left
dxL = max([0, VFLObj.xi_01]);
xL = -VFLObj.xData(2);
LL = VFLObj.Lg;
xbL = xL + 0.5*dxL;
nL = -dxL/LL;

% Variables
kL = VFLObj.hl_col/LL;
kR = VFRObj.hl_col/LR;

kT = kL*kR/(kL+kR);
xT = xbL-xbR;
mT = nL-nR;

%yc = TFObj.yc_l; 
yc = VFRObj.alpha_l*LR; 
  
Fl = force(yc,xT,mT,kT,0*VFRObj.NonLinMode*VFRObj.etal_col);
% ------------------------------------
  
end


function fcol = force(y,x,m,k,n)
  fcol = (1/4)*k*y*(2*x+m*y)*(2+n*(2*x^2+2*m*x*y+(m*y)^2));
end

