function CalcuParam(TFObj,VFLObj,VFRObj)
% Function to calculate Left and Right Vocal Fold parameters

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

if(dxR == 0 && dxL == 0)
    TFObj.yc_u = min(LR,LL)*(xbL>=xbR);
else
    yc = (xbL-xbR)/(nR-nL);
    TFObj.yc_u = max(0,yc);
end
    
% Alpha
VFRObj.alpha_u = max(0,min(TFObj.yc_u/LR,1));
VFLObj.alpha_u = max(0,min(TFObj.yc_u/LL,1));
% VFLObj.alpha_u = alpha_u(VFLObj);
% VFRObj.alpha_u = alpha_u(VFRObj);


% Area 
AtrR = ((1-VFRObj.alpha_u)^2)*LR*dxR/2;
AtrL = ((1-VFLObj.alpha_u)^2)*LL*dxL/2;
Arec = max(0,-(xbL-xbR))*(LL+LR)/2;

VFRObj.au = AtrR+AtrL+Arec;
VFLObj.au = AtrR+AtrL+Arec;
% VFRObj.au = area_u(VFRObj);
% VFLObj.au = area_u(VFRObj);

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

if(dxR == 0 && dxL == 0)
    TFObj.yc_l = min(LR,LL)*(xbL>=xbR);
else
    yc = (xbL-xbR)/(nR-nL);
    TFObj.yc_l = max(0,yc);
end
    
% Alpha
VFRObj.alpha_l = max(0,min(TFObj.yc_l/LR,1));
VFLObj.alpha_l = max(0,min(TFObj.yc_l/LL,1));
% VFLObj.alpha_l = alpha_l(VFLObj);
% VFRObj.alpha_l = alpha_l(VFRObj);

% Area 
AtrR = ((1-VFRObj.alpha_l)^2)*LR*dxR/2;
AtrL = ((1-VFLObj.alpha_l)^2)*LL*dxL/2;
Arec = max(0,-(xbL-xbR))*(LL+LR)/2;

VFRObj.al = AtrR+AtrL+Arec;
VFLObj.al = AtrR+AtrL+Arec;
% VFRObj.al = area_l(VFRObj);
% VFLObj.al = area_l(VFRObj);

end

function apl = alpha_u(TBCMobj) % CHECK!
    % Function for computing alpha_u the portion of the upper mass in
    % contact (0<=alpha_u<=1).
      dx_u = max([0, TBCMobj.xi_02]); % Posterior position of the upper mass in the TBCM.
      x0u = dx_u/2; % Rest posiotion for the upper mass
      xu = TBCMobj.xData(1); %- 0.5*deltau0;
      
      if dx_u == 0
          apl = 1*(xu<=0);
      else
          apl = max([0, min([-(xu - x0u)/dx_u,1])]);
      end
end
    
function alp = alpha_l(TBCMobj) % CHECK!
    % Function for computing alpha_l the portion of the lower mass in
    % contact (0<=alpha_l<=1). 
      dx_l = max([0, TBCMobj.xi_01]);% Posterior position of the lower mass in the TBCM.
      x0l = dx_l/2; % Rest posiotion for the lower mass
      xl = TBCMobj.xData(2); % - 0.5*deltal0;
      
      if dx_l == 0
          alp = 1*(xl<=0);
      else
          alp = max([0, min([-(xl - x0l)/dx_l,1])]);
      end
end
    
function au = area_u(TBCMobj)
    % Function for computing au the glottal area for the upper portion
      dx_u = max([0, TBCMobj.xi_02]);
      xu = TBCMobj.xData(1);
      xh = max([0,xu-dx_u/2]);
      h = TBCMobj.Lg;
      au = 2*max([0, (1-TBCMobj.alpha_u)*h*(xh+(1-TBCMobj.alpha_u)*dx_u/2)]);
      
end
    
function al = area_l(TBCMobj)
    % Function for computing al the glottal area for the lower portion
      xl = TBCMobj.xData(2);
      dx_l = max([0, TBCMobj.xi_01]);
      xh = max([0,xl-dx_l/2]);
      h = TBCMobj.Lg;
      al = 2*max([0, (1-TBCMobj.alpha_l)*h*(xh+(1-TBCMobj.alpha_l)*dx_l/2)]);
    end
    