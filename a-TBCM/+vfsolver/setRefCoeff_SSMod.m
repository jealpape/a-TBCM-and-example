function [A_VTnew] = setRefCoeff_SSMod(A_VTold, RefCoeff_old, RefCoeff_new, whichTract)
    % Definition of simulation parameters
    
    if nargin==3
        whichTract="supraglottal";
    end
    
    %% State space model for supraglottal tract
%     re=1;
    % Getting supraglottal tract area function
%     sup_Areafun=vfsolver.getTract(supratract_type);
    [L_SSMod,~]=size(A_VTold);
    % Attenuations factors
    if isstring(whichTract)&&strcmp(whichTract,"supraglottal")
        L_SSMod=(L_SSMod-2)/2;
    elseif isstring(whichTract)&&strcmp(whichTract,"subglottal")
        L_SSMod=(L_SSMod)/2;
    else
        error('An incorrect tract selection has been introduced! Please, select one of the valid options: "supraglottal" or "subglottal".')
    end
    A_VTnew = A_VTold;
    A_VTnew(1,1) = A_VTnew(1,1)*RefCoeff_new/RefCoeff_old;
    A_VTnew(L_SSMod+2,1) = A_VTnew(L_SSMod+2,1)*RefCoeff_new/RefCoeff_old;
    A_VTnew(L_SSMod+1,2) = A_VTnew(L_SSMod+1,2)*RefCoeff_new/RefCoeff_old;
    A_VTnew(L_SSMod+1,L_SSMod+1) = A_VTnew(L_SSMod+1,L_SSMod+1)*RefCoeff_new/RefCoeff_old;
    
    
end