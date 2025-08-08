function [A_sup,Gamma_sup,H_sup,A_sub,Gamma_sub,H_sub] = getTract_SSMod(fsampling,sup_Areafun,sub_Areafun,re,rs)
    % Definition of simulation parameters
    rho = 1.146; % [kg m^-3] Density of the air
    c = 350; % [m/s] sound velocity
    fs = fsampling;
    Delta_z = c/(2*fs); % lenght of each tube section [m]
    %% State space model for supraglottal tract
    L_sup=length(sup_Areafun);
    % Attenuations factors
    A_attsup  = diag(1 - (3.8e-3./sqrt(sup_Areafun))*Delta_z);
    % Reflections coefficients
    r_coefsup = (sup_Areafun(1:L_sup-1) - sup_Areafun(2:L_sup))./(sup_Areafun(1:L_sup-1) + sup_Areafun(2:L_sup));
    % Radiation parameters
    I = (2*fs/c)*(8/3)*sqrt(sup_Areafun(end)/pi^3);
    R = 128/(9*pi);
    b1 = + I - R + R*I;
    b2 = + I + R + R*I;
    c1 = + I - R - R*I;
    c2 = - I - R + R*I;
    lambda1 = (b2+c2)*A_attsup(end,end)/b2;
    lambda2 = (-b1+c1)*A_attsup(end,end)/b2;
    lambda3 = b1/b2;
    gamma1 = c2/b2;
    gamma2 = c1/b2*A_attsup(end,end);
    gamma3 = b1/b2;
    K_Ug = rho*c/sup_Areafun(1);
    Mat1_aux=zeros(L_sup,L_sup); Mat1_aux(end,end)=1;
    Mat2_aux=zeros(L_sup,1); Mat2_aux(end)=1;
    % State space model matrices
    A_sup = zeros(2*L_sup+2,2*L_sup+2);
    A_sup = [gamma3*Mat1_aux+diag(1-r_coefsup,1)*A_attsup, ...
             diag([r_coefsup,gamma1])*A_attsup,...
             gamma2*Mat2_aux, 0*Mat2_aux; ... % End row block 1
             diag([re,-r_coefsup(1:end)])*A_attsup,...
             diag(1+r_coefsup(1:end),-1)*A_attsup,...
             0*Mat2_aux, 0*Mat2_aux; ...  % End row block 2
             0*Mat2_aux',Mat2_aux.',0,0; ...  % End row block 3
             0*Mat2_aux',lambda1*Mat2_aux.',lambda2,lambda3];  % End row block 4
    SelVect = ones(2*L_sup+2,1);
    SelVect([L_sup 2*L_sup+1 2*L_sup+2]) = 0*SelVect([L_sup 2*L_sup+1 2*L_sup+2]); % Work!
    SelMat1 = diag(SelVect); SelMat2 = eye(2*L_sup+2) - diag(SelVect);
    
    Gamma_sup = zeros(2*L_sup+2,1); Gamma_sup(L_sup+1) = K_Ug;
    Gamma_sup = A_sup*Gamma_sup+Gamma_sup;
    Gamma_sup=sparse(Gamma_sup);
    
    A_sup = (SelMat1*A_sup^2+SelMat2*A_sup);
    A_sup = sparse(A_sup);
    H_sup = zeros(1,2*L_sup+2); H_sup(end) = 1; H_sup = sparse(H_sup);

    %% State space model for subglottal tract
    L_sub=length(sub_Areafun);
    % Attenuations factors
    A_attsub  = diag(1 - (11.2e-3./sqrt(sub_Areafun))*Delta_z);
    % Reflections coefficients
    r_end = -0.0;
    r_coefsub = (sub_Areafun(1:L_sub-1) - sub_Areafun(2:L_sub))./(sub_Areafun(1:L_sub-1) + sub_Areafun(2:L_sub));
    K_Ug = -rho*c/sub_Areafun(1);
    Mat1_aux=zeros(L_sub,L_sub); Mat1_aux(end,end)=1;
    Mat2_aux=zeros(L_sub,1); Mat2_aux(end)=1;
    % State space model matrices
    A_sub = zeros(2*L_sub,2*L_sub);
    A_sub = [diag(1-r_coefsub,1)*A_attsub, ...
             diag([r_coefsub,r_end])*A_attsub; ... % End row block 1
             diag([rs,-r_coefsub(1:end)])*A_attsub,...
             diag(1+r_coefsub(1:end),-1)*A_attsub];  % End row block 2
    SelVect = ones(2*L_sub,1);
%     SelVect([1 L_sub+1]) = 0*SelVect([1 L_sub+1]);
    SelMat1 = diag(SelVect); SelMat2 = eye(2*L_sub) - diag(SelVect);
    
    Gamma_sub = zeros(2*L_sub,2); Gamma_sub(L_sub+1,1) = K_Ug; Gamma_sub(L_sub,2) = 1.0;
    Gamma_sub = A_sub*Gamma_sub+Gamma_sub;
    Gamma_sub=sparse(Gamma_sub);
    
    A_sub = (SelMat1*A_sub^2+SelMat2*A_sub);
    A_sub = sparse(A_sub);
    H_sub = zeros(1,2*L_sub); H_sub(end) = 1; H_sub = sparse(H_sub);


end