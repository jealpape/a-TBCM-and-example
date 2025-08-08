%%
% getStateSpaceWRAModel: Function for computing the internal matrices for
% the state space representation of the half-time Wave Reflection
% Analogue (WRA) solver taking into account the attenuation effects due to
% wave propagation in the tubulets sections. The WRA model is based on
% Sec. 6.4 in [1], and the state space formulation is based on Galindo's
% ideas [2]. 
%
% Structure: getStateSpaceWRAModel(VTobj)
%            getStateSpaceWRAModel(VTobj,re_n)
%            [A_ss, Gamma_ss] = getStateSpaceWRAModel(...)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% re_n: is the supraglottal reflection coefficient (=1 by default),
% A_ss: system matrix for the state space WRA model,
% Gamma_ss: input matrix for the state space WRA model.
%
% References:
% [1] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation, 1st
%     editio. National Center for Voice and Speech, 2006. 
% [2] G. E. Galindo, Bayesian Estimation of a Subject-Specific Model of
%     Voice Production for the Clinical Assessment of Vocal Function, Ph.D.
%     dissertation, Universidad Técnica Federico Santa María, Valparaíso,
%     2017.
%
% Coded by Gabriel Alzamendi, January 2020.
function varargout = getStateSpaceWRAModel(VTobj,varargin)
    re_n = 1;
    if (nargin == 2)
      re_n = varargin{1};
    end

    % Definition of simulation parameters
    rho = VTobj.RHO_AIR; % [kg m^-3] Density of the air
    c = VTobj.C_AIR; % [m/s] speed of sound
    fs = VTobj.fs;
    Delta_z = VTobj.Delta_z; % c/(2*fs); % lenght of each tube section [m]
    N_AreaSection=VTobj.N_AreaSection;
    
    % Attenuations factors
    A_att  = diag(1 - 2*(3.8e-3./sqrt(VTobj.AreaFunction))*Delta_z); 
    
    % Reflections coefficients
    r_coef = (VTobj.AreaFunction(1:N_AreaSection-1) - VTobj.AreaFunction(2:N_AreaSection)) ./ ...
             (VTobj.AreaFunction(1:N_AreaSection-1) + VTobj.AreaFunction(2:N_AreaSection));
         
    % Radiation parameters
    I = (2*fs/c)*(8/3)*sqrt(VTobj.AreaFunction(end)/pi^3);
    R = 128/(9*pi^2);
    b1 = + I - R + R*I;
    b2 = + I + R + R*I;
    c1 = + I - R - R*I;
    c2 = - I - R + R*I;
    lambda1 = (b2+c2)*A_att(end)/b2;
    lambda2 = (-b1+c1)*A_att(end)/b2;
    lambda3 = b1/b2;
    gamma1 = c2/b2*A_att(end);
    gamma2 = c1/b2*A_att(end);
    gamma3 = b1/b2;
    
    % Input acoustic impedance
    Z_Ug = rho*c/VTobj.AreaFunction(1);
    
    % Model matrices for the state space representation of WRA
    Mat1_aux=zeros(N_AreaSection,N_AreaSection); Mat1_aux(end,end)=1;
    Mat2_aux=zeros(N_AreaSection,1); Mat2_aux(end)=1;
    
%     A_mat = zeros(2*N_AreaSection+2,2*N_AreaSection+2);
    A_mat = [gamma3*Mat1_aux+diag(1-r_coef,1)*A_att, ...
             diag([r_coef;gamma1])*A_att,...
             gamma2*Mat2_aux, 0*Mat2_aux; ... % End row block 1
             diag([re_n;-r_coef(1:end)])*A_att,...
             diag(1+r_coef(1:end),-1)*A_att,...
             0*Mat2_aux, 0*Mat2_aux; ...  % End row block 2
             0*Mat2_aux',Mat2_aux.',0,0; ...  % End row block 3
             0*Mat2_aux',lambda1*Mat2_aux.',lambda2,lambda3];  % End row block 4
    SelVect = ones(2*N_AreaSection+2,1);
    SelVect([N_AreaSection 2*N_AreaSection+1 2*N_AreaSection+2]) = 0*SelVect([N_AreaSection 2*N_AreaSection+1 2*N_AreaSection+2]); % Work!
    SelMat1 = diag(SelVect); SelMat2 = eye(2*N_AreaSection+2) - diag(SelVect);
    
    Gamma_mat = zeros(2*N_AreaSection+2,1); Gamma_mat(N_AreaSection+1) = Z_Ug;
    Gamma_mat = A_mat*Gamma_mat+Gamma_mat;
    Gamma_mat=sparse(Gamma_mat);
    
    A_mat = (SelMat1*A_mat^2+SelMat2*A_mat);
    A_mat = sparse(A_mat);

    %% Storing resulting matrices
    VTobj.A_ss = A_mat;
    VTobj.Gamma_ss = Gamma_mat;
    VTobj.SSWRAvarOK = true;
    
    if nargout == 2
      varargout{1} =  A_mat;
      varargout{2} =  Gamma_mat;
    elseif (nargout==1)||(nargout>2)
      error('It is requested an incorrect number of output varaibles!')  
    end
end
