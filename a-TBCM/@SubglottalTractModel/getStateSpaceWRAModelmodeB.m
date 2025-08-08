%%
% getStateSpaceWRAModelmodeA: Function for computing the internal matrices 
% for the state space representation of the half-time Wave Reflection
% Analogue (WRA) solver taking into account the attenuation effects due to
% wave propagation in the tubulets sections. The WRA model is based on
% Sec. 6.4 in [1], and the state space formulation is based on Galindo's
% ideas [2]. 
%
% Structure: getStateSpaceWRAModelmodeB(SGTObj)
%            getStateSpaceWRAModelmodeB(SGTObj,rs_n)
%            [A_ss, Gamma_ss] = getStateSpaceWRAModelmodeB(...)
%
% where
%
% SGTObj: is an object from VocalTractModel (handle) class,
% rs_n: is the supraglottal reflection coefficient (=1 by default),
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
function varargout = getStateSpaceWRAModelmodeB(SGTObj,varargin)
    rs_n = 1;
    if (nargin == 2)
      if isnumeric(varargin{1})&&(abs(varargin{1})<=1)
        rs_n = varargin{1};
      else
        error('Incorrect ''rs_n'' subglottal reflection coefficient. Correct value mast fullf rs_n need to be |rs_n|<=1.')
      end
    end

    % Definition of simulation parameters
    rho = SGTObj.RHO_AIR; % [kg m^-3] Density of the air
    c = SGTObj.C_AIR; % [m/s] speed of sound
    fs = SGTObj.fs;
    Delta_z = SGTObj.Delta_z; % c/(2*fs); % lenght of each tube section [m]
    N_AreaSection=SGTObj.N_AreaSection;
    
    % Attenuations factors
    A_att  = diag(1 - (11.2e-3./sqrt(SGTObj.AreaFunction))*Delta_z); 
    
    % Reflections coefficients
    r_end = SGTObj.r_end;
    r_coef = (SGTObj.AreaFunction(1:N_AreaSection-1) - SGTObj.AreaFunction(2:N_AreaSection)) ./ ...
             (SGTObj.AreaFunction(1:N_AreaSection-1) + SGTObj.AreaFunction(2:N_AreaSection));
    
    % Input acoustic impedance
    Z_Ug = -rho*c/SGTObj.AreaFunction(1);
    
    % Model matrices for the state space representation of WRA
%     Mat1_aux=zeros(N_AreaSection,N_AreaSection); Mat1_aux(end,end)=1;
%     Mat2_aux=zeros(N_AreaSection,1); Mat2_aux(end)=1;
    
%     A_mat = zeros(2*N_AreaSection+2,2*N_AreaSection+2);
    A_mat = [diag(1-r_coef,1)*A_att, ...
             diag([r_coef;r_end])*A_att; ... % End row block 1
             diag([rs_n;-r_coef(1:end)])*A_att,...
             diag(1+r_coef(1:end),-1)*A_att]; % End row block 2
    SelVect = ones(2*N_AreaSection,1);
%     SelVect([N_AreaSection]) = 0*SelVect([N_AreaSection]); % Work!
    SelMat1 = diag(SelVect); SelMat2 = eye(2*N_AreaSection) - diag(SelVect);
    
    Gamma_mat = zeros(2*N_AreaSection,2); Gamma_mat(N_AreaSection+1,1) = Z_Ug;
    Gamma_mat(N_AreaSection,2) = 1;
    Gamma_mat = A_mat*Gamma_mat+Gamma_mat;
%     Gamma_mat=sparse(Gamma_mat);
    
    A_mat = (SelMat1*A_mat^2+SelMat2*A_mat);
%     A_mat = sparse(A_mat);

    %% Storing resulting matrices
    SGTObj.A_ss = A_mat;
    SGTObj.Gamma_ss = Gamma_mat;
    SGTObj.SSWRAvarOK = true;
    
    if nargout == 2
      varargout{1} =  A_mat;
      varargout{2} =  Gamma_mat;
    elseif (nargout==1)||(nargout>2)
      error('It is requested an incorrect number of output varaibles!')  
    end
end
