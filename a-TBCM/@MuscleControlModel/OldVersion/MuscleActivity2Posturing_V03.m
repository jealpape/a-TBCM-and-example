% G. A. Alzamendi, October 2018
% Code for the static posturing of the Larynx and the vocal folds as a
% result of the muscular activation of LCA, IA, PCA, CT and TA muscles. 
%
% Inputs:
% a_MAL: vector containing the LCA,IA,PCA,CT,TA normalized muscle activity level (ranging between 0 and 1)
%
% Outputs:
% VFPostur: Struct with the positional variables describing the glottal
% process, the cricoid-arytenoid joint (CAJ), and vocal folds strain
% (translational, rotational, adductory and total).

function VFPostur=MuscleActivity2Posturing_V03(a_MAL,LarynxParam,varargin)%,LarynxData) 
    VFPostur=struct;
    %% Physiollogical parameters definition
    % Vocal fold length at rest
    a_MALv7=[a_MAL; zeros(2,1)];
    L0=LarynxParam.L0; % [m] Vocal fold length
    if nargin==2
        ColisionLim=false;
    elseif (nargin==3)&&(strcmp(varargin{1},'ColisionLim'))
        if varargin{2}
            ColisionLim=true;
        else
            ColisionLim=false;
        end
    else
        error('Error on calling ''MuscleActivity2Posturing'' function!')
    end
            
        
    %% Iterative approach for solving vocal fold posturing based on 
    % (Titze and Hunter, 2007; Titze, 2006; Titze and Story, 2002)
    
    %% Initialization procedure
    % Initialization of internal variables
    xi_a=0; % CAJ horizontal displacement
    psi_a=0; % CAJ vertical displacement
    theta_a=0; % CAJ rotation angle, positive for counter-clockwise (adductori) rotation
    kappa=7e-3; % Stiffness for CAJ rotation
    k_x=60; % Stiffness for CAJ horizontal displacement
    k_y=74; % Stiffness for CAJ vertical displacement
    k_t=1500; % [N/m] Translational stiffness
    k_r=0.082; % [N m/rad] Rotational stiffness
    F_musc=zeros(7,1); % Force vector for intrinsic muscles of the larynx
    Torque_boundary=0; % Auxilary torque for restricting adductori displacement of glottal process
%     sigma_interstress=zeros(7,1); % Active internal stress
    E_young_new=zeros(7,1); % Young's modulus
    
    % Direction cosines. Order: LCA,IA,PCA,CT,TA, VF ligament, VF mucosa
    alpha = [LarynxParam.LCA.alpha; LarynxParam.IA.alpha; LarynxParam.PCA.alpha; ...
             LarynxParam.CT.alpha; LarynxParam.TA.alpha; ...
             LarynxParam.Lig.alpha; LarynxParam.Muc.alpha];
    beta = [LarynxParam.LCA.beta; LarynxParam.IA.beta; LarynxParam.PCA.beta; ...
            LarynxParam.CT.beta; LarynxParam.TA.beta; ...
            LarynxParam.Lig.beta; LarynxParam.Muc.beta];
    % Directional moment arms. Order: LCA,IA,PCA,CT,TA, VF ligament, VF mucosa
    gamma = [LarynxParam.LCA.gamma; LarynxParam.IA.gamma; LarynxParam.PCA.gamma; ...
             LarynxParam.CT.gamma; LarynxParam.TA.gamma; ...
             LarynxParam.Lig.gamma; LarynxParam.Muc.gamma];
    % CAJ rotation angle at rest
    theta0_rest=acos((LarynxParam.xbar_02-LarynxParam.x_CAJ)/LarynxParam.R_CA);
    
    % Muscular internal variables. Order: LCA,IA,PCA,CT,TA, VF ligament, VF mucosa
    CroSecArea = [LarynxParam.LCA.CrossArea; LarynxParam.IA.CrossArea; LarynxParam.PCA.CrossArea; ...
                  LarynxParam.CT.CrossArea; LarynxParam.TA.CrossArea; ...
                  LarynxParam.Lig.CrossArea; LarynxParam.Muc.CrossArea]; % Cross Sectional Area
    epsilon1 = [LarynxParam.LCA.epsilon1; LarynxParam.IA.epsilon1; LarynxParam.PCA.epsilon1; ...
                LarynxParam.CT.epsilon1; LarynxParam.TA.epsilon1; ...
                LarynxParam.Lig.epsilon1; LarynxParam.Muc.epsilon1]; % Strain at zero stress
    epsilon2 = [LarynxParam.LCA.epsilon2; LarynxParam.IA.epsilon2; LarynxParam.PCA.epsilon2; ...
                LarynxParam.CT.epsilon2; LarynxParam.TA.epsilon2; ...
                LarynxParam.Lig.epsilon2; LarynxParam.Muc.epsilon2]; % Strain at exponential stress
    B = [LarynxParam.LCA.B; LarynxParam.IA.B; LarynxParam.PCA.B; ...
         LarynxParam.CT.B; LarynxParam.TA.B; ...
         LarynxParam.Lig.B; LarynxParam.Muc.B]; % Exponential stress constant
    sigma0 = [LarynxParam.LCA.sigma0; LarynxParam.IA.sigma0; LarynxParam.PCA.sigma0; ...
              LarynxParam.CT.sigma0; LarynxParam.TA.sigma0; ...
              LarynxParam.Lig.sigma0; LarynxParam.Muc.sigma0]; % [Pa] Passive stress at resting length
    sigma2 = [LarynxParam.LCA.sigma2; LarynxParam.IA.sigma2; LarynxParam.PCA.sigma2; ...
              LarynxParam.CT.sigma2; LarynxParam.TA.sigma2; ...
              LarynxParam.Lig.sigma2; LarynxParam.Muc.sigma2]; % [Pa] Scaling of exponential stress
    b = [LarynxParam.LCA.b; LarynxParam.IA.b; LarynxParam.PCA.b; ...
         LarynxParam.CT.b; LarynxParam.TA.b; ...
         0; 0]; % coeff. for active stress-strain
    epsilonm = [LarynxParam.LCA.epsilonm; LarynxParam.IA.epsilonm; LarynxParam.PCA.epsilonm; ...
                LarynxParam.CT.epsilonm; LarynxParam.TA.epsilonm; ...
                0;0]; % coeff. for active stress-strain
    sigmam = [LarynxParam.LCA.sigmam; LarynxParam.IA.sigmam; LarynxParam.PCA.sigmam; ...
              LarynxParam.CT.sigmam; LarynxParam.TA.sigmam; ...
              0; 0]; % [Pa] Scaling of exponential stress
    
    % Strain variable initialization
    epsilon_translational=0; % Translational strain
    epsilon_rotation=0; % Rotational strain
    epsilon_adductory=0; % Adductory strain
    epsilon_total=epsilon_translational+epsilon_rotation+epsilon_adductory; % Total strain
    epsilon_CT=-L0/LarynxParam.CT.Length*...
                (epsilon_rotation*LarynxParam.CTJ.w_CT/LarynxParam.CTJ.h_TA...
                +epsilon_translational/LarynxParam.CTJ.cos_phi); % CT strain
    
    % Initialization of method hyperparameters
    iter_Max = 100;
    zerovec7 = zeros(7,1);
    
    %% Iterative computation of vocal fold posturing
    for cont_iter=1:iter_Max
        if cont_iter<=1
            lambda_conv=1;
        elseif cont_iter<=5
            lambda_conv=0.9;
        else
            lambda_conv=0.7;
        end
        % Storing total strain variable computed in the previous iteration
        epsilon_total_prev=epsilon_total;
        % ''Effective'' stiffnesses
        if abs(xi_a)<3.0e-3 % [m]
            k_x_new=60*(1+2e5*xi_a^2); % [N/m] translational x displacements stiffness
        else
            k_x_new=60*(1+2e5*xi_a^2+2e13*(abs(xi_a)-3e-3)^4); % [N/m] translational x displacements stiffness
%             k_x_new=60*(1+2e5*xi_a^2+2e10*(abs(xi_a)-3e-3)^4); % [N/m] translational x displacements stiffness
        end
        k_x=lambda_conv*k_x+(1-lambda_conv)*k_x_new;
        
        if abs(psi_a)<3.0e-3 % [m]
            k_y_new=200*(1+2e5*psi_a^2); % [N/m] translational y displacements stiffness
        else
            k_y_new=200*(1+2e5*psi_a^2+2e13*(abs(psi_a)-3e-3)^4); % [N/m] translational y displacements stiffness
%             k_y_new=200*(1+2e5*psi_a^2+2e10*(abs(psi_a)-3e-3)^4); % [N/m] translational y displacements stiffness
        end
        k_y=lambda_conv*k_y+(1-lambda_conv)*k_y_new;
        
%         if abs(theta_a)<0.3 % [rad]
            kappa_new=5e-3*(1+10*theta_a^2); % [N m/rad] rotation stiffness
%         else
%             kappa_new=5e-3*(1+10*theta_a^2+1e4*(abs(theta_a)-0.3)^4); % [N m/rad] rotation stiffness
% %             kappa_new=5e-3*(1+10*theta_a^2+1e2*(abs(theta_a)-0.3)^4); % [N m/rad] rotation stiffness
%         end
        kappa=lambda_conv*kappa+(1-lambda_conv)*kappa_new;
        
        % Translational stiffness
        if abs(L0*epsilon_translational)<=1.5e-3
            k_t_new=500*(1+1e6*L0^2*epsilon_translational^2)*...
                        (1+30/pi*L0/LarynxParam.CTJ.h_TA*abs(epsilon_rotation)); % [N/m]
        else
            k_t_new=500*(1+1e6*L0^2*epsilon_translational^2+2e14*(abs(L0*epsilon_translational)-1.5e-3)^4)*...
                        (1+30/pi*L0/LarynxParam.CTJ.h_TA*abs(epsilon_rotation)); % [N/m]
%             k_t_new=500*(1+1e6*L0^2*epsilon_translational^2+2e10*(abs(L0*epsilon_translational)-1.5e-3)^4)*...
%                         (1+3/pi*L0/LarynxParam.CTJ.h_TA*abs(epsilon_rotation)); % [N/m]
        end
        k_t=lambda_conv*k_t+(1-lambda_conv)*k_t_new;
        
        % Rotational stiffness
        theta_CTJ=L0/LarynxParam.CTJ.h_TA*epsilon_rotation; % Rotation angle of the CTJ
        if abs(theta_CTJ)<=0.2
            k_r_new= 0.05*(1+1500*L0*abs(epsilon_translational))*...
                        (1+40*theta_CTJ^2);% [N m/rad]
        else
            k_r_new= 0.05*(1+1500*L0*abs(epsilon_translational))*...
                        (1+40*theta_CTJ^2+2e5*(abs(theta_CTJ)-0.2)^4);% [N m/rad]
% %             k_r_new= 0.05*(1+1500*L0*abs(epsilon_translational))*...
% %                         (1+40*theta_CTJ^2+2e3*(abs(theta_CTJ)-0.2)^4);% [N m/rad]
        end
        k_r=lambda_conv*k_r+(1-lambda_conv)*k_r_new;

        % Young's modulus computation for LCA, IA, PCA, CT, TA, VF ligament, VF mucosa
        epsilony=[zeros(3,1); epsilon_CT; epsilon_total; zeros(2,1)*epsilon_total];
        % It is applied the rules for computing Young's modulus
        % reported in (Titze, 2006), Eqs. (2.62) and (2.63) -Pag. 77-.
        Ind_epsilon = epsilony<=epsilon2;
        % For epsilony <= epsilon2
        E_young_new(Ind_epsilon)=-sigma0(Ind_epsilon)./epsilon1(Ind_epsilon); % [Pa]
        % For epsilony > epsilon2
        E_young_new(~Ind_epsilon)=-sigma0(~Ind_epsilon)./epsilon1(~Ind_epsilon) + ...
                    B(~Ind_epsilon).*sigma2(~Ind_epsilon).* ...
                    (exp(B(~Ind_epsilon).*(epsilony(~Ind_epsilon)-epsilon2(~Ind_epsilon)))-1); % [Pa]
        E_young = E_young_new;
        % Active internal stress estimation for LCA,IA,PCA,CT,TA muscles
        sigma_interstress_new=sigmam.*a_MALv7.*max([zerovec7,1-b.*(epsilony-epsilonm).^2],[],2); % Eq. (12)
%         sigma_interstress=lambda_conv*sigma_interstress+(1-lambda_conv)*sigma_interstress_new;
        sigma_interstress = sigma_interstress_new;
        % Force computation corresponding to LCA,IA,PCA,CT,TA muscles
        F_musc_new=CroSecArea.*(sigma_interstress+E_young.*epsilony); % Eq. (11)
        F_musc = lambda_conv*F_musc+(1-lambda_conv)*F_musc_new;
        
        % Motion of the arytenoid cartilage around CAJ
        xi_a_new=1/k_x*sum(alpha.*F_musc); % [m] Motion on x direction
        xi_a=lambda_conv*xi_a+(1-lambda_conv)*xi_a_new; % [m] Motion on x direction
        
        psi_a_new=1/k_y*sum(beta.*F_musc); % [m] Motion on y direction
        psi_a=lambda_conv*psi_a+(1-lambda_conv)*psi_a_new; % [m] Motion on y direction
        theta_a_new=1/kappa*(sum(gamma.*F_musc)+Torque_boundary); % [rad] Rotation around CAJ
        theta_a_proof=lambda_conv*theta_a+(1-lambda_conv)*theta_a_new; % [rad] Rotation around CAJ
        
        % Checking boundary condition of glottal closure
        if ColisionLim
            xi_02=LarynxParam.x_CAJ-(LarynxParam.x_CAJ-LarynxParam.xbar_02)*cos(theta_a_proof)...
                                    +LarynxParam.y_CAJ*sin(theta_a_proof)+xi_a;
            if xi_02<0
                ang_aux=acos((LarynxParam.x_CAJ+xi_a)/LarynxParam.R_CA);
                theta_a_closure=pi-theta0_rest-ang_aux;
                theta_a_closure=0.99*theta_a_closure;
                Torque_boundary=kappa*theta_a_closure-sum(gamma.*F_musc);
                theta_a=lambda_conv*theta_a+(1-lambda_conv)*theta_a_closure;
            else
                theta_a=theta_a_proof;
            end
        else
            theta_a=theta_a_proof;
        end
        
        % Vocal fold strains
        epsilon_translational_new= 1/(k_t*L0)*(F_musc(4)*LarynxParam.CTJ.cos_phi-sum(F_musc(5:7))); % Translational strain of the CTJ
        epsilon_translational=lambda_conv*epsilon_translational+(1-lambda_conv)*epsilon_translational_new;
        
        epsilon_rotation_new=LarynxParam.CTJ.h_TA/(k_r*L0)*...
                        (LarynxParam.CTJ.w_CT*F_musc(4)-LarynxParam.CTJ.h_TA*sum(F_musc(5:7))); % Rotational strain of the CTJ
        epsilon_rotation=lambda_conv*epsilon_rotation+(1-lambda_conv)*epsilon_rotation_new;
        
        epsilon_adductory_new=-1/L0*(LarynxParam.y_CAJ*(1-cos(theta_a))-(LarynxParam.x_CAJ-LarynxParam.xbar_02)*sin(theta_a)+psi_a);
        epsilon_adductory=lambda_conv*epsilon_adductory+(1-lambda_conv)*epsilon_adductory_new;
        
        epsilon_total=epsilon_translational+epsilon_rotation+epsilon_adductory; % Total Vocal Fold strain

        epsilon_CT_new=-L0/LarynxParam.CT.Length*...
                        (epsilon_rotation*LarynxParam.CTJ.w_CT/LarynxParam.CTJ.h_TA...
                        +epsilon_translational/LarynxParam.CTJ.cos_phi);
%             epsilon_CT=epsilon_CT_new;
        epsilon_CT=lambda_conv*epsilon_CT+(1-lambda_conv)*epsilon_CT_new;



        delta_eps=abs(epsilon_total-epsilon_total_prev);
%         fprintf('Iter=% 2u,\t epsilon_tot=%1.1e,\t epsilon_tot_prev=%1.1e,\t delta_eps=%1.1e\n',cont_iter,epsilon_total,epsilon_total_prev,delta_eps);
        if cont_iter>1&&delta_eps<1e-6%&&xi_02>0
            break;
        end
    end
    VFPostur.epsilon_translational=epsilon_translational;
    VFPostur.epsilon_rotation=epsilon_rotation;
    VFPostur.epsilon_adductory=epsilon_adductory;
    VFPostur.epsilon_total=epsilon_total;
    VFPostur.epsilon_CT=epsilon_CT;
    VFPostur.xi_a=xi_a;
    VFPostur.psi_a=psi_a;
    VFPostur.theta_a=theta_a;
    VFPostur.xi_02=LarynxParam.x_CAJ-(LarynxParam.x_CAJ-LarynxParam.xbar_02)*cos(theta_a)...
                    +LarynxParam.y_CAJ*sin(theta_a)+xi_a;
    VFPostur.psi_02=LarynxParam.y_CAJ-LarynxParam.y_CAJ*cos(theta_a)...
                    -(LarynxParam.x_CAJ-LarynxParam.xbar_02)*sin(theta_a)...
                    +psi_a-L0*(epsilon_rotation+epsilon_translational);
end
