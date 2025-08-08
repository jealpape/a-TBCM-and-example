function SimulateWRA(fsampling,AreaFun,Ug,re,B_in,F_in)
    % Definition of simulation parameters
    rho = 1.146; % [kg m^-3] Density of the air
    c = 350; % [m/s] sound velocity
    fs = fsampling;
    Delta_z = c/(2*fs); % lenght of each tube section [m]
    L_AreaFun=length(AreaFun);
    % Attenuations factors
    A_att  = 1 - (3.8e-3./sqrt(AreaFun))*Delta_z; A_att = A_att';
    % Reflections coefficients
    r_coef = (AreaFun(1:L_AreaFun-1) - AreaFun(2:L_AreaFun))./(AreaFun(1:L_AreaFun-1) + AreaFun(2:L_AreaFun));
    r_coef = r_coef';
    % Radiation parameters
    I = (2*fs/c)*(8/3)*sqrt(AreaFun(end)/pi^3);
    R = 128/(9*pi);
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
    K_Ug = rho*c/AreaFun(1);
    % Initialization
    Resp = zeros(size(Ug));
    B1 = zeros(size(Ug));
    F1 = zeros(size(Ug));
    N_sim=length(Ug);
    if nargin==4
        B_k = zeros(L_AreaFun,1); % Backward pressure waves
        F_k = zeros(L_AreaFun,1); % Forware pressure waves
    elseif nargin==6
        B_k = B_in;
        F_k = F_in;
    end
            
    F_pass = F_k;
    for contSim = 1:N_sim
        
        F_L_aux = F_pass(end); % Forware pressure at L section and previous iteration
        B_pass=B_k;
        F_pass=F_k;
        
        %% Sound wave propagation
        % Forward pressure wave at odd junctions and half sample computation
        F_k(1) = K_Ug*Ug(contSim) + re*A_att(1)*B_k(1);
        F_k(3:2:end) = (1+r_coef(2:2:end-1)).*A_att(2:2:end-1).*F_k(2:2:end-1) ...
                                - r_coef(2:2:end-1).*A_att(3:2:end).*B_k(3:2:end);
%         for cont_pos = 1:2:L_AreaFun
%             if cont_pos == 1
%                 F_k(cont_pos) = K_Ug*Ug(contSim) + re*A_att(1)*B_k(cont_pos);
%             else
%                 F_k(cont_pos) = (1+r_coef(cont_pos-1))*A_att(cont_pos-1)*F_k(cont_pos-1) ...
%                                 - r_coef(cont_pos-1)*A_att(cont_pos)*B_k(cont_pos);
%             end
%         end
        % Backward pressure wave at even junctions and half sample computation
        B_k(2:2:end-1) = r_coef(2:2:end-1).*A_att(2:2:end-1).*F_k(2:2:end-1) ...
                            + (1-r_coef(2:2:end-1)).*A_att(3:2:end-1).*B_k(3:2:end-1);
        B_k(end) = gamma2*F_L_aux+gamma1*F_k(end)...
                                      + gamma3*B_k(end);
%         for cont_pos = 2:2:L_AreaFun
%             if cont_pos == L_AreaFun
%                 B_k(cont_pos) = gamma2*F_L_aux+gamma1*F_k(cont_pos)...
%                                       + gamma3*B_k(cont_pos);
%             else
%                 B_k(cont_pos) = r_coef(cont_pos)*A_att(cont_pos)*F_k(cont_pos) ...
%                             + (1-r_coef(cont_pos))*A_att(cont_pos+1)*B_k(cont_pos+1);
%             end
%         end
        
        % Forward pressure wave at even junctions and integer sample computation
        F_k(2:2:end) = (1+r_coef(1:2:end)).*A_att(1:2:end-1).*F_k(1:2:end-1) ...
                                - r_coef(1:2:end).*A_att(2:2:end).*B_k(2:2:end);
%         for cont_pos = 2:2:L_AreaFun
%             F_k(cont_pos) = (1+r_coef(cont_pos-1))*A_att(cont_pos-1)*F_k(cont_pos-1) ...
%                                 - r_coef(cont_pos-1)*A_att(cont_pos)*B_k(cont_pos);
%         end
        % Backward pressure wave at odd junctions and integer sample computation
        B_k(1:2:end-1) = r_coef(1:2:end).*A_att(1:2:end-1).*F_k(1:2:end-1) ...
                            + (1-r_coef(1:2:end)).*A_att(2:2:end).*B_k(2:2:end);
%         for cont_pos = 1:2:L_AreaFun
%             B_k(cont_pos) = r_coef(cont_pos)*A_att(cont_pos)*F_k(cont_pos) ...
%                             + (1-r_coef(cont_pos))*A_att(cont_pos+1)*B_k(cont_pos+1);
%         end
        
        %% Output: radiated pressure
        B1(contSim) = B_k(1);
        F1(contSim) = F_k(1);
        if contSim == 1
            Resp(contSim) = lambda2*F_L_aux + lambda1*F_pass(end) + lambda3*0;
        else
            Resp(contSim) = lambda2*F_L_aux + lambda1*F_pass(end) + lambda3*Resp(contSim-1);
        end
        
    end

end