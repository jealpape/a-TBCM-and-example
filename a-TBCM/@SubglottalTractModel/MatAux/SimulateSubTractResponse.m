function [Resp,B1,F1,B_k,F_k] = SimulateSubTractResponse(fsampling,AreaFun,Ug,rs,B_in,F_in)
%     if rem(length(AreaFun),2)==1
%         AreaFun=AreaFun(1:end-1);
%     end
    % Definition of simulation parameters
    rho = 1.146; % [kg m^-3] Density of the air
    c = 350; % [m/s] sound velocity
    fs = fsampling;
    Delta_z = c/(2*fs); % lenght of each tube section [m]
    L_AreaFun=length(AreaFun);
    % Attenuations factors 
    A_att  = 1 - (11.2e-3./sqrt(AreaFun))*Delta_z; A_att=A_att';
    % Reflections coefficients
    r_end = -0.5;
    r_coef = (AreaFun(1:L_AreaFun-1) - AreaFun(2:L_AreaFun))./(AreaFun(1:L_AreaFun-1) + AreaFun(2:L_AreaFun));
    r_coef=r_coef';
    K_Ug = -rho*c/AreaFun(1);
    % Initialization
    Resp = zeros(size(Ug));
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
        F_k(1) = K_Ug*Ug(contSim) + rs*A_att(1)*B_k(1);
        F_k(3:2:L_AreaFun) = (1+r_coef(2:2:L_AreaFun-1)).*A_att(2:2:L_AreaFun-1).*F_k(2:2:L_AreaFun-1) ...
                                - r_coef(2:2:L_AreaFun-1).*A_att(3:2:L_AreaFun).*B_k(3:2:L_AreaFun);
%         for cont_pos = 1:2:L_AreaFun
%             if cont_pos == 1
%                 F_k(cont_pos) = K_Ug*Ug(contSim) + rs*A_att(1)*B_k(cont_pos);
%             else
%                 F_k(cont_pos) = (1+r_coef(cont_pos-1))*A_att(cont_pos-1)*F_k(cont_pos-1) ...
%                                 - r_coef(cont_pos-1)*A_att(cont_pos)*B_k(cont_pos);
%             end
%         end
        % Backward pressure wave at even junctions and half sample computation
        B_k(2:2:L_AreaFun-1) = r_coef(2:2:L_AreaFun-1).*A_att(2:2:L_AreaFun-1).*F_k(2:2:L_AreaFun-1) ...
                            + (1-r_coef(2:2:L_AreaFun-1)).*A_att(3:2:L_AreaFun).*B_k(3:2:L_AreaFun);
        if rem(L_AreaFun,2)==0
            B_k(L_AreaFun) = r_end*A_att(L_AreaFun)*F_k(L_AreaFun);
        end
%         for cont_pos = 2:2:L_AreaFun
%             if cont_pos == L_AreaFun
%                 B_k(cont_pos) = r_end*A_att(cont_pos)*F_k(cont_pos);
%             else
%                 B_k(cont_pos) = r_coef(cont_pos)*A_att(cont_pos)*F_k(cont_pos) ...
%                             + (1-r_coef(cont_pos))*A_att(cont_pos+1)*B_k(cont_pos+1);
%             end
%         end
        
        % Forward pressure wave at even junctions and integer sample computation
        F_k(2:2:L_AreaFun) = (1+r_coef(1:2:L_AreaFun-1)).*A_att(1:2:L_AreaFun-1).*F_k(1:2:L_AreaFun-1) ...
                                - r_coef(1:2:L_AreaFun-1).*A_att(2:2:L_AreaFun).*B_k(2:2:L_AreaFun);
%         for cont_pos = 2:2:L_AreaFun
%             F_k(cont_pos) = (1+r_coef(cont_pos-1))*A_att(cont_pos-1)*F_k(cont_pos-1) ...
%                                 - r_coef(cont_pos-1)*A_att(cont_pos)*B_k(cont_pos);
%         end
        % Backward pressure wave at odd junctions and integer sample computation
        B_k(1:2:L_AreaFun-1) = r_coef(1:2:L_AreaFun-1).*A_att(1:2:L_AreaFun-1).*F_k(1:2:L_AreaFun-1) ...
                            + (1-r_coef(1:2:L_AreaFun-1)).*A_att(2:2:L_AreaFun).*B_k(2:2:L_AreaFun);
        if rem(L_AreaFun,2)==1
            B_k(L_AreaFun) = r_end*A_att(L_AreaFun)*F_k(L_AreaFun);
        end
%         for cont_pos = 1:2:L_AreaFun
%             if cont_pos == L_AreaFun
%                 B_k(cont_pos) = r_end*A_att(cont_pos)*F_k(cont_pos);
%             else
%                 B_k(cont_pos) = r_coef(cont_pos)*A_att(cont_pos)*F_k(cont_pos) ...
%                             + (1-r_coef(cont_pos))*A_att(cont_pos+1)*B_k(cont_pos+1);
%             end
%         end
        
        %% Output: Lump-end forward pressure wave
        B1(contSim) = B_k(1);
        F1(contSim) = F_k(1);
        Resp(contSim) = F_k(end);
        
    end

end