function F = modelTBCMforces(x,Pe,Ps,constants)
  
  %% Inputs
    % Constants
      Lg = constants.Lg; % Length of the glottis [m]
      T_u = constants.T_u; % Thickness of the upper portion of the glottis [m]
      T_l = constants.T_l; % Thickness of the lower portion of the glottis [m]
      
      k_u = constants.k_u;  % Upper spring constant [kg s^{-2}] or [N m^{-1}]
      k_l = constants.k_l;  % Lower spring constant [kg s^{-2}] or [N m^{-1}]
      k_b = constants.k_b;  % Body spring constant [kg s^{-2}] or [N m^{-1}]
      k_c = constants.k_c;  % Coupling spring constant [kg s^{-2}] or [N m^{-1}]
      h_uCol = constants.h_uCol;  % Upper collision spring constant [kg s^{-2}] or [N m^{-1}]
      h_lCol = constants.h_lCol;  % Lower collision spring constant [kg s^{-2}] or [N m^{-1}]

      eta_u = constants.eta_u; % Non-linear upper spring constant [m^{-2}]
      eta_l = constants.eta_l; % Non-linear lower spring constant [m^{-2}]
      eta_b = constants.eta_b; % Non-linear body spring constant [m^{-2}]
      eta_uCol = constants.eta_uCol; % Non-linear upper collision spring constant [m^{-2}]
      eta_lCol = constants.eta_lCol; % Non-linear lower collision spring constant [m^{-2}]
    
      zeta_u = constants.zeta_u; % Upper damping ratio [-]
      zeta_l = constants.zeta_l; % Lower damping ratio [-]
      zeta_b = constants.zeta_b; % Body damping ratio [-]
      zeta_uCol = constants.zeta_uCol; % Upper additional damping ratio during collision [-]
      zeta_lCol = constants.zeta_lCol; % Lower additional damping ratio during collision [-]
      
      m_u = constants.m_u; % Mass of the upper block [kg]
      m_l = constants.m_l; % Mass of the lower block [kg]
      m_b = constants.m_b; % Mass of the body block [kg]

      x_u0 = constants.x_u0; % Rest position of the upper mass [m]
      x_l0 = constants.x_l0; % Rest position of the lower mass [m]
      x_b0 = constants.x_b0; % Rest position of the body mass [m]

      x_uCol = constants.x_uCol; % Collision point of the upper mass [m]
      x_lCol = constants.x_lCol; % Collision point of the lower mass [m]
      
      x_du0 = constants.x_du0; % Posterior displacement of the upper mass [m]
      x_dl0 = constants.x_dl0; % Posterior displacement of the lower mass [m]

      Delta_t = constants.Delta_t; % time step [s]

    % Variables
      x_u = x(1); % Position of the upper mass [m]
      x_l = x(2); % Position of the lower mass [m]
      x_b = x(3); % Position of the body mass [m]
      
      v_u = x(4); % Velocity of the upper mass [m s^{-1}]
      v_l = x(5); % Velocity of the lower mass [m s^{-1}]
      v_b = x(6); % Velocity of the body mass [m s^{-1}]
      
    % Contact point: point where the contours hits the collision plane
      contanctPoint_u = min(Lg,max(0,(x_uCol - x_u)*Lg/x_du0));
      contanctPoint_l = min(Lg,max(0,(x_lCol - x_l)*Lg/x_dl0));

      alpha_u = contanctPoint_u/Lg;
      alpha_l = contanctPoint_l/Lg;      
      
  %% Paper Equations
    % Spring Equations
      F_ku = -k_u*(((x_u-x_u0)-(x_b-x_b0)) + eta_u*((x_u-x_u0)-(x_b-x_b0))^3); % Spring force for upper mass [N]
      F_kl = -k_l*(((x_l-x_l0)-(x_b-x_b0)) + eta_l*((x_l-x_l0)-(x_b-x_b0))^3); % Spring force for lower mass [N]
      F_kb = -k_b*((x_b-x_b0) + eta_b*(x_b-x_b0)^3);                           % Spring force for body mass [N]

      F_kc = -k_c*((x_l-x_l0)-(x_u-x_u0));                                    % Coupling Spring force [N]

    % Collision spring forces [N]
      F_kuCol = -h_uCol*alpha_u*(((0.5*alpha_u*x_du0) + (x_u-x_uCol)) + eta_uCol*(0.25*(alpha_u*x_du0)^3 + (alpha_u*x_du0)^2*(x_u-x_uCol)^1 + 1.5*(alpha_u*x_du0)^1*(x_u-x_uCol)^2 + (x_u-x_uCol)^3)); 
      F_klCol = -h_lCol*alpha_l*(((0.5*alpha_l*x_dl0) + (x_l-x_lCol)) + eta_lCol*(0.25*(alpha_l*x_dl0)^3 + (alpha_l*x_dl0)^2*(x_l-x_lCol)^1 + 1.5*(alpha_l*x_dl0)^1*(x_l-x_lCol)^2 + (x_l-x_lCol)^3));
      
    % Damping Equations
      F_du = -2*zeta_u*sqrt(m_u*k_u)*(v_u-v_b);  % Damping force of upper mass [N]
      F_dl = -2*zeta_l*sqrt(m_l*k_l)*(v_l-v_b);  % Damping force of lower mass [N]
      F_db = -2*zeta_b*sqrt(m_b*k_b)*(v_b);      % Damping force of body mass [N]
      
    % Collision Damping forces [N]
      F_duCol = -2*zeta_uCol*sqrt(m_u*k_u)*(v_u-v_b)*alpha_u;
      F_dlCol = -2*zeta_lCol*sqrt(m_l*k_l)*(v_l-v_b)*alpha_l;
      
    % Pressure Equations
      z_ant = 0;  % Anterior position
      z_pos = Lg; % Posterior position
      
      if (x_du0 - x_dl0) == 0
        z_cro = 0;
      else
        z_cro = min(Lg,max(0,(x_l - x_u)*Lg/(x_du0 - x_dl0))); % Crossing point of the upper mass with the lower mass
      end
      
      if x_du0 == 0
        z_uClo = Lg;
      else
        z_uClo = min(Lg,max(0,(x_uCol - x_u)*Lg/x_du0)); % Crossing point of the upper mass with the collision border
      end
      
      if x_dl0 == 0
        z_lClo = Lg;
      else
        z_lClo = min(Lg,max(0,(x_lCol - x_l)*Lg/x_dl0)); % Crossing point of the lower mass with the collision border  
      end

      % Sorting of the different points
        if (z_cro <= z_lClo) && (z_lClo <= z_uClo)        % C < L < U
          z_sort = [z_ant,z_cro,z_lClo,z_uClo,z_pos];
        elseif (z_cro <= z_uClo) && (z_uClo <= z_lClo)    % C < U < L
          z_sort = [z_ant,z_cro,z_uClo,z_lClo,z_pos];
        elseif (z_lClo <= z_cro) && (z_cro <= z_uClo)     % L < C < U
          z_sort = [z_ant,z_lClo,z_cro,z_uClo,z_pos];
        elseif (z_uClo <= z_cro) && (z_cro <= z_lClo)     % U < C < L
          z_sort = [z_ant,z_uClo,z_cro,z_lClo,z_pos];
        elseif (z_lClo <= z_uClo) && (z_uClo <= z_cro)    % L < U < C
          z_sort = [z_ant,z_lClo,z_uClo,z_cro,z_pos];
        elseif (z_uClo <= z_lClo) && (z_lClo <= z_cro)    % U < L < C
          z_sort = [z_ant,z_uClo,z_lClo,z_cro,z_pos];
        end

      % Initialization of variables 
        R_el = 0; % Resistance of the lower pressure 
        L_ou = 0; % Length of open area of the upper mass
        L_ol = 0; % Length of open area of the lower mass
        a_u = 0;  % Open area of the upper mass
        a_l = 0;  % Open area of the lower mass 
        a_m = 0;  % Minimum area
        per_u = 0; % perimeter of the open portion in the upper mass
        per_l = 0; % perimeter of the open portion on the lower mass
        per_m = 0; % minimum perimeter (for hydraulic diameter calculation)
        per_prev_u = 0; % previous border (displacement) to account for overlaping sections on upper area
        per_prev_l = 0; % previous border (displacement) to account for overlaping sections on lower area
        per_prev_m = 0; % previous border (displacement) to account for overlaping sections on minimal areas

      % Integral Auxiliar Values
        A = x_u;
        B = x_du0/Lg;
        C = x_l;
        D = x_dl0/Lg;

      auxZ2 = z_sort(1);
      
      % Check for each section
        for k = 1:4 
          auxZ1 = auxZ2;
          auxZ2 = z_sort(k+1);
          if auxZ1 ~= auxZ2 % Check if the limits points of the region are different. (just to speed it up a little bit).
            % Border displacements
            x_uz1 = max(0,x_u+auxZ1*(x_du0/Lg)); % upper anterior displacement
            x_uz2 = max(0,x_u+auxZ2*(x_du0/Lg)); % upper posterior displacement
            x_lz1 = max(0,x_l+auxZ1*(x_dl0/Lg)); % lower anterior displacement
            x_lz2 = max(0,x_l+auxZ2*(x_dl0/Lg)); % lower posterior displacement
            
            % Mean Displacements
            x_uz = (x_uz1 + x_uz2)/2;
            x_lz = (x_lz1 + x_lz2)/2;

            % Area
            a_uz = (x_uz1+x_uz2)*abs(auxZ2- auxZ1);
            a_lz = (x_lz1+x_lz2)*abs(auxZ2- auxZ1);
            
            % Perimeter
            per_uz = 2*(x_uz1 + x_uz2) + 2*sqrt((x_uz2-x_uz1)^2 + (auxZ2-auxZ1)^2);
            per_lz = 2*(x_lz1 + x_lz2) + 2*sqrt((x_lz2-x_lz1)^2 + (auxZ2-auxZ1)^2);

            if x_uz > 0     % Is the upper area Open?
              L_ou = L_ou + auxZ2 - auxZ1; % Increase the open length
              a_u = a_u + a_uz;
              per_u = per_u + per_uz - per_prev_u;
              per_prev_u = x_uz1;
            end

            if x_lz > 0  % Is the lower area open?
              L_ol = L_ol + auxZ2 - auxZ1; % Increase the open length
              a_l = a_l + a_lz;
              per_l = per_l + per_lz - per_prev_l;
              per_prev_l = x_lz1;
            end

            % Lower Pressure Resistance
            if x_lz > 0  % Is the lower area open?
              if x_uz > 0     % Is the upper area Open?
                if x_lz > x_uz  % Is convergent?
                  a_m = a_m + a_uz; % The upper area is the minimum area;
                  per_m = per_m + per_uz - per_prev_m;
                  per_prev_m = x_uz1;
                  
                  % Avoid numerical problems on limits Z1 and Z2
                    Iz = (auxZ2 - auxZ1)*(x_uz/x_lz)^2;
%                   if     (B == 0) && (D == 0)
%                     Iz = ((A/C)^2)*(auxZ2 - auxZ1);
%                   elseif (B == 0) && (D ~= 0)
%                     Iz = (-A^2/((C+D*auxZ2)*D))-(-A^2/((C+D*auxZ1)*D));
%                   elseif (B ~= 0) && (D == 0)
%                     Iz = ((A + B*auxZ2)^3/(3*B*C^2))-((A + B*auxZ1)^3/(3*B*C^2));
%                   elseif (B ~= 0) && (D ~= 0)
%                     Iz = ((B/D)^2*auxZ2 - (B*C-A*D)^2/((C+D*auxZ2)*D^3) + (2*B/D^3)*(A*D - B*C)*log(C+D*auxZ2))-((B/D)^2*auxZ1 - (B*C-A*D)^2/((C+D*auxZ1)*D^3) + (2*B/D^3)*(A*D - B*C)*log(C+D*auxZ1)); % [B ~= 0, D ~= 0]
%                   end

                  R_el = R_el + Iz;
                else % No, it's divergent
                  a_m = a_m + a_lz; % The lower area is the minimum area;
                  per_m = per_m + per_lz - per_prev_m;
                  per_prev_m = x_lz1;

                  R_el = R_el + 1*(auxZ2 - auxZ1);              
                end
              else % No, Upper area is closed
                R_el = R_el + 0*(auxZ2 - auxZ1);
              end
            else % No, Lower area is closed
              R_el = R_el + 0*(auxZ2 - auxZ1);
            end
          end
        end
        
      % force due to the pressure
        F_eu = T_u*Pe*L_ou;
        F_el = T_l*Ps*L_ol - T_l*(Ps - Pe)*R_el;

    % Total Force Equations
      F_u = F_ku + F_du - F_kc + F_eu + F_kuCol + F_duCol;         % Total Force in the Upper Mass [N]
      F_l = F_kl + F_dl + F_kc + F_el + F_klCol + F_dlCol;         % Total Force in the Lower Mass [N]
      F_b = F_kb + F_db - (F_ku + F_du + F_duCol + F_kl + F_dl + F_dlCol);  % Total Force in the Body Mass [N] (Typo in the original paper)
 
    % Output forces
    F.u = F_u;
    F.l = F_l;
    F.b = F_b;

    F.ku = F_ku;
    F.kl = F_kl;
    F.kb = F_kb;
    F.kc = F_kc;

    F.du = F_du;
    F.dl = F_dl;
    F.db = F_db;

    F.kuCol = F_kuCol;
    F.klCol = F_klCol;
    
    F.duCol = F_duCol;
    F.dlCol = F_dlCol;
    
    F.eu = F_eu;
    F.el = F_el;    
end