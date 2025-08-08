function x = modelBCM(x,Pe,Ps,constants)
  
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
      
      Delta_t = constants.Delta_t; % time step [s]

    % Variables
      x_u = x(1); % Position of the upper mass [m]
      x_l = x(2); % Position of the lower mass [m]
      x_b = x(3); % Position of the body mass [m]
      
      v_u = x(4); % Velocity of the upper mass [m s^{-1}]
      v_l = x(5); % Velocity of the lower mass [m s^{-1}]
      v_b = x(6); % Velocity of the body mass [m s^{-1}]
      
  %% Paper Equations
    % Spring Equations
      F_ku = -k_u*(((x_u-x_u0)-(x_b-x_b0))+ eta_u*((x_u-x_u0)-(x_b-x_b0))^3); % Spring force for upper mass [N]
      F_kl = -k_l*(((x_l-x_l0)-(x_b-x_b0))+ eta_l*((x_l-x_l0)-(x_b-x_b0))^3); % Spring force for lower mass [N]
      F_kb = -k_b*((x_b-x_b0)+ eta_b*(x_b-x_b0)^3);                           % Spring force for body mass [N]

      F_kc = -k_c*((x_l-x_l0)-(x_u-x_u0));                                    % Coupling Spring force [N]

      if (x_u <= x_uCol); F_kuCol = -h_uCol*((x_u-x_uCol) + eta_uCol*(x_u-x_uCol)^3); else F_kuCol = 0; end % Collision force for the upper mass [N] 
      if (x_l <= x_lCol); F_klCol = -h_lCol*((x_l-x_lCol) + eta_lCol*(x_l-x_lCol)^3); else F_klCol = 0; end % Collision force for the lower mass [N]

    % Damping Equations
      d_u = 2*zeta_u*sqrt(m_u*k_u);   % Damping coefficient for upper mass [kg s^{-1}]
      d_l = 2*zeta_l*sqrt(m_l*k_l);   % Damping coefficient for lower mass [kg s^{-1}]
      d_b = 2*zeta_b*sqrt(m_b*k_b);   % Damping coefficient for body mass [kg s^{-1}]

      F_du = -d_u*(v_u-v_b);  % Damping force of upper mass [N]
      F_dl = -d_l*(v_l-v_b);  % Damping force of lower mass [N]
      F_db = -d_b*(v_b);      % Damping force of body mass [N]
      
      if (x_u <= 0); F_duCol = -2*zeta_uCol*sqrt(m_u*k_u)*(v_u-v_b); else F_duCol = 0; end
      if (x_l <= 0); F_dlCol = -2*zeta_lCol*sqrt(m_l*k_l)*(v_l-v_b); else F_dlCol = 0; end

    % Pressure Equations
      au = max(0,2*x_u*Lg); % Upper glottal area [m]
      al = max(0,2*x_l*Lg); % Lower glottal area [m]

      if (au == 0) && (al == 0)     % Both masses colliding
        P_u = 0;                        % Upper mass pressure [N m^{-2}]
        P_l = 0;                        % Lower mass pressure [N m^{-2}]
      elseif au == 0                 % Upper mass colliding
        P_u = 0;                        % Upper mass pressure [N m^{-2}]
        P_l = Ps;                      % Lower mass pressure [N m^{-2}]
      elseif al == 0                 % Lower mass colliding
        P_u = Pe;                      % Upper mass pressure [N m^{-2}]
        P_l = 0;                        % Lower mass pressure [N m^{-2}]
      elseif al > au                % Bernoulli
        P_u = Ps-(Ps-Pe)*(au/al);% Pe;                      % Upper mass pressure [N m^{-2}]
        P_l = Ps-(Ps-Pe)*(au/al)^2;% Lower mass pressure [N m^{-2}]
      elseif al <= au               % Jet
        P_u = Pe;                      % Upper mass pressure [N m^{-2}]
        P_l = Pe;                      % Lower mass pressure [N m^{-2}]
      end
      
      F_eu = P_u*Lg*T_u;  % Upper pressure force [N]
      F_el = P_l*Lg*T_l;  % Lower pressure force [N]

    % Total Force Equations
      F_u = F_ku + F_du - F_kc + F_eu + F_kuCol + F_duCol;         % Total Force in the Upper Mass [N] 
      F_l = F_kl + F_dl + F_kc + F_el + F_klCol + F_dlCol;         % Total Force in the Lower Mass [N]
      F_b = F_kb + F_db - (F_ku + F_du + F_duCol + F_kl + F_dl + F_dlCol);  % Total Force in the Body Mass [N] (Typo in the original paper)
      
    % Acceleration
      a_u = F_u/m_u;
      a_l = F_l/m_l;
      a_b = F_b/m_b;

  %% Output
    if strcmpi(constants.solver,'TTS')
      x(1) = x_u + Delta_t*v_u + 0.5*Delta_t^2*a_u;      % Position of upper mass [m]
      x(2) = x_l + Delta_t*v_l + 0.5*Delta_t^2*a_l;      % Position of lower mass [m]
      x(3) = x_b + Delta_t*v_b + 0.5*Delta_t^2*a_b;      % Position of body mass [m]
      x(4) = v_u + Delta_t*a_u;      % Velocity of upper mass [m s^{-1}]
      x(5) = v_l + Delta_t*a_l;      % Velocity of lower mass [m s^{-1}]
      x(6) = v_b + Delta_t*a_b;      % Velocity of body mass [m s^{-1}]
      x(7) = a_u; 
      x(8) = a_l; 
      x(9) = a_b; 
    else
      x(1) = v_u;      % Position of upper mass [m]
      x(2) = v_l;      % Position of lower mass [m]
      x(3) = v_b;      % Position of body mass [m]
      x(4) = a_u;      % Velocity of upper mass [m s^{-1}]
      x(5) = a_l;      % Velocity of lower mass [m s^{-1}]
      x(6) = a_b;      % Velocity of body mass [m s^{-1}]
    end
end