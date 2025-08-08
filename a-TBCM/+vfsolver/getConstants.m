function constants = getConstants(Config,light)
  if ~exist('light','var')
    light = false;
    % Light version for the particle filtering process
  end
  switch upper(Config.constants)
    case 'GALINDO2016'
      constants.enviroment.rho = 1.14;        % Density of the air 1.14 [kg m^-3]
      constants.enviroment.c = 350;           % Sound Velocity in [m s]
      constants.enviroment.mu = 18.36922e-6;  % Air Viscosity [Pa s]
      constants.enviroment.Re = 1600;         % Reynolds Number Threshold - Titze 2006 "The myoelastic theory of phonation" P.263

      constants.simulation.fs = Config.fs; % Samplign Frequency
      constants.simulation.T = Config.SimTime; % Simulation Time
      constants.simulation.Delta_t = 1/constants.simulation.fs; % Simulation step
 
      if light
        constants.simulation.K = 1;
      else
        constants.simulation.K = round(Config.fs*Config.SimTime); % Simulation Steps
      end
      
      % Parameters
      if light
        constants.parameter.ps = Config.ps(1);
        constants.parameter.pe = Config.pe(1);
      else
        constants.parameter.ps = distributeInput(Config.ps,constants.simulation.K);
        constants.parameter.pe = distributeInput(Config.pe,constants.simulation.K);
      
        % Introductory step for pressures of 10[ms]
        if (constants.simulation.K >= 0.02*constants.simulation.fs)
          Hann_wind = hann(ceil(0.02*constants.simulation.fs))';
          Hann_wind = Hann_wind(1:ceil(end/2));
          constants.parameter.ps(1:length(Hann_wind)) = constants.parameter.ps(1:length(Hann_wind)).*Hann_wind;
          constants.parameter.pe(1:length(Hann_wind)) = constants.parameter.pe(1:length(Hann_wind)).*Hann_wind; 
        end
      end
      
      % Arytenoid Configuration
      if light
        if strcmpi(Config.kModel,'BCM')
          Rot = 0;
          Dis = 0;
        else
          Rot = Config.rotation(1);
          Dis = Config.displacement(1);
        end
      else
        if strcmpi(Config.kModel,'BCM')
          Rot = distributeInput(0,constants.simulation.K);
          Dis = distributeInput(0,constants.simulation.K);
        else
          Rot = distributeInput(Config.rotation,constants.simulation.K);
          Dis = distributeInput(Config.displacement,constants.simulation.K);
        end
      end
      
      
      if strcmpi(Config.sex,'male') % kim et al - AORL - 2004 - Comparison of Humman, Canine, And ovine laryngeal dimensions
        La =  17.51e-3;
      else
        La =  12.65e-3;
      end
      gam = 0.2; %Ratio open to the flow channel
      
      La = La*Config.proportion; % Proportion adjustement

      if light
        Act = Config.Act(1);
        Ata = Config.Ata(1);
        Alc = Config.Alc(1);
      else
        Act = distributeInput(Config.Act,constants.simulation.K);
        Ata = distributeInput(Config.Ata,constants.simulation.K);
        Alc = distributeInput(Config.Alc,constants.simulation.K);
      end
        
      % Model parameter
        constants.model.fixed.La = La;                   % Aritenoid Length [m]
        constants.model.fixed.ratio = gam;               % Aritenoid ratio for the exposed section [-]
        for k = 1:constants.simulation.K
          constants.model.variable(k).Act = Act(k);
          constants.model.variable(k).Ata = Ata(k);
          constants.model.variable(k).Alc = Alc(k);

          constants.model.variable(k).rotation = Rot(k);            % Aritenoid rotation [º]
          constants.model.variable(k).displacement = Dis(k);        % Aritenoid displacement [m]
          
          rules = vfsolver.getRules(Act(k),Ata(k),Alc(k),Config.sex,Config.proportion);

          constants.model.variable(k).Lg  = rules.L;             % Length of the glottis [m]
          constants.model.variable(k).T_u = rules.T2;            % Thickness of the upper portion of the glottis [m]
          constants.model.variable(k).T_l = rules.T1;            % Thickness of the lower portion of the glottis [m]
          constants.model.variable(k).k_u = rules.k2;           % Upper spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).k_l = rules.k1;           % Lower spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).k_b = rules.kb;           % Body spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).k_c = rules.kc;           % Coupling spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).h_uCol = 3*rules.k2;       % Upper collision spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).h_lCol = 3*rules.k1;       % Lower collision spring constant [kg s^{-2}] or [N m^{-1}]
          constants.model.variable(k).eta_u = 100*1e4;          % Non-linear upper spring constant [m^{-2}]
          constants.model.variable(k).eta_l = 100*1e4;          % Non-linear lower spring constant [m^{-2}]
          constants.model.variable(k).eta_b = 100*1e4;          % Non-linear body spring constant [m^{-2}]
          constants.model.variable(k).eta_uCol = 500*1e4;        % Non-linear upper collision spring constant [m^{-2}]
          constants.model.variable(k).eta_lCol = 500*1e4;        % Non-linear lower collision spring constant [m^{-2}]
          constants.model.variable(k).zeta_u = 0.6;             % Upper damping ratio [-]
          constants.model.variable(k).zeta_l = 0.1;             % Lower damping ratio [-]
          constants.model.variable(k).zeta_b = 0.1;             % Body damping ratio [-]
          constants.model.variable(k).zeta_uCol = 1;            % Upper additional damping ratio during collision [-]
          constants.model.variable(k).zeta_lCol = 1;            % Lower additional damping ratio during collision [-]
          constants.model.variable(k).m_u = rules.m2;           % Mass of the upper block [kg]
          constants.model.variable(k).m_l = rules.m1;           % Mass of the lower block [kg]
          constants.model.variable(k).m_b = rules.mb;           % Mass of the body block [kg]
          constants.model.variable(k).x_u0 = rules.x_02;         % Rest position of the upper mass [m]
          constants.model.variable(k).x_l0 = rules.x_01;         % Rest position of the lower mass [m]
          constants.model.variable(k).x_b0 = 0.3*1e-2;           % Rest position of the body mass [m]
          constants.model.variable(k).x_uCol = 0;                % Collision medial displacement for the upper mass [m]
          constants.model.variable(k).x_lCol = 0;                % Collision medial displacement for the lower mass [m]

          constants.model.variable(k).PGD = 2*(Dis(k) + La*sind(Rot(k))); % Model PGD [m]
          constants.model.variable(k).PGO = gam*La*cosd(Rot(k))*(constants.model.variable(k).PGD - gam*La*sind(Rot(k)));  % Model PGO [m^2]
          constants.model.variable(k).MGO = (Dis(k) + La*sind(Rot(k)))*constants.model.variable(k).Lg; % Model MGO [m^2];
          constants.model.variable(k).x_du0 = constants.model.variable(k).PGD/2; % upper triangular displacement for the upper mass
          constants.model.variable(k).x_dl0 = constants.model.variable(k).PGD/2; % lower triangular displacement for the upper mass
          
          if (constants.model.variable(k).Lg <= 0) || (constants.model.variable(k).T_u <= 0) || (constants.model.variable(k).T_l <= 0)
            error('vfsolver:getConstant:Dude, what a fuck!');
          end      
        end
      
      % WRA
        constants.wra.delta_z = constants.enviroment.c/(2*constants.simulation.fs);  % lenght of each tube section [m]
        
        if ~light
          c = constants.enviroment.c;
          fs = constants.simulation.fs;
          Delta_z = c/(2*fs); % lenght of each tube section [m]
        end
        
        % Supraglottal
          TRACT = vfsolver.getTract(Config.supTract); % Glottal tract

          TRACT(1) = TRACT(1)*Config.compression; % compression of the supreaglottal tube <<<<<<----------- Check this.
          
          if ~light
            N = length(TRACT); % Number of sections in the tract

            DIST = Delta_z:Delta_z:N*Delta_z; % Position of each segment [m]

            % Attenuations
            a_plus  = 1 - (3.8e-3./sqrt(TRACT(2:N)))*Delta_z;
            a_minus = 1 - (3.8e-3./sqrt(TRACT(1:N-1)))*Delta_z;

            % Reflections
            r = (TRACT(1:N-1) - TRACT(2:N))./(TRACT(1:N-1) + TRACT(2:N));

            % Radiation
            I = (2*fs/c)*(8/3)*sqrt(TRACT(end)/pi^3);
            R = 128/(9*pi);
            b1 = + I - R + R*I;
            b2 = + I + R + R*I;
            c1 = + I - R - R*I;
            c2 = - I - R + R*I;      

            % Matrices
            S_U = [zeros(N-2,1) eye(N-2,N-2); zeros(1,N-1)];
            S_L = [zeros(1,N-1); eye(N-2,N-2) zeros(N-2,1)];
            M_U = [1; zeros(N-2,1)];
            M_L = [zeros(N-2,1); 1];

            A_matrix = [...
              ((eye(N-1) - diag(r))*diag(a_plus)*S_U)   (diag(r)*diag(a_minus)*S_L);...
              (-diag(r)*diag(a_plus)*S_U)               ((eye(N-1) + diag(r))*diag(a_minus)*S_L)];

            B_matrix = [...
              ((eye(N-1) - diag(r))*diag(a_plus)*M_L)   (diag(r)*diag(a_minus)*M_U);...
              (-diag(r)*diag(a_plus)*M_L)               ((eye(N-1) + diag(r))*diag(a_minus)*M_U)];

            B_matrix = (A_matrix + eye(2*(N-1)))*B_matrix;
            A_matrix = A_matrix^2;

            C_matrix = zeros(3,2*(N-1));
            C_matrix(1,end) = (c2/b2)*a_plus(end);
            C_matrix(2,end) = ((c2+b2)/b2)*a_plus(end);
            C_matrix(3,end) = a_plus(end);

            D_matrix = zeros(3,3);
            D_matrix(1,1) = (b1/b2);
            D_matrix(2,2) = (b1/b2);

            D_matrix(1,3) = (c1/b2);
            D_matrix(2,3) = ((c1-b1)/b2);
          end
        
          % Storage
          constants.wra.supraglottal.area = TRACT;% Glottal tract [m^2]
          if ~light
            constants.wra.supraglottal.N = N; % Number of sections in the tract
            constants.wra.supraglottal.dist = DIST; % Position of each segment [m]
            constants.wra.supraglottal.A_matrix = A_matrix; % SSM A matrix
            constants.wra.supraglottal.B_matrix = B_matrix; % SSM B matrix
            constants.wra.supraglottal.C_matrix = C_matrix; % SSM C matrix
            constants.wra.supraglottal.D_matrix = D_matrix; % SSM D matrix (multiplied by Y vector)
          end

        % Subglottal
          TRACT = vfsolver.getTract(Config.subTract); % Glottal tract
          
          if ~light
            Delta_z = constants.enviroment.c/(2*constants.simulation.fs); % lenght of each tube section [m]
            N = length(TRACT); % Number of sections in the tract

            DIST = Delta_z:Delta_z:N*Delta_z; % Position of each segment [m]

            % Attenuations
            a_plus  = 1 - (11.2e-3./sqrt(TRACT(2:N)))*Delta_z;
            a_minus = 1 - (11.2e-3./sqrt(TRACT(1:N-1)))*Delta_z;

            % Reflections
            r = (TRACT(1:N-1) - TRACT(2:N))./(TRACT(1:N-1) + TRACT(2:N));

            % Matrices
            S_U = [zeros(N-2,1) eye(N-2,N-2); zeros(1,N-1)];
            S_L = [zeros(1,N-1); eye(N-2,N-2) zeros(N-2,1)];
            M_U = [1; zeros(N-2,1)];
            M_L = [zeros(N-2,1); 1];

            A_matrix = [...
              ((eye(N-1) - diag(r))*diag(a_plus)*S_U)   (diag(r)*diag(a_minus)*S_L);...
              (-diag(r)*diag(a_plus)*S_U)               ((eye(N-1) + diag(r))*diag(a_minus)*S_L)];

            B_matrix = [...
              ((eye(N-1) - diag(r))*diag(a_plus)*M_L)   (diag(r)*diag(a_minus)*M_U);...
              (-diag(r)*diag(a_plus)*M_L)               ((eye(N-1) + diag(r))*diag(a_minus)*M_U)];

            B_matrix = (A_matrix + eye(2*(N-1)))*B_matrix;
            A_matrix = A_matrix^2;

            % no reflection, no radiation
            C_matrix = zeros(3,2*(N-1));
            D_matrix = zeros(3,3);
          end
        
          % Storage
          constants.wra.subglottal.area = TRACT;% Glottal tract [m^2]
          if ~light
            constants.wra.subglottal.N = N; % Number of sections in the tract
            constants.wra.subglottal.dist = -DIST; % Position of each segment [m]
            constants.wra.subglottal.A_matrix = A_matrix; % SSM A matrix
            constants.wra.subglottal.B_matrix = B_matrix; % SSM B matrix
            constants.wra.subglottal.C_matrix = C_matrix; % SSM C matrix
            constants.wra.subglottal.D_matrix = D_matrix; % SSM D matrix (multiplied by Y vector)
          end

    otherwise
      error('vfsolver:getConstant:constants',['The constant type ''' Config.constants ''' does not exist.']);
  end
end

function out = distributeInput(parameter,K)
  N = length(parameter);
  if N==1
    parameter = parameter*[1 1];
    N = 2;
  end
  
  out = zeros(1,K);
  step = round(linspace(1,K,N));
  
  for n = 1:N-1
    out(step(n):step(n+1)) = linspace(parameter(n),parameter(n+1),step(n+1)-step(n)+1);
  end
end