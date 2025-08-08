function out = solver(varargin)
  % System Configuration
  if (nargin == 0)
    config = vfsolver.setConfig;
    warning('vfsolver:default','Using default model configuration, see setConfig.m for more info');
  elseif (nargin > 1)
    config = vfsolver.setConfig(varargin{:});
  else
    config = varargin{1};
  end
  
  constants = vfsolver.getConstants(config);
  
  % Disociation of constants
    % Kinetics
    constKin.variable = rmfield(constants.model.variable,{'PGD','PGO','MGO','displacement','rotation','Act','Ata','Alc'});
    constKin.fixed.solver = config.kSolver;
    constKin.fixed.Delta_t = constants.simulation.Delta_t;

    % Hidraulic
    constFlow.fixed.Ae = constants.wra.supraglottal.area(1);
    constFlow.fixed.As = constants.wra.subglottal.area(1);
    constFlow.fixed.rho = constants.enviroment.rho;
    constFlow.fixed.mu = constants.enviroment.mu;
    constFlow.fixed.c = constants.enviroment.c;
    constFlow.fixed.solver = config.fSolver;
    for k = 1:constants.simulation.K
      constFlow.variable(k).T = constants.model.variable(k).T_u + constants.model.variable(k).T_l;
      constFlow.variable(k).L = constants.model.variable(k).Lg;
      constFlow.variable(k).PGO = constants.model.variable(k).PGO;
      constFlow.variable(k).PGD = constants.model.variable(k).PGD;
    end
    
%     constFlowCell = [fields(constFlow) struct2cell(constFlow)]';
    
    % Noise
    constNoise.fixed.rho = constants.enviroment.rho;
    constNoise.fixed.mu = constants.enviroment.mu;
    constNoise.fixed.Re = constants.enviroment.Re;
    constNoise.fixed.on = strcmpi(config.turbulence,'on');
    for k = 1:constants.simulation.K
      constNoise.variable(k).PGO = constants.model.variable(k).PGO;
    end
%     constNoiseCell = [fields(constNoise) struct2cell(constNoise)]';
    
    
  % State
    % Pressure
      P.sup = zeros(2*(constants.wra.supraglottal.N-1),1);
      P.sup_in = zeros(2,1);
      P.sup_out = zeros(3,1);

      P.sub = zeros(2*(constants.wra.subglottal.N-1),1);
      P.sub_in = zeros(2,1);
      P.sub_out = zeros(3,1);
      
      Pe_plus = 0;
      Ps_plus = 0;
      Pe_minus = 0;
      Ps_minus = 0;
    
    % Kinematics
      x = zeros(9,1);
      x(1) = constKin.variable(1).x_u0;   % Position
      x(2) = constKin.variable(2).x_l0;   % Position
      x(3) = constKin.variable(3).x_b0;   % Position
      x(4) = 0; % Velocity
      x(5) = 0; % Velocity
      x(6) = 0; % Velocity
      x(7) = 0; % Acceleration
      x(8) = 0; % Acceleration
      x(9) = 0; % Acceleration
      
    % Flow
      F.u_g = 0;
      F.u_n = 0;
      F.u = 0;
  
  % Output
    out.x = zeros(9,constants.simulation.K);
    out.Amin = zeros(1,constants.simulation.K);
    out.Pr = zeros(1,constants.simulation.K);
    out.Ug = zeros(1,constants.simulation.K);
    
  % Time step simulation
  for k = 1:constants.simulation.K
    % Pressures
      Pe = constants.parameter.pe(k) + Pe_plus + Pe_minus;
      Ps = constants.parameter.ps(k) + Ps_plus + Ps_minus;
      
    % Kinematics 
      % model solver x = solveModel(Solver,Model,x,Pe,Ps,constKin)
        x = vfsolver.solveModel(config.kSolver,config.kModel,x,Pe,Ps,structcat(constKin.fixed,constKin.variable(k)));
        F = vfsolver.modelTBCMforces(x,Pe,Ps,structcat(constKin.fixed,constKin.variable(k)));
        
      % Model noise
        if sum(config.noise)
          x = x + normrnd(0,config.noise(:));
        end
      % Minimal area and perimeter
        [Amin,Pmin] = vfsolver.getArea(x,config.kModel,constKin.variable(k).Lg,constKin.variable(k).x_du0,constKin.variable(k).x_dl0);
    
    % Flow
      % Glottal Flow Ug = soveFlow(Ps_plus,Pe_minus,Am,varargin)  - ['Ae',0,'As',0,'rho',1.14,'mu',18.36922e-6,'c',350,'L',1,'T',1,'PGO',0,'PGD',0]
        F.u_g = vfsolver.solveFlow(Ps_plus + constants.parameter.ps(k)/2,Pe_minus + constants.parameter.pe(k)/2,Amin,structcat(constFlow.fixed,constFlow.variable(k)));
      % Turbulent Flow Un = solveNoise(Qg,Am,Pm,noise,Ap,rho,mu,Re_c) - ['PGO',0,'rho',1.14,'mu',18.36922e-6,'Re',1200]
        F.u_n = vfsolver.solveNoise(F.u_g,Amin,Pmin,structcat(constNoise.fixed,constNoise.variable(k)));
      % Total flow
        F.u = F.u_g + constNoise.fixed.on*F.u_n;
        
    % Acoustics
      % Flow to acoustic pressures [Pe_plus,Ps_minus] = vfsolver.flow2pressure(Ug,Ps_plus,Pe_minus,Am,varargin)  - ['Ae',0,'As',0,'rho',1.14,'c',350,'PGO',0]
        [Pe_plus,Ps_minus] = vfsolver.flow2pressure(F.u, Ps_plus, Pe_minus, Amin, structcat(constFlow.fixed,constFlow.variable(k)));
        %   P_input = [Pn-;P0+]
        P.sup_in(2) = Pe_plus;
        P.sub_in(2) = Ps_minus;
          
      % Supraglottal Pressures  [P_vector,P_input,P_output] = solveWRA(A_matrix,B_matrix,C_matrix,D_matrix,P_vector,P_input,P_output)
        [P.sup,P.sup_in,P.sup_out] = vfsolver.solveWRA(constants.wra.supraglottal.A_matrix,constants.wra.supraglottal.B_matrix,constants.wra.supraglottal.C_matrix,constants.wra.supraglottal.D_matrix,P.sup,P.sup_in,P.sup_out);
      % Subglottal Pressures  [P_vector,P_input,P_output] = solveWRA(A_matrix,B_matrix,C_matrix,D_matrix,P_vector,P_input,P_output)
        [P.sub,P.sub_in,P.sub_out] = vfsolver.solveWRA(constants.wra.subglottal.A_matrix,constants.wra.subglottal.B_matrix,constants.wra.subglottal.C_matrix,constants.wra.subglottal.D_matrix,P.sub,P.sub_in,P.sub_out);
      % Pressures feedback
        if strcmpi(config.interaction,'on')
          Ps_plus = P.sub(1);
          Pe_minus = P.sup(1);
        else
          Ps_plus = 0;
          Pe_minus = 0; 
        end
        
    % Oral Flow
        F.u_o = (P.sup(end) - P.sup(end/2))/(constFlow.fixed.rho*constFlow.fixed.c/constants.wra.supraglottal.area(end)); % Rude approximation of oral flow
        
    % Storing outputs
      out.x(:,k) = x;
      out.Amin(k) = Amin;
      out.Pr(k) = P.sup_out(2);
      out.Ug(k) = F.u;
      out.Uo(k) = F.u_o;
      
      out.pe(k) = constants.parameter.pe(k);      
      out.pe_plus(k) = Pe_plus;      
      out.pe_minus(k) = Pe_minus;  
      
      out.ps(k) = constants.parameter.ps(k);      
      out.ps_plus(k) = Ps_plus;      
      out.ps_minus(k) = Ps_minus;
      
      out.F_kCol(k) = F.kuCol + F.klCol;
  end
end


function out = structcat(varargin)  
  for i=1:nargin
    auxSTR = varargin{i};
    if ~isstruct(auxSTR)
      error('structcat: The input arguments must be structures!');
    end
    
    auxNames = fieldnames(auxSTR);
    for n = 1:length(auxNames)
      out.(auxNames{n}) = auxSTR.(auxNames{n});
    end
  end
end
