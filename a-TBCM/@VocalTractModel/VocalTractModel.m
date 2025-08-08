classdef VocalTractModel < handle
% Handle class for modeling acoustic propagation through the vocal tract
% and voice signal production. Different methods are implemented for
% setting and handle the physiological vocal tract area functions, and for
% simulating the modulation effects of the vocal tract into the glottal
% source.  
%
% Coded by Gabriel Alzamendi, January 2020.

  properties (Constant, Hidden)
   C_AIR = 350; % [m/s] speed of sound
   RHO_AIR = 1.146; % [kg m^-3]  air density
   SOLVER_LIST = {'WRA','StateSpaceWRA'}; % List of available tract solvers
   
   NOTRACTERRORMSG = 'It is required an area function describing the vocal tract!'
  end
  
  properties (SetAccess = protected)
    % Dynamic state variable
    Model = 'VocalTract';
    xData = 0; % State variable descring backward and forward acoustic waves 
    n_IterCont = 0; % Simulation time index
    % General object description
    Delta_z = 0; % [m] Lenght of the tract sections 
    AreaFunction = []; % [m^2] Area function describing vocal tract configuration
    N_AreaSection = 0; % Number of sections in the area function (an even number is required!)
    solver = 'WRA'; % Solver ofr the vocal tract simulation
    sex = ''; % Sex information
    
    % Internal variables for StateSpaceWRA method
    A_ss = 0; % State matrix
    Gamma_ss = 0; % Imput matrix
    SSWRAvarOK = false; % boolean variable informing whether the internal 
                        % variables for StateSpaceWRA method were computed 
                        % correctily (= true) or not (= false).
    ComputeSSWRAModAlways = false; % boolean variable informing whether the 
                        % internal variables for StateSpaceWRA method need 
                        % to be computed for each iteration (= true) 
                        % or not (= false).
    
    % Simulation parameters
    fs = 0; % [Hz] Sampling frequency for the numerical simulation
    Ts = 0; % [s] Sampling period for the numerical simulation
    SimParamOK = false;  % If 'true' simulation parameters are set,
                         % otherwise simulation parameters are missing
  end
  
  properties (Dependent)
    % Acoustic phenomena in the tract
    Pressure % Pressure distribution throughout the tract. 
    Airflow % Air flow (volume velocity) distribution throughout the tract.
  end
    
  methods
    % Class constructor
    function VTobj = VocalTractModel(varargin)
      if nargin ==0
        VTobj.sex = 'male';
      elseif (nargin == 1)&&(ischar(varargin{1}))
        if (strcmpi(varargin{1},'male'))
          VTobj.sex = 'male';
        elseif (strcmpi(varargin{1},'female'))
          VTobj.sex = 'female';
        else
          error('Acceptable ''Sex'' options are ''male'' or ''female'' ')
        end
      end
    end
    
    function SetSolver(VTobj,varargin)
    % Function for setting the solver for the computational simulation of
    % the vocal tract model.
      if nargin == 1
        VTobj.solver = 'WRA'; % Default setting
      elseif (nargin == 2)&&(ischar(varargin{1}))&&any(strcmp(varargin{1},VTobj.SOLVER_LIST))
        VTobj.solver = varargin{1}; % solver setting
      else
        SolverOpts = '';
        for OptsAux = VTobj.SOLVER_LIST
          SolverOpts = strcat(SolverOpts, ' ''', char(OptsAux), ''', ');
        end
        SolverOpts = SolverOpts(1:end-1);
        error('The acceptable solver options are:%s. \n',SolverOpts);
      end    
    end
        
    function InitModel(VTobj)
    % Function for initializing the dynamic state and simulation time index
    % prior to run the simulation of the vocal tract model
      switch VTobj.solver
        case 'WRA'
          VTobj.InitWRA;
        case 'StateSpaceWRA'
          VTobj.InitStateSpaceWRA;
        otherwise
          error('Incorrect model solver!')
      end
      
      VTobj.n_IterCont = 0; % Simulation time index
    end
        
    function varargout = Simulate(VTobj,Ug_n,varargin)
    % Function for simulating acoustic progation in the vocal tract
    %
    % Structure: SimulateWRA(VTobj,Ug_n)
    %            SimulateWRA(VTobj,Ug_n,re_n)
    %            PressureWaves = SimulateWRA(...)
    %
    % where
    %
    % VTobj: is an object from VocalTractModel (handle) class,
    % Ug_n: is the air flow (volume velocity) value for the n-th instant,
    % re_n: is the supraglottal reflection coefficient (=1 by default),
    % PressureWaves: is vector gathering the backward (Bn) and forward (Fn)
    %                acoustic wave components and the radiated acoustic
    %                pressure Pout
    %                PressureWaves = [B1 B2 ... BL F1 F2 ... FL FL_p Pout],
    %                with FL_p an auxiliari variable.
      re_n = 1;
      if (nargin == 3)&&isnumeric(varargin{1})&&(abs(varargin{1})<=1)
        re_n = varargin{1};
      end
      
      switch VTobj.solver
        case 'WRA'
          PressureWaves = VTobj.SimulateWRA(Ug_n,re_n);
        case 'StateSpaceWRA'
          if (VTobj.ComputeSSWRAModAlways)||(nargin == 3)
            VTobj.getStateSpaceWRAModel(re_n);
          end
          PressureWaves = VTobj.SimulateStateSpaceWRA(Ug_n); % Not implemented
        otherwise
          error('Incorrect model solver!')
      end
      
      % Increase time acumulator
      VTobj.n_IterCont = VTobj.n_IterCont+1; % Simulation time index
      
      % Output variables
      if nargout == 1
        varargout{1} =  PressureWaves;
      elseif (nargout>1)
        error('It is requested more output varaibles than allowed!')  
      end
      
    end
    
    function Pressure = get.Pressure(VTobj) % CHECK!
    % Function for computing Pressure distribution in the vocal tract
      switch VTobj.solver
        case 'WRA'
          Pressure = VTobj.PressureWRA;
        case 'StateSpaceWRA'
          Pressure = VTobj.PressureWRA;
        otherwise
          error('Incorrect model solver!')
      end
    
    end
    
    function Airflow = get.Airflow(VTobj) % CHECK!
    % Function for computing Pressure distribution in the vocal tract
      switch VTobj.solver
        case 'WRA'
          Airflow = VTobj.AirflowWRA;
        case 'StateSpaceWRA'
          Airflow = VTobj.AirflowWRA;
        otherwise
          error('Incorrect model solver!')
      end
    
    end
    
    % Functions defined on separate files
    varargout = getSimpleVocalTract(VTobj,varargin)
    
    varargout = getMaleVocalTract_Story1996(VTobj,varargin)
    
    varargout = getMaleVocalTract_Story2008(VTobj,varargin)
    
    setSimulationParameter(VTobj,Param)
    
    plotTract(VTobj)
    
  end
    
  methods (Access = protected, Hidden)
        
    %% Methods for the WRA solver
    function InitWRA(VTobj)
    % Function for initializing the dynamic internal state for WRA (Wave
    % Reflection analogue) method.
    
      % Checking the area function describing the tract  
      if isempty(VTobj.AreaFunction)
        error(VTobj.NOTRACTERRORMSG)  
      end
      
      % Initializing the internal state
      VTobj.xData = zeros(2*VTobj.N_AreaSection+2,1);
    end
        
    function Pressure = PressureWRA(VTobj)
    % Function for computing Pressure distribution in the vocal tract
    % according to the WRA model
      Pressure = [VTobj.xData(1:VTobj.N_AreaSection) + ...
                 VTobj.xData(VTobj.N_AreaSection+1:2*VTobj.N_AreaSection);...
                 VTobj.xData(2*VTobj.N_AreaSection+2)];
    end
        
    function Airflow = AirflowWRA(VTobj)
    % Function for computing airflow distribution in the vocal tract
    % according to the WRA model
      crho = VTobj.C_AIR * VTobj.RHO_AIR;
      Airflow = VTobj.AreaFunction/crho .* ...
               (VTobj.xData(VTobj.N_AreaSection+1:2*VTobj.N_AreaSection) - ...
                VTobj.xData(1:VTobj.N_AreaSection) );
    end
    
    
    % Functions defined on separate files
    varargout = SimulateWRA(VTobj,Ug,re)
    
    %% Methods for the StateSpaceWRA solver
    function InitStateSpaceWRA(VTobj,varargin)
    % Function for initializing the dynamic internal state for state space
    % WRA (Wave Reflection analogue) method.
    
      re_n = 1;
      if (nargin == 2)
        re_n = varargin{1};
      end
    
      % Checking the area function describing the tract  
      if isempty(VTobj.AreaFunction)
        error(VTobj.NOTRACTERRORMSG)  
      end
      
      % Computing StateSpaceWRA internal variables
      VTobj.getStateSpaceWRAModel(re_n);
      
      % Initializing the internal state
      VTobj.xData = zeros(2*VTobj.N_AreaSection+2,1);
    end
    
    % Functions defined on separate files
    varargout = getStateSpaceWRAModel(VTobj,varargin)
    
    function varargout = SimulateStateSpaceWRA(VTobj,Ug_n)
    % Function for simulating the acoustic propagation in the state space
    % formulation of the WRA method.
      PressureWaves = VTobj.A_ss*VTobj.xData + VTobj.Gamma_ss*Ug_n;
      VTobj.xData = PressureWaves;
      
      % Output variables
      if nargout == 1
        varargout{1} =  PressureWaves;
      elseif (nargout>1)
        error('It is requested more output varaibles than allowed!')  
      end
    end    
      
  end
end