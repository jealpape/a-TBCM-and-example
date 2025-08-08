classdef SubglottalTractModel < handle
% Handle class for modeling acoustic propagation through the subglottal
% tract. Different methods are implemented for setting and handle the
% physiological subglottal tract area functions, and for simulating the
% acoustic effects of the subglottal tract into the glottal source.
% Differente strategies are implemented for dealing with the hidrostatic
% pressure in the trachea/lung end.
%
% Coded by Gabriel Alzamendi, January 2020.

  properties (Constant, Hidden)
   C_AIR = 350; % [m/s] speed of sound
   RHO_AIR = 1.146; % [kg m^-3]  air density
   SOLVER_LIST = {'WRA','StateSpaceWRA'}; % List of available tract solvers
   
   NOTRACTERRORMSG = 'It is required an area function describing the subglottal tract!'
  end
  
  properties (SetAccess = protected)
    % Dynamic state variable
    Model = 'SubglottalTract';
    xData = 0; % State variable descring backward and forward acoustic waves 
    n_IterCont = 0; % Simulation time index
    % General object description
    Delta_z = 0; % [m] Lenght of the tract sections 
    AreaFunction = []; % [m^2] Area function describing vocal tract configuration
    N_AreaSection = 0; % Number of sections in the area function (an even number is required!)
    solver = 'WRA'; % Solver ofr the vocal tract simulation
    sex = ''; % Sex information
    r_end = -0.2; % Reflection coefficient for the traqueal junt
    
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
    function SGTObj = SubglottalTractModel(varargin)
      if nargin ==0
        SGTObj.sex = 'male';
      elseif (nargin == 1)&&(ischar(varargin{1}))
        if (strcmpi(varargin{1},'male'))
          SGTObj.sex = 'male';
        elseif (strcmpi(varargin{1},'female'))
          SGTObj.sex = 'female';
        else
          error('Acceptable ''Sex'' options are ''male'' or ''female'' ')
        end
      end
    end
    
    function SetSolver(SGTObj,varargin)
    % Function for setting the solver for the computational simulation of
    % the vocal tract model.
      if nargin == 1
        SGTObj.solver = 'WRA'; % Default setting
      elseif (nargin == 2)&&(ischar(varargin{1}))&&any(strcmp(varargin{1},SGTObj.SOLVER_LIST))
        SGTObj.solver = varargin{1}; % solver setting
      else
        SolverOpts = '';
        for OptsAux = SGTObj.SOLVER_LIST
          SolverOpts = strcat(SolverOpts, ' ''', char(OptsAux), ''', ');
        end
        SolverOpts = SolverOpts(1:end-1);
        error('The acceptable solver options are:%s. \n',SolverOpts);
      end    
    end
        
    function InitModel(SGTObj)
    % Function for initializing the dynamic state and simulation time index
    % prior to run the simulation of the vocal tract model
%       switch SGTObj.solver
%         case 'WRA'
%           SGTObj.InitWRA;
%         case 'StateSpaceWRA'
%           SGTObj.InitStateSpaceWRA;
%         otherwise
%           error('Incorrect model solver!')
%       end
      
      % Checking the area function describing the tract  
      if isempty(SGTObj.AreaFunction)
        error(SGTObj.NOTRACTERRORMSG)  
      end
      
      % Initializing the internal state
      SGTObj.xData = zeros(2*SGTObj.N_AreaSection,1);
      
      SGTObj.n_IterCont = 0; % Simulation time index
    end
        
    function Pressure = get.Pressure(SGTObj) % CHECK!
    % Function for computing Pressure distribution in the vocal tract
      switch SGTObj.solver
        case 'WRA'
          Pressure = SGTObj.PressureWRA;
        case 'StateSpaceWRA'
          Pressure = SGTObj.PressureWRA;
        otherwise
          error('Incorrect model solver!')
      end
    
    end
    
    function Airflow = get.Airflow(SGTObj) % CHECK!
    % Function for computing Pressure distribution in the vocal tract
      switch SGTObj.solver
        case 'WRA'
          Airflow = SGTObj.AirflowWRA;
        case 'StateSpaceWRA'
          Airflow = SGTObj.AirflowWRA;
        otherwise
          error('Incorrect model solver!')
      end
    
    end
    
    % Functions defined on separate files
    varargout = getSubglottalTract(SGTObj,varargin)
%     
%     varargout = getMaleVocalTract_Story1996(SGTObj,varargin)
%     
    setSimulationParameter(SGTObj,Param)
    
    varargout = Simulate(SGTObj,Ug_n,varargin)
    
    plotTract(SGTObj)
    
  end
    
  methods (Access = protected, Hidden)
        
    %% Methods for the WRA solver        
    function Pressure = PressureWRA(SGTObj)
    % Function for computing Pressure distribution in the vocal tract
    % according to the WRA model
      Pressure = SGTObj.xData(1:SGTObj.N_AreaSection) + ...
                 SGTObj.xData(SGTObj.N_AreaSection+1:2*SGTObj.N_AreaSection);
    end
        
    function Airflow = AirflowWRA(SGTObj)
    % Function for computing airflow distribution in the vocal tract
    % according to the WRA model
      crho = SGTObj.C_AIR * SGTObj.RHO_AIR;
      Airflow = SGTObj.AreaFunction/crho .* ...
               (SGTObj.xData(1:SGTObj.N_AreaSection) - ...
                SGTObj.xData(SGTObj.N_AreaSection+1:2*SGTObj.N_AreaSection) );
    end
    
    
    % Functions defined on separate files
    varargout = SimulateWRAmodeA(SGTObj,Ug,re)
    
    varargout = SimulateWRAmodeB(SGTObj,Ug_n,PL_flag,PL_n,varargin)
    
    %% Methods for the StateSpaceWRA solver
    
    function varargout = SimulateStateSpaceWRAmodeA(SGTObj,Ug_n)
    % Function for simulating the acoustic propagation in the state space
    % formulation of the WRA method. Only glottal airflow Ug_n is assumed
    % as input variable.
      PressureWaves = SGTObj.A_ss*SGTObj.xData + SGTObj.Gamma_ss*Ug_n;
      SGTObj.xData = PressureWaves;
      
      % Output variables
      if nargout == 1
        varargout{1} =  PressureWaves;
      elseif (nargout>1)
        error('It is requested more output varaibles than allowed!')  
      end
    end
    
    function varargout = SimulateStateSpaceWRAmodeB(SGTObj,Ug_n,PL_n)
    % Function for simulating the acoustic propagation in the state space
    % formulation of the WRA method. Both lung/tracheal pressure and
    % glottal airflow Ug_n are considered as input variables.
      PressureWaves = SGTObj.A_ss*SGTObj.xData + SGTObj.Gamma_ss*[Ug_n;PL_n];
      SGTObj.xData = PressureWaves;
      
      % Output variables
      if nargout == 1
        varargout{1} =  PressureWaves;
      elseif (nargout>1)
        error('It is requested more output varaibles than allowed!')  
      end
    end    
    
    % Functions defined on separate files
    varargout = getStateSpaceWRAModelmodeA(SGTObj,varargin)
      
  end
end