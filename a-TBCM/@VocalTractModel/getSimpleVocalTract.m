%%
% getSimpleVocalTract: Function to get the idealized vocal tract
% configurations and other relevant simulation parameters, and to
% configure the vocal tract model accordingly.
%
% Structure: getSimpleVocalTract(VTobj,VTConf)
%            Data = getSimpleVocalTract(VTobj,VTConf)
%            [AreaFun,Fsamp] = getSimpleVocalTract(VTobj,VTConf)
%            [AreaFun,Fsamp,Deltaz] = getSimpleVocalTract(VTobj,VTConf)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% VTConf: string for the selected idealized vocal tract configuration. 
%         Available options = {'neutral','/a/','/i/','uniform'}.
% Data: Struct gathering the relevan parameter for the selected vocal tract
%       configuration.
% AreaFun: Area function for the selected vocal tract configuration.
% Fsamp: Sampling frequency.
% Deltaz: length if the vocal tract sections.
%
% References:
% [1] M. Zañartu, “Influence of Acoustic Loading on the Flow-Induced
%     Oscillations of Single Mass Models of the Human Larynx,” M.S. Thesis,
%     School of Electrical and Computer Engineering, Purdue University, May
%     2006.
% [2] Story, Brad H. 1995.  Physiologically Based Speech Simulation with an
%     Enhanced Wave-Reflection Model ofthe Vocal Tract. Iowa City, IA:
%     University of Iowa dissertation. 
%
% Coded by Gabriel Alzamendi, January 2020.
% Based on original code by Matías Zañartu and Gabriel Galindo
function varargout = getSimpleVocalTract(VTobj,varargin)
  VOCALTRACTOPTIONS = {'neutral','/a/','/i/','/u/','uniform','no-interaction'};
  N_sections = 44;
  SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the idealized
                     %  vocal tract using WRA.
  LenghtSec = 3.977e-3; % [m] Length of vocal tract sections
  if (nargin == 1)
    tract_type = 'neutral';
  elseif (nargin == 2)&&any(strcmp(varargin{1},VOCALTRACTOPTIONS))
    tract_type = char(varargin{1});
  else
    SolverOpts = '';
    for OptsAux = VOCALTRACTOPTIONS
      SolverOpts = strcat(SolverOpts, ' ''', char(OptsAux), ''', ');
    end
    SolverOpts = SolverOpts(1:end-1);
    error('Vocal tract configuration not available! The acceptable options are:%s. \n',SolverOpts);  
  end
  % Available idealized vocal tract configurations
  switch tract_type
    case 'neutral' 
      section = 1e-4*[0.5*ones(10,1); 4*ones(N_sections - 10,1)];
    case 'uniform' 
      section = 1e-4*[4*ones(N_sections,1)];
    case 'no-interaction' 
      section = 1e-4*[30*ones(N_sections,1)];
    case '/a/' 
      section = 1e-4*[0.5*ones(N_sections/2,1); 3*ones(N_sections/2,1)];
    case '/i/' 
      section = 1e-4*[3*ones(N_sections/2,1); 0.5*ones(N_sections/2,1)];
    case '/u/' 
      section = 1e-4*[3*ones(N_sections-2,1); 0.5*ones(2,1)];
  end
  
  % Setting model parameters
  VTobj.AreaFunction = section;
  VTobj.N_AreaSection = N_sections;
  VTobj.Delta_z = LenghtSec;
  VTobj.sex = 'male';
  VTobj.setSimulationParameter(SampFreq);
  
  % Function Output
  switch nargout
    case 1
      Info.AreaFunction = section;
      Info.N_AreaSection = N_sections;
      Info.Delta_z = LenghtSec;
      Info.sex = 'male';
      Info.fs = SampFreq;
      Info.Ts = 1/SampFreq;
      varargout{1} = Info;
    case 2
      varargout{1} = section;
      varargout{2} = SampFreq;
    case 3
      varargout{1} = section;
      varargout{2} = SampFreq;
      varargout{3} = LenghtSec;
  end
  
end