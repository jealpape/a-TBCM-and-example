%%
% getSubglottalTract: Function to get the subglottal tract
% configurations reported in the literature [1,2,3] and the relevant 
% simulation parameters, and to configure the vocal tract model 
% accordingly.
%
% Structure: getSubglottalTract(VTobj,VTConf)
%            Data = getSubglottalTract(VTobj,VTConf)
%            [AreaFun,Fsamp] = getSubglottalTract(...)
%            [AreaFun,Fsamp,Deltaz] = getSubglottalTract(...)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% VTConf: string for the selected idealized vocal tract configuration.
%         Default option is 'Story_smooth'.
%         Available options = {'sub1','sub2','sub3', 'sub4', 'sub5', ...
%                              'short weibel', 'short weibel 2', ...
%                              'Story_simple', 'Story_raw','Story_smooth'}.
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
% [3] E. R. Weibel, Morphometry of the Human Lung, First Edition, New York,
%     Springer, 1963. 
%
% Coded by Gabriel Alzamendi, January 2020.
% Based on original code by Matías Zañartu and Gabriel Galindo
function varargout = getSubglottalTract(VTobj,varargin)
  VOCALTRACTOPTIONS = {'sub1','sub2','sub3', 'sub4', 'sub5', ...
                       'short weibel', 'short weibel 2', ...
                       'Story_simple', 'Story_raw', 'Story_smooth', ...
                       'Takemoto_weibel','no-interaction'};
  
  if (nargin == 1)
    tract_type = 'Story_smooth';
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
      
    case 'sub1' % Subglotal tract shape 1, from Weibel
      section = 1e-4*[
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2,...
        2.2	2.0 2.0	2.0	2.0	2.0	2.0	1.8	1.8	1.8,...
        2.7	2.7	2.7	3.6	3.6	3.6	4.9	4.9	6.7	6.7,...
        9.8	9.8	13.8 13.8 23.1 30.7 ]; % 39.1];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.
      
    case 'sub2' % Subglotal tract shape 2, from Story
      section = 1e-4*[
        1.0	1.3	1.5	1.8	2.0	2.0	2.4	2.7	3.0	3.2,...
        3.3	3.4	3.5	3.4	3.3	3.1	2.9	2.7	2.6	2.5,...
        2.4	2.3	2.3	2.2	2.1	1.8	1.6	2.0	3.5	4.0];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

    case 'sub3' % Modified Subglotal tract shape 2, from Story
      section = 1e-4*[
        1.0	1.3	1.6	1.8	2.1	2.4	2.6	2.8	2.9	3.0,...
        3.1	3.2	3.2	3.2	3.1	3.1	3.0	2.8	2.7	2.6,...
        2.4	2.3	2.1	2.0	2.0	1.9	2.0	2.1	2.2	2.5];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

    case 'sub4' % WEIBEL DESIGN
      section = 1e-4*[
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33,...
        2.33	2.33	2.13	2.13	2.13	2.13	2.13	2.00	2.00	2.48,...
        2.48	2.48	3.11	3.11	3.11	3.96	3.96	5.10	5.10	6.95,...
        6.95	9.56	13.40	19.60	28.80	44.50	69.40	113   338   4688];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.


    case 'sub5' % Modified Subglotal tract shape with boundary correction
      section = 1e-4*[
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
        2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2,...
        2.2	2.0 2.0	2.0	2.0	2.0	2.0	1.8	1.8	1.8,...
        2.7	2.7	2.7	3.6	3.6	3.6	4.9	4.9	6.7	6.7,...
        9.8 1.1];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

%     case 'sub6' % Like sub5 but with no boundary correction
%       section = 1e-4*[
%         2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
%         2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
%         2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7	2.7,...
%         2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2	2.2,...
%         2.2	2.0 2.0	2.0	2.0	2.0	2.0	1.8	1.8	1.8,...
%         2.7	2.7	2.7	3.6	3.6	3.6	4.9	4.9	6.7	6.7,...
%         9.8];
      
    case 'short weibel' % Weibel Design
      section = 1e-4*[
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33,...
        ];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

    case 'short weibel 2' % Weibel Design
      section = 1e-4*[
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
%         2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33,...
        ];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.
    
    case 'Story_simple'
      section = 1e-4*[  
        0.30 0.33 0.40 0.52 0.66 0.82 0.97 1.10 1.22 1.31 ...
        1.37 1.42 1.45 1.47 1.48 1.49 1.49 1.50 1.50 1.50 ...
        1.50 1.50 1.50 1.50 1.50 1.50 1.50 1.50 1.50 1.50 ...
        1.50 1.50 1.50 1.50 ];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.
    
    case 'Story_raw'
      section = 1e-4*[  
        0.30 0.33 0.40 0.52 1.03 1.47 1.71 1.74 2.77 2.94 ...
        3.03 2.48 2.49 3.19 3.51 3.30 3.03 3.62 3.64 3.20 ...
        2.87 2.85 2.58 2.52 2.56 2.43 2.48 2.74 2.31 2.12 ...
        1.90 2.06 2.45 3.20 4.00 5.00 6.00 8.00 ];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.
    
    case 'Story_smooth'
      section = 1e-4*[  
        0.31 0.34 0.41 0.67 1.01 1.37 1.73 2.06 2.36 2.62 ...
        2.82 2.99 3.11 3.19 3.23 3.25 3.24 3.22 3.18 3.12 ...
        3.05 2.97 2.87 2.76 2.64 2.51 2.38 2.26 2.15 2.10 ...
        2.11 2.24 2.53 3.05 4.00 5.00 6.00 8.00 ];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

    case 'Takemoto_weibel'  % Use fs=C/(2*.25e-2)
%       len = 0.25e-2;
%       sex = 'Unknown';
%       study = 'Unknown';
%       type = 'Unknown';
      section = 1e-4*[
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.54,...
        2.54	2.54	2.54	2.54	2.54	2.54	2.54	2.33	2.33	2.33,...
        2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33	2.33,...
        2.33	2.33	2.33];
      LenghtSec = 0.25e-2; % [m] Length of tract sections
      SampFreq = 70e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.
    
    case 'no-interaction'
      section = 1e-4*[30*ones(1,38)];
      LenghtSec = 3.977e-3; % [m] Length of tract sections
      SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the 
                         % subglottal using WRA.

                         
      section = section(end:-1:1); % reverse the tract
  end
  
  section = section(:);
  N_sections = length(section);
  
  % Setting model parameters
  VTobj.AreaFunction = section;
  VTobj.N_AreaSection = N_sections;
  VTobj.Delta_z = LenghtSec;
  VTobj.sex = '';
  VTobj.setSimulationParameter(SampFreq);
  
  % Function Output
  switch nargout
    case 1
      Info.AreaFunction = section;
      Info.N_AreaSection = N_sections;
      Info.Delta_z = LenghtSec;
      Info.sex = '';
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