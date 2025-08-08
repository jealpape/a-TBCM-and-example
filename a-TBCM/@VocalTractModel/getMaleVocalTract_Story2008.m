%%
% getMaleVocalTract_Story2008: Function to get the male vocal tract
% configurations reported in [1], and the relevant simulation parameters,
% and to configure the vocal tract model accordingly.
%
% Structure: getMaleVocalTract_Story1996(VTobj,VTConf)
%            Data = getMaleVocalTract_Story1996(VTobj,VTConf)
%            [AreaFun,Fsamp] = getMaleVocalTract_Story1996(VTobj,VTConf)
%            [AreaFun,Fsamp,Deltaz] = getMaleVocalTract_Story1996(VTobj,VTConf)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% VTConf: string for the selected idealized vocal tract configuration. 
%         Available options = {'/i/','/I/','/e/','/E/','/AE/','/V/','/A/','/O/',
%                              '/o/', '/U/','/u/' }.
% Data: Struct gathering the relevan parameter for the selected vocal tract
%       configuration.
% AreaFun: Area function for the selected vocal tract configuration.
% Fsamp: Sampling frequency.
% Deltaz: length if the vocal tract sections.
%
% References:
% [1] B. H. Story, “Comparison of magnetic resonance imaging-based vocal
%     tract area functions obtained from the same speaker in 1994 and
%     2002,” J. Acoust. Soc. Am., vol. 123, no. 1, pp. 327–335, Jan. 2008.
%
% Coded by Gabriel Alzamendi, January 2020.
% Based on original code by Matías Zañartu and Gabriel Galindo
function varargout = getMaleVocalTract_Story2008(VTobj,varargin)
  VOCALTRACTOPTIONS = {'/i/','/I/','/e/','/E/','/AE/','/V/','/A/', ...
                       '/O/','/o/','/U/','/u/' };
  SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the idealized
                     %  vocal tract using WRA.
  
  if (nargin == 1)
    tract_type = '/A/';
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
    case '/i/'  % Close front unrounded vowel: IPA number 301
      section = 10^-4*[
          0.51 0.59 0.62 0.72 1.24 2.30 3.30 3.59 3.22 2.86 ...
          3.00 3.61 4.39 4.95 5.17 5.16 5.18 5.26 5.20 5.02 ...
          4.71 4.13 3.43 2.83 2.32 1.83 1.46 1.23 1.08 0.94 ...
          0.80 0.67 0.55 0.46 0.40 0.36 0.35 0.35 0.38 0.51 ...
          0.74 0.92 0.96 0.91 ];
      LenghtSec = 0.384e-2; % [m] Length of vocal tract sections
        
    case '/I/'  % Near-close near-front unrounded vowel: IPA number 319
      section = 10^-4*[
         0.28 0.21 0.21 0.30 0.47 0.71 1.12 1.48 1.35 1.05 ...
         0.92 0.92 1.19 1.94 2.83 3.31 3.48 3.60 3.64 3.49 ...
         3.20 2.90 2.59 2.21 1.87 1.54 1.20 0.92 0.74 0.59 ...
         0.52 0.54 0.59 0.65 0.71 0.67 0.61 0.57 0.50 0.48 ...
         0.54 0.73 0.93 0.82 ];
      LenghtSec = 0.376e-2; % [m] Length of vocal tract sections
      
    case '/e/'  % Open-mid front unrounded vowel: IPA number 303
      section = 10^-4*[
         0.29 0.26 0.30 0.40 0.55 0.74 0.99 1.09 0.90 0.69 ...
         0.77 1.31 2.13 2.74 3.03 3.23 3.33 3.27 3.09 2.84 ...
         2.66 2.46 2.14 1.79 1.44 1.17 1.00 0.88 0.80 0.81 ...
         0.85 0.84 0.86 0.96 1.18 1.35 1.48 1.62 1.49 1.29 ...
         1.24 1.17 1.04 0.95 ];
      LenghtSec = 0.386e-2; % [m] Length of vocal tract sections
      
    case '/E/'  % Open-mid front unrounded vowel: IPA number 303
      section = 10^-4*[
         0.37 0.30 0.24 0.23 0.29 0.41 0.58 0.82 0.97 0.82 ...
         0.62 0.54 0.54 0.68 1.09 1.62 2.03 2.35 2.51 2.39 ...
         2.22 2.13 2.00 1.78 1.58 1.43 1.31 1.23 1.24 1.38 ...
         1.61 1.82 1.96 2.01 2.00 1.95 1.77 1.48 1.30 1.21 ...
         1.10 0.99 0.91 0.78 ];
      LenghtSec = 0.393e-2; % [m] Length of vocal tract sections
      
    case '/AE/'  % Near-open front unrounded vowel: IPA number 325
      section = 10^-4*[
         0.31 0.21 0.18 0.23 0.33 0.50 0.78 0.96 0.85 0.63 ...
         0.46 0.36 0.33 0.46 0.73 1.00 1.30 1.66 1.97 2.06 ...
         2.03 2.01 1.89 1.66 1.49 1.42 1.37 1.34 1.41 1.58 ...
         1.82 2.19 2.63 2.97 3.17 3.40 3.56 3.57 3.58 3.44 ...
         3.15 3.38 3.99 4.17 ];
      LenghtSec = 0.366e-2; % [m] Length of vocal tract sections
      
    case '/V/'  % Open-mid back unrounded vowel: IPA number 314
      section = 10^-4*[
         0.23 0.34 0.47 0.60 0.77 1.06 1.26 1.09 0.80 0.65 ...
         0.56 0.47 0.37 0.24 0.17 0.18 0.23 0.27 0.28 0.29 ...
         0.33 0.52 0.97 1.50 1.91 2.23 2.65 3.29 4.13 5.00 ...
         5.77 6.33 6.61 6.63 6.45 6.04 5.39 4.42 3.29 2.37 ...
         1.74 1.36 1.17 0.99 ];
      LenghtSec = 0.390e-2; % [m] Length of vocal tract sections
      
    case '/A/'  % Open back unrounded vowel: IPA number 305
      section = 10^-4*[
         0.56 0.62 0.66 0.78 0.97 1.16 1.12 0.82 0.55 0.45 ...
         0.37 0.29 0.21 0.15 0.16 0.25 0.34 0.43 0.54 0.61 ...
         0.67 0.98 1.76 2.75 3.52 4.08 4.74 5.61 6.60 7.61 ...
         8.48 9.06 9.29 9.26 9.06 8.64 7.91 6.98 6.02 5.13 ...
         4.55 4.52 4.71 4.72 ];
      LenghtSec = 0.388e-2; % [m] Length of vocal tract sections
      
    case '/O/'  % Open-mid back rounded vowel: IPA number 306
      section = 10^-4*[
         0.27 0.43 0.54 0.67 0.83 0.92 0.89 0.73 0.55 0.44 ...
         0.37 0.28 0.18 0.14 0.15 0.14 0.13 0.14 0.18 0.20 ...
         0.23 0.49 1.03 1.58 2.06 2.61 3.35 4.34 5.51 6.70 ...
         7.75 8.63 9.29 9.59 9.42 8.78 7.82 6.50 4.95 3.47 ...
         2.15 1.38 1.11 0.90 ];
      LenghtSec = 0.395e-2; % [m] Length of vocal tract sections
      
    case '/o/'  % Close-mid back rounded vowel: IPA number 307
      section = 10^-4*[
         0.38 0.45 0.57 0.77 1.31 1.92 1.74 1.11 0.75 0.59 ...
         0.57 0.68 0.73 0.67 0.58 0.49 0.44 0.42 0.49 0.53 ...
         0.38 0.30 0.45 0.61 0.71 0.79 0.86 1.01 1.41 2.09 ...
         3.00 4.10 5.16 6.22 7.34 8.15 8.61 8.37 6.76 4.37 ...
         2.30 1.06 0.58 0.47 ];
      LenghtSec = 0.417e-2; % [m] Length of vocal tract sections
      
    case '/U/'  % Near-close near-back rounded vowel: IPA number 321
      section = 10^-4*[
         0.37 0.38 0.49 0.62 0.85 1.28 1.62 1.47 1.04 0.81 ...
         1.03 1.44 1.49 1.28 1.06 0.85 0.69 0.56 0.42 0.27 ...
         0.27 0.41 0.51 0.49 0.47 0.50 0.58 0.80 1.19 1.62 ...
         2.27 3.24 4.16 5.00 5.70 6.11 6.21 6.29 6.24 4.91 ...
         2.61 1.09 0.63 0.59 ];
      LenghtSec = 0.440e-2; % [m] Length of vocal tract sections
      
    case '/u/'  % Close back rounded vowel: IPA number 308
      section = 10^-4*[
         0.54 0.61 0.66 0.75 1.13 1.99 2.83 2.90 2.52 2.40 ...
         2.83 3.56 3.99 3.89 3.50 3.04 2.64 2.44 2.31 2.07 ...
         1.80 1.52 1.14 0.74 0.42 0.22 0.14 0.20 0.47 0.89 ...
         1.15 1.42 2.17 3.04 3.69 4.70 5.74 5.41 3.82 2.34 ...
         1.35 0.65 0.29 0.16 ];
      LenghtSec = 0.445e-2; % [m] Length of vocal tract sections
      
  end
  section = section(:);
  N_sections = length(section);
  
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