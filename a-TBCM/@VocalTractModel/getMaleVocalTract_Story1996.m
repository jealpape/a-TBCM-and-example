%%
% getMaleVocalTract_Story1996: Function to get the male vocal tract
% configurations reported in [1], and the relevant simulation parameters,
% and to % configure the vocal tract model accordingly.
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
%         Available options = {'/i/','/I/','/E/','/AE/','/V/','/A/','/O/',
%                              '/o/', '/U/','/u/','/3/' }.
% Data: Struct gathering the relevan parameter for the selected vocal tract
%       configuration.
% AreaFun: Area function for the selected vocal tract configuration.
% Fsamp: Sampling frequency.
% Deltaz: length if the vocal tract sections.
%
% References:
% [1] B. H. Story, I. R. Titze, and E. A. Hoffman, “Vocal tract area
%     functions from magnetic resonance imaging,” J Acoust Soc Am, vol.
%     100, no. 1, pp. 537–554, Jul. 1996. 
%
% Coded by Gabriel Alzamendi, January 2020.
% Based on original code by Matías Zañartu and Gabriel Galindo
function varargout = getMaleVocalTract_Story1996(VTobj,varargin)
  VOCALTRACTOPTIONS = {'/i/','/I/','/E/','/AE/','/V/','/A/','/O/', ...
                       '/o/','/U/','/u/','/3/' };
  SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the idealized
                     %  vocal tract using WRA.
  LenghtSec = 3.977e-3; % [m] Length of vocal tract sections
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
         0.33 0.30 0.36 0.34 0.68 0.50 2.43 3.15 2.66 2.49 ...
         3.39 3.80 3.78 4.35 4.50 4.43 4.68 4.52 4.15 4.09 ...
         3.51 2.95 2.03 1.66 1.38 1.05 0.60 0.35 0.32 0.12 ...
         0.10 0.16 0.25 0.24 0.38 0.28 0.36 0.65 1.58 2.05 ...
         2.01 1.58 ];
        
    case '/I/'  % Near-close near-front unrounded vowel: IPA number 319
      section = 10^-4*[
        0.20 0.17 0.18 0.18 0.10 1.08 1.66 1.64 1.19 0.92 ,...
        1.13 2.48 2.76 2.97 3.43 3.32 3.48 3.96 3.79 3.88,...
        3.47 2.98 2.62 2.37 1.99 1.90 1.70 1.44 1.45 1.06,...
        0.87 0.75 1.06 1.29 1.78 1.83 1.70 1.97 1.92 1.62,...
        1.36 1.18];
      
    case '/E/'  % Open-mid front unrounded vowel: IPA number 303
      section = 10^-4*[
        0.21 0.13 0.16 0.14 0.06 0.78 1.25 1.24 0.99 0.72,...
        0.73 1.06 1.77 1.97 2.46 2.70 2.92 3.03 2.84 2.84,...
        2.83 2.36 2.14 2.00 1.78 1.81 1.79 1.50 1.37 1.36,...
        1.43 1.83 2.08 2.59 2.54 2.11 2.34 2.74 2.19 1.60,...
        ];
      
    case '/AE/'  % Near-open front unrounded vowel: IPA number 325
      section = 10^-4*[
        0.22 0.26 0.26 0.16 0.13 0.21 0.83 1.50 1.35 0.99,...
        0.69 1.35 2.32 2.13 1.94 2.17 2.85 3.26 3.73 3.80,...
        3.69 3.87 3.68 3.20 3.26 3.29 3.19 3.23 3.23 3.40,...
        3.78 3.84 3.98 4.41 4.56 4.79 4.39 4.42 4.23 4.56,...
        4.31 3.94];
      
    case '/V/'  % Open-mid back unrounded vowel: IPA number 314
      section = 10^-4*[
        0.33 0.28 0.23 0.15 0.17 0.33 0.39 1.02 1.22 1.14,...
        0.82 0.76 0.66 0.80 0.72 0.66 1.08 0.91 1.09 1.06,...
        1.09 1.17 1.39 1.55 1.89 2.17 2.46 2.65 3.13 3.81,...
        4.30 4.57 4.94 5.58 5.79 5.51 5.49 4.69 4.50 3.21,...
        2.79 2.11 1.98 1.17];
      
    case '/A/'  % Open back unrounded vowel: IPA number 305
      section = 10^-4*[
         0.45 0.20 0.26 0.21 0.32 0.30 0.33 1.05 1.12 0.85 ...
         0.63 0.39 0.26 0.28 0.23 0.32 0.29 0.28 0.40 0.66 ...
         1.20 1.05 1.62 2.09 2.56 2.78 2.86 3.02 3.75 4.60 ...
         5.09 6.02 6.55 6.29 6.27 5.94 5.28 4.70 3.87 4.13 ...
         4.25 4.27 4.69 5.03 ];
      
    case '/O/'  % Open-mid back rounded vowel: IPA number 306
      section = 10^-4*[
        0.61 0.28 0.19 0.10 0.07 0.30 0.18 1.13 1.42 1.21,...
        0.69 0.51 0.43 0.66 0.57 0.32 0.43 0.45 0.53 0.60,...
        0.77 0.65 0.58 0.94 2.02 2.50 2.41 2.62 3.29 4.34,...
        4.78 5.24 6.07 7.08 6.81 6.20 5.89 5.04 4.29 2.49,...
        1.84 1.33 1.19 0.88];
      
    case '/o/'  % Close-mid back rounded vowel: IPA number 307
      section = 10^-4*[
        0.18 0.17 0.23 0.28 0.59 1.46 1.60 1.11 0.82 1.01,...
        2.72 2.71 1.96 1.92 1.70 1.66 1.52 1.28 1.44 1.28,...
        0.89 1.25 1.38 1.09 0.71 0.46 0.39 0.32 0.57 1.06,...
        1.38 2.29 2.99 3.74 4.39 5.38 7.25 7.00 4.57 2.75,...
        1.48 0.68 0.39 0.14];
      
    case '/U/'  % Near-close near-back rounded vowel: IPA number 321
      section = 10^-4*[
        0.32 0.39 0.39 0.43 0.56 1.46 2.20 2.06 1.58 1.11,...
        1.11 1.26 1.30 0.98 0.93 0.83 0.61 0.97 0.75 0.93,...
        0.53 0.65 0.95 0.99 1.07 1.39 1.47 1.79 2.34 2.68,...
        3.36 3.98 4.74 5.48 5.69 5.57 4.99 4.48 3.07 1.67,...
        1.13 0.64 0.15 0.22];
      
    case '/u/'  % Close back rounded vowel: IPA number 308
      section = 10^-4*[
        0.40 0.38 0.28 0.43 0.55 1.72 2.91 2.88 2.37 2.10,...
        3.63 5.86 5.63 5.43 4.80 4.56 4.29 3.63 3.37 3.16,...
        3.31 3.22 2.33 2.07 2.07 1.52 0.74 0.23 0.15 0.22,...
        0.22 0.37 0.60 0.76 0.86 1.82 2.35 2.55 3.73 5.47,...
        4.46 2.39 1.10 0.77 0.41 0.86];
      
    case '/3/'  % Open-mid central unrounded vowel: IPA number 326
      section = 10^-4*[
        0.41 0.38 0.40 0.29 0.13 0.53 1.58 1.56 1.22 1.19,...
        1.00 0.77 0.92 1.19 1.27 1.35 1.48 1.56 1.61 1.87,...
        2.10 2.01 2.62 2.96 3.07 3.11 2.77 2.67 2.47 2.34,...
        2.25 1.90 1.32 0.76 0.44 0.45 0.92 2.05 3.25 3.63,...
        3.59 3.07 2.25 1.20];
      
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