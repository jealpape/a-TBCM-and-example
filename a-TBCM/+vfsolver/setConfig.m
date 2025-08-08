function out = setConfig(varargin)
%   out = inputParser;
  out = struct();
  
  out = getOneArgument(out,'noise',0,@isnumeric,varargin{:}); % state noise [x_u x_l x_b v_u v_l v_b]';
  out = getOneArgument(out,'stochastic',0,@islogical,varargin{:}); % state noise [x_u x_l x_b v_u v_l v_b]';
  
%   out.addOptional('output','basic',@(x)any(validatestring(x,{'basic','medium','full'})));
%   out.addOptional('fs',70000,@isnumeric); % [s^-1]  
%   % Kinematics
%   out.addOptional('kSolver','tts',@(x)any(validatestring(x,{'off','tts','ode4'})));  
%   out.addOptional('kModel','tbcm',@(x)any(validatestring(x,{'bcm','bcm+','tbcm'})));  
%   % Flow  
%   out.addOptional('fSolver','galindo2016',@(x)any(validatestring(x,{'off','titze2002','titze2008','lucero2015','galindo2016'}))); % incompleto
%   out.addOptional('fModel','titze',@(x)any(validatestring(x,{'titze'}))); % incompleto
%   out.addOptional('turbulence','on',@(x)any(validatestring(x,{'on','off'})));
%   % Pressures
%   out.addOptional('interaction','on',@(x)any(validatestring(x,{'on','off'})));  
%   out.addOptional('ps',726,@isnumeric); % [Pa]
%   out.addOptional('pe',0,@isnumeric); % [Pa]
%   out.addOptional('SimTime',0.2,@isnumeric); %[s]
%   % Acoustics
%   out.addOptional('supTract','Takemoto_e',@isstr); % Supra-glottal Tract
%   out.addOptional('subTract','Takemoto_weibel',@isstr); % Supra-glottal Tract
%   out.addOptional('compression',1,@isnumeric); % Supra-glottal Tract
%   % Configuration
%   out.addOptional('sex','male',@(x)any(validatestring(x,{'male','female'})));
%   out.addOptional('Act',0.04,@isnumeric); % Activation of the Cricothyroid [-]
%   out.addOptional('Ata',0.18,@isnumeric); % Activation of the Thyroarytenoid [-]
%   out.addOptional('Alc',0.5,@isnumeric); % Activation of the Lateral Cricoarytenoid [-]
%   out.addOptional('rotation',1,@isnumeric); % Rotation of the arythenoids [º]
%   out.addOptional('displacement',0,@isnumeric); % displacement of the arythenoids [m]
%   out.addOptional('constants','galindo2016',@isstr); % set of constants to use
% 
%   out.addOptional('noise',0,@isnumeric); % state noise [x_u x_l x_b v_u v_l v_b]';
%   out.addOptional('stochastic',0,@islogical); % state noise [x_u x_l x_b v_u v_l v_b]';
%   
%   out.parse(varargin{:});
  
  if out.stochastic
    auxPos = find(strcmpi(varargin,'act'));
    if isempty(auxPos)
      varargin = [varargin {'Act' min(1,max(0,out.Act + normrnd(0,0.03)))}];
    else
      varargin{auxPos+1} = min(1,max(0,out.Act + normrnd(0,0.03)));
    end
    
    auxPos = find(strcmpi(varargin,'ata'));
    if isempty(auxPos)
      varargin = [varargin {'Ata' min(1,max(0,out.Ata + normrnd(0,0.22)))}];
    else
      varargin{auxPos+1} = min(1,max(0,out.Ata + normrnd(0,0.22)));
    end
    
    auxPos = find(strcmpi(varargin,'ps'));
    if isempty(auxPos)
      varargin = [varargin {'ps' max(0,out.ps + normrnd(0,0.99))}];
    else
      varargin{auxPos+1} = max(0,out.ps + normrnd(0,0.99));
    end
    
    auxPos = find(strcmpi(varargin,'rotation'));
    if isempty(auxPos)
      varargin = [varargin {'rotation' out.rotation + normrnd(0,0.4)}];
    else
      varargin{auxPos+1} = out.rotation + normrnd(0,0.4);
    end
  end
  
  % Simulation
  out = getOneArgument(out,'output','basic',@(x)any(validatestring(x,{'basic','medium','full'})),varargin{:});
  out = getOneArgument(out,'fs',70000,@isnumeric,varargin{:}); % [s^-1]  
  % Kinematics
  out = getOneArgument(out,'kSolver','tts',@(x)any(validatestring(x,{'off','tts','ode4'})),varargin{:});  
  out = getOneArgument(out,'kModel','tbcm',@(x)any(validatestring(x,{'bcm','bcm+','tbcm'})),varargin{:});  
  % Flow  
  out = getOneArgument(out,'fSolver','galindo2016',@(x)any(validatestring(x,{'off','titze2002','titze2008','lucero2015','galindo2016'})),varargin{:}); % incompleto
  out = getOneArgument(out,'fModel','titze',@(x)any(validatestring(x,{'titze'})),varargin{:}); % incompleto
  out = getOneArgument(out,'turbulence','on',@(x)any(validatestring(x,{'on','off'})),varargin{:});
  % Pressures
  out = getOneArgument(out,'interaction','on',@(x)any(validatestring(x,{'on','off'})),varargin{:});  
  out = getOneArgument(out,'ps',726,@isnumeric,varargin{:}); % [Pa]
  out = getOneArgument(out,'pe',0,@isnumeric,varargin{:}); % [Pa]
  out = getOneArgument(out,'SimTime',0.2,@isnumeric,varargin{:}); %[s]
  % Acoustics
  out = getOneArgument(out,'supTract','Takemoto_e',@isstr,varargin{:}); % Supra-glottal Tract
  out = getOneArgument(out,'subTract','Takemoto_weibel',@isstr,varargin{:}); % Supra-glottal Tract
  out = getOneArgument(out,'compression',1,@isnumeric,varargin{:}); % Supra-glottal Tract
  % Configuration
  out = getOneArgument(out,'sex','male',@(x)any(validatestring(x,{'male','female'})),varargin{:});
  out = getOneArgument(out,'Act',0.04,@isnumeric,varargin{:}); % Activation of the Cricothyroid [-]
  out = getOneArgument(out,'Ata',0.18,@isnumeric,varargin{:}); % Activation of the Thyroarytenoid [-]
  out = getOneArgument(out,'Alc',0.5,@isnumeric,varargin{:}); % Activation of the Lateral Cricoarytenoid [-]
  out = getOneArgument(out,'rotation',1,@isnumeric,varargin{:}); % Rotation of the arythenoids [º]
  out = getOneArgument(out,'displacement',0,@isnumeric,varargin{:}); % displacement of the arythenoids [m]
  out = getOneArgument(out,'proportion',1,@isnumeric,varargin{:}); % Proportion of the relation between default anatomical parameters and the measured ones.
  out = getOneArgument(out,'constants','galindo2016',@isstr,varargin{:}); % set of constants to use
end

function structure = getOneArgument(structure,name,defaultValue,validation,varargin)
  auxPos = find(strcmpi(varargin,name));

  if ~isempty(auxPos)
    givenValue = varargin{auxPos+1};
  else
    givenValue = [];
  end
  
  if ~isempty(givenValue)
    if validation(givenValue)
      structure.(name) = givenValue;
    else
      warning_string = [name ':' num2str(givenValue) ' is not a valid value, default ''' defaultValue ''' is used'];
      warning('setConfig:argument',warning_string);
      structure.(name) = defaultValue;
    end
  else
    structure.(name) = defaultValue;
  end
end