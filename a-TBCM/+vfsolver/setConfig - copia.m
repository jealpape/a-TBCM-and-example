function out = setConfig(varargin)
  out = inputParser;

  % Simulation
  out.addOptional('output','basic',@(x)any(validatestring(x,{'basic','medium','full'})));
  out.addOptional('fs',70000,@isnumeric); % [s^-1]  
  % Kinematics
  out.addOptional('kSolver','tts',@(x)any(validatestring(x,{'off','tts','ode4'})));  
  out.addOptional('kModel','tbcm',@(x)any(validatestring(x,{'bcm','bcm+','tbcm'})));  
  % Flow  
  out.addOptional('fSolver','galindo2016',@(x)any(validatestring(x,{'off','titze2002','titze2008','lucero2015','galindo2016'}))); % incompleto
  out.addOptional('fModel','titze',@(x)any(validatestring(x,{'titze'}))); % incompleto
  out.addOptional('turbulence','on',@(x)any(validatestring(x,{'on','off'})));
  % Pressures
  out.addOptional('interaction','on',@(x)any(validatestring(x,{'on','off'})));  
  out.addOptional('ps',726,@isnumeric); % [Pa]
  out.addOptional('pe',0,@isnumeric); % [Pa]
  out.addOptional('SimTime',0.2,@isnumeric); %[s]
  % Acoustics
  out.addOptional('supTract','Takemoto_e',@isstr); % Supra-glottal Tract
  out.addOptional('subTract','Takemoto_weibel',@isstr); % Supra-glottal Tract
  out.addOptional('compression',1,@isnumeric); % Supra-glottal Tract
  % Configuration
  out.addOptional('sex','male',@(x)any(validatestring(x,{'male','female'})));
  out.addOptional('Act',0.04,@isnumeric); % Activation of the Cricothyroid [-]
  out.addOptional('Ata',0.18,@isnumeric); % Activation of the Thyroarytenoid [-]
  out.addOptional('Alc',0.5,@isnumeric); % Activation of the Lateral Cricoarytenoid [-]
  out.addOptional('rotation',1,@isnumeric); % Rotation of the arythenoids [º]
  out.addOptional('displacement',0,@isnumeric); % displacement of the arythenoids [m]
  out.addOptional('constants','galindo2016',@isstr); % set of constants to use

  out.addOptional('noise',0,@isnumeric); % state noise [x_u x_l x_b v_u v_l v_b]';
  out.addOptional('stochastic',0,@islogical); % state noise [x_u x_l x_b v_u v_l v_b]';
  
  out.parse(varargin{:});
  
  if out.Results.stochastic
    auxPos = find(strcmpi(varargin,'act'));
    if isempty(auxPos)
      varargin = [varargin {'Act' min(1,max(0,out.Results.Act + normrnd(0,0.03)))}];
    else
      varargin{auxPos+1} = min(1,max(0,out.Results.Act + normrnd(0,0.03)));
    end
    
    auxPos = find(strcmpi(varargin,'ata'));
    if isempty(auxPos)
      varargin = [varargin {'Ata' min(1,max(0,out.Results.Ata + normrnd(0,0.22)))}];
    else
      varargin{auxPos+1} = min(1,max(0,out.Results.Ata + normrnd(0,0.22)));
    end
    
    auxPos = find(strcmpi(varargin,'ps'));
    if isempty(auxPos)
      varargin = [varargin {'ps' max(0,out.Results.ps + normrnd(0,0.99))}];
    else
      varargin{auxPos+1} = max(0,out.Results.ps + normrnd(0,0.99));
    end
    
    auxPos = find(strcmpi(varargin,'rotation'));
    if isempty(auxPos)
      varargin = [varargin {'rotation' out.Results.rotation + normrnd(0,0.4)}];
    else
      varargin{auxPos+1} = out.Results.rotation + normrnd(0,0.4);
    end
    
    out.parse(varargin{:});
  end
  
end