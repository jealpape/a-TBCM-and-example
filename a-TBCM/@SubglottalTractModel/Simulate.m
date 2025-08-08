function varargout = Simulate(SGTObj,Ug_n,varargin)
% Function for simulating acoustic progation in the subglottal tract
%
% Structure: Simulate(SGTObj,Ug_n)
%            Simulate(SGTObj,Ug_n,rs_n)
%            Simulate(SGTObj,Ug_n,'PL',PL_n)
%            Simulate(SGTObj,Ug_n,'PL',PL_n,rs_n)
%            PressureWaves = Simulate(...)
%
% where
%
% SGTObj: is an object from VocalTractModel (handle) class,
% Ug_n: is the air flow (volume velocity) value for the n-th instant,
% 'PL': flag informing that the lung/tracheal pressure is assigned,
% PL_n: lung/tracheal pressure level,
% rs_n: is the supraglottal reflection coefficient (=1 by default),
% PressureWaves: is vector gathering the backward (Bn) and forward (Fn)
%                acoustic wave components and the radiated acoustic
%                pressure Pout
%                PressureWaves = [B1 B2 ... BL F1 F2 ... FL].
%
  rs_n = 1;
  if (nargin == 3)||(nargin == 5)
    if isnumeric(varargin{end})&&(abs(varargin{end})<=1)
      rs_n = varargin{end};
    else
      error('Incorrect ''rs_n'' subglottal reflection coefficient. Correct value must fullfil rs_n need to be |rs_n|<=1.')
    end
  end

  % Simulation using mode A
  if (nargin >= 2)&&(nargin <= 3) 
    switch SGTObj.solver
      case 'WRA'
      PressureWaves = SGTObj.SimulateWRAmodeA(Ug_n,rs_n);
    case 'StateSpaceWRA'
      if (SGTObj.ComputeSSWRAModAlways)||(SGTObj.n_IterCont == 0)||(nargin == 3)
        SGTObj.getStateSpaceWRAModelmodeA(rs_n);
      end
      PressureWaves = SGTObj.SimulateStateSpaceWRAmodeA(Ug_n); % Not implemented
    otherwise
      error('Incorrect model solver!')
    end
  % Simulation using mode B  
  elseif (nargin >= 4)&&(nargin <= 5) 
    
    if ~strcmp(varargin{1},'PL')||~isnumeric(varargin{2})
      error('Error in the input parameters! See function description.')
    end
    
    switch SGTObj.solver
      case 'WRA'
      PressureWaves = SGTObj.SimulateWRAmodeB(Ug_n,varargin{1},varargin{2},rs_n);
    case 'StateSpaceWRA'
      if (SGTObj.ComputeSSWRAModAlways)||(SGTObj.n_IterCont == 0)||(nargin == 5)
        SGTObj.getStateSpaceWRAModelmodeB(rs_n);
      end
      PressureWaves = SGTObj.SimulateStateSpaceWRAmodeB(Ug_n,varargin{2});
    otherwise
      error('Incorrect model solver!')
    end
  % Incorrect number of input parameters   
  else
    error('Error in the input parameters! See function description.')
  end

  % Increase time acumulator
  SGTObj.n_IterCont = SGTObj.n_IterCont+1; % Simulation time index

  % Output variables
  if nargout == 1
    varargout{1} =  PressureWaves;
  elseif (nargout>1)
    error('It is requested more output varaibles than allowed!')  
  end

end