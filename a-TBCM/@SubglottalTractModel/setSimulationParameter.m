%%
% setSimulationParameter: Function for setting the simulation parameters.
%
% Structure 1: setSimulationParameter(SGTObj,fs)
% where
% SGTObj: is an object from SubglottalTractModel (handle) class,
% fs: is the sampling frequency in Hertz (real numbre, >=1).
%
% Structure 2: setSimulationParameter(SGTObj,Ts)
% where
% BCMObj: is an object from BodyCoverModel (handle) class,
% Ts: is the sampling period in seconds (real numbre, <1).
%
% Coded by Gabriel Alzamendi, Januari 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function setSimulationParameter(SGTObj,Param)
  if isnumeric(Param)&&(Param>0)
    if (Param>=1) % Param is sampling frequency in Hertz
      SGTObj.fs = Param;
      SGTObj.Ts = 1/Param;
    else % Param is sampling period in seconds
      SGTObj.Ts = Param;
      SGTObj.fs = 1/Param;
    end
    SGTObj.SimParamOK = true;
  else
    error('Incorrect input arguments! See function structures.')
  end
end