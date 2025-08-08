%%
% setSimulationParameter: Function for setting the simulation parameters.
%
% Structure 1: setSimulationParameter(VTobj,fs)
% where
% VTobj: is an object from SubglottalTractModel (handle) class,
% fs: is the sampling frequency in Hertz (real numbre, >=1).
%
% Structure 2: setSimulationParameter(VTobj,Ts)
% where
% BCMObj: is an object from BodyCoverModel (handle) class,
% Ts: is the sampling period in seconds (real numbre, <1).
%
% Coded by Gabriel Alzamendi, Januari 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function setSimulationParameter(VTobj,Param)
  if isnumeric(Param)&&(Param>0)
    if (Param>=1) % Param is sampling frequency in Hertz
      VTobj.fs = Param;
      VTobj.Ts = 1/Param;
    else % Param is sampling period in seconds
      VTobj.Ts = Param;
      VTobj.fs = 1/Param;
    end
    VTobj.SimParamOK = true;
  else
    error('Incorrect input arguments! See function structures.')
  end
end