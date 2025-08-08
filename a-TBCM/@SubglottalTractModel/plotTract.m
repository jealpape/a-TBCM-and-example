%%
% getSimpleVocalTract: Function to get the idealized vocal tract
% configurations and other relevant simulation parameters, and to
% configure the vocal tract model accordingly.
%
% Coded by Gabriel Alzamendi, January 2020.
function plotTract(VTobj)
  % Checking the area function describing the tract  
  if isempty(VTobj.AreaFunction)
    error(VTobj.NOTRACTERRORMSG)  
  end

  ind_data = VTobj.Delta_z*(0:VTobj.N_AreaSection);
  areadata = [VTobj.AreaFunction; VTobj.AreaFunction(end)]*1e4;
  
  stairs(ind_data,areadata)
  title('Subglottal tract configuration')
  ylabel('Area function [cm^2]'); ylim([0 ceil(1.1*max(VTobj.AreaFunction*1e4))]);
  xlabel('Distance from the glottis [m]')

end
