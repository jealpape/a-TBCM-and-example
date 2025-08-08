function plotCAJointSpringConstants
% Function for plotting the translational and rotational 
% stiffnesses in the CA joint.
  N_val = 50;  
  disp = linspace(0,4,N_val)*1e-3;
  rot = linspace(0,0.5,N_val);
  k_x = zeros(N_val,1);
  k_y = zeros(N_val,1);
  kappa = zeros(N_val,1);
  for cont_n = 1:N_val
    [k_x(cont_n), k_y(cont_n)] = MuscleControlModel.CalcCAJoinTranslationalStiffness(disp(cont_n),disp(cont_n));
    kappa(cont_n) = MuscleControlModel.CalcCAJoinRotationalStiffness(rot(cont_n));
  end
  figure
  ax1 = subplot(1,2,1);
  hold on
  l1 = plot(disp*1e3,k_x); l1.LineWidth = 1.5;
  l2 = plot(disp*1e3,k_y); l2.LineWidth = 1.5;
  l3 = plot(disp*1e3,82*ones(size(disp)),'--'); l3.Color = l1.Color;  l3.LineWidth = 1.5;
  l4 = plot(disp*1e3,253*ones(size(disp)),'--'); l4.Color = l2.Color; l4.LineWidth = 1.5;
  hold off
  ax1.Box = 'on'; ax1.XGrid = 'on'; ax1.YGrid = 'on';
  ax1.XLim = [0 4]; ax1.XTick = [0:1:4];
  ax1.YLim = [0 1500]; ax1.YTick = [0:500:1500];
  ax1.XLabel.String = 'Translational displacement in mm'; ax1.XLabel.Interpreter = 'latex';
  ax1.YLabel.String = 'Translational stiffness in N/m'; ax1.YLabel.Interpreter = 'latex';
  tex1 = text(3.2,350,'$k_x$'); tex1.Interpreter = 'latex';
  tex2 = text(2.6,600,'$k_y$'); tex2.Interpreter = 'latex';

  ax2 = subplot(1,2,2);
  hold on
  l5 = plot(rot,kappa);  l5.LineWidth = 1.5;
  l6 = plot(rot,0.007*ones(size(rot)),'--'); l6.Color = l5.Color; l6.LineWidth = 1.5;
  hold off
  ax2.Box = 'on'; ax2.XGrid = 'on'; ax2.YGrid = 'on';
  ax2.XLim = [0 0.5]; ax2.XTick = [0:0.1:0.5];
  ax2.YLim = [0 0.03]; ax2.YTick = [0:0.01:0.03];
  ax2.XLabel.String = 'Rotational angle in rad'; ax2.XLabel.Interpreter = 'latex';
  ax2.YLabel.String = 'Rotational stiffness in Nm/rad'; ax2.YLabel.Interpreter = 'latex';
  tex1 = text(0.35,0.017,'$\kappa$'); tex1.Interpreter = 'latex';
end
    