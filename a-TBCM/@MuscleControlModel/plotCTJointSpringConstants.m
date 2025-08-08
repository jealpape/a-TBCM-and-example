
    
    function plotCTJointSpringConstants
    % Function for plotting the translational and rotational 
    % stiffnesses in the CA joint.
      N_val = 50;  
      disp = linspace(0,2,N_val)*1e-3;
      rot = 0:0.02:0.2; N_rot = length(rot); % linspace(0,0.3,N_val);
      k_t = zeros(N_val,N_rot);
      for cont_n = 1:N_val
        for cont_r = 1:N_rot  
          k_t(cont_n,cont_r) = MuscleControlModel.CalcCTJoinTranslationalStiffness(disp(cont_n),rot(cont_r)*MuscleControlModel.h_TA);
        end
      end
      figure
      ax1 = subplot(1,2,1);
      hold on
      l1 = plot(disp*1e3,k_t,'b'); % l1.LineWidth = 1.5;
%       l2 = plot(disp*1e3,k_y); l2.LineWidth = 1.5;
%       l3 = plot(disp*1e3,82*ones(size(disp)),'--'); l3.Color = l1.Color;  l3.LineWidth = 1.5;
%       l4 = plot(disp*1e3,253*ones(size(disp)),'--'); l4.Color = l2.Color; l4.LineWidth = 1.5;
      hold off
      ax1.Box = 'on'; ax1.XGrid = 'on'; ax1.YGrid = 'on';
      ax1.XLim = [0 2]; ax1.XTick = [0:0.5:2];
      ax1.YLim = [0 6000]; ax1.YTick = [0:1000:6000];
      ax1.XLabel.String = 'Translational displacement $\Delta L_t$ in mm'; ax1.XLabel.Interpreter = 'latex';
      ax1.YLabel.String = 'Translational stiffness in N/m'; ax1.YLabel.Interpreter = 'latex';
      tex1 = text(1.1,800,'$\theta = 0.0$ rad'); tex1.Interpreter = 'latex';
      tex2 = text(0.6,4000,'$\theta = 0.2$ rad'); tex2.Interpreter = 'latex';
      
      disp = (0:0.1:1)*1e-3; N_disp = length(disp);
      rot = linspace(0,0.3,N_val); % linspace(0,0.3,N_val);
      k_r = zeros(N_val,N_disp);
      for cont_n = 1:N_val
        for cont_d = 1:N_disp  
          k_r(cont_n,cont_d) = MuscleControlModel.CalcCTJoinRotationalStiffness(disp(cont_d),rot(cont_n)*MuscleControlModel.h_TA);
        end
      end
      ax2 = subplot(1,2,2);
      hold on
      l5 = plot(rot,k_r,'b');  % l5.LineWidth = 1.5;
%       l6 = plot(rot,0.007*ones(size(rot)),'--'); l6.Color = l5.Color; l6.LineWidth = 1.5;
      hold off
      ax2.Box = 'on'; ax2.XGrid = 'on'; ax2.YGrid = 'on';
      ax2.XLim = [0 0.3]; ax2.XTick = [0:0.1:0.3];
      ax2.YLim = [0 0.6]; ax2.YTick = [0:0.1:0.6];
      ax2.XLabel.String = 'Rotational angle $\theta$ in rad'; ax2.XLabel.Interpreter = 'latex';
      ax2.YLabel.String = 'Rotational stiffness $k_r$ in Nm/rad'; ax2.YLabel.Interpreter = 'latex';
      tex3 = text(0.15,0.08,'$\Delta L_t = 0$ mm'); tex3.Interpreter = 'latex';
      tex4 = text(0.06,0.3,'$\Delta L_t = 1$ mm'); tex4.Interpreter = 'latex';
    end