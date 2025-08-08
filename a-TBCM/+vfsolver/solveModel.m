function x = solveModel(Solver,Model,x,Pe,Ps,constKin)
  %x = [C.x_u C.x_l C.x_b C.v_u C.v_l C.v_b]';
  switch upper(Solver)
    case 'TTS'
      switch upper(Model)
        case {'BCM','BCM+'}
          x = vfsolver.modelBCM(x,Pe,Ps,constKin);
        case 'TBCM'               
          x = vfsolver.modelTBCM(x,Pe,Ps,constKin);
      end
    case 'ODE4'
      switch upper(Model)
        case {'BCM','BCM+'}
          x = vfsolver.ode4(@(t,x)vfsolver.modelBCM(x,Pe,Ps,constKin),(1:3)*constKin.Delta_t,x);
          x = x(2,:)';
          aux_x = vfsolver.modelBCM(x,Pe,Ps,constKin);
        case 'TBCM'               
          x = vfsolver.ode4(@(t,x)vfsolver.modelTBCM(x,Pe,Ps,constKin),(1:3)*constKin.Delta_t,x);
          x = x(2,:)';
          aux_x = vfsolver.modelTBCM(x,Pe,Ps,constKin);
      end
      x(7) = aux_x(4);
      x(8) = aux_x(5);
      x(9) = aux_x(6);
  end
end