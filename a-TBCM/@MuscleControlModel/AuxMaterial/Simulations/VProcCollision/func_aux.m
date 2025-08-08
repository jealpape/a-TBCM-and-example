function f_val = func_aux(theta_a, Variables)  
  xCAJ = Variables.xCAJ;
  yCAJ = Variables.yCAJ;
  x__02 = Variables.x__02;
  xi_a = Variables.xi_a;
  
  x02 = xCAJ - (xCAJ-x__02)*cos(theta_a) + yCAJ*sin(theta_a) + xi_a;
  f_val = x02;

end