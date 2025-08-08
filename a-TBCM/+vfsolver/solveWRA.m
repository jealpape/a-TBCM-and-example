function [P_vector,P_input,P_output] = solveWRA(A_matrix,B_matrix,C_matrix,D_matrix,P_vector,P_input,P_output)
%   P_vector = [P-,P+]'
%   P_input = [Pn-;P0+]
%   P_output = [Pn-; Pradiated; P(n-1)+]

  % Output
  P_output = C_matrix*P_vector + D_matrix*P_output;
  % State
  P_vector = A_matrix*P_vector + B_matrix*P_input;
  % Input
  P_input = [0 0; 0 1]*P_input + [1 0 0;0 0 0]*P_output;
end