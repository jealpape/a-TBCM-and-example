%%
% Simulate: Function for produce the one-step simulation of the body 
% cover model of the vocal folds.
%
% Structure: SimulatePosture(MCObj,a_Act)
% BCMObj: is an object from MuscleControlModel (handle) class,
% a_Act: is the intrinsic muscle activation vector 
%        a_Act=[a_LCA, a_IA, a_PCA, a_CT, a_TA]^T.
%
% References:
% [1] I. R. Titze and E. J. Hunter, “A two-dimensional biomechanical model
%     of vocal fold posturing,” J. Acoust. Soc. Am., vol. 121, no. 4, pp.
%     2254–2260, Apr. 2007. 
% [2] I. R. Titze, The Myoelastic Aerodynamic Theory of Phonation, 1st
%     editio. National Center for Voice and Speech, 2006.
%
% Coded by Gabriel Alzamendi, February 2020.

function SimulatePosture(MCObj,a_Act)
  % Strain variables
  coef = 0.7;
  theta_a = MCObj.xData_CAJ(7);
  eps_LCA = coef * -(MCObj.LarMuscObj(1).gamma/MCObj.LarMuscObj(1).L0)*theta_a;
  eps_IA = coef * -(MCObj.LarMuscObj(2).gamma/MCObj.LarMuscObj(2).L0)*theta_a;
  eps_PCA = coef * -(MCObj.LarMuscObj(3).gamma/MCObj.LarMuscObj(3).L0)*theta_a;
  
  eps_VF = MCObj.getStrainVF;
  eps_CT = MCObj.strain_CT;

  % Simulate the dynamic of intrinsic muscle
  a_ALL = [a_Act(:); 0; 0]; % Activation from 5 intrinsic muscle + passive Ligament + passive Mucosa
  MuscStrain = [eps_LCA; eps_IA; eps_PCA; eps_CT; eps_VF; eps_VF; eps_VF];
  dMuscStrain = zeros(7,1);
  for cont_m = MCObj.LarMusIndex
    MCObj.LarMuscObj(cont_m).SimulateStress(a_ALL(cont_m), ...
                                MuscStrain(cont_m), dMuscStrain(cont_m));  
  end
  
  % Simulate the Cricoarytenoid (CAJ) Joint
  MCObj.SimulateCAJ;
  
  % Simulate the CricoTHYROID (CtJ) Joint
  MCObj.SimulateCTJ;
  
  % Incrementing time simulation index
  MCObj.n_IterCont = MCObj.n_IterCont+1;
  
  % Saving strain results for time delay
  MCObj.deps_VF = (MCObj.getStrainVF-eps_VF)/MCObj.Ts;
  MCObj.deps_CT = (MCObj.strain_CT-eps_CT)/MCObj.Ts;

end