%%
% CollisionForces: Function implementing the computation of the 
% collision forces produced during closure phase in accordance to the  
% triangular body-cover model of the vocal folds.
%
% Structure: [Fu_col,Fl_col] = TotalCollisionForces(TBCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fu_col: is the collision force in the upper mass,
% Fl_col: is the collision force in the lower mass.
%
% Coded by Jesús Parra, March 2021.

function P_col = CollisionPressure(TBCMobj)
  
  [Fu,Fl] = CollisionForces(TBCMobj);  
  F = Fu + Fl;
  if(F == 0)
      P_col = 0;
  else
      P_col = F/TBCMobj.acont;
  end  
   
end

