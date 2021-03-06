function [rhsv, rhsp] = AkuRHS1D_HDG(vel,p, timelocal)
%Evaluates RHS flux in 1D wave equation

Globals1D;

%material properties 
rho = ones(Np,K); c = ones(Np,K);

 p_dash = (1./rho).*p;
 v_dash = rho.*(c.^2).*vel;
 p_c = c.*p;
 v_c = c.*vel;
 
%HDG flux
dp = zeros(Nfp*Nfaces,K); 
dp(:) = ((p_dash(vmapM)+p_dash(vmapP))/2).*nx(:) + ((v_c(vmapM)-v_c(vmapP))/2).*nx(:).*nx(:);
du = zeros(Nfp*Nfaces,K); 
du(:) = ((v_dash(vmapM)+vel(vmapP))/2).*nx(:) + (p_c(vmapM)-p_c(vmapP));

%Boundary condition %Dirichlet condition
%do nothing with velocity
dp (mapI) = 0;
dp (mapO) = 0;

du (mapI) = (v_dash(vmapI)).*nx(mapI) + (p_c(vmapI));
du (mapO) = (v_dash(vmapO)).*nx(mapO) + (p_c(vmapO));
     

% compute right hand sides of the semi-discrete PDE
rhsv = rx.*(M_inv*S_trans*p_dash) - LIFT*(Fscale.*dp);
rhsp = rx.*(M_inv*S_trans*v_dash) - LIFT*(Fscale.*du);

% rhsv = rx.*(G * p_dash) - LIFT*(Fscale.*dp);
% rhsp = rx.*(G * v_dash) - LIFT*(Fscale.*du);


return
