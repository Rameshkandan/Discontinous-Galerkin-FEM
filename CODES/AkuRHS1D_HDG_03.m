clearfunction [rhsv, rhsp] = AkuRHS1D_HDG_03(vel,p, timelocal,bc,que)
%Evaluates RHS flux in 1D wave equation

Globals1D;
 

 if que == 3
 % homeogenous material condition
 rho = 1.2*ones(Np,K); c = 340*ones(Np,K);
 else
 % heterogeneous material layers 
 
 rho1 = [0.16*ones(1,K/3),1.2*ones(1,K/3),0.16*ones(1,K/3)];
 rho = ones(Np,1)* rho1;
 
 c1 = [1000*ones(1,K/3),340*ones(1,K/3),1000*ones(1,K/3)];
 c = ones(Np,1)* c1;
 
 end
 
 p_dash = (1./rho).*p;
 v_dash = rho.*(c.^2).*vel;
 p_c = c.*p;
 v_c = c.*vel;
 
%HDG flux
dp = zeros(Nfp*Nfaces,K); 
dp(:) = ((p_dash(vmapM)+p_dash(vmapP))/2).*nx(:) + ((v_c(vmapM)-v_c(vmapP))/2).*nx(:).*nx(:);
du = zeros(Nfp*Nfaces,K); 
du(:) = ((v_dash(vmapM)+vel(vmapP))/2).*nx(:) + (p_c(vmapM)-p_c(vmapP));

if bc == 'ABSORBING'
%Boundary condition %Absorbing boundary conditions)

  dp (mapI) = v_c(vmapI).*nx(mapI).*nx(mapI); 
  dp (mapO) = v_c(vmapO).*nx(mapO).*nx(mapO); 

  du (mapI) = (v_dash(vmapI)).*nx(mapI) + (p_c(vmapI) - (v_dash(vmapI).*nx(mapI)));
  du (mapO) = (v_dash(vmapO)).*nx(mapO) + (p_c(vmapO) - (v_dash(vmapO).*nx(mapO))) ;

else
 
%Boundary condition %Dirichlet condition
%do nothing with velocity
dp (mapI) = 0;
dp (mapO) = 0;

du (mapI) = (v_dash(vmapI)).*nx(mapI) + (p_c(vmapI));
du (mapO) = (v_dash(vmapO)).*nx(mapO) + (p_c(vmapO));
     
end

% compute right hand sides of the semi-discrete PDE
rhsv = rx.*(M_inv*S_trans*p_dash) - LIFT*(Fscale.*dp);
rhsp = rx.*(M_inv*S_trans*v_dash) - LIFT*(Fscale.*du);



return
