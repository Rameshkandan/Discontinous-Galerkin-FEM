function [error_inf_norm_LF,error_inf_norm_HDG,t] = Aku1D_LF_HDG(vel_HDG,p_HDG,vel,p,FinalTime)

% Evaluation of v and P with RK method 

Globals1D;
time = 0;

% Runge-Kutta residual storage for both the vectors Velocity and Pressure  
resv_LF = zeros(Np,K); resp_LF = zeros(Np,K); 
resv_HDG = zeros(Np,K); resp_HDG = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)- x(2,:)));
CFL=0.3; dt   = CFL*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

t = 0:dt:FinalTime;

% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;        
        [rhsv, rhsp] = AkuRHS1D(vel,p, timelocal);
        [rhsv1, rhsp1] = AkuRHS1D_HDG(vel_HDG,p_HDG, timelocal);
        
        %separate computations for velocity and pressure for LF and HDG Flux
        resv_LF = rk4a(INTRK)*resv_LF + dt*rhsv;
        resp_LF = rk4a(INTRK)*resp_LF + dt*rhsp;
        
        resv_HDG = rk4a(INTRK)*resv_HDG + dt*rhsv1;
        resp_HDG = rk4a(INTRK)*resp_HDG + dt*rhsp1;
        
        %result values
        vel = vel+rk4b(INTRK)*resv_LF;
        p = p+rk4b(INTRK)*resp_LF;
        
       
        vel_HDG = vel_HDG+rk4b(INTRK)*resv_LF;
        p_HDG = p_HDG+rk4b(INTRK)*resp_LF;
        

    end;

    time = time+dt;

end;

%Error estimation with Infinity norm
p_exact = zeros(Np*K,1);
p_exact = sin(pi*x(:)).*sin(pi*FinalTime);
error_inf_norm_LF = norm((p_exact - p(:)),inf);
error_inf_norm_HDG = norm((p_exact - p_HDG(:)),inf);
return
