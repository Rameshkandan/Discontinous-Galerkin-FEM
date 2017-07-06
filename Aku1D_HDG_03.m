function [vel,p,p_ans,t] = Aku1D_HDG_03(vel,p, FinalTime,bc,que)

% Evaluation of v and p using RK scheme 

Globals1D;
time = 0;

% Runge-Kutta residual storage for both the vectors Velocity and Pressure  
resv = zeros(Np,K); resp = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)- x(2,:)));
dt = 0.000001;
Nsteps = floor(FinalTime/dt); 

t = 0:dt:FinalTime;

p_ans = zeros(Np*K, Nsteps+1);
p_ans_filter = zeros(Np*K, Nsteps+1);

%applying filter  
 F = Filter1D(N,floor(N/2),32); 


p_ans(:,1) = p(:);

% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;        

         [rhsv, rhsp] = AkuRHS1D_HDG_03(vel,p, timelocal,bc,que);
         
        %separate computations for velocity and pressure
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resp = rk4a(INTRK)*resp + dt*rhsp;

   
        %Applying filters
          resv = F*resv;
         resp = F*resp; 
         
        vel = vel+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;

          vel = F*vel;
          p = F*p;

  
    end;
    % Increment time
    time = time+dt;
    
    
    p_ans(:,tstep+1) = p(:);  

%  enable to see the wave motion    
    plot(x,p); title('pressure'); hold on;
    plot(x(:),p_ans(:,1));
    hold off; ylim([-1 1]);xlim([0 1]);
    title('Pressure Approximation');
     xlabel('X');
    pause(0.0001)
    
    
end;
return
