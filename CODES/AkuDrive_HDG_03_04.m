clear all;
warning off;
% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation 
N = 10; 

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,1.0,15);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
vel = zeros(Np,K); p = exp(-((x-0.5).^2)/(0.02)^2);

% Solve Problem
FinalTime = 0.003;

% ABSORBING = absorbing boundary conditions
% DIRICHLET = dirichlet boundary conditions

% last parameter for the material properties for different questions
% que =  3 for question 3 enables homogeneous material properties
% que = 4 for heterogeneous layers
bound_Cond = 'ABSORBING';
que=4;
[vel,p,p_ans,t] = Aku1D_HDG_03(vel,p,FinalTime,bound_Cond,que);

% time vector
t_ref = 0:0.0006:0.003;

%plotting graph at different time intervals
% q=1;
% for i = 1:size(t,2)
%     for j = 1:size(t_ref,2)
%      if t(i) == t_ref(j)
% %        figure(q);
%          plot(x(:),p_ans(:,i));
%          hold on
% %         title(['Pressure at t=' num2str(t(i))])
% %         xlabel('X');
% %         ylabel('Pressure P');
% %        % if q==1
% %             ylim([-1 1]);
% %         %else
% %           %  ylim([-0.5 0.5]);
% %         %end
% %             xlim([0 1]);
% %         grid on
% %         hold off
% %         q = q+1;
%      end
%     end
% end
% %
% legend('t=0','t = 0.0006','t = 0.0012','t = 0.0016','t = 0.0024','t =0.003'); %plot(t, p_ans(55,:))
% xlabel('X');
% ylabel('Pressure P');
% if bound_Cond == 'ABSORBING' 
%     if que == 3
%     title('Plot with absorbing boundary conditions in homogenous material N = 10, K = 15');
%  else  
% title('Plot with absorbing boundary conditions in heterogenous material N = 10, K = 15');
%     end    
% else
%  if   que == 3
% title('Plot with dirichlet boundary conditions in homogenous material N = 10, K = 15');
% else
% title('Plot with dirichlet boundary conditions in heterogenous material N = 10, K = 15');
%  end
% end
