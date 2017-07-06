
% Driver script for solving the 1D acoustic equations
Globals1D;

N11 = [1,2,3, 4 ];
K11 = [5,10,20,40,80];
for i = 1:size(K11,2)
     for j = 1:size(N11,2)

 N = N11(j);
% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,1.0,K11(i));

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
p = zeros(Np,K); vel = cos(pi*x);


% Solve Problem
FinalTime = 1;
[error_inf_norm_LF,error_inf_norm_HDG,t] = Aku1D_LF_HDG(vel,p,vel,p,FinalTime);

%max rror estimation
error_max_LF(i,j) = error_inf_norm_LF;
error_max_HDG(i,j) = error_inf_norm_HDG;
    end
end


%convergence plot : semilog plot : error Vs N
figure(1)
for i = 1:size(K11,2)
warning off
subplot(2,1,1)
semilogy(N11,error_max_LF(i,:));
title('Error convergence plot with LF');
legend('K=5','K=10','K=20','K=40','K=80'); 
hold on

warning off
subplot(2,1,2)
semilogy(N11,error_max_HDG(i,:))
title('Error convergence plot with HDG');
legend('K=5','K=10','K=20','K=40','K=80'); 
hold on
end
ylabel('||Error||_I_n_f_i_n_i_t_y');
xlabel('Degree of the polynomial N');

%convergence plot : loglog plot : error Vs K
figure(2)
for i = 1:size(N11,2)
warning off
subplot(2,1,1)
loglog(K11,error_max_LF(:,i));
title('Error convergence plot with LF');
legend('N=1','N=2','N=3','N=4'); 
warning off
hold on

warning off
subplot(2,1,2)
loglog(K11,error_max_HDG(:,i));
title('Error convergence plot with HDG');
legend('N=1','N=2','N=3','N=4');
warning off
hold on
end
ylabel('||Error||_I_n_f_i_n_i_t_y');
xlabel('No of elements K');


disp('ERROR CONVERGENCE LF FLUX');
construct_Table(K11,N11,error_max_LF);
disp('ERROR CONVERGENCE HDG FLUX');
construct_Table(K11,N11,error_max_HDG);



