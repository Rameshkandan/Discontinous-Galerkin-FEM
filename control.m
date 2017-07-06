clear all
K = 20;
Tf = 10;
periodic = 1;
N = 11;
Cr = 0.025/N;

N1 = N+1;
y=zeros(2*K,1);
for k=1:K
    y(2*k-1) = (k-1)/K;
    y(2*k) = k/K;
end
%x = x.^0.6;

[Mk,Sk,xunit] = getUnitMatrices(N);

for k=1:K
    x((N1*k-N):N1*k) = y(2*k-1)+(y(2*k)-y(2*k-1))*xunit;
end

h = zeros(K,1);
for k=1:K
    h(k) = x(N1*k)-x(N1*k-N);
end

dt = Cr * min(h);
NT = round(Tf/dt);

disp(['Number of elements: ' num2str(K) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

% matrices
alpha = 0;                                        % 0 = upwind, 1 = central
Fk = 0.5*[1 1; -1 -1] + 0.5*(1-alpha)*[1 -1; -1 1];


M = sparse(N1*K,N1*K);
Minv = sparse(N1*K,N1*K);
S = sparse(N1*K,N1*K);
F = sparse(N1*K,N1*K);
F(1,1) = Fk(2,2);
F(N1*K,N1*K) = Fk(1,1);
for k=1:K
    M((N1*k-N):N1*k,(N1*k-N):N1*k) = h(k)*Mk;
    Minv((N1*k-N):N1*k,(N1*k-N):N1*k) = inv(h(k)*Mk);
    S((N1*k-N):N1*k,(N1*k-N):N1*k) = Sk;
    if k<K
        F(N1*k:(N1*k+1),N1*k:(N1*k+1)) = Fk;
    end
end
Fbound = zeros(N1*K,2);
if periodic == 0
    Fbound(1,1) = Fk(2,1);
    Fbound(N1*K,2) = Fk(1,2);
else
    F(1,N1*K) = Fk(2,1);
    F(N1*K,1) = Fk(1,2);
end

v=zeros(N1*K,NT+1);
% for i=1:length(x)
%     if x(i) < 0.4
%         v(i) = 0;
%     elseif x(i) < 0.5
%         v(i) = 10*(x(i)-0.4);
%     elseif x(i) < 0.6
%         v(i) = 10*(0.6-x(i));
%     else
%         v(i) = 0;
%     end
% end
%v(:,1) = sin(2*pi*(exp(2*x)-1)/(exp(2)-1));
v(:,1) = exp(sin(2*pi*x));
%v(:,1) = 2*x - 2*(x>0.5).*ones(size(x));
for n=1:NT
    ubound = [sin(2*pi*(0-n*dt)); sin(2*pi*(1-n*dt))];
    k1 = Minv*((S'-F)*v(:,n) - Fbound*ubound);
    
    ubound = [sin(2*pi*(0-(n+0.5)*dt)); sin(2*pi*(1-(n+0.5)*dt))];
    k2 = Minv*((S'-F)*(v(:,n)+0.5*dt*k1) - Fbound*ubound);
    
    k3 = Minv*((S'-F)*(v(:,n)+0.5*dt*k2) - Fbound*ubound);
    
    ubound = [sin(2*pi*(0-(n+1)*dt)); sin(2*pi*(1-(n+1)*dt))];
    k4 = Minv*((S'-F)*(v(:,n)+dt*k3) - Fbound*ubound);
    
    v(:,n+1) = v(:,n) + dt/6*(k1+2*k2+2*k3+k4);
end

figure(1)
plot(x,v(:,1),'b',x,v(:,end),'r-')
xlabel('x')
ylabel('u_h(x)')
title(['N=' num2str(N) ', K=' num2str(K) ', dt = ' num2str(dt)])
legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],'Location','NorthEast')

% compute maximum error, assuming that we went back to the initial state
u = v(:,1);
err = v(:,end) - u;
disp(['Error in maximum norm: ' num2str(norm(err,inf))])

% for n=1:10:N
%     plot(x,v(:,1),'b',x,v(:,n+1),'r-')
%     xlabel('x')
%     ylabel('u_h(x)')
%     legend('u_h(x,0)',['u_h(x,' num2str((n+1)*dt) ')'],'Location','NorthEast')
%     axis([0 1 -1.2 1.2])
%     pause(0.05)
% end