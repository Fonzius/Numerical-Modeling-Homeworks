
%% Point 3
% For μ = 10, 1, 0.1 plot the finite difference solution at time
% 0.1,0.5,2.3

clc
close all
clear all

NT = 1000;

time = linspace(0,2.5,NT);
t_index = zeros(3,1);
[~,t_index(1)] = min(abs(time-0.1));
[~,t_index(2)] = min(abs(time-0.5));
[~,t_index(3)] = min(abs(time-2.3));

mu = [10,1,0.1];
NX = 1000;

h = 1/NX;
dt = 2.5/NT;

x = linspace(0,1,NX);
figure()
% title(["Results for  h=",num2str(h)," and dt=",num2str(dt)]);
for ii = 1:3
    [u_ex,u] = Burger1D(mu(ii),NX,NT);
    subplot(3,2,2*ii-1)
    plot(x,u_ex(:,t_index(1)));
    hold on;
    plot(x,u_ex(:,t_index(2)));
    hold on;
    plot(x,u_ex(:,t_index(3)));
    legend("t=0.1s","t=0.5s","t=2.3s");
    xlabel(['Exact Solution for \mu = ',num2str(mu(ii))]);
    
    subplot(3,2,2*ii)
    plot(x,u(:,t_index(1)));
    hold on;
    plot(x,u(:,t_index(2)));
    hold on;
    plot(x,u(:,t_index(3)));
    legend("t=0.1s","t=0.5s","t=2.3s");
    xlabel(['Numerical Solution for \mu = ',num2str(mu(ii))]);
end



%% Compute the L2 and L∞ norm of the error at T = 2.5 
% for different values of the space h and time step ∆t

clc
close all
clear all

mu = [0.1,10];

NX = [250,300,350,400,450,500,550,600,650,700,850,900,950,1000,1250,1500,1750,2000,2500,3000,4000];
NT = [500,2000];

for m = 1:length(mu)
    L2   = zeros(length(NX),length(NT));
    Linf = zeros(length(NX),length(NT));
    for n = 1:length(NT)
        for i = 1:length(NX)
            [~,~,L2(i,n),Linf(i,n)] = Burger1D(mu(m),NX(i),NT(n));
            fprintf("NT: %i, NX: %i \n",[NT(n),NX(i)]);
        end
    end
    
    h = 1./NX;
    dt = 2.5./NT;
    figure(m)
    for n = 1:length(NT)
        subplot(length(NT),2,2*n-1);
        loglog(L2(:,n));
        xlabel(['L2 error for dt=',num2str(dt(n)),', \mu=',num2str(mu(m))]);
        
        subplot(length(NT),2,2*n);
        loglog(Linf(:,n));
        xlabel(['L∞ error for dt=',num2str(dt(n)),', \mu=',num2str(mu(m))]);
    end
end

    
    



