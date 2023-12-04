clc
clear all
close all
% Intervals
T = 2.5;
I = [0 1];

mu = [0.1,1,10];
dt = 0.0005;
dx = 0.0005;
    
NT_ex = round(2.5/dt);
NX_ex = round(1/dx);
x_ax_ex = linspace(0,1,NX_ex);

NT = [500, 2000];
NX = [500, 2000];

for m = 1:length(mu)
    figure(m)
    %compute reference solution
    fprintf('computing for mu = %.2f \n', mu(m));
    u_ex = computeBurgerSol(mu(m),NX_ex,NT_ex); 
    subplot(3,2,[1,2])
    plot(x_ax_ex,u_ex(:,1));
    hold on;
    plot(x_ax_ex,u_ex(:,2));
    hold on;
    plot(x_ax_ex,u_ex(:,3));
    title('Reference Solution: h = dt = 0.0005');

    %compute solutions for other NX,NT
    for t = 1:length(NT)
        for x = 1:length(NX)
            fprintf("NT: %i, NX: %i \n",[NT(t),NX(x)]);
            u = computeBurgerSol(mu(m),NX(x),NT(t));
            h = 1/NX(x);
            dt_ = 2.5/NT(t);
            x_ax = linspace(0,1,NX(x));

            subplot(3,2,3+(2*(t-1))+(x-1))
            plot(x_ax,u(:,1));
            hold on;
            plot(x_ax,u(:,2));
            hold on;
            plot(x_ax,u(:,3));
            title(['h = ',num2str(h), ' dt = ', num2str(dt_)]);
        end
    end

end


%% Compute the L2 and L∞ norm of the error at T = 2.5 
% for different values of the space h and time step ∆t

clc
close all
clear all

mu = [0.1,1,10];
dt = 0.0005;
dx = 0.0005;
    
NT_ex = round(2.5/dt);
NX_ex = round(1/dx);

mu = [0.1,1,10];

NX = [50,100,200,400,500,1000,2000];
NT = [500,2000];

for m = 1:length(mu)

    u_ex = computeBurgerSol(mu(m),NX_ex,NT_ex);
    
    L2   = zeros(length(NX),length(NT));
    Linf = zeros(length(NX),length(NT));
    for n = 1:length(NT)
        for i = 1:length(NX)
            u = computeBurgerSol(mu(m),NX(i),NT(n));
            dx = 1/NX(i);
            [L2(i,n),Linf(i,n)] = computeErrors(u_ex,u,dx);
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

    
    












