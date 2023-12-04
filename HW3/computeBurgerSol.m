function [u] = computeBurgerSol(mu,NX,NT)
% NX = 1000;
% NT = 1000;
% 
% 
dx = 1/NX;
dt = 2.5/NT;
% 
% mu = 0.1;

lambda = dt/dx;

x_axis = linspace(0,1,NX);
t_axis = linspace(0,2.5,NT)';

t_index = zeros(3,1);
[~,t_index(1)] = min(abs(t_axis-0.1));
[~,t_index(2)] = min(abs(t_axis-0.5));
[~,t_index(3)] = min(abs(t_axis-2.3));

% Initial conditions
SOL = zeros(NX+2,NT+1);

% exact solution for computing the error
SOL(2:end-1,1) = 1/4 .* cos(x_axis.*pi);
SOL(end-1,1) = 0;

for ii = 2:NX
    SOL(ii,2) = SOL(ii,1) - lambda/2 .* SOL(ii,1).*(SOL(ii+1,1) - SOL(ii-1,1))...
                + lambda*mu/dx .*(SOL(ii+1,1) - 2.*SOL(ii,1) + SOL(ii-1,1));
end
SOL(2,2) = 0;

for n = 2:NT
    
    f = zeros(NX+2,1);
    A = zeros(NX+2,NX+2);
    A(1,1) = 1;
    A(1,2) = -1;
    A(end,end) = 1;
    A(end-1,end-1) = 1;
    
    alpha = 1 + (4*mu)/(3*dx^2) .* dt;

    for ii = 2:NX
        A(ii,ii) = alpha;
    end
    
    gamma = -2/3 .* dt .* (mu/(dx^2) + (2*SOL(1,n)-SOL(1,n-1))./(2.*dx));
    A(2,1) = gamma;

    for ii = 2:NX-1
        beta  = -2/3 .* dt .* (mu/(dx^2) - (2*SOL(ii,n)-SOL(ii,n-1))./(2.*dx));
        gamma = -2/3 .* dt .* (mu/(dx^2) + (2*SOL(ii,n)-SOL(ii,n-1))./(2.*dx));

        A(ii,ii+1) = beta;
        A(ii+1,ii) = gamma;

        f(ii) = 4/3 .* SOL(ii,n) - 1/3 .* SOL(ii,n-1);
    end

    f(end-1) = 0;
    f(end) = 0;
    A(end,end) = 1;
    
    SOL(:,n+1) = A\f;    
    %constant extrapolation
    % SOL(1,n+1) = SOL(2,n+1);
    % SOL(end,n+1) = SOL(end-1,n+1);

end

for n = 1:NT
    SOL(:,n) = SOL(:,n) - 1/4 .* exp(-mu.*t_axis(n));
end

SOL(1,:) = [];
SOL(end,:) = [];

u(:,1) = SOL(:,t_index(1));
u(:,2) = SOL(:,t_index(2));
u(:,3) = SOL(:,t_index(3));
end

