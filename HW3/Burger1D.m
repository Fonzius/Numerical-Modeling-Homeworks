% Solution and error of the Burger's equation
%
%    u_t + u u_x = 0
%
% x \in I=[a,b]  and  t \in [0,T] with the
% initial data  u(x,0) = u_0(x) with different methods
%   Lax-Friedrichs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,SOL,L2,L_inf] = Burger1D(mu,NX,NT)

% Intervals
T = 2.5;
I = [0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE EXACT SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_axis = linspace(I(1),I(2),NX);
t_axis = linspace(0,T,NT)';

c0 = trapz(x_axis, exp(-1/(2*pi*mu)*(1-cos(pi.*x_axis))));

num = zeros(NX,NT);
denum = zeros(NX,NT);

denum(:,:) = c0;

for n = 1:NT
    cn = 2* trapz(x_axis, exp(-1/(2*pi*mu).*(1-cos(pi.*x_axis))).*cos(n.*pi*x_axis));
    num = num + (cn.*n.*exp(-n^2 .* pi^2 .* mu .* t_axis) .* sin(n.*pi.*x_axis))';
    denum = denum + (cn.*exp(-n^2 .* pi^2 .* mu.* t_axis) .* cos(n.*pi.*x_axis))';
end

% exact solution for computing the error
u = 2.*mu.*pi.*(num./denum);

dt = T/NT;
dx = (I(2)-I(1))/NX;
lambda = dt/dx;

% Initial conditions
SOL = zeros(NX+2,NT+1);


SOL(2:end-1,1) = sin(x_axis.*pi);
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
    A(2,2) = 1;
    A(end,end) = 1;
    A(end-1,end-1) = 1;

    alpha = 1 + (4*mu)/(3*dx^2) .* dt;

    for ii = 3:NX
        A(ii,ii) = alpha;
    end

    for ii = 3:NX-1
        beta  = -2/3 .* dt .* (mu/(dx^2) - (2*SOL(ii,n)-SOL(ii,n-1))./(2.*dx));
        gamma = -2/3 .* dt .* (mu/(dx^2) + (2*SOL(ii,n)-SOL(ii,n-1))./(2.*dx));

        A(ii,ii+1) = beta;
        A(ii+1,ii) = gamma;

        f(ii) = 4/3 .* SOL(ii,n) - 1/3 .* SOL(ii,n-1);
    end
    
    SOL(:,n+1) = A\f;
    
    %constant extrapolation
    % SOL(1,n+1) = SOL(2,n+1);
    % SOL(end,n+1) = SOL(end-1,n+1);

end

SOL(1,:) = [];
SOL(end,:) = [];

L2 = sqrt(dx.*sum(abs(SOL(:,NT)-u(:,NT)).^2));
L_inf = max(abs(SOL(:,NT)-u(:,NT)));











