function [L2,Linf] = computeErrors(u_ex,u,dx)
dim_ex = length(u_ex);
dim = length(u);
step = floor(dim_ex/dim);
u_ex = u_ex(1:step:end, :);

L2 = sqrt(dx.*sum(abs(u_ex(:,end)-u(:,end)).^2));
Linf = max(abs(u_ex(:,end)-u(:,end)));

end
