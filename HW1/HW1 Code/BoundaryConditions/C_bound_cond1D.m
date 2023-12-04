function [A,b,u_g] = C_bound_cond1D(A,b,femregion,Dati)
%% [A,b,u_g] = C_bound_cond1D(A,b,femregion,Dati)
%==========================================================================
% Assign Dirchlet boundary conditions
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          femregion   : (struct)  see C_create_femregion.m
%          Dati        : (struct)  see C_dati.m

%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          u_g         : (sparse(ndof,1) real) evaluation of Dirichlet conditions 
%

%fprintf('============================================================\n')
%fprintf('Assign Dirichlet boundary conditions ... \n');
%fprintf('============================================================\n')


ndof = length(b);
u_g = sparse(ndof,1);

if(strcmp(Dati.bc,'D'))
    boundary_points = femregion.boundary_points;
elseif(strcmp(Dati.bc,'M') || strcmp(Dati.bc,'M_impulse') || strcmp(Dati.bc,'R_A'))
    boundary_points = femregion.boundary_points(1);
end

x = femregion.dof(boundary_points,1);
t = Dati.t;

if(strcmp(Dati.bc,'D') || strcmp(Dati.bc,'M'))
    u_g(boundary_points) = eval(Dati.exact_sol); % Compute the lifting operator ug
elseif(strcmp(Dati.bc,'M_impulse') || strcmp(Dati.bc,'R_A'))
    if( abs(t-0.2) <= 0.1)
        u_g(boundary_points) = eval(Dati.exact_sol); % Compute the lifting operator ug
    end
end

x_g = sparse(ndof,1);
A_0 = A;

b_0 = b - A*u_g; % modify the load vector --> F(v) = F(v) - a(ug,v)

% Reduce the system A in order to solve the pb with 
% homogeneous Dirichlet conditions 
for k = 1:length(boundary_points)
    A_0(boundary_points(k),:) = 0;
    A_0(:,boundary_points(k)) = 0;
    A_0(boundary_points(k),boundary_points(k)) = 1;
    b_0(boundary_points(k)) = 0;
end

b = b_0;
A = A_0;
