%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2]
%                'bc'                 % boundary conditions DD-Dirichlet,
%                                     % NN-Neumann, PP-peridic,
%                                     % AA-absorbing
%                'c2'                 % c^2 wave speed
%                'T'                  % final time
%                'dt'                 % time step 
%                'u0',                % Initial condition u               
%                'v0',                % Initial condition  du/dt        
%                'exact_sol',         % set the exact solution
%                'force',             % set the forcing term
%                'grad_exact_1',      % set the first componenet of the gradient of the exact solution
%                'fem',               % set finite element space
%                'nqn_1D',            % number of quadrature nodes for integrals over lines
%                'refinement_vector', % set the level of refinement for the grid
%                'visual_graph',      % if you want to display the graphical results ['Y','N']
%                'print_out',         % if you want to print out the results ['Y','N']
%                'plot_errors'        % you want to print the computed errors ['Y','N']
% 
%========================================================================================================
%
%  REFRACTION
%
% Test 1 - wave propagation with analytucal solution
% Plane wave with dirichlet condition - verification purpose
%
%*************************************************************************
%
% Test 2 - wave propagation rightward - heterogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = (sin(2*pi*t) - x.*0.5*sin(2*pi*t)) t<= 0.5
% u(L,t) = 0 t>=0
%
%          | 2 x <= 1;
% c^2(x) = |   
%          | 1 x > 1;
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test2',6)
%
%
%==========================================================================
%
%  GIBBS - PHENOMENON
%
% Test 3 - Triangular wave - homogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = u(L,t) = 0 t>=0
%
% c^2(x) = 1; 
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test3',6)
%
%*************************************************************************
%
% Test 4 - Square wave - homogeneous bar 
% Null initial conditions - null force
% Dirichlet conditions 
% u(0,t) = u(L,t) = 0 t>=0
%
% c^2(x) = 1; 
%
% Example of usage:
% [errors,solutions,femregion,Dati] = C_main1D('Test4',6)
%
%
%==========================================================================



function [Dati]=C_dati(test)

if strcmp(test,'Test1')==1
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds       
               'bc',           'DD',...                          
               ... % boundary conditions                      
               'c2',             '1+0.*x', ...
               ... % wave speed ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0.*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*sin(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'sin(2*pi*x).*sin(2*pi*t)',...      
               ... % Definition of exact solution
               'grad_exact',     '2*pi*cos(2*pi*x).*sin(2*pi*t)',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P2',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );

elseif strcmp(test,'Test2')==1
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,2],...                          
               ... % Domain bounds       
               'bc',           'DD',...                          
               ... % boundary conditions                      
               'c2',             '1.*(x<=1) + 2.*(x>1)', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ...  %damping coefficient
               'T',               3, ...
               ... % Final time ...
               'dt',            0.005, ...
               ... % Time step
               'u0',      '0.*x', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '(sin(2*pi*t) - x.*0.5*sin(2*pi*t)).*(t<=0.5) + 0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P3',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );
           
                      

elseif strcmp(test,'Test3')==1
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-5,5],...                          
               ... % Domain bounds       
               'bc',           'DD',...                          
               ... % boundary conditions                      
               'c2',             '1 + 0.*x', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ... % wave speed ...
               'T',               2, ...
               ... % Final time ...
               'dt',            0.01, ...
               ... % Time step
               'u0',      '(x+1).*(x>=-1).*(x<0) + (1-x).*(x<=1).*(x>=0)', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P4',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );  
  
           
elseif strcmp(test,'Test4')==1
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-5,5],...                          
               ... % Domain bounds       
               'bc',           'DD',...                          
               ... % boundary conditions                      
               'c2',             '1 + 0.*x', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ... % wave speed ...
               'T',               2, ...
               ... % Final time ...
               'dt',            0.005, ...
               ... % Time step
               'u0',      '1.*(x>=-1).*(x<0) + 1.*(x<=1).*(x>=0)', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     '0.*t',...   
               ... % c2du/dx(L,t) 
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );  
             
elseif strcmp(test,'HW2ex2')==1
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-1,1],...                          
               ... % Domain bounds       
               'bc',           'M',...    %Mixed                      
               ... % boundary conditions                      
               'c2',             '1 + 0.*x', ...
               ... % wave speed ...
               'kk',             '0.*x', ...
               ... % wave speed ...
               'T',               0.3, ...
               ... % Final time ...
               'dt',            0.0005, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*sin(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'sin(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution
               'grad_exact',     '2*pi*cos(2*pi*x)*sin(2*pi*t)',...   
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',     '0.*t',...   
               ... % c2du/dx(0,t) 
               'neumann2',     'sin(2*pi*t)',...   
               ... % c2du/dx(L,t) 
               'fem',              'P5',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [4,5,6,7],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );  
end



