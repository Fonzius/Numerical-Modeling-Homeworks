%=======================================================================================================
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================

function [Dati]=C_dati(test)

if test=="Test1"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,3],...                          
               ... % Domain bounds   
               'bc',              'R', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing 
               'c2',               1, ...
               ... % Diffusive term ...
               'T',               0.5, ...
               ... % Final time ...
               'dt',            0.01, ...
               ... % Time step
               'u0',      '0.*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*sin(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'sin(2*pi*x).*sin(2*pi*t)',...      
               ... % Definition of exact solution or lifting function Rg
               'grad_exact',     '2*pi*cos(2*pi*x).*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '2*pi*sin(2*pi*t)',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '2*pi*sin(2*pi*t)',...
               ... %  c2du/dx(L,t) = g2(t)               
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

elseif test=="Test2"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-1,1],...                          
               ... % Domain bounds   
               'bc',              'P', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing 
               'c2',               10^2/4^2, ...
               ... % Diffusive term ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.01, ...
               ... % Time step
               'u0',      'cos(4*pi*x)', ...
               ... % Initial condition u               
               'v0',      '10*pi*sin(4*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(4*pi*x-10*pi*t)',...      
               ... % Definition of exact solution or lifting function Rg
               'grad_exact',     '-4*pi*sin(4*pi*x-10*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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

elseif test=="Test3"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-5,5],...                          
               ... % Domain bounds   
               'bc',              'A', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing 
               'c2',               2, ...
               ... % Diffusive term ...
               'T',               5, ...
               ... % Final time ...
               'dt',            0.01, ...
               ... % Time step
               'u0',      'exp(-x.^2)', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution or lifting function Rg
               'grad_exact',     '0.*x.*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
           
elseif test=="HW1_3Sa"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed
               'c2',               1, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.0005, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors
               'shape',            '1 + 0.*x'...
               ...% Section function
               );         
           
elseif test=="HW1_3Sb"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing 
               'c2',               1, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '(8*x + 4).*(2*pi*sin(2*pi*x)*sin(2*pi*t))',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '(1+2*x).^2'...
               ...% Section function
               );        

elseif test=="HW1_3Sc"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing 
               'c2',               1, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '-(15*(pi^2)/2).*cos(5*pi*x).*sin(2*pi*x).*sin(2*pi*t)',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '1 - 3*sin(5*pi*x)/4'...
               ...% Section function
               ); 

elseif test=="HW1_4Sa"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M_impulse', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '1 + 0.*x'...
               ...% Section function
               ); 

elseif test=="HW1_4Sb"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M_impulse', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '(1+2*x).^2'...
               ...% Section function
               ); 

elseif test=="HW1_4Sc"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M_impulse', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '1 - 3*sin(5*pi*x)/4'...
               ...% Section function
               ); 

elseif test=="HW1_5Sa"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'R_A', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               ... % R_A = Robinson and Absorbing
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '1 + 0.*x'...
               ...% Section function
               ); 

elseif test=="HW1_5Sb"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'R_A', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               ... % R_A = Robinson and Absorbing
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '(1+2*x).^2'...
               ...% Section function
               ); 

elseif test=="HW1_5Sc"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'R_A', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed,
               ... % M_impulse = Mixed with impulsive excitation
               ... % R_A = Robinson and Absorbing
               'c2',               4, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.001, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '0*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.5.*(1 + cos(pi*(t - 0.2)/0.1)).*(1 + 0.*x)',...      
               ... % Lifting function 
               'grad_exact',      '0*x*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors 
               'shape',            '1 - 3*sin(5*pi*x)/4'...
               ...% Section function
               ); 

elseif test=="HW1_6Sa"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed
               'c2',               343^2, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.5, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors
               'shape',            '1 + 0.*x'...
               ...% Section function
               );         
elseif test=="HW1_6Sb"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed
               'c2',               343^2, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.5, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors
               'shape',            '(1+2*x).^2'...
               ...% Section function
               );         
elseif test=="HW1_6Sc"
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds   
               'bc',              'M', ...
               ... % D = Dirichlet, N = Neumann, R = Robin, 
               ... % P = Periodic, A = Absorbing, M = Mixed
               'c2',               343^2, ...
               ... % Wave velocity ...
               'T',               1, ...
               ... % Final time ...
               'dt',            0.5, ...
               ... % Time step
               'u0',      '0*x', ...
               ... % Initial condition u               
               'v0',      '2*pi*cos(2*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(2*pi*x)*sin(2*pi*t)',...      
               ... % Definition of exact solution 
               'grad_exact',     '-2*pi*sin(2*pi*x)*sin(2*pi*t)',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'neumann1',        '0.*x.*t',...
               ... % -c2du/dx(0,t) = g1(t)
               'neumann2',        '0.*x.*t',...
               ... %  c2du/dx(L,t) = g2(t)               
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
               'plot_errors',       'Y', ...
               ...% Compute Errors
               'shape',            '1 - 3*sin(5*pi*x)/4'...
               ...% Section function
               );         
end

