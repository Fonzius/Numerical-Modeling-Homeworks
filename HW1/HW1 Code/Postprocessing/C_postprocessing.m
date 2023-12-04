function [solutions]=C_postprocessing(Dati,femregion,uh)
%% [solutions]=C_postprocessing(Dati,femregion,uh)
%==========================================================================
% POST PROCESSING OF THE SOLUTION
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%
%    OUTPUT:
%          solutions   : (struct) containg solution vector uh and
%                        analytical solution u_ex
%

fprintf('============================================================\n')
fprintf('Post-processing the solution ... \n');
fprintf('============================================================\n')


%==========================================================================
% EVALUATION OF THE EXACT SOLUTION
%==========================================================================

[u_ex] = C_eval_exact_sol(femregion,Dati.exact_sol,Dati.T);
if(Dati.name == "HW1_4Sa" || Dati.name == "HW1_4Sb" || Dati.name == "HW1_4Sc" ...
   || Dati.name == "HW1_5Sa" || Dati.name == "HW1_5Sb" || Dati.name == "HW1_5Sc")
    u_ex = u_ex - 1;
end

%==========================================================================
% PLOT SOLUTION
%==========================================================================

if(Dati.visual_graph == 'Y')
    C_pointwise_sol(femregion, uh, u_ex, Dati);
end

%==========================================================================
% SAVE SOLUTIONS
%==========================================================================

solutions = struct('u_ex',u_ex,'uh',uh);
