function [errors,solutions,femregion,Dati] = C_main1D_5b(TestName,nRef)
%==========================================================================
% Solution of the Wave Equation with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [errors,solutions,femregion,Dati] = C_main1D('Test1',3)



addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati_5b(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[M_nbc,A_nbc] = C_matrix1D_lf(Dati,femregion);

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
Dati.t = 0;
[b_nbc] = C_rhs1D_lf(Dati,femregion);

%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================
x = femregion.coord;
u0 = eval(Dati.u0);
v0 = eval(Dati.v0);

%% First step of leapfrog ...
% 1) Compute the rhs
b_nbc = (M_nbc - 0.5*Dati.dt^2*A_nbc)*u0 + Dati.dt*M_nbc*v0 + 0.5*Dati.dt^2 *b_nbc;

% 2) Dirichlet boundary conditions
[M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
u1 = M\b;
u1 = u1 + u_g;

% 5) Plot the obtained solution --> C_snapshot_1D
[u1] = C_snapshot_1D(femregion, u1, Dati);

%% End of first step leapfrog

fprintf('============================================================\n')
fprintf('Starting time-loop ... \n');
fprintf('============================================================\n')



for t = Dati.dt : Dati.dt : Dati.T - Dati.dt
    
    fprintf('time = %5.3e \n',t);
    
    %==========================================================================
    % BUILD FINITE ELEMENTS RHS a time t
    %==========================================================================
    Dati.t = t;
    [b_nbc] = C_rhs1D_lf(Dati,femregion);
    
    % Repeat steps 1) to 5) for the general time step
    b_nbc = (2*M_nbc - Dati.dt^2*A_nbc)*u1 - M_nbc*u0 + Dati.dt^2 *b_nbc;
    [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
    u2 = M\b;
    u2 = u2 + u_g;
    [u2] = C_snapshot_1D(femregion, u2, Dati);
    
    % Put a pause between one step and the other to see the plot
    %pause(0.01);
    
    % Update the solution
    u0 = u1;
    u1 = u2;
    %hold on
end

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing_lf(Dati,femregion,u2);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors_lf(Dati,femregion,solutions);
end


