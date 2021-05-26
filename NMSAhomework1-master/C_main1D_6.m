function [errors,solutions,femregion,Dati] = C_main1D_6(TestName,nRef)
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
%    [errors,solutions,femregion,Dati] = C_main1D_4('Test1',3)
 


addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati_6(TestName);
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
% BUILD FINITE ELEMENT MATRICES and RHS
%==========================================================================

[M_nbc,A_nbc,f_nbc] = C_matrix1D(Dati,femregion);

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
[A,b,ug] = C_bound_cond1D(A_nbc,f_nbc,femregion,Dati);

%==========================================================================
% SOLVE THE LINEAR SYSTEM
%==========================================================================
n = Dati.n;
w = Dati.omega;

uh = (-A +(n^2)*(w^2)*M_nbc)\b;

%==========================================================================
% ASSIGN DIRICHLET BOUNDARY CONDITIONS -- through the lifting ug
%==========================================================================

uh = uh + ug; 

T = Dati.T;
dt = Dati.dt;
t = 0:dt:T;
omega = Dati.omega;
for i = 1:length(t)
    Dati.t = t(i);
    uh_plot = uh.*cos(omega*t(i));
    [uh_plot] = C_snapshot_1D(femregion, uh_plot, Dati);
end

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,uh);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors(Dati,femregion,solutions);
end



