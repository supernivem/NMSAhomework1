%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2]
%                'c2'                 % c^2 wave speed
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

function [Dati]=C_dati(test)
L=1;
if test=='Test1'
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,L],...                          
               ... % Domain bounds       
               'c2',               1, ...
               ... % Diffusive term ...
               'omega',               1, ...
               ... % Angular velocity   
               'n',               1, ...
               ... % Index of refraction 
               'g1',              '0',...
               ... % Boundary condition in x = 0
               'g2',              '0',... 
               ... % Boundary condition in x = L  
               'exact_sol',       'sin(2*pi*x)',...      
               ... % Definition of exact solution
               'grad_exact',     '2*pi*cos(2*pi*x)',...    
               ... % du/dx 
               'force',           '1*1*sin(2*pi*x)-4*pi^2*sin(2*pi*x)',...  
               ... % Forcing term %d^2u/dx^2+n^2w^2u
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [2,3,4,5],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );

end


