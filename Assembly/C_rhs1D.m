function [f] = C_rhs1D(Dati,femregion)
%% [f] = C_matrix1D(Dati,femregion)
%==========================================================================
% Assembly of the rhs f at time t
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          f           : (sparse(ndof,1) real) rhs vector


addpath FESpace
addpath Assembly

% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_1D);


f = sparse(ndof,1);     % Global Load vector


for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %==============================================
    % FORCING TERM --RHS
    %==============================================

    % Local load vector
    [load] = C_loc_rhs1D(Dati.force,dphiq,BJ,w_1D,pphys_1D,nln,Dati.t);    

    % Assembly phase for the load vector
    f(iglo) = f(iglo) + load;

end
