%% PROBLEM PARAMETERS
% geometric dimensions of the domain:
l_x = 1;
l_y = 1;

% subdivision of the domain
n_el_x = 20;
n_el_y = 20;


% pde parameters:
mu = 1e-3; % diffusion parameter:
w_f = @(x) [ones(size(x, 1), 1), -ones(size(x, 1), 1)]; % convection
rhs_f = @(x) ones(size(x, 1), 1); 
dbc_f = @(x) zeros(size(x, 1));


%% INITIALIZE ENVIRONMENT
% elementary integrals on the reference triangle:
ints = elementary_integrals();

% build a square mesh:
mesh = build_square_mesh(l_x, l_y, n_el_x, n_el_y);


%% MATRIX ASSEMBLY
rhs = assemble_p1_convection_rhs(mesh, ints, rhs_f, 0.5, w_f);
mat = assemble_p1_convection_matrix(mesh, ints, w_f, mu, 0.5);

%inflow_nodes = input_flow_surface_nodes(mesh, w_f);
[mat, rhs] = set_boundary_conditions(mesh, mat, rhs, dbc_f, mesh.boundary_nodes);


%% LINEAR PROBLEM
u = linsolve(mat, rhs);
