%% PROBLEM PARAMETERS
% geometric dimensions of the domain:
l_x = 1;
l_y = 3/10;

% subdivision of the domain:
n_el_x = 20;
n_el_y = 6;

% initial rotation of the domain:
theta = 0;


% pde parameters:
mu = 0; % diffusion parameter:
delta = 0.5; % coefficient of the GLS stabilisation
tau = 0.0125; % timestep of the time discretisation
w_f = @(x) 0.5 * [ones(size(x, 1), 1), zeros(size(x, 1), 1)]; % convection
rhs_f = @(x) zeros(size(x, 1), 1); 
dbc_f = @(x) ones(size(x, 1));
u_init_f = @(x) 0 * ones(size(x, 1), 1);
u_exact_f = @(x, t) x(:, 1) < 0.5 * t;

%% INITIALIZE ENVIRONMENT
% elementary integrals on the reference triangle:
ints = functional.elementary_integrals_p1();

% build a square mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, theta);
inflow_nodes = geometry.input_flow_surface_nodes(mesh, w_f);


%% MATRIX ASSEMBLY
mat = convection2d.assemble_p1_convection_evolutive_matrix(mesh, ints, w_f, mu, delta, tau);
mat = convection2d.set_mat_boundary_conditions(mesh, mat, dbc_f, inflow_nodes);


%% INITIAL CONDITION
u = u_init_f(mesh.nodes);
u = convection2d.set_rhs_boundary_conditions(mesh, u, dbc_f, inflow_nodes);

sol = {u};
%% TIME EVOLUTION
for k = 1:80
  printf('time step #%i\n', k); fflush(stdout);
  rhs = convection2d.assemble_p1_convection_evolutive_rhs(mesh, ints, rhs_f, delta, w_f, u, tau);
  rhs = convection2d.set_rhs_boundary_conditions(mesh, rhs, dbc_f, inflow_nodes);

  %% LINEAR PROBLEM
  u = linsolve(mat, rhs);
  sol{end + 1} = u;
end
