clear all;

% geometrical parameters of the domain:
l_x = 1;
l_y = 1;


% subdivisions of the domain:
n_el_x = 2;
n_el_y = 2;


% pde parameters:
rho = 1;
laminar_viscosity = 1;
smagorinsky_coefficient = 1;
smagorinsky_caracteristic_length = max(l_x, l_y);
force_f = @(x) [x(:, 2) > 0.25, zeros(size(x, 1), 1)];


% build the mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% problem sizes:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_u_dof = n_nodes + n_elems;
n_p_dof = n_nodes;
n_dof = 2 * n_u_dof + n_p_dof;

p_range = 2 * n_u_dof + (1:n_nodes);
u_x_range = 1:n_u_dof;
u_y_range = n_u_dof + (1:n_u_dof);
% build the degree of freedom mapping between the elements and the linear system:
dof_map = [mesh.elements, (n_nodes + 1:n_u_dof)'];

% get the elementary integrals:
ints = functional.elementary_integrals_p1bulle();

% initialise the solution with a guess to start the newton iteration:
x = ones(n_dof, 1);

u_x = x(1:n_u_dof);
u_y = x(n_u_dof + (1:n_u_dof));
p = x(2 * n_u_dof + (1:n_nodes));

mat = navierstokes2d.assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, ...
							  u_x, u_y, ...
							  rho, laminar_viscosity, ...
							  smagorinsky_coefficient, ...
							  smagorinsky_caracteristic_length);

for n = 1:1
  h = 1.e-8;
  delta = (1:n_u_dof == n)';
  d_rhs = (navierstokes2d.assemble_p1bullep1_stationary_rhs(mesh, dof_map, ints, ...
							    u_x + delta * h, u_y, p, ...
							    force_f, ...
							    rho, ...
							    laminar_viscosity, ...
							    smagorinsky_coefficient, ...
							    smagorinsky_caracteristic_length) ...
	   - navierstokes2d.assemble_p1bullep1_stationary_rhs(mesh, dof_map, ints, ...
							      u_x, u_y, p, ...
							      force_f, ...
							      rho, ...
							      laminar_viscosity, ...
							      smagorinsky_coefficient, ...
							      smagorinsky_caracteristic_length) ) ...
	  / h;

  
end
  full([d_rhs(u_x_range),  mat(u_x_range, n)])
