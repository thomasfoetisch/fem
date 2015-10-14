% geometrical parameters of the domain:
l_x = 1;
l_y = 1;


% subdivisions of the domain:
n_el_x = 1;
n_el_y = 1;


% pde parameters:
ctx = struct();
ctx.rho = 1;
ctx.laminar_viscosity = 1;
ctx.smagorinsky_coefficient = 1;
ctx.smagorinsky_caracteristic_length = max(l_x, l_y);
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
[ints, basis] = functional.elementary_integrals_p1bulle();

% initialise the solution with a guess to start the newton iteration:
x = navierstokes2d.build_newton_initial_solution(mesh);

u_x = x(1:n_u_dof);
u_y = x(n_u_dof + (1:n_u_dof));
p = x(2 * n_u_dof + (1:n_nodes));


% numerical approx of the derivatine:
n = 1
h = 1.e-8;
delta = (1:n_u_dof == n)';

d_rhs = (navierstokes2d.quadrature(@(x) navierstokes2d.viscosity_term_rhs(x, mesh, dof_map, ctx, u_x + h*delta, u_y, 1, basis, 1, n)) - ...
	navierstokes2d.quadrature(@(x) navierstokes2d.viscosity_term_rhs(x, mesh, dof_map, ctx, u_x, u_y, 1, basis, 1, n))) / h
mat = navierstokes2d.quadrature(@(x) navierstokes2d.viscosity_term_1(x, mesh, dof_map, ctx, u_x, u_y, 1, basis, 1, 1, 1, n) + ...
				navierstokes2d.viscosity_term_2(x, mesh, dof_map, ctx, u_x, u_y, 1, basis, 1, 1, 1, n))
