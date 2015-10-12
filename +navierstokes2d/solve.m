function [u_x, u_y, p] = solve(mesh, force_f, rho, laminar_viscosity, smagorinsky_coefficient, smagorinsky_caracteristic_length)

% problem sizes:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_u_dof = n_nodes + n_elems;
n_p_dof = n_nodes;
n_dof = 2 * n_u_dof + n_p_dof;

n_newton_iterations = 10;

% build the degree of freedom mapping between the elements and the linear system:
dof_map = [mesh.elements, (n_nodes + 1:n_u_dof)'];


% get the elementary integrals:
ints = functional.elementary_integrals_p1bulle();


% initialise the solution with a guess to start the newton iteration:
x = navierstokes2d.build_newton_initial_solution(mesh);


% newton iteration:
for k = 1:n_newton_iterations

  u_x = x(1:n_u_dof);
  u_y = x(n_u_dof + (1:n_u_dof));
  p = x(2 * n_u_dof + (1:n_nodes));

  mat = navierstokes2d.assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, ...
							    u_x, u_y, ...
							    rho, ...
							    laminar_viscosity, ...
							    smagorinsky_coefficient, ...
							    smagorinsky_caracteristic_length);

  rhs = navierstokes2d.assemble_p1bullep1_stationary_rhs(mesh, dof_map, ints, ...
							 u_x, u_y, p, ...
							 force_f, ...
							 rho, ...
							 laminar_viscosity, ...
							 smagorinsky_coefficient, ...
							 smagorinsky_caracteristic_length);
  [mat, rhs] = navierstokes2d.set_boundary_conditions(mesh, mat, rhs);
  
  delta_x = mat \ rhs;
  printf('| delta_x | = %f\n', sqrt(sum(delta_x.^2)));
  x = x - delta_x(1:(2 * n_u_dof + n_p_dof));
end

u_x = u_x(1:n_nodes);
u_y = u_y(1:n_nodes);
