function [u_x, u_y, p] = solve(mesh, force_f, laminar_viscosity, smagorinsky_coefficient, smagorinsky_caracteristic_length)

% problem sizes:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_u_dof = n_nodes + n_elems;
n_p_dof = n_nodes;
n_dof = 2 * n_u_dof + n_p_dof;

n_newton_iterations = 0;

% build the degree of freedom mapping between the elements and the linear system:
dof_map = [mesh.elements, (n_nodes + 1:n_u_dof)'];


% get the elementary integrals:
ints = functional.elementary_integrals_p1bulle();


% initialise the solution with a guess to start the newton iteration:
x = navierstokes2d.build_newton_initial_solution(mesh);


% newton iteration:
for k = 1:n_newton_iterations
  mat = navierstokes2d.assemble_p1bullep1_stationary_matrix();
  rhs = navierstokes2d.assemble_p1bullep1_stationary_rhs();

  [mat, rhs] = navierstokes2d.set_boundary_conditions(mat, rhs);

  delta_x = mat \ rhs;
  x = x - delta_x;
end


%u_x = x(1:n_nodes);
%u_y = x((n_u_dof + 1):(n_u_dof + n_nodes));
%p = x((2*n_u_dof + 1):(2*n_u_dof + n_nodes));

% dummy output for successful execution:
u_x = zeros(n_nodes, 1);
u_y = zeros(n_nodes, 1);
p = zeros(n_nodes, 1);
