function x = build_newton_initial_solution(mesh)
  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  x = zeros(n_dof, 1);

  x(1:n_nodes) = mesh.nodes(:, 1);
  x(n_u_dof + (1:n_nodes)) = mesh.nodes(:, 2);
  
  x(n_nodes + (1:n_elems)) = 1;
  x(n_u_dof + n_nodes + (1:n_elems)) = 1;
