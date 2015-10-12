function [mat, rhs] = set_boundary_conditions(mesh, mat, rhs)

  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  constraint = spalloc(2 * length(mesh.boundary_nodes) + 1, size(rhs, 1), 0);
  for i = 1:2
    for n = 1:length(mesh.boundary_nodes)
      constraint((i - 1) * length(mesh.boundary_nodes) + n, (i - 1) * n_u_dof + mesh.boundary_nodes(n)) = 1;
    end
  end
  constraint(end, end) = 1;

  mat = [mat, constraint'; constraint, zeros(size(constraint, 1))];
  rhs = [rhs; zeros(size(constraint, 1), 1)];
