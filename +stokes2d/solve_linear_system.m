function [u_x, u_y, p, x] = solve_linear_system(mesh, mat, rhs)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % solve and extract the components:
  x = mat \ rhs;
  u_x = x(1:n_nodes);
  u_y = x((n_u_dof + 1):(n_u_dof + n_nodes));
  p = x((2 * n_u_dof + 1): 2 * n_u_dof + n_nodes);
    
