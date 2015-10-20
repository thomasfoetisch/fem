function [u_x, u_y, p, x] = solve_linear_system(mesh, mat, rhs, rot_inv)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % solve and extract the components:
  x = mat \ rhs;
  x = full(rot_inv * x);

  % split each component:
  u_x = x(1:n_u_dof);
  u_y = x((n_u_dof + 1):(n_u_dof + n_u_dof));
  p = x((2 * n_u_dof + 1): 2 * n_u_dof + n_nodes);
    
