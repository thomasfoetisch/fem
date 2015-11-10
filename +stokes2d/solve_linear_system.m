function [u_x, u_y, p, x] = solve_linear_system(mesh, mat, rhs, rot_inv, lifting)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % solve and extract the components:
  %x = mat \ rhs;
  [L, U] = ilu(mat, struct('type', 'ilutp', 'droptol', 1e-1, 'udiag', 1));
  x = gmres(mat, rhs, 1000, 1e-10, 1000, L, U);
  x(1:(2*n_u_dof)) = x(1:(2*n_u_dof)) + lifting;
  x = full(rot_inv * x);

  % split each component:
  u_x = x(1:n_u_dof);
  u_y = x((n_u_dof + 1):(n_u_dof + n_u_dof));
  p = x((2 * n_u_dof + 1): 2 * n_u_dof + n_nodes);
    
