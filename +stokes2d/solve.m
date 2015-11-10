function [u_x, u_y, p] = solve(mesh, mu, force_f, boundary_conditions)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % elementary integrals:
  ints = functional.elementary_integrals_p1bulle();

  % build the degree of freedom map to take the bubbles into account:
  dof_map = [mesh.elements, ((n_nodes + 1):n_u_dof)'];

  % build the lifting:
  [lifting_u_x, lifting_u_y] = stokes2d.build_lifting(mesh, boundary_conditions);

  % build the linear system:
  mat = stokes2d.assemble_p1bullep1_matrix(mesh, dof_map, ints, mu);

  % build the linear system rhs:
  rhs = stokes2d.assemble_p1bullep1_rhs(mesh, dof_map, ints, [lifting_u_x, lifting_u_y], force_f);

  % set the boundary condition for the velocity, and remove the kernel for the pressure:
  [mat, rhs, rot_inv] = stokes2d.set_boundary_conditions(mesh, mat, rhs, boundary_conditions);

  % solve the system and extract the different variables:
  [u_x, u_y, p] = stokes2d.solve_linear_system(mesh, mat, rhs, rot_inv, [lifting_u_x; lifting_u_y]);
