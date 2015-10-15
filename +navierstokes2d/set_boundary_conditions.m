function [mat, rhs] = set_boundary_conditions(mesh, mat, rhs)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;


  % normal homogeneous dirichlet conditions:
  


  % homogeneous dirichlet condition for the velocity:
  for i = 1:2
    for n = 1:length(mesh.boundary_nodes)
      mat(mesh.boundary_nodes(n), :) = 0;
      mat(mesh.boundary_nodes(n), mesh.boundary_nodes(n)) = 1;
      rhs(mesh.boundary_nodes(n)) = 0;

      mat(n_u_dof + mesh.boundary_nodes(n), :) = 0;
      mat(n_u_dof + mesh.boundary_nodes(n), n_u_dof + mesh.boundary_nodes(n)) = 1;
      rhs(n_u_dof + mesh.boundary_nodes(n)) = 0;
    end
  end

  % pin the pressure in one point:
  mat(end, :) = 0;
  mat(end, end) = 1;
  rhs(end) = 0;
