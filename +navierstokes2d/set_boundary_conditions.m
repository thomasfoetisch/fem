function [mat, rhs, T_inv] = set_boundary_conditions(mesh, mat, rhs)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;


  % normal homogeneous dirichlet conditions:
  T = spalloc(size(mat, 1), size(mat, 1), 0);
  T_inv = spalloc(size(mat, 1), size(mat, 1), 0);

  T(sub2ind(size(T), 1:n_dof, 1:n_dof)) = 1;
  T_inv(sub2ind(size(T), 1:n_dof, 1:n_dof)) = 1;
  
  for node = 1:size(mesh.boundary_nodes, 1);
    local_t = [mesh.normals_nodes(node, :)', mesh.tangents_nodes(node, :)'];
    local_t_inv = inv(local_t);

    node_id = mesh.boundary_nodes(node);
    
    T(node_id, node_id) = local_t(1, 1);
    T(node_id, n_u_dof + node_id) = local_t(1, 2);
    T(n_u_dof + node_id, node_id) = local_t(2, 1);
    T(n_u_dof + node_id, n_u_dof + node_id) = local_t(2, 2);

    T_inv(node_id, node_id) = local_t_inv(1, 1);
    T_inv(node_id, n_u_dof + node_id) = local_t_inv(1, 2);
    T_inv(n_u_dof + node_id, node_id) = local_t_inv(2, 1);
    T_inv(n_u_dof + node_id, n_u_dof + node_id) = local_t_inv(2, 2);
  end

  % change variables on the degrees of freedom in the matrix:
  mat = T_inv*mat*T;

  % the rhs need to be rotated too in order to have the basis function for the solution space and the test space syncronized.
  rhs = T_inv*rhs;

  
  % homogeneous dirichlet condition for the velocity:
  for n = 1:length(mesh.boundary_nodes)
    % normal component:
    mat(mesh.boundary_nodes(n), :) = 0;
    mat(mesh.boundary_nodes(n), mesh.boundary_nodes(n)) = 1;
    rhs(mesh.boundary_nodes(n)) = 0;

    % tangent component:
    %mat(n_u_dof + mesh.boundary_nodes(n), :) = 0;
    %mat(n_u_dof + mesh.boundary_nodes(n), n_u_dof + mesh.boundary_nodes(n)) = 1;
    %rhs(n_u_dof + mesh.boundary_nodes(n)) = 0;
  end

  % homogeneous dirichlet condition on the corners;
  for n = 1:length(mesh.corner_nodes)
   % normal component:
   mat(mesh.corner_nodes(n), :) = 0;
   mat(mesh.corner_nodes(n), mesh.corner_nodes(n)) = 1;
   rhs(mesh.corner_nodes(n)) = 0;

   % tangent component:
   mat(n_u_dof + mesh.corner_nodes(n), :) = 0;
   mat(n_u_dof + mesh.corner_nodes(n), n_u_dof + mesh.corner_nodes(n)) = 1;
   rhs(n_u_dof + mesh.corner_nodes(n)) = 0;
  end

  % pin the pressure in one point:
  mat(end, :) = 0;
  mat(end, end) = 1;
  rhs(end) = 0;
