function [rot, rot_inv] = build_dof_rotation_matrix(mesh, nodes, new_basis)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  % normal homogeneous dirichlet conditions:
  rot = spalloc(n_dof, n_dof, 0);
  rot_inv = spalloc(n_dof, n_dof, 0);

  % initialize with the identity matrix:
  for i = 1:size(rot, 1)
    rot(i, i) = 1;
    rot_inv(i, i) = 1;
  end
  
  % loop over all the selected boundary nodes:
  for node = 1:length(nodes);
    local_rot = new_basis(:, :, node);
    local_rot_inv = inv(local_rot);

    node_id = mesh.boundary_nodes(nodes(node));

    % assemble the rotation matrix and the related inverse rotation:
    rot(node_id, node_id) = local_rot(1, 1);
    rot(node_id, n_u_dof + node_id) = local_rot(1, 2);
    rot(n_u_dof + node_id, node_id) = local_rot(2, 1);
    rot(n_u_dof + node_id, n_u_dof + node_id) = local_rot(2, 2);

    rot_inv(node_id, node_id) = local_rot_inv(1, 1);
    rot_inv(node_id, n_u_dof + node_id) = local_rot_inv(1, 2);
    rot_inv(n_u_dof + node_id, node_id) = local_rot_inv(2, 1);
    rot_inv(n_u_dof + node_id, n_u_dof + node_id) = local_rot_inv(2, 2);
  end
  
