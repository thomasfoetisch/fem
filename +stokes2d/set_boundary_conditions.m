function [mat, rhs, rot_inv] = set_boundary_conditions(mesh, mat, rhs, boundary_conditions)
  % modify the matrix and the rhs of the stationary stokes 2d problem to enforce the essential boundary conditions.

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;


  nodes_to_rotate = unique(sort([boundary_conditions.u_n{1}; boundary_conditions.u_t{1}]));
  % build the new basis for each boundary nodes and associated rotation matrix:
  normal_tangent_basis = cat(2,
			     reshape(permute(mesh.normals_nodes, [2, 1]), 2, 1, []),
			     reshape(permute(mesh.tangents_nodes, [2, 1]), 2, 1, []));
  [rot, rot_inv] = stokes2d.build_dof_rotation_matrix(mesh, nodes_to_rotate, normal_tangent_basis(:, :, nodes_to_rotate));

  % change variables on the degrees of freedom in the matrix and rhs:
  mat = rot_inv * mat * rot;
  rhs = rot_inv * rhs;


  % 1. inhomogeneous dirichlet conditions on: u_x:
  if not(isempty(boundary_conditions.u_x))
     nodes = mesh.boundary_nodes(boundary_conditions.u_x{1});
     if is_function_handle(boundary_conditions.u_x{2})
       values = boundary_conditions.u_x{2}(mesh.nodes(nodes, :));
     else
       values = boundary_conditions.u_x{2};
     end
     
    for n = 1:length(nodes)
      mat(nodes(n), :) = 0;
      mat(nodes(n), nodes(n)) = 1;
      rhs(nodes(n)) = values(n);
    end
  end

  % 2. inhomogeneous dirichlet conditions on: u_y:
  if not(isempty(boundary_conditions.u_y))
     nodes = mesh.boundary_nodes(boundary_conditions.u_y{1});

     if is_function_handle(boundary_conditions.u_y{2})
       values = boundary_conditions.u_y{2}(mesh.nodes(nodes, :));
     else
       values = boundary_conditions.u_y{2};
     end
     
    for n = 1:length(nodes)
      mat(n_u_dof + nodes(n), :) = 0;
      mat(n_u_dof + nodes(n), n_u_dof + nodes(n)) = 1;
      rhs(n_u_dof + nodes(n)) = values(n);
    end
  end

  % 3. inhomogeneous dirichlet conditions on: u_n:
  if not(isempty(boundary_conditions.u_n))
     nodes = mesh.boundary_nodes(boundary_conditions.u_n{1});
     
     if is_function_handle(boundary_conditions.u_n{2})
       values = boundary_conditions.u_n{2}(mesh.nodes(nodes, :));
     else
       values = boundary_conditions.u_n{2};
     end

    for n = 1:length(nodes)
      mat(nodes(n), :) = 0;
      mat(nodes(n), nodes(n)) = 1;
      rhs(nodes(n)) = values(n);
    end
  end

  % 4. inhomogeneous dirichlet conditions on: u_t:
  if not(isempty(boundary_conditions.u_t))
     nodes = mesh.boundary_nodes(boundary_conditions.u_t{1});
     
     if is_function_handle(boundary_conditions.u_t{2})
       values = boundary_conditions.u_t{2}(mesh.nodes(nodes, :));
     else
       values = boundary_conditions.u_t{2};
     end

    for n = 1:length(nodes)
      mat(n_u_dof + nodes(n), :) = 0;
      mat(n_u_dof + nodes(n), n_u_dof + nodes(n)) = 1;
      rhs(n_u_dof + nodes(n)) = values(n);
    end
  end
  
  % pin the pressure in one point:
  mat(end, :) = 0;
  mat(end, end) = 1;
  rhs(end) = 0;


    
