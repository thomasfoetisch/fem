function [lifting_u_x, lifting_u_y] = build_lifting(mesh, boundary_conditions)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;


  % allocate memory for the lifting:
  lifting_u_x = zeros(n_u_dof, 1);
  lifting_u_y = zeros(n_u_dof, 1);
  

  % 1. inhomogeneous dirichlet conditions on: u_x:
  if not(isempty(boundary_conditions.u_x))
    nodes = mesh.boundary_nodes(boundary_conditions.u_x{1});
    if is_function_handle(boundary_conditions.u_x{2})
      values = boundary_conditions.u_x{2}(mesh.nodes(nodes, :));
    else
      values = boundary_conditions.u_x{2};
    end
    
    lifting_u_x(nodes) = values;
  end

  % 2. inhomogeneous dirichlet conditions on: u_y:
  if not(isempty(boundary_conditions.u_y))
    nodes = mesh.boundary_nodes(boundary_conditions.u_y{1});

    if is_function_handle(boundary_conditions.u_y{2})
      values = boundary_conditions.u_y{2}(mesh.nodes(nodes, :));
    else
      values = boundary_conditions.u_y{2};
    end

    lifting_u_y(nodes) = values;
  end

  % 3. inhomogeneous dirichlet conditions on: u_n:
  if not(isempty(boundary_conditions.u_n))
    nodes = mesh.boundary_nodes(boundary_conditions.u_n{1});
    
    if is_function_handle(boundary_conditions.u_n{2})
      values = boundary_conditions.u_n{2}(mesh.nodes(nodes, :));
    else
      values = boundary_conditions.u_n{2};
    end

    lifting_u_x(nodes) = values;
  end

  % 4. inhomogeneous dirichlet conditions on: u_t:
  if not(isempty(boundary_conditions.u_t))
    nodes = mesh.boundary_nodes(boundary_conditions.u_t{1});
    
    if is_function_handle(boundary_conditions.u_t{2})
      values = boundary_conditions.u_t{2}(mesh.nodes(nodes, :));
    else
      values = boundary_conditions.u_t{2};
    end

    lifting_u_y(nodes) = values;
  end

