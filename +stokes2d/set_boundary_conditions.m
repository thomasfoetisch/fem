function [mat, rhs] = set_boundary_conditions(mesh, dof_map, ints, mat, rhs, pressure_constraint_type)

  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % set the dirichlet boundary conditions for the velocity:
  for sol_dim = 1:2
    for n = 1:length(mesh.boundary_nodes)
      mat((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), :) = 0;
      mat((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), ...
	  (sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n)) = 1;
      rhs((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), 1) = 0;
    end
  end


  % process the pressure constraint method:
  if strcmpi(pressure_constraint_type, 'mean')
    mean_pressure = 1;
  elseif strcmpi(pressure_constraint_type, 'node')
      mean_pressure = 0;
  else
      error('Unknown pressure contraint method.');
  end


  % set the dirichlet boundary conditions for the pressure:
  if mean_pressure
    mean_p = spalloc(1, n_dof + 1, n_nodes + 1);
    for el = 1:n_elems
      for i = 1:3
	mean_p(1, 2 * n_u_dof + dof_map(el, i)) = mean_p(1, 2 * n_u_dof + dof_map(el, i)) + mesh.jac(el) * ints.phi(i);
      end
    end
    mean_p(1, end) = 0;
    mat = [mat, mean_p(1:end-1)'; mean_p];
    rhs = [rhs; 0];
  else
    mat(end, :) = 0;
    mat(end, end) = 1;
  end



    
