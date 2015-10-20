function rhs = assemble_p1bullep1_stationary_rhs(mesh, dof_map, ints, basis,...
                                                 u_x, u_y, pressure, ...
						 force_f, ...
						 ctx)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  
  % matrix allocation:
  rhs = spalloc(n_dof, 1, 0);
  

  % p1 approximation of the force field:
  force_nodes = force_f(mesh.nodes);

  % helpers:
  u = {u_x, u_y};
  kronecker = eye(2);

  for el = 1:n_elems
    printf('rhs element %i\n', el);

    dofs = dof_map(el, :);
    
    % assemble F_{i, k}:
    for i = 1:2
      for k = 1:4
	contrib = 0;

	% convection:
	for q = 1:4
	  for p = 1:4
	    for j = 1:2
	      for d = 1:2
		contrib = contrib + ctx.rho * mesh.jac(el) * u{i}(dofs(q)) * u{j}(dofs(p)) * mesh.jmt(j, d, el) * ints.phiphidphi(p, k, d, q);
	      end
	    end
	  end
	end

	% viscosity:
	if 1
	  contrib = contrib + navierstokes2d.quadrature5(@(x) navierstokes2d.viscosity_term_rhs(x, mesh, dof_map, ctx, u_x, u_y, el, basis, i, k));
	end

	% pressure:
	for q = 1:2
	  for p = 1:3
	    contrib = contrib - mesh.jac(el) * mesh.jmt(i, q, el) * pressure(dofs(p)) * ints.phidphi(p, q, k);
	  end
	end
	
	% force:
	for q = 1:3
	  contrib = contrib - mesh.jac(el) * force_nodes(dofs(q), i) * ints.phiphi(q, k);
	end

	rhs((i - 1) * n_u_dof + dof_map(el, k), 1) = ...
	rhs((i - 1) * n_u_dof + dof_map(el, k), 1) + contrib;
      end
    end


    % assemble G_{k}:
    for k = 1:3
      contrib = 0;

      % divergence:
      for q = 1:2
	for j = 1:4
	  for p = 1:2
	    contrib = contrib - mesh.jac(el) * u{q}(dofs(j)) * mesh.jmt(q, p, el) * ints.phidphi(k, p, j);
	  end
	end
      end

      rhs(2 * n_u_dof + dof_map(el, k), 1) = ...
      rhs(2 * n_u_dof + dof_map(el, k), 1) + contrib;
    end
  end

  edge_dof_map = mesh.boundary_edges;

  n_edges = size(mesh.boundary_edges, 1);

  for ed = 1:n_edges
    for i = 1:2
      for k = 1:2
	%contrib = navierstokes2d.quadrature9_1d(@(x) navierstokes2d.boundary_term_rhs(x, mesh, u_x, u_y, pressure, ed, basis, i, k));
	  
	%rhs((i - 1) * n_u_dof + edge_dof_map(ed, k), 1) = ...
	%rhs((i - 1) * n_u_dof + edge_dof_map(ed, k), 1) + contrib; 
      end
    end
  end
