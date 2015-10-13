function mat = assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, basis, ...
						    u_x, u_y, ...
						    ctx)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  
  % matrix allocation:
  mat = spalloc(n_dof, n_dof, 0);
  
  % helpers:
  u = {u_x, u_y};
  kronecker = eye(2);

  for el = 1:n_elems
    printf('mat element %i\n', el);

    dofs = dof_map(el, :);


    elem_u = zeros(2, 4, 2, 4); % i, k, m, n
    for i = 1:2
      for k = 1:4
	for m = 1:2
	  for n = 1:4

	    contrib = 0;
	    for p = 1:2
	      for q = 1:4
		contrib = contrib + ctx.rho * mesh.jac(el) ...
				    * mesh.jmt(m, p, el) ...
				    * u{i}(dofs(q)) ...
				    * ints.phiphidphi(n, k, p, q);
		for j = 1:2
		  contrib = contrib + ctx.rho * mesh.jac(el) ...
				      * mesh.jmt(j, p, el) ...
				      * u{j}(dofs(q)) ...
				      * kronecker(i, m) ...
				      * ints.phiphidphi(q, k, p, n);
		end
	      end 
	    end

	    % viscosity .1
	    if 1
	      contrib = contrib + navierstokes2d.quadrature(@(x) navierstokes2d.viscosity_term_1(x, mesh, dof_map, ctx, u_x, u_y, el, basis, i, k, m, n));
	    end % enable/disable
	    
	    if 1
	       contrib = contrib + navierstokes2d.quadrature(@(x) navierstokes2d.viscosity_term_2(x, mesh, dof_map, ctx, u_x, u_y, el, basis, i, k, m, n));
	    end % enable/disable

	    elem_u(i, k, m, n) = contrib;
	    mat((i - 1) * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) = ...
	    mat((i - 1) * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) + contrib;
	  end 
	end
      end
    end

    elem_p = zeros(2, 4, 3); % i, k, m
    for i = 1:2
      for k = 1:4
	for m = 1:3
	    
	    contrib = 0;
	    for p = 1:2
	      contrib = contrib - mesh.jac(el) * mesh.jmt(i, p, el) * ints.phidphi(m, p, k);
	    end
	    
	    elem_p(i, k, m) = contrib;
	    mat((i - 1) * n_u_dof + dof_map(el, k), 2 * n_u_dof + dof_map(el, m)) = ...
	    mat((i - 1) * n_u_dof + dof_map(el, k), 2 * n_u_dof + dof_map(el, m)) + contrib;
	end
      end
    end


    elem_div = zeros(3, 2, 4); % k, m, n
    for k = 1:3
      for m = 1:2
	for n = 1:4

	    contrib = 0;
	    for p = 1:2
		contrib = contrib - mesh.jac(el) * mesh.jmt(m, p, el) * ints.phidphi(k, p, n);
	    end

	    elem_div(i, k, m) = contrib;
	    mat(2 * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) = ...
	    mat(2 * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) + contrib;
	end
      end
    end

  end
