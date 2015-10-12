function mat = assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, ...
						    u_x, u_y, ...
						    rho, ...
						    laminar_viscosity, ...
						    smagorinsky_coefficient, ...
						    smagorinsky_caracteristic_length)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  
  % matrix allocation:
  mat = spalloc(n_dof, n_dof, 0);
  

  % necessary fields:
  u = {u_x, u_y};
  grad_u_x = navierstokes2d.build_grad_u(mesh, ints, u_x);
  grad_u_y = navierstokes2d.build_grad_u(mesh, ints, u_y);

  epsilon = navierstokes2d.build_epsilon(mesh, grad_u_x, grad_u_y);

  mu = navierstokes2d.build_mu(mesh, epsilon, ...
			       laminar_viscosity, ...
			       smagorinsky_coefficient, ...
			       smagorinsky_caracteristic_length);

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
		contrib = contrib + rho * mesh.jac(el) * mesh.jmt(m, p, el) * u{i}(dofs(q)) * ints.phiphidphi(n, k, p, q);
		for j = 1:2
		  contrib = contrib + rho * mesh.jac(el) * mesh.jmt(j, p, el) * u{j}(dofs(q)) *  kronecker(i, m) * ints.phiphidphi(q, k, p, n);
		end
	      end 
	    end

	    for p = 1:2
	      for q = 1:2
		for r = 1:2
		  for s = 1:2
		    contrib = contrib + rho * 4 * smagorinsky_coefficient * smagorinsky_caracteristic_length^2 ...
					* mesh.jac(el) * mesh.jmt(p, q, el) * mesh.jmt(s, r, el) * epsilon(m, p, el) * epsilon(s, i, el) * ints.dphidphi(q, n, r, k);
		  end
		end
	      end
	    end

	    for p = 1:2
	      for q = 1:2
		contrib = contrib + rho * mu(el) * mesh.jac(el) * mesh.jmt(i, p, el) * mesh.jmt(m, q, el) * ints.dphidphi(p, n, q, k);
		for r = 1:2
		  for s = 1:2
		    contrib = contrib + rho * mu(el) * mesh.jac(el) * kronecker(i, m) * mesh.jmt(p, s, el) * mesh.jmt(p, r, el) * ints.dphidphi(r, n, s, k);
		  end
		end
	      end
	    end

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
		contrib = contrib + mesh.jac(el) * mesh.jmt(m, p, el) * ints.phidphi(k, p, n);
	    end

	    elem_div(i, k, m) = contrib;
	    mat(2 * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) = ...
	    mat(2 * n_u_dof + dof_map(el, k), (m - 1) * n_u_dof + dof_map(el, n)) + contrib;
	end
      end
    end

  end
