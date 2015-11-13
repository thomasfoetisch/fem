function rhs = assemble_p1bullep1_rhs(mesh, dof_map, ints, lifting, force_f)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = n_u_dof * 2 + n_nodes;

  % allocate some space
  rhs = spalloc(n_dof, 1, 2 * n_u_dof);
  f_dof = force_f(mesh.nodes);
  
  % build the rhs:
  for el = 1:n_elems
    printf('rhs element %i / %i\n', el, n_elems);

    for sol_dim = 1:2
      for i = 1:4
	for k = 1:3 % no bubble for the force field
	  rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	  = rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	    + mesh.jac(el) * f_dof(dof_map(el, k), sol_dim) * ints.phiphi(i, k);
	end
      end
    end


    % lifting contribution to the velocity:
    for i = 1:4
      for lift_dim = 1:2
	contrib = 0;
	for m = 1:2
	  for p = 1:2
	    for q = 1:2
	      for k = 1:4
		contrib = contrib - mesh.jac(el) * lifting(dof_map(el, k), lift_dim) * mesh.jmt(m, p, el) * mesh.jmt(m, q, el) * ints.dphidphi(p, k, q, i);
	      end
	    end
	  end
	end

	rhs((lift_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	= rhs((lift_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	  + contrib;
      end
    end

    % lifting contribution to the pressure equations:
    for i = 1:3
      contrib = 0;
      for k = 1:4
	for j = 1:2
	  for p = 1:2
	    contrib = contrib + mesh.jac(el) * mesh.jmt(j, p, el) * lifting(dof_map(el, k), j) * ints.phidphi(i, p, k);
	  end
	end
      end
      rhs(2*n_u_dof + dof_map(el, i), 1) ...
      = rhs(2*n_u_dof + dof_map(el, i), 1) ...
	+ contrib;
    end
  end
