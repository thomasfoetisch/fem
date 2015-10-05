function rhs = assemble_p1bullep1_stokes_rhs(mesh, dof_map, ints, force_f)

  % system sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % allocate some space
  rhs = spalloc(n_dof, 1, 2 * n_u_dof);
  f_dof = force_f(mesh.nodes);
  
  % build the rhs:
  for el = 1:n_elems
    for sol_dim = 1:2
      for i = 1:4
	for k = 1:3 % no bubble for the force field
	  rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ... 
	  = rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	    + mesh.jac(el) * f_dof(dof_map(el, k), sol_dim) * ints.phiphi(i, k);
	end
      end
    end
  end
