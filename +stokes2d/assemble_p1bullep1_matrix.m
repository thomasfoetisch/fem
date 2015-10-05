function mat = assemble_p1bullep1_matrix(mesh, dof_map, ints, mu)

  % build the matrix:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_dof = n_u_dof * 2 + n_nodes;

  % allocate the matrix:
  mat = spalloc(n_dof, n_dof, 0);

  for el = 1:n_elems
    % elliptic term:
    elem_u = zeros(4, 2, 4, 2);
    for sol_dim = 1:2
      for test_dim = 1:2
	if sol_dim == test_dim
	  for i = 1:4
	    for j = 1:4
	      for k = 1:2
		for p = 1:2
		  for q = 1:2
		    elem_u(i, sol_dim, j, test_dim) ...
		    = elem_u(i, sol_dim, j, test_dim) ...
		      + mu * mesh.jac(el) ...
			* mesh.jmt(k, p, el) ...
			* mesh.jmt(k, q, el) ...
			* ints.dphidphi(p, i, q, j);
		  end
		end
	      end
	    end
	  end
	end
      end
    end

    % pressure term:
    elem_p = zeros(3, 4, 2);
    for test_dim = 1:2 
      for i = 1:3
	for j = 1:4
	  for p = 1:2
	    elem_p(i, j, test_dim) ...
	    = elem_p(i, j, test_dim) ...
	      - mesh.jac(el) ...
		* mesh.jmt(test_dim, p, el) ...
		* ints.phidphi(i, p, j);
	  end
	end
      end
    end

    % incompressibility constraint:
    elem_div = zeros(4, 2, 3);
    for sol_dim = 1:2
      for j = 1:3
	for i = 1:4
	  for p = 1:2
	    elem_div(i, sol_dim, j) ...
	    = elem_div(i, sol_dim, j) ...
	      + mesh.jac(el) ...
		* mesh.jmt(sol_dim, p, el) ...
		* ints.phidphi(j, p, i);
	  end
	end
      end
    end


    % accumulate everything in the linear system:
    for i = 1:4
      for j = 1:4
	for sol_dim = 1:2
	  for test_dim = 1:2
	    if test_dim == sol_dim
	      mat((test_dim - 1) * n_u_dof + dof_map(el, j),...
		  (sol_dim - 1) * n_u_dof + dof_map(el, i)) ...
	      = mat((test_dim - 1) * n_u_dof + dof_map(el, j), ...
		    (sol_dim - 1) * n_u_dof + dof_map(el, i)) ...
		+ elem_u(i, sol_dim, j, test_dim);
	    end
	  end
	end
      end
    end

    for i = 1:3
      for j = 1:4
	for test_dim = 1:2
	  mat((test_dim - 1) * n_u_dof + dof_map(el, j),
	      2 * n_u_dof + dof_map(el, i)) ...
	  = mat((test_dim - 1) * n_u_dof + dof_map(el, j),
		2 * n_u_dof + dof_map(el, i)) ...
	    + elem_p(i, j, test_dim);
	end
      end
    end

    for i = 1:4
      for sol_dim = 1:2
	for j = 1:3
	  mat(2 * n_u_dof + dof_map(el, j),
	      (sol_dim - 1) * n_u_dof +  dof_map(el, i)) ...
	  = mat(2 * n_u_dof + dof_map(el, j),
		(sol_dim - 1) * n_u_dof +  dof_map(el, i)) ...
	    + elem_div(i, sol_dim, j);
	end
      end
    end
  end
