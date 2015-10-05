function rhs = assemble_p1_gls_rhs(mesh, ints, rhs_f, delta, w_f)

  % assemble the right hand side:
  rhs = zeros(size(mesh.nodes, 1), 1);
  f_nodes = rhs_f(mesh.nodes);

  % evaluate the convection field on the nodes:
  w_nodes = w_f(mesh.nodes);

  for el = 1:size(mesh.elements, 1)
    elem = zeros(3, 1);
      
    if delta > 0
      h_k = sqrt(2 * mesh.jac(el)); %% i'm not sure if it's right...
      barycenter = sum([mesh.nodes(mesh.elements(el, :), 1),
			mesh.nodes(mesh.elements(el, :), 2)], 1) / 3;
      w_barycenter = w_f(barycenter);
      w_inv = 1 / sqrt(sum(w_barycenter.^2));

      coefficient = delta * h_k * w_inv * mesh.jac(el);

      for j = 1:3
	contrib = 0;
	for n = 1:3
	  for d = 1:2
	    for m = 1:3
	      for p = 1:2
		  contrib = contrib + f_nodes(mesh.elements(el, n)) ...
				      * w_nodes(mesh.elements(el, m), d) ...
				      * mesh.jmt(d, p, el) ...
				      * ints.phiphidphi(n, m, p, j);
	      end
	    end
	  end
	end
	elem(j) = elem(j) + coefficient * contrib;
      end
    end

    % compute the source term contribution:
    for j = 1:3
      contrib = 0;
      for p = 1:3
	contrib = contrib ...
		  + f_nodes(mesh.elements(el, p)) ...
		    * mesh.jac(el) ...
		    * ints.phiphi(j, p);
      end
      elem(j) = elem(j) + contrib;
    end

    % accumulate the element contribution into the rhs:
    for j = 1:3
      rhs(mesh.elements(el, j)) = rhs(mesh.elements(el, j)) + elem(j);
    end
  end
end
