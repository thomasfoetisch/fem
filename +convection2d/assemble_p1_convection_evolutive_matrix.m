function mat = assemble_p1_convection_evolutive_matrix(mesh, ints, w_f, mu, delta, tau)

  % initialise the matrix:
  mat = spalloc(size(mesh.adj, 1), ...
  		size(mesh.adj, 1), ...
  		sum(mesh.adj(:)));


  % evaluate the convection field on the nodes:
  w_nodes = w_f(mesh.nodes);


  % loop over all the elements:U
  for el = 1:size(mesh.elements, 1)
    elem = zeros(3, 3);

    % add the time dependance term:
    for i = 1:3
      for j = 1:3
	contrib = 0;

	contrib = contrib ...
		  + mesh.jac(el) * ints.phiphi(i, j) / tau;

	elem(i, j) = elem(i, j) + contrib;
      end
    end

    % build the convection term:
    for i = 1:3
      for j = 1:3
	contrib = 0;
	for d = 1:2
	  for p = 1:3
	    for n = 1:2
	      contrib = contrib ...
			+ mesh.jac(el) ...
			  * w_nodes(mesh.elements(el, p), d) ...
			  * mesh.jmt(d, n, el) ...
			  * ints.phiphidphi(p, i, n, j);
	    end
	  end
	end

	% accumulate the contribution:
	elem(i, j) = elem(i, j) + contrib;
      end
    end

    % build the diffusion term:
    for i = 1:3
      for j = 1:3
	contrib = 0;
	for p = 1:2
	  for m = 1:2
	    for n = 1:2

	      contrib = contrib + ... 
			mu ...
			* mesh.jac(el) ...
			* mesh.jmt(p, m, el) ...
			* mesh.jmt(p, n, el) ...
			* ints.dphidphi(m, i, n, j);
	    end
	  end
	end
	elem(i, j) = elem(i, j) + contrib;
      end
    end

    % build the GLS stabilisation term if needed:
    if delta > 0
      h_k = mesh.h(el);
      barycenter = sum([mesh.nodes(mesh.elements(el, :), 1),
			mesh.nodes(mesh.elements(el, :), 2)], 1) / 3;
      w_barycenter = w_f(barycenter);
      w_inv = 1 / sqrt(sum(w_barycenter.^2));

      coefficient = delta * h_k * w_inv * mesh.jac(el);

      for i = 1:3
	for j = 1:3
	  contrib = 0;
	  for n = 1:3
	    for d = 1:2
	      for m = 1:3
		for e = 1:2
		  for p = 1:2
		    for q = 1:2
		      contrib = contrib + w_nodes(mesh.elements(el, n), d) ...
					  * w_nodes(mesh.elements(el, m), e) ...
					  * mesh.jmt(d, p, el) ...
					  * mesh.jmt(e, q, el) ...
					  * ints.phiphidphidphi(n, m, p, i, q, j);
		    end
		  end
		end
	      end
	    end
	  end
	  elem(i, j) = elem(i, j) + coefficient * contrib;
	end
      end
    end

    % accumulate the element contribution to the matrix:
    for i = 1:3
      for j = 1:3
    	mat(mesh.elements(el, i), ...
	    mesh.elements(el, j)) = mat(mesh.elements(el, i), ...
					mesh.elements(el, j)) ...
				    + elem(i, j);
      end
    end
  end
end
