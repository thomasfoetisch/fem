function mat = assemble_p1_laplacian_matrix(mesh, ints)
  
  % assemble the matrix:
  mat = spalloc(size(mesh.adj, 1), ...
  		size(mesh.adj, 1), ...
  		sum(mesh.adj(:)));

  for el = 1:size(mesh.elements, 1)
    for i = 1:3
      for j = 1:3
	contrib = 0;
	for m = 1:2
	  for n = 1:2
	    for p = 1:2
	      contrib = contrib + mesh.jac(el) ...
				  * mesh.jmt(p, n, el) ...
	                          * mesh.jmt(p, m, el) ...
				  * ints.dphidphi(n, i, m, j);
	    end
	  end
	end
	mat(mesh.elements(el, i), ...
	    mesh.elements(el, j)) = mat(mesh.elements(el, i), ...
					mesh.elements(el, j)) ...
				    + contrib;
      end
    end
  end
end
