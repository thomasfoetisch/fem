function rhs = assemble_p1_rhs(mesh, ints, rhs_f)

  % assemble the right hand side:
  rhs = zeros(size(mesh.nodes, 1), 1);
  f_nodes = rhs_f(mesh.nodes);

  for e = 1:size(mesh.elements, 1)
    elem = zeros(3, 1);
      
    % compute the source term contribution:
    for j = 1:3
      contrib = 0;
      for p = 1:3
	contrib = contrib ...
		  + f_nodes(mesh.elements(e, p)) ...
		    * mesh.jac(e) ...
		    * ints.phiphi(j, p);
      end
      elem(j) = elem(j) + contrib;
    end

    % accumulate the element contribution into the rhs:
    for j = 1:3
      rhs(mesh.elements(e, j)) = rhs(mesh.elements(e, j)) + elem(j);
    end
  end
end
