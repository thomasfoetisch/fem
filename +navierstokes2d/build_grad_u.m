function gu = build_grad_u(mesh, ints, u)

  gu = zeros(size(mesh.elements, 1), 2);
  for el = 1:size(mesh.elements, 1)
    for d = 1:2
      for k = 1:3
	for p = 1:2
          gu(el, d) = gu(el, d) + mesh.jmt(d, p, el) ...
				  * ints.gradphi_p1(p, k) ...
				  * u(mesh.elements(el, k));
	end
      end
    end
  end

    
