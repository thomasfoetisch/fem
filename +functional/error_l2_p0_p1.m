function sp = error_l2_p0_p1(mesh, f_p0, g_p1)

  % midedge quadrature exact for p2
  sp = 0;
  for e = 1:size(mesh.elements, 1)
    for k = 1:3
      kk = mod(k, 3) + 1;
      
      midpoint = (mesh.nodes(mesh.elements(e, k),:) 
		  + mesh.nodes(mesh.elements(e, kk), :))/2;
      
      sp = sp + mesh.jac(e)/6 * (f_p0(e)
				 - (g_p1(mesh.elements(e, k)) + g_p1(mesh.elements(e, kk)))/2)^2;
      
    end
  end
  sp = sqrt(sp);
end
