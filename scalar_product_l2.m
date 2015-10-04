function sp = scalar_product_l2(mesh, u, v)

  % midedge quadrature exact for p2
  sp = 0;
  for e = 1:size(mesh.elements, 1)
    for k = 1:3
      kk = mod(k, 3) + 1;
      sp = sp + mesh.jac(e)/6 * (u(mesh.elements(e, k)) + u(mesh.elements(e, kk)))/2 ...
      * (v(mesh.elements(e, k)) + v(mesh.elements(e, kk)))/2;
    end
  end
end
