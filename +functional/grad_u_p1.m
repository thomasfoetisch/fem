function du = grad_u_p1(mesh, ints, field_p1)
  du = zeros(size(mesh.elements, 1), 2);

  for el = 1:size(mesh.elements, 1)
    du(el, :) = (mesh.jmt(:, :, el) * ints.dphi) * field_p1(mesh.elements(el, :));
  end

end
