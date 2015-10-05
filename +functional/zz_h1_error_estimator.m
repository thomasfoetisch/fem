function err = zz_h1_error_estimator(mesh, ints, u)
  du = grad_u_p1(mesh, ints, u);
  gu = zeros(size(mesh.nodes));

  gu(:, 1) = build_clement_interpolant(mesh, du(:, 1));
  gu(:, 2) = build_clement_interpolant(mesh, du(:, 2));

  err = sqrt(error_l2_p0_p1(mesh, du(:, 1), gu(:, 1))^2 ...
	     + error_l2_p0_p1(mesh, du(:, 2), gu(:, 2))^2);
end
