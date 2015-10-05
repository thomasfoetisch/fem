function zzgrad = zz_grad_u(mesh, gradu)
  zzgrad = zeros(size(mesh.nodes));

  zzgrad(:, 1) = build_clement_interpolant(mesh, gradu(:, 1));
  zzgrad(:, 2) = build_clement_interpolant(mesh, gradu(:, 2));
end
