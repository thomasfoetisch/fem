
f = @(x) sin(sum(x.^2, 2) * 2);
error = [];

for i = 3:7
  printf('i = %i\n', i);


  mesh = build_square_mesh(1, 1, 2^i, 2^i, 0);
  u = f(mesh.nodes);
  

  error(end + 1) = zz_h1_error_estimator(mesh, ints, u);
end
