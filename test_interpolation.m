clear all;
setenv('GNUTERM', 'qt')

mesh1 = geometry.build_square_mesh(1, 1, 1, 1, 0);
mesh2 = geometry.build_square_mesh(0.9, 0.9, 20, 20, 0);


function_f = @(x) sin(2 * pi * sqrt(x(:, 1).^2 + x(:, 2).^2));
field1 = function_f(mesh1.nodes);


edge_adj1 = geometry.build_elements_edge_adjacency(mesh1);
bc1 = geometry.build_barycentric_coordinates(mesh1);
tic();
[field2, average_path_length] = geometry.interpolate_on_mesh(mesh1, edge_adj1, bc1, field1, mesh2.nodes);
time = toc();

printf('Interpolation of %i points took %f, with an average path length of %f\n', size(mesh2.nodes, 1), time, average_path_length);

vis = 0;
if vis
  figure(1);
  trisurf(mesh1.elements, mesh1.nodes(:, 1), mesh1.nodes(:, 2), field1);
  figure(2);
  trisurf(mesh2.elements, mesh2.nodes(:, 1), mesh2.nodes(:, 2), field2);
end
