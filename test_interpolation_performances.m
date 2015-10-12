clear all;
more off;
setenv('GNUTERM', 'qt')

function_f = @(x) sin(2 * pi * sqrt(x(:, 1).^2 + x(:, 2).^2));


n_elems = [];
interp_time = [];
interp_time_per_point = [];
average_path_length = [];
for n_elem = 5:40
  printf('Elements per side: %i\n', n_elem);

  mesh1 = geometry.build_square_mesh(1, 1, n_elem, n_elem, 0);
  mesh2 = geometry.build_square_mesh(0.9, 0.9, 50, 50, 0);

  field1 = function_f(mesh1.nodes);

  edge_adj1 = geometry.build_elements_edge_adjacency(mesh1);
  bc1 = geometry.build_barycentric_coordinates(mesh1);
  tic();
  [field2, apl, vpl] = geometry.interpolate_on_mesh(mesh1, edge_adj1, bc1, field1, mesh2.nodes);
  time = toc();

  n_elems(end + 1) = size(mesh1.elements, 1);
  interp_time(end + 1) = time;
  interp_time_per_point(end + 1) = time / size(mesh2.nodes, 1);
  average_path_length(end + 1) = apl;
  variance_path_length(end + 1) = vpl;
end


vis = 1;
if vis
   figure(1); cla;
   plot(n_elems, interp_time);
   figure(2); cla;
   plot(n_elems, interp_time_per_point);
   figure(3); cla;
   errorbar(n_elems, average_path_length, variance_path_length);
end
