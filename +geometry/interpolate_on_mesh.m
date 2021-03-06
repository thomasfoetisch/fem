function [values, average_path_length, variance_path_length] = interpolate_on_mesh(mesh, edge_adj, baryc, field_p1, points);

  values = zeros(size(points, 1), size(field_p1, 2));
  hint = 1;
  path_lengths = zeros(size(points, 1), 1);
  for n = 1:size(points, 1)
    %[hint, coordinates] = geometry.find_element_at_point_naive(mesh, edge_adj, baryc, points(n, :));
    [hint, coordinates, path_length] = geometry.find_element_at_point(mesh, edge_adj, baryc, points(n, :), hint);
    values(n, :) = coordinates' * field_p1(mesh.elements(hint, :), :);

    % if it is a p1-bubble field:
    if size(field_p1, 1) > size(mesh.nodes, 1)
       bubble = @(x, y) 27 * x * y * (1 - x - y);
       values(n, :) = values(n, :) + field_p1(size(mesh.nodes, 1) + hint, :) * bubble(coordinates(2), coordinates(3));
    end
    
    path_lengths(n) = path_length;
  end
  average_path_length = mean(path_lengths);
  variance_path_length = var(path_lengths);
  
