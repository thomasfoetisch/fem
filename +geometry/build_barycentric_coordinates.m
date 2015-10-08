function coordinates = build_barycentric_coordinates(mesh)
  coordinates = zeros(3, 3, size(mesh.elements, 1));

  for el = 1:size(mesh.elements, 1)
      m = [mesh.nodes(mesh.elements(el, 1), :)', ...
           mesh.nodes(mesh.elements(el, 2), :)', ...
           mesh.nodes(mesh.elements(el, 3), :)'];
    coordinates(:, :, el) = inv([m; [1 1 1]]);
  end

