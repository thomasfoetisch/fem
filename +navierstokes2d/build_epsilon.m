function epsilon = build_epsilon(mesh, gu_x, gu_y)

  gu = {gu_x, gu_y};

  epsilon = zeros(2, 2, size(mesh.elements, 1));
  for el = 1:size(mesh.elements, 1)
    for i = 1:2
      for j = 1:2
	epsilon(i, j, el) = 0.5 * (gu{i}(el, j) + gu{j}(el, i));
      end
    end
  end
