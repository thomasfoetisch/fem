function mesh = build_square_mesh(l_x, l_y, n_el_x, n_el_y, theta)

  % build the nodes:
  nodes = zeros((n_el_x + 1), (n_el_y + 1), 2);

  h_x = l_x / n_el_x;
  h_y = l_y / n_el_y;

  for i = 1:(n_el_x + 1)
    for j = 1:(n_el_y + 1)
      nodes(i, j, 1) = (i-1) * h_x - l_x/2;
      nodes(i, j, 2) = (j-1) * h_y - l_y/2;
    end
  end

  nodes = reshape(nodes, [(n_el_x + 1) * (n_el_y + 1), 2]);
  nodes_indices = reshape(int16(1:((n_el_x + 1) * (n_el_y + 1))),
			  [n_el_x + 1, n_el_y + 1]);

  rot = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
  nodes = (rot(theta)*nodes')';

  % build elements:
  elements = int16(zeros(3, 2, n_el_x ,n_el_y));
  for i = 1:n_el_x
    for j = 1:n_el_y
      % first triangle:
      elements(1, 1, i, j) = nodes_indices(i,     j);
      elements(2, 1, i, j) = nodes_indices(i + 1, j);
      elements(3, 1, i, j) = nodes_indices(i,     j + 1);

      % second triangle:
      elements(1, 2, i, j) = nodes_indices(i + 1, j + 1);
      elements(2, 2, i, j) = nodes_indices(i,     j + 1);
      elements(3, 2, i, j) = nodes_indices(i + 1, j);
    end
  end

  elements = reshape(elements, [3, n_el_x * n_el_y * 2]);
  elements = transpose(elements);


  % random permutation of the node numbering:
  perm = randperm(size(nodes, 1));
  perm = 1:size(nodes, 1);
  [~, inv_perm] = sort(perm);

  nodes = nodes(perm, :);
  elements = reshape(inv_perm(reshape(elements, 1, [])), size(elements));

  for e = 1:size(elements, 1)
    %elements(e, :) = elements(e, randperm(3));
  end
  elements = int16(elements);


  % build the adjacency matrix:
  adj = int8(zeros(size(nodes, 1), size(nodes, 1)));
  for i = 1:size(elements, 1)
    for e = 1:3
      adj(elements(i, e), elements(i, mod(e + 1, 3) + 1)) = 1;
      adj(elements(i, mod(e + 1, 3) + 1), elements(i, e)) = 1;
    end
  end
  for n = 1:size(nodes, 1)
    adj(n, n) = 1;
  end


  % build the jacobian matrix:
  jmt = zeros(2, 2, size(elements, 1));
  jac = zeros(size(elements, 1), 1);
  for e = 1:size(elements, 1)
    origin = nodes(elements(e, 1), :)';
    m = [nodes(elements(e, 2), :)' - origin, nodes(elements(e, 3), :)' - origin];
    jmt(:, :, e) = inv(m)';
    jac(e) = abs(det(m));
  end


  % build the boundary node list:
  edges = int16(zeros(3*size(elements, 1), 3));
  edges_number = 0;
  for e = 1:size(elements)
    for k = 1:3
      kk = mod(k, 3) + 1;
      edge = sort(elements(e, [k, kk]));
      loc = find(edges(:, 1) == edge(1) & edges(:, 2) == edge(2));
      if isempty(loc)
	edges(edges_number + 1, :) = [edge, 1];
	edges_number = edges_number + 1;
      else
	edges(loc, 3) = edges(loc, 3) + 1;
      end
    end
  end

  boundary_nodes = sort(unique(edges(find(edges(:, 3) == 1), 1:2)));
  inner_nodes = setdiff(1:size(nodes, 1), boundary_nodes);

  % build the barycenters of the elements:
  barycenters = zeros(size(elements, 1), 2);
  for e = 1:size(elements, 1)
    barycenters(e, :) = sum([nodes(elements(e, :), 1),  nodes(elements(e, :), 2)], 1)/3;
  end


  % build the diameter of the elements:
  diameters = zeros(size(elements, 1), 1);
  dist = @(x) sqrt(sum(x.^2, 2));
  for el = 1:size(elements, 1)
    diameters(el) = max(dist([diff(nodes(elements(el, [1, 2]), :));
			      diff(nodes(elements(el, [1, 3]), :));
			      diff(nodes(elements(el, [2, 3]), :))]));
  end

  
  % build the node-element adjacency list:
  node_el_adj = cell(1, size(nodes, 1));
  for node = 1:size(nodes, 1)
    node_el_adj(node) = find(sum(elements == node, 2));
  end

  % return the result in a struct:
  mesh = struct('nodes', nodes,
		'elements', elements,
		'adj', adj,
		'boundary_nodes', boundary_nodes,
		'inner_nodes', inner_nodes,
		'jac', jac,
		'jmt', jmt,
		'barycenters', barycenters,
		'h', diameters);

  mesh.node_el_adj = node_el_adj;

end
