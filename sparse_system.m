function [nodes, elements, u, err_l2] = sparse_system(h)

  % parameters:
  l_x = 1;
  l_y = 1;

  n_el_x = round(l_x / h);
  n_el_y = round(l_y / h);

  rhs_f = @(x) (- 4) * ones(size(x, 1), 1);
  u_exact_f = @(x) x(:, 1).^2 + x(:, 2).^2;

  % elementary integrals
  one = 1/2;
  phi = [1/6, 1/6, 1/6];
  dphi = [-1, 1, 0,
          -1, 0, 1];
  phiphi = [1/12, 1/24, 1/24; 1/24, 1/12, 1/24; 1/24, 1/24, 1/12]; % brin d'acier haha

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
  [~, inv_perm] = sort(perm);

  nodes = nodes(perm, :);
  elements = reshape(inv_perm(reshape(elements, 1, [])), size(elements));

  for e = 1:size(elements, 1)
    elements(e, :) = elements(e, randperm(3));
  end
  elements = int16(elements);


  % build the barycenters of the elements:
  barycenters = zeros(size(elements, 1), 2);
  for e = 1:size(elements, 1)
    barycenters(e, :) = sum([nodes(elements(e, :), 1),  nodes(elements(e, :), 2)], 1)/3;
  end


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


  % assemble the matrix:
  mat = spalloc(size(nodes, 1),
  	      size(nodes, 1),
  	      sum(adj(:)));

  for e = 1:size(elements, 1)
    for i = 1:3
      for j = 1:3
	contrib = jac(e)/2 * (jmt(:, :, e) * dphi(:, i))' * (jmt(:, :, e) * dphi(:, j));
	mat(elements(e, i), elements(e, j)) = mat(elements(e, i), elements(e, j)) + contrib;
      end
    end
  end


  % assemble the right hand side:
  rhs = zeros(size(nodes, 1), 1);
  f_nodes = rhs_f(nodes);

  for e = 1:size(elements, 1)
    for j = 1:3
      for p = 1:3
	rhs(elements(e, j)) = rhs(elements(e, j)) + f_nodes(elements(e, p)) * jac(e) * phiphi(j, p);
      end
    end
  end


  % build the boundary node list:
  edge_adj = zeros(size(elements, 1), size(elements, 1));
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


  % set boundary conditions
  g_nodes = u_exact_f(nodes);

  [i, j, v] = find(mat);
  for bcn = 1:size(boundary_nodes, 1)
    idx = find(i == boundary_nodes(bcn));
    mat(i(idx), j(idx)) = 0;
    %mat(j(idx), i(idx)) = 0; % if we wanna kill the columns, we need to modify the rhs
    mat(boundary_nodes(bcn), boundary_nodes(bcn)) = 1;
    rhs(boundary_nodes(bcn)) = g_nodes(boundary_nodes(bcn));
  end


  % solve the linear system:
  u = linsolve(mat, rhs);
  u_exact = g_nodes;


  % compute the L_2 error:
  err_sqr = 0;
  mid_edges = zeros(size(elements, 1) * 3, 2);
  for e = 1:size(elements, 1)
    for k = 1:3
      kk = mod(k, 3) + 1;
      mid_edge = sum([nodes(elements(e, [k, kk]), 1), nodes(elements(e, [k, kk]), 2)], 1) / 2;
      err_sqr = err_sqr + jac(e)/6 * (g(mid_edge) - 1/2 * (u(elements(e, k)) + u(elements(e, kk))))^2;
      mid_edges(3 * (e - 1) + k, :) = mid_edge;
    end
  end
  
  
  % return the L_2 error:
  err_l2 = sqrt(err_sqr);

end
