function values = viscosity_term_1(points, mesh, dof_map, ctx, u_x, u_y, element, basis, i, k, m, n)
  n_points = size(points, 1);

  u = {u_x(dof_map(element, :)), u_y(dof_map(element, :))};

  epsilon = zeros(2, 2, n_points);
  epsilon_norm = zeros(n_points, 1);
  %for d = 1:n_points
  for i_ = 1:2
    for j_ = 1:2
      for p = 1:4
	for r = 1:2
	  epsilon(i_, j_, :) = epsilon(i_, j_, :) + ...
			       reshape(1/2 * (u{i_}(p) * mesh.jmt(j_, r, element) * basis.dphi{r, p}(points) ...
					      + u{j_}(p) * mesh.jmt(i_, r, element) * basis.dphi{r, p}(points)), 1, 1, []);
	end
      end
    end
  end
  %end

  for i_ = 1:n_points
    epsilon_norm(i_) = sqrt(sum(epsilon(:, :, i_)(:).^2, 1));
  end

  values = zeros(n_points, 1);
  
  %for d = 1:n_points
    for p = 1:2
      for q = 1:2
	for r = 1:2
	  for s = 1:2
	    values(:) = values(:) ...
			+ 2 * ctx.rho * ctx.smagorinsky_coefficient * ctx.smagorinsky_caracteristic_length^2 ...
			  * mesh.jac(element) ...
			  .* squeeze(epsilon(m, r, :)) .* squeeze(epsilon(i, s, :)) ...
			  .* mesh.jmt(r, p, element) * mesh.jmt(s, q, element) ...
			  .* basis.dphi{p, n}(points) ...
			  .* basis.dphi{q, k}(points);
	  end
	end
      end
    end
  %end

  non_null_components = find(abs(epsilon_norm) > 1.e-10);
  values(non_null_components) = values(non_null_components) ./ epsilon_norm(non_null_components);
