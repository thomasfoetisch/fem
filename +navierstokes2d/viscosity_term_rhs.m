function values = viscosity_term_rhs(points, mesh, dof_map, ctx, u_x, u_y, element, basis, i, k)

  n_points = size(points, 1);
  u = {u_x(dof_map(element, :)), u_y(dof_map(element, :))};

  epsilon = zeros(2, 2, n_points);
  epsilon_norm = zeros(n_points, 1);
  for d = 1:n_points
  for i_ = 1:2
    for j_ = 1:2
      for p = 1:4
	for r = 1:2
	  epsilon(i_, j_, d) = epsilon(i_, j_, d) + ...
			       1/2 * (u{i_}(p) * mesh.jmt(j_, r, element) * basis.dphi{r, p}(points(d, :)) ...
				      + u{j_}(p) * mesh.jmt(i_, r, element) * basis.dphi{r, p}(points(d, :)));
	end
      end
    end
  end
  end

  for i_ = 1:n_points
    epsilon_norm(i_) = sqrt(sum(epsilon(:, :, i_)(:).^2, 1));
  end

  mu = ctx.laminar_viscosity(ones(n_points, 1)) + ctx.smagorinsky_coefficient * ctx.rho * ctx.smagorinsky_caracteristic_length * epsilon_norm;

  values = zeros(n_points, 1);
  for d = 1:n_points
    for p = 1:2
      for q = 1:2
	values(d) = values(d) + 2 * mesh.jac(element) * mu(d) * epsilon(q, i, d) * mesh.jmt(q, p, element) * basis.dphi{p, k}(points(d, :));
      end
    end
  end
