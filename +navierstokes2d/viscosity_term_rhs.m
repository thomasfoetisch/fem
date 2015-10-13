function values = viscosity_term_rhs(points, mesh, dof_map, ctx, u_x, u_y, element, basis, i, k)

  n_points = size(points, 1);
  dof = mesh.elements(element, :);
  u = {u_x(dof_map(element, :)), u_y(dof_map(element, :))};

  epsilon = zeros(2, 2, n_points);
  for i = 1:2
    for j = 1:2
      for p = 1:4
	for r = 1:2
	  for s = 1:2
	    epsilon(i, j, :) = 1/2 * (u{i}(p) * mesh.jmt(j, r, element) * basis.dphi{r, p}(points) ...
				      + u{j}(p) * mesh.jmt(i, s, element) * basis.dphi{s, p}(points));
	  end
	end
      end
    end
  end
  epsilon_norm = sqrt(sum(reshape(epsilon, [4, n_points]).^2, 1));

  mu = ctx.laminar_viscosity + ctx.smagorinsky_coefficient * ctx.rho * ctx.smagorinsky_caracteristic_length * epsilon_norm;


  values = zeros(n_points, 1);
  for d = 1:n_points
    for p = 1:2
      for q = 1:2
	values(d) = values(d) + ctx.rho * 2 * mesh.jac(element) * mu(d) * epsilon(q, i, d) * mesh.jmt(q, p, element) * basis.dphi{p, k}(points(k, :));
      end
    end
  end
