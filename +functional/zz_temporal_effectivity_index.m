function eff = zz_temporal_effectivity_index(mesh, ints, sol, w, u_exact_f, tau, time)
  err = error_l2(mesh, sol{end}, @(x) u_exact_f(x, time));

  est = zeros(length(sol), 1);

  for k = 1:length(sol)
    gradu = grad_u_p1(mesh, ints, sol{k});
    zzgradu = zz_grad_u(mesh, gradu);

    est(k) = est(k) + sum(mesh.h .* residual_error(mesh, grad_u_p1(mesh, ints,
							     sol{k}), w)...
			  .* sqrt(error_l2_p0_p1_per_element(mesh,
							     gradu(:, 1),
							     zzgradu(:, 1)).^2
				  + error_l2_p0_p1_per_element(mesh,
							       gradu(:, 2),
							       zzgradu(:, 2)).^2));

  end
  
  time_integral = 0;
  for k = 1:(length(est) - 1)
      time_integral = time_integral + tau * (est(k) + est(k + 1)) / 2;
  end

  eff = sqrt(time_integral) / err;
end


function result = residual_error(mesh, gradu, w_nodes)
  result = zeros(size(mesh.elements, 1), 1);
  for el = 1:size(mesh.elements, 1)
    for k = 1:3
      kk = mod(k, 3) + 1;

      result(el) = result(el) + mesh.jac(el)/6 ...
				* (gradu(el, :) ...
				   * (w_nodes(mesh.elements(el, k), :)' ...
				      + w_nodes(mesh.elements(el, kk), :)')/2)^2;
    end
    result(el) = sqrt(result(el));
  end
end
