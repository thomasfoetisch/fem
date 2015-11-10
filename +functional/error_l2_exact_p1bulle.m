function global_integral = error_l2_exact_p1bulle(mesh, f_exact, f_p1bulle)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_dof = n_nodes + n_elems;

  dof_map = [mesh.elements, ((n_nodes + 1):n_dof)'];
  [~, basis] = functional.elementary_integrals_p1bulle();
  quadrature = functional.get_quadrature_tri(9);

  global_integral = 0;
  for el = 1:size(mesh.elements, 1)
    a = mesh.nodes(mesh.elements(el, 1), :);
    T = [mesh.nodes(mesh.elements(el, 2), :)' - a', ...
	 mesh.nodes(mesh.elements(el, 3), :)' - a']';

    local_integral = 0;
    for q = 1:size(quadrature.nodes, 1)
      field_value = 0;
      for i = 1:4
	field_value = field_value + f_p1bulle(dof_map(el, i)) * basis.phi{i}(quadrature.nodes(q, :));
      end
      local_integral = local_integral ...
		       + quadrature.weights(q) ...
			 * (f_exact(a + quadrature.nodes(q, :) * T) ...
			    - field_value)^2;
    end
    global_integral = global_integral + mesh.jac(el) * 1/4 * local_integral;
  end
