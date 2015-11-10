mesh = geometry.build_square_mesh(1, 1, 40, 40, 0);

% system sizes:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_dof = n_nodes + n_elems;

[~, basis] = functional.elementary_integrals_p1bulle();

quadrature = functional.get_quadrature_tri(9);

u_x_exact = @(x) - cos(pi * x(:, 1)) .* sin(pi * x(:, 2));
u_y_exact = @(x)   sin(pi * x(:, 1)) .* cos(pi * x(:, 2));

dof_map = [mesh.elements, ((n_nodes + 1):n_dof)'];


f_analytic = u_x_exact;
field = zeros(n_dof, 1);
field(1:n_nodes) = f_analytic(mesh.nodes);

global_integral = 0;
for el = 1:size(mesh.elements, 1)
  a = mesh.nodes(mesh.elements(el, 1), :);
  T = [mesh.nodes(mesh.elements(el, 2), :)' - a', ...
       mesh.nodes(mesh.elements(el, 3), :)' - a']';

  local_integral = 0;
  for q = 1:size(quadrature.nodes, 1)
    field_value = 0;
    for i = 1:4
      field_value = field_value + field(dof_map(el, i)) * basis.phi{i}(quadrature.nodes(q, :));
    end
    local_integral = local_integral ...
		       + quadrature.weights(q) ...
			 * (f_analytic(a + quadrature.nodes(q, :) * T) ...
			    - 0*field_value);
  end
    global_integral = global_integral + mesh.jac(el) * 1/4 * local_integral;
end
