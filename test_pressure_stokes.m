v = zeros(2 * n_u_dof + n_nodes, 1);
p_nodes = mesh.nodes(:, 1);
for el = 1:n_elems
  for i = 1:4
    for p = 1:3
      for k = 1:2
	v(dof_map(el, i), 1) = v(dof_map(el, i), 1) - mesh.jac(el) * p_nodes(dof_map(el, p)) * mesh.jmt(1, k, el) * ints.phidphi(p, k, i);
	v(n_u_dof + dof_map(el, i), 1) = v(n_u_dof + dof_map(el, i), 1) - mesh.jac(el) * p_nodes(dof_map(el, p)) * mesh.jmt(2, k, el) * ints.phidphi(p, k, i);
      end
    end
  end
end

for sol_dim = 1:2
  for n = 1:length(mesh.boundary_nodes)
    v((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), 1) = 0;
  end
end

[v, full(rhs)];



x_exact = [zeros(2*n_u_dof, 1); p_nodes];
