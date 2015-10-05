function [mat, rhs] = set_boundary_conditions(mesh, mat, rhs, dbc_f, dbc_nodes)

  [i, j, v] = find(mat);
  for bcn = 1:size(dbc_nodes, 1)
    idx = find(i == dbc_nodes(bcn));
    mat(i(idx), j(idx)) = 0;
    %mat(j(idx), i(idx)) = 0; % if we wanna kill the columns, we need to modify the rhs
    mat(dbc_nodes(bcn), ...
	dbc_nodes(bcn)) = 1;
    rhs(dbc_nodes(bcn)) = dbc_f(mesh.nodes(dbc_nodes(bcn), :));
  end
end
    
