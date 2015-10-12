function x = build_newton_initial_solution(mesh)
  x = zeros(2*(size(mesh.nodes, 1) + size(mesh.elements, 1)) + size(mesh.nodes, 1), 1);
