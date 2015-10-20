function x = build_newton_initial_solution(mesh, force_f, ctx, boundary_conditions)

  [u_x, u_y, p] = stokes2d.solve(mesh, ctx.laminar_viscosity, force_f, boundary_conditions);
  x = [u_x; u_y; p];
