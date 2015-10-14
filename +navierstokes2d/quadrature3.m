function value = quadrature(function_f)
  points = [0, 0; 1, 0; 0, 1; 0.5, 0; 0, 0.5; 0.5, 0.5; 1/3, 1/3];
  weights = [3, 3, 3, 8, 8, 8, 27] / 60;

  value = 1/2 * weights * function_f(points);
return;


