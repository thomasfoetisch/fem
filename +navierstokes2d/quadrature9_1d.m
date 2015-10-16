function value = quadrature9_1d(func_f)

  points = [sqrt(3/7 - 2/7 * sqrt(6/5)), -sqrt(3/7 - 2/7 * sqrt(6/5)), sqrt(3/7 + 2/7 * sqrt(6/5)), -sqrt(3/7 + 2/7 * sqrt(6/5))];
  weights = [(18 + sqrt(30))/36, (18 + sqrt(30))/36, (18 - sqrt(30))/36, (18 - sqrt(30))/36];

  value = 1/2 * weights * func_f(points)';
