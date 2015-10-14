function value = quadrature5(function_f)

  s15 = sqrt(15);
  a = (6 - s15) / 21;
  b = (6 + s15) / 21;
  c = (9 - 2 * s15) / 21;
  d = (9 + 2 * s15) / 21;

  points = [1/3, 1/3,
	    a, a,
	    a, d,
	    d, a,
	    b, b,
	    d, c,
	    c, d];

  w1 = (155 - s15) / 1200;
  w2 = (155 + s15) / 1200;
  weights = [0.225, w1, w1, w1, w2, w2, w2];

  value = 1/2 * weights * function_f(points);


