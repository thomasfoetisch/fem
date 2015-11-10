function quadrature = get_quadrature_tri(order)

  pw = [0.0451890097844 0.0451890097844 0.0519871420646
	0.0451890097844 0.9096219804312 0.0519871420646
	0.9096219804312 0.0451890097844 0.0519871420646
	0.7475124727339 0.0304243617288 0.0707034101784
	0.2220631655373 0.0304243617288 0.0707034101784
	0.7475124727339 0.2220631655373 0.0707034101784
	0.2220631655373 0.7475124727339 0.0707034101784
	0.0304243617288 0.7475124727339 0.0707034101784
	0.0304243617288 0.2220631655373 0.0707034101784
	0.1369912012649 0.2182900709714 0.0909390760952
	0.6447187277637 0.2182900709714 0.0909390760952
	0.1369912012649 0.6447187277637 0.0909390760952
	0.2182900709714 0.6447187277637 0.0909390760952
	0.2182900709714 0.1369912012649 0.0909390760952
	0.6447187277637 0.1369912012649 0.0909390760952
	0.0369603304334 0.4815198347833 0.1032344051380
	0.4815198347833 0.0369603304334 0.1032344051380
	0.4815198347833 0.4815198347833 0.1032344051380
	0.4036039798179 0.1927920403641 0.1881601469167
	0.4036039798179 0.4036039798179 0.1881601469167
	0.1927920403641 0.4036039798179 0.1881601469167];

  quadrature = struct();
  quadrature.nodes = pw(:, [1 2]);
  quadrature.weights = pw(:, 3);
