nullModel <- stats::formula(expr ~ 1 + (1 | strain))
interceptModel <- stats::formula(expr ~ 1 + covariate + (1 | strain))
slopeModel <- stats::formula(expr ~ 1 + covariate + (1 + covariate | strain))
