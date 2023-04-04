# exponential autocorrelation function
rho_exp = function(H,gamma_2)
{
	exp(-H/gamma_2)
}

# spherical spatial autocorrelation function
rho_sph = function(H,gamma_2)
{
	(1 - 1.5*H/gamma_2 + 0.5*(H/gamma_2)^3)*(H/gamma_2 < 1)
}

# Gaussian spatial autocorrelation function
rho_gau = function(H,gamma_2)
{
	exp(-(H/gamma_2)^2)
}

# rational quadratic spatial autocorrelation function
rho_rat = function(H,gamma_2)
{
	1/(1+(H/gamma_2)^2)
}

# AR1 time series model
rho_AR1 = function(H,gamma_2)
{
	gamma_2^H
}

# spatial CAR model
rho_CAR = function(W, M, gamma_2)
{
	solve(diag(dim(W)[1]) - gamma_2*W, M)
}

# spatial SAR model
rho_SAR = function(W, M, gamma_2)
{
	solve((diag(dim(W)[1]) - gamma_2*W) %*% (diag(dim(W)[1]) - gamma_2*t(W)))
}

# uniform kernel function
kern_unif = function(H, gamma_2)
{
	(H/gamma_2 < 1)*1
}

# Epanechnikov kernel function
kern_epi = function(H, gamma_2)
{
	(1 - (H/gamma_2)^2)*(H/gamma_2 < 1)
}

# Epanechnikov kernel function
kern_quad = function(H, gamma_2)
{
	(1 - (H/gamma_2)^2)^2*(H/gamma_2 < 1)
}

# Epanechnikov kernel function
kern_tri = function(H, gamma_2)
{
	(1 - (H/gamma_2)^2)^3*(H/gamma_2 < 1)
}

# Gaussian kernel function
kern_gau = function(H, gamma_2)
{
	exp(-(H/gamma_2)^2)
}
