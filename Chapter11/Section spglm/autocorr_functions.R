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
