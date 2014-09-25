function [PDF_PHI_inv] = PDF_PHI_inv(x)

	% inverse of probability density function of standard normal distribution with mu = 0 and sigma = 1

	PDF_PHI_inv = (-2*log(x*sqrt(2*pi)))^(1/2);

end
