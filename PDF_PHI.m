function [PDF_PHI] = PDF_PHI(x)

	% probability density function of standard normal distribution with mu = 0 and sigma = 1

	PDF_PHI = (exp((-1/2)*x^2)) / sqrt(2*pi);

end

