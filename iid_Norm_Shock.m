% This code is based Makoto Nakajima Notes, http://www.compmacro.com/makoto/comp.php
% We want to approximate an iid process (shock)
% with z(i) and p(i) where z is realization and p is the probability of that realization

% create random variable we want to approximate

% this is still work in progress somehow and there might be some bugs in the code

% this code approximates mean of eps (random variable)

pd = makedist('Normal');

time = 2000;

eps = zeros(1,time);

for t = 1:time

	eps(:,t) = random(pd);

end

% 1. Set n, number of realizations of z(i)

n = 2000;

% 2. Set upper and lower bounds of z(i)

lambda = 3;

mu     = 0.001; % since we know that mean of Normal Distribution is 0

sigma_z = 1;	% variance of Normal Distribution is 1

z_up  = mu + lambda * sigma_z;

z_low = mu - lambda * sigma_z;

% 3. Make sure that z(i)'s are evenly spaced

z = zeros(n,1);

for i = 1:n

	z(i,1) = z_low + 2*(1/(n-1))*lambda*sigma_z*(i-1);

end

% 4. Construct midpoints

m = zeros(n-1,1);

for t = 1:n-1

	m(t,1) = (z(t+1,1) + z(t,1))/2;

end

% 5. Compute the probabilities for corresponding z

% First probability

p = zeros(n,1);

p(1,1) = normcdf((m(1) - mu)/sigma_z);

% Last one

p(n,1) = 1 - normcdf((m(n-1)-mu)/sigma_z);

% probabilities in between

for i = 2:n-1

	p(i,1) = normcdf((m(i,1)-mu)/sigma_z) - normcdf((m(i-1,1)-mu)/sigma_z);

end

% calculate the mean of the approximation process

mean_approx = sum(z.*p);

% compute variance

for k = 1:length(z)

	var_approx = sum(p(k,1)*(z(k,1)- mean_approx)^2);

end

% check the approximation of mean

if abs(mean(eps)-mean_approx) < 10e-3

	disp 'Mean Approximated Correctly'

else

	disp 'Mean Approximated Incorectly' 

end


