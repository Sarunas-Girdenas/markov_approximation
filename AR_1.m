% This code is based Makoto Nakajima Notes
% We want to approximate an iid process (shock)
% with z(i) and p(i) where z is realization and p is the probability of that realization

% very good resource of AR(1) approximations is http://people.su.se/~mflod/code.htm

% this code approximates mean of AR(1) process r(t) = (1-rho)*mu + rho*r(t-1) + eps(t+1) 

pd = makedist('Normal');

time = 100;

eps = zeros(1,time);

for t = 1:time + 1

	eps(:,t) = random(pd);

end

% 1. Set n, number of realizations of z(i)

n      = 100;

rho    = 0.9;

lambda = 3;

mu      = mean(eps);	% since we know that mean of Normal Distribution is 0

sigma_z = std(eps)/(sqrt(1-rho^2));	% variance of Normal Distribution is 1

% AR(1) process we wish to approximate

r = zeros(1,time);

r(1) = 1;

for k = 2:time

	r(k) = (1-rho)*mu + rho*r(k-1) + eps(k+1); 

end

% 2. Set the upper and lower bounds of z

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

p = zeros(n,n); % p(i,j)

% initialize probability, j=1

for i = 1:n

	p(i,1) = normcdf((m(1) - (1-rho)*mu - rho*z(i))/sigma_z);

	p(i,n) = 1 - normcdf((m(n-1) - (1-rho)*mu - rho*z(i))/sigma_z);

end

% probabilities in between


for i = 1:n

	for j = 2:n-1

		p(i,j) = normcdf((m(j,1) - (1-rho)*mu - rho*z(i,1))/sigma_z) - normcdf((m(j-1,1) - (1-rho)*mu - rho*z(i,1))/sigma_z);

	end

end
