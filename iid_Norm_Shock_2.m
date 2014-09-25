% This code is based Makoto Nakajima Notes
% We want to approximate an iid process (shock)
% with z(i) and p(i) where z is realization and p is the probability of that realization
% in this code we are using Ada & Cooper (2003) method for approximation

% create random variable we want to approximate
% This code approximates mean of eps

pd = makedist('Normal');

time = 1000;

eps = zeros(1,time);

for t = 1:time

	eps(:,t) = random(pd);

end

% 1. Set n, number of realizations of z(i)

n      = 1000;

mu     = mean(eps); % mean of normal distribution

sigma_z = 1;	% variance of normal distribution

% 2. Construct m(i) using the inverse CDF of Normal distribution of N(0,1)

m = zeros(n,1);

for k = 1:n

	m(k,1) = (norminv(k/n))*sigma_z + mu;

end

% 3. Compute the first and the last z: z(1) and z(n)

z = zeros(n,1);

z(1) = mu - sigma_z*n*PDF_PHI((m(1)-mu)/sigma_z);

z(n) = mu + sigma_z*n*PDF_PHI((m(n-1)-mu)/sigma_z);

% 4. Construct rest of the z(i)'s

for i = 2:n-1

	z(i) = mu - sigma_z*n*(PDF_PHI((m(i)-mu)/sigma_z)-PDF_PHI((m(i-1)-mu)/sigma_z));

end

% 5. Construct p

p = zeros(n,1);

for i = 1:n

	p(i,1) = 1/n;

end

% calculate the mean of the approximation process

mean_approx = sum(z.*p);

% compute variance

for k = 1:length(z)

	var_approx = sum(p(k,1)*(z(k,1)- mean_approx)^2);

end

% check the approximation of mean

if abs(mean(eps) - mean_approx) < 10e-3

	disp 'Mean Approximated Correctly'

else

	disp 'Mean Approximated Incorectly' 

end


