
# This is the script that approximate iid random shock
# This is Python version of the code 'iid_Norm_Shock.m'

from scipy.stats import norm
import random
import numpy as np

# create random variable we want to approximate

time = 2000

eps = np.zeros([1,time])

for t in range(1,time):

	eps[:,t] = random.random()

# 1. Set the number of realizations of z(i)

n = 2000

# 2. Set the upper and lower bounds of z(i)

lambdaa = 3

mu = np.mean(eps)

sigma_z = 1;

z_up = mu + lambdaa * sigma_z

z_low = mu - lambdaa * sigma_z	

# 3. Make sure that z(i)'s are evenly spaced

z = np.zeros([n,1])

for i in range(0,n):

	z[i] = z_low + 2*(1.0/(n-1))*lambdaa*sigma_z*(i-1)

# 4. Construct midpoints

m = np.zeros([n-1,1])

for t in range(0,n-1):

	m[t] = (z[t+1] + z[t])/2.0

# 5. Compute probabilities for corresponding z(i)

p = np.zeros([n,1])

# First probability

p[0] = norm.cdf((m[0] - mu)/sigma_z)

# Last probability

p[n-1] = 1 - norm.cdf((m[n-2] - mu)/sigma_z) # since python counts from 0 

# Probabilities in between

for i in range(1,n-2):

	p[i] = norm.cdf((m[i] - mu)/sigma_z) - norm.cdf((m[i-1] - mu)/sigma_z)

# Calculate the mean of the approximation process

mean_approx = np.sum(z*p)

if np.abs(np.mean(eps) - mean_approx) < 10e-3:

	print 'Mean Approximated Correctly'

else:

	print 'Mean Approximated Incorectly'

