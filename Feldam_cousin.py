import gammapy.stats
import numpy as np
from scipy import stats

x_bins = np.arange(0, 50)
mu_bins = np.linspace(0, 15, 15 / 0.005 + 1, endpoint=True)
matrix = [stats.poisson(mu + 3.0).pmf(x_bins) for mu in mu_bins]
acceptance_intervals = gstats.fc_construct_acceptance_intervals_pdfs(matrix, 0.9)
LowerLimitNum, UpperLimitNum, _ = gstats.fc_get_limits(mu_bins, x_bins, acceptance_intervals)

print(gstats.fc_find_limit(1, UpperLimitNum, mu_bins))
