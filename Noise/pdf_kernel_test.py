"""Created on Tue Nov 26 13:52:43 2019 @author: dadhikar """

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

np.random.seed(0)
gaussian1 = -6 + 3 * np.random.randn(1700)
gaussian2 = 4 + 1.5 * np.random.randn(300)
gaussian_mixture = np.hstack([gaussian1, gaussian2])

df = pd.DataFrame(gaussian_mixture, columns=['data'])

# non-parametric pdf
nparam_density = stats.kde.gaussian_kde(df.values.ravel())
x = np.linspace(-20, 20, 200)
nparam_density = nparam_density(x)

# parametric fit: assume normal distribution
loc_param, scale_param = stats.norm.fit(df)
param_density = stats.norm.pdf(x, loc=loc_param, scale=scale_param)

fig, ax = plt.subplots(figsize=(10, 6))

# ax.hist(df.values, bins=30, normed=True)
ax.plot(x, nparam_density, 'r-', label=r'non-parametric density')
ax.plot(x, param_density, 'k--', label='parametric density')
ax.set_ylim([0, 0.15])
ax.legend(loc='best')
plt.show()
