import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import pickle
from matplotlib.gridspec import GridSpec
from scipy.stats import lognorm

plt.close('all')
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,2,1)

s = 0.954
mean, var, skew, kurt = lognorm.stats(s, moments='mvsk')
x = np.linspace(lognorm.ppf(0.01, s),
                lognorm.ppf(0.99, s), 100)
# ax.plot(x, lognorm.pdf(x, s),
#        'r-', lw=5, alpha=0.6, label='lognorm pdf')
r = lognorm.rvs(s, size=1000)
ax.hist(r, density=True, bins='auto', histtype='stepfilled', alpha=0.2)
ax.set_xlim([x[0], x[-1]])
ax.legend(loc='best', frameon=False)

# ax.scatter(ESP['df']['datetime'],ESP['df']['PB_quantity'])

ax2 = fig.add_subplot(1,2,2)
ax2.hist(np.log(r), density=True, bins='auto', histtype='stepfilled', alpha=0.2)

# ax.set_yscale('log')
plt.show(block=False)
plt.pause(0.1)
