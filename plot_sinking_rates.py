import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

size = [2, 3, 3, 4, 5, 17, 4, 11, 1, 4]
size_err = [1, 1, 2, 2, 4, 11, 2, 7, 0, 5]
sink = [53, 60, 52, 95, 78, 95, np.nan, 43, 78, 60]
sink_err = [17, 26, 17, 121, 17, 26, np.nan, 26, 60, 17]

fig = plt.figure(figsize=(8,8))
ax = fig.gca()

ax.plot(size,sink,marker='o',linestyle='none')
ax.set_xlabel('Longest dimension (mm)')
ax.set_ylabel('Sinking rates (m day-1)')

plt.show(block=False)
plt.pause(0.1)
# plt.close()


