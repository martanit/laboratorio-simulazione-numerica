import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys

file_name = sys.argv[1]
print(file_name)
prova = np.loadtxt(file_name)
plt.hist(prova, 1000)
plt.show()

