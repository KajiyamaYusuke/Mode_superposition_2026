import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../output/modeforce.dat")

time = data[:, 0]
fi0 = data[:, 1]
fi1 = data[:, 2]

plt.figure(figsize=(10,4))
plt.plot(time, fi0, label='Mode 1 Force', color='blue')
plt.plot(time, fi1, label='Mode 2 Force', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Modal Force [N]')
plt.title('Modal Force over time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()