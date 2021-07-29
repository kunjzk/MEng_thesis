import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd

with open("control.txt", "rb") as fp:   # Unpickling
        c = pickle.load(fp)

with open("xaxis.txt", "rb") as fp1:   # Unpickling
        x = pickle.load(fp1)

Y_bs_var = pd.DataFrame(c).rolling(30).std().dropna().values
Y_bs_var_mean = np.mean(Y_bs_var)
print(Y_bs_var_mean)
# plt.plot(x, Y_bs_var_mean)
# plt.show()
