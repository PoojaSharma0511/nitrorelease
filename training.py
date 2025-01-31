import numpy as np
from scipy.linalg import cholesky, cho_solve, solve_triangular,cho_factor
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic
from sklearn.gaussian_process.kernels import RBF

X = np.zeros((196,2))
y = np.zeros(196)

#loading the data
data = np.loadtxt("3NP_train_196_6march.out")
X_orig=data[:,:2]
y_orig=data[:,2]

#normalize the y
y_avg=sum(y_orig)/len(y)
y = y_orig - y_avg
cons = np.sqrt(len(y)/np.sum(y*y))
y = cons*y
print(cons,y_avg)

X_train = X_orig
y_train = y

#training the model
kernel = 1.1**2 * RBF(length_scale=0.1, length_scale_bounds=(1e-2, 10.0))
gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
gaussian_process.fit(X_train, y_train)
print("kernel = ",gaussian_process.kernel_)   # this will print the final kernel with optimized hyperparameters which are used in the prediction program.


