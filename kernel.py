import numpy as np
from scipy.linalg import cholesky, cho_solve, solve_triangular,cho_factor
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic
from sklearn.gaussian_process.kernels import RBF

X = np.zeros((196,2))
y = np.zeros(196)

data_train = np.loadtxt("3NP_train_196_6march.out")
X_train = data_train[:,:2]
y_train = data_train[:, 2]

#load test data set
data_test = np.loadtxt("test_nitrophenol_625pts.out")
x_test=data_test[:,:2]
y_test=data_test[:,2]

n=len(X)
print(n)
kern=np.zeros((n,n))
kk=np.zeros(n)
p1=0.643**2
p2=0.102
for i in range(n):
  for j in range(n):
    distsq=np.sum((X_train[i,:]-X_train[j,:])**2/(p2*p2))
    kern[i,j]=np.exp(-0.5*distsq)
    
#generating kernel and alpha matrix to be used for prediction at the time of dynamics
fil = open('kernel_3NP_196_6march.out', 'w')
for i in range(n):
  print(X_train[i,0],X_train[i,1],file=fil)
fil.close()

kern_sav=kern

kern[np.diag_indices_from(kern)] += 1e-10
L=cholesky(kern,lower=True,check_finite=False)
alpha = cho_solve((L,True),y_train,check_finite=False)

fil = open('alpha_3NP_196_6march.out', 'w')
for i in range(len(alpha)):
  print(alpha[i],file=fil)
fil.close()


##testing

m=len(x_test)
y_pred1=np.zeros(m)
for i in range(m):
  for j in range(n):
    distsq=np.sum((x_test[i,:]-X_train[j,:])**2/(p2*p2))
    kk[j]=np.exp(-0.5*distsq)
  y_pred1[i]=np.dot(kk,alpha)

y_pred1 = y_pred1/cons + y_avg


fil = open('test_NP_51x51_196_23feb_predicted.out', 'w')
for i in range(m):
  print(x_test[i,0],x_test[i,1],y_pred1[i], file=fil)
fil.close()

