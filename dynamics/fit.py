import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the decaying exponential function
def decay_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Load data from file (replace 'data_file.txt' with your file name)
data = np.loadtxt('pop.out')

# Assuming the data has two columns (x and y)
x_data = data[:, 0]
y_data = data[:, 1]

# Perform the curve fitting
popt, pcov = curve_fit(decay_func, x_data, y_data)

# Get the fitted parameters
a_fit, b_fit, c_fit = popt

# Generate points for the fitted curve
x_fit = np.linspace(min(x_data), max(x_data), 100)
y_fit = decay_func(x_fit, a_fit, b_fit, c_fit)

# Plot the original data and the fitted curve
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_fit, y_fit, 'r-', label='Fitted Curve')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Decaying Exponential Curve Fitting')
plt.legend()
plt.grid(True)
plt.show()

# Print the fitted parameters
print("Fitted Parameters:")
print("a =", a_fit)
print("b =", b_fit)
print("c =", c_fit)

f=open('fitted.txt','w')
print(a_fit,b_fit,c_fit, file=f)
f.close()
