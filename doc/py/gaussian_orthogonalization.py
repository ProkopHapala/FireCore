import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu)/sigma)**2) / (sigma * np.sqrt(2*np.pi))

# Parameters
x = np.linspace(-5, 15, 1000)
x0, sigma1 = 0.0, 1.0
x1, sigma2 = 3.0, 1.5

# Create two gaussians
phi1 = gaussian(x, x0, sigma1)
phi2 = gaussian(x, x1, sigma2)

# Orthogonalize phi2 against phi1
phi2_ortho = phi2 - np.dot(phi2, phi1)*phi1/np.dot(phi1, phi1)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x, phi1,       'r-', label='Original φ₁ (red)')
plt.plot(x, phi2,       'b-', label='Original φ₂ (blue)')
plt.plot(x, phi2_ortho, 'g-', label='Orthogonalized φ₂ (green)')
plt.title('Gaussian Orthogonalization')
plt.legend()
plt.grid(True)
plt.savefig('gaussian_orthogonalization.png')
plt.show()
