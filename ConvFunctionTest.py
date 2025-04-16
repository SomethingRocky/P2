import numpy as np
import pandas as pd  # Import pandas for DataFrame manipulation
from main import heatConvCoeff

# Define test parameters
k = 200.0  # Thermal conductivity (W/m·K)
density = 2700.0  # Density of the material (kg/m³)
dx = 0.01  # Grid spacing in x-direction (m)
dy = 0.01  # Grid spacing in y-direction (m)
dz = 0.01  # Grid spacing in z-direction (m)
area = 1.0  # Surface area (m²)

# Calculate h to ensure coeff = 2
target_coeff = 2
h = target_coeff * (density * dx * dy * dz * k) / area

# Other parameters
matrixSize = 5*3*3  # Total matrix size (3x3x5 grid)
matrix2dSize = 9  # 2D matrix size (3x3 grid)
zConfines = (0, 2)  # Confined layers in the z-direction

# Call the function
coeff_matrix = heatConvCoeff(h, k, density, area, matrixSize, matrix2dSize, zConfines, dx, dy, dz)

# Convert the coefficient matrix to a DataFrame and replace 0 with blank spaces
coeff_df = pd.DataFrame(coeff_matrix).replace(0, "")

# Print the DataFrame
print("Heat Convection Coefficient DataFrame:")
print(coeff_df)

