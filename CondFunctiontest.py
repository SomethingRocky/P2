import numpy as np
import pandas as pd
from main import *

x = xCondCoeff(2,
               3*3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

y = yCondCoeff(2,
               3*3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

z = zCondCoeff(2,
               3*3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

# Replace zeroes with blank spaces
x_df = pd.DataFrame(x).replace(0, "")
y_df = pd.DataFrame(y).replace(0, "")
z_df = pd.DataFrame(z).replace(0, "")
T_df = pd.DataFrame(x+y+z).replace(0, "")

print("\ndT/dx")
print(x_df)
print("\ndT/dy")
print(y_df)
print("\ndT/dz")
print(z_df)
print("\ndT/dt")
print(T_df)
