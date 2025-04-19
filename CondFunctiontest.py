import numpy as np
import pandas as pd
from main import *
import time  # Import the time module

# Start timing
#start_time = time.time()

x = dTdx(2,
               3*3*3,
               3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

y = dTdy(2,
               3*3*3,
               3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

z = dTdz(2,
               3*3*3,
               3*3,
               {"x": 3, "y": 3, "z": 3},
               (0,2))

# End timing
#end_time = time.time()

# Replace zeroes with blank spaces
x_df = pd.DataFrame(x).replace(0, "")
y_df = pd.DataFrame(y).replace(0, "")
z_df = pd.DataFrame(z).replace(0, "")
T_df = pd.DataFrame(x+y+z).replace(0, "")

print("\ndT/dx")
#print(x_df)
print("\ndT/dy")
#print(y_df)
print("\ndT/dz")
#print(z_df)
print("\ndT/dt")
#print(T_df)


# Print the elapsed time
#print(f"\nExecution time: {end_time - start_time:.4f} seconds")
