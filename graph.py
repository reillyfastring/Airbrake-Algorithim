import pandas as pd
import matplotlib.pyplot as plt

# Read the file
# The header is space-separated while the data rows are comma-separated.
with open("Data.txt", "r") as f:
    header_line = f.readline().strip()
    col_names = header_line.split()  # Splits the header by whitespace

# Read the remaining data using the correct delimiter (comma)
df = pd.read_csv("Data.txt", skiprows=1, header=None, names=col_names, sep=",")

# Plot 1: Regular and Adjusted Position vs. Time
plt.figure(figsize=(10, 6))
plt.scatter(df['Time'], df['Position'], color='blue', label='Regular Position')
plt.scatter(df['Time'], df['AdjustedPosition'], color='red', label='Adjusted Position')
plt.xlabel("Time")
plt.ylabel("Position")
plt.title("Regular vs. Adjusted Position vs. Time")
plt.legend()
plt.grid(True)
plt.show()

# Plot 2: Regular and Adjusted Velocity vs. Time
plt.figure(figsize=(10, 6))
plt.scatter(df['Time'], df['Velocity'], color='blue', label='Regular Velocity')
plt.scatter(df['Time'], df['AdjustedVelocity'], color='red', label='Adjusted Velocity')
plt.xlabel("Time")
plt.ylabel("Velocity")
plt.title("Regular vs. Adjusted Velocity vs. Time")
plt.legend()
plt.grid(True)
plt.show()

# Plot 3: Regular and Adjusted Angle vs. Time
plt.figure(figsize=(10, 6))
plt.scatter(df['Time'], df['Angle'], color='blue', label='Regular Angle')
plt.scatter(df['Time'], df['AdjustedAngle'], color='red', label='Adjusted Angle')
plt.xlabel("Time")
plt.ylabel("Angle")
plt.title("Regular vs. Adjusted Angle vs. Time")
plt.legend()
plt.grid(True)
plt.show()

# Plot 4: DeploymentAngle vs. Time
plt.figure(figsize=(10, 6))
plt.scatter(df['Time'], df['DeploymentAngle'], color='green', label='Deployment Angle')
plt.xlabel("Time")
plt.ylabel("Deployment Angle")
plt.title("Deployment Angle vs. Time")
plt.legend()
plt.grid(True)
plt.show()
