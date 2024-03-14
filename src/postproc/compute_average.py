import os

import numpy as np


# Function to calculate average and standard deviation for each column
def calculate_stats(data):
    averages = np.mean(data, axis=0)
    std_devs = np.std(data, axis=0)
    return averages, std_devs


# Read data from the file
file_path = os.environ.get("FILE")  # Replace with the actual path to your file
with open(file_path, "r") as file:
    lines = file.readlines()

# Parse data from the lines
data = []
for line in lines:
    values = [float(val) for val in line.split()]
    data.append(values)

# Convert data to a NumPy array for easy calculations
data = np.array(data)

# Calculate averages and standard deviations
averages, std_devs = calculate_stats(data)

# Print the results
print(
    averages[0],
    std_devs[0],
    averages[1],
    std_devs[1],
    averages[2],
    std_devs[2],
    averages[3],
    std_devs[3],
    averages[4],
    std_devs[4],
)
