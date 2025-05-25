  # Enable all CPU cores
from tqdm import tqdm
import pandas as pd
import re

# Read the file content
def parse_to_datafame(file_path):
  with open(file_path, "r") as file:
      content = file.read()
  # Define a regex pattern to extract the test case and its details
  pattern = r"\s+structure\s+(?P<structure>.+)\s+nat\s+(?P<nat>\d+)\s+nprim\s+(?P<nprim>\d+)\s+maxl\s+(?P<maxl>\d+)\s+max_nnl\s+(?P<max_nnl>\d+)\s+mean_nnl\s+(?P<mean_nnl>\d+)\s+max_nsh\s+(?P<max_nsh>\d+)\s+cpu_time\s+(?P<cpu_time>[\d.]+)\s+gpu_gb_in\s+(?P<gpu_gb_in>[\d.]+)\s+gpu_gb_out\s+(?P<gpu_gb_out>[\d.]+)\s+gpu_gb_total\s+(?P<gpu_gb_total>[\d.]+)\s+gpu_between_atoms\s+(?P<gpu_between_atoms>[\d.]+)\s+gpu_in_atoms\s+(?P<gpu_in_atoms>[\d.]+)\s+gpu_time\s+(?P<gpu_time>[\d.]+)\s+gpu_transfer_time\s+(?P<gpu_transfer_time>[\d.]+)"

  # Use re.finditer to extract all matches
  matches = re.finditer(pattern, content)

  # Create a list of dictionaries to store the extracted data
  data = []
  for match in matches:
      data.append(match.groupdict())

  # Convert the list of dictionaries into a pandas DataFrame
  df = pd.DataFrame(data)

# Convert numeric columns to appropriate data types
  numeric_columns = [
    'nat', 'nprim', 'maxl', 'max_nnl', 'mean_nnl', 'max_nsh',
    'cpu_time', 'gpu_gb_in', 'gpu_gb_out', 'gpu_gb_total',
    'gpu_between_atoms', 'gpu_in_atoms', 'gpu_time', 'gpu_transfer_time'
  ]
  # # Convert the numeric columns to float
  df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric)

  # Display the DataFrame
  df = df.sort_values('nat')
  return df

import os
import subprocess

def get_gpu_name():
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True, text=True, check=True
        )
        gpu_names = result.stdout.strip().split('\n')
        if gpu_names:
            print("Available GPUs:")
            for idx, name in enumerate(gpu_names):
                print(f"{idx}: {name}")
            idx = 0
            return gpu_names[idx].replace(' ', '_')
        else:
            raise RuntimeError("No GPUs found.")
    except Exception as e:
        print(f"Could not get GPU name: {e}")
        return input("Enter GPU name: ").replace(' ', '_')

if __name__ == "__main__":
    option = get_gpu_name()
    os.makedirs(f"notes/capture/{option}", exist_ok=True)
    os.makedirs(f"notes/data/{option}", exist_ok=True)

    with open(f"notes/capture/{option}/gpu.txt", "w") as outfile:
        # Also write full nvidia-smi's and lscpu's output to `notes/capture/hardware.txt`
        with open(f"notes/capture/{option}/hardware.txt", "w") as hwfile:
            subprocess.run(["nvidia-smi"], stdout=hwfile, stderr=subprocess.STDOUT)
            hwfile.write("\n\n")
            subprocess.run(["lscpu"], stdout=hwfile, stderr=subprocess.STDOUT)

        # Write build results
        subprocess.run(
            ["meson", "test", "-C", "build/", "hamiltonian", "--verbose", "-t", "0"],
            stdout=outfile,
            stderr=subprocess.STDOUT
        )
    df = parse_to_datafame(f"notes/capture/{option}/gpu.txt")
    df.to_csv(f"notes/data/{option}/gpu.csv", index=False)
    print(df)