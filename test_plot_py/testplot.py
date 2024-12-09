import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Print current working directory
print("Current working directory:", os.getcwd())

# Function to process a CSV file and generate required outputs
def process_csv(file_path):
    print(f"Processing {file_path}")
    # Load the dataset
    df = pd.read_csv(file_path)

    # Ensure required columns exist
    required_columns = {'k', 't', 'SPSS(K)', 'TIME_SELECTING_KMERS', 'TIME_SPSS_CONSTRUCTION', 'TIME_BUILD_FMI'}
    if not required_columns.issubset(df.columns):
        print(f"Skipping {file_path}: Missing required columns.")
        return

    # Generate a table for running times for different k values
    time_table = df.groupby('k')[['TIME_SELECTING_KMERS', 'TIME_SPSS_CONSTRUCTION', 'TIME_BUILD_FMI']].mean()
    time_table = time_table.round(2)  # Round to 2 decimal places
    time_table.reset_index(inplace=True)

    # Save the running time table
    output_file = file_path.replace('.csv', '_running_times.csv')
    time_table.to_csv(output_file, index=False)
    print(f"  - Running time table saved as: {output_file}")

# Generate a combined plot for a specific k across all datasets
def combined_plot(csv_files, k, output_name):
    plt.figure()
    for file_path in csv_files:
        # Load the dataset
        df = pd.read_csv(file_path)

        # Ensure required columns exist
        if 'k' not in df.columns or 't' not in df.columns or 'SPSS(K)' not in df.columns:
            print(f"Skipping {file_path}: Missing required columns.")
            continue

        # Filter data for the specified k
        df_k = df[df['k'] == k]
        if df_k.empty:
            print(f"No data for k={k} in {file_path}")
            continue

        plt.plot(df_k['t'], df_k['SPSS(K)'], marker='o', label=file_path)

    plt.title(f"SPSS(K) vs t for k={k} (All Datasets)")
    plt.xlabel("t (Threshold)")
    plt.ylabel("SPSS(K)")
    plt.legend()
    plt.grid()
    combined_plot_path = output_name
    plt.savefig(combined_plot_path)
    plt.close()
    print(f"Combined plot saved as: {combined_plot_path}")

# Detect and display all CSV files in the current directory
csv_files = glob.glob("*.csv")
print("Detected CSV files:", csv_files)

# Main processing for multiple CSV files
if not csv_files:
    print("No CSV files found in the current directory.")
else:
    for csv_file in csv_files:
        process_csv(csv_file)

    # Generate combined plots for k=21 and k=31
    combined_plot(csv_files, k=21, output_name="combined_spss_k21_plot.png")
    combined_plot(csv_files, k=31, output_name="combined_spss_k31_plot.png")

# Test write permissions
try:
    with open("test_file.txt", "w") as f:
        f.write("Test output")
    print("Write permission test successful. 'test_file.txt' created.")
    os.remove("test_file.txt")  # Cleanup test file
except Exception as e:
    print("Write permission test failed:", e)
