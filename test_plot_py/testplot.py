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

    # Filter data for k=21 and plot SPSS(K) vs t
    df_k21 = df[df['k'] == 21]
    plt.figure()
    plt.plot(df_k21['t'], df_k21['SPSS(K)'], marker='o', label=f"{file_path}")
    plt.title(f"SPSS(K) vs t for k=21 ({file_path})")
    plt.xlabel("t (Threshold)")
    plt.ylabel("SPSS(K)")
    plt.legend()
    plt.grid()
    plot_path = file_path.replace('.csv', '_spss_k21_plot.png')
    plt.savefig(plot_path)
    plt.close()
    print(f"  - Plot saved as: {plot_path}")

    # Generate a table for running times for different k values
    time_table = df.groupby('k')[['TIME_SELECTING_KMERS', 'TIME_SPSS_CONSTRUCTION', 'TIME_BUILD_FMI']].mean()
    time_table = time_table.round(2)  # Round to 2 decimal places
    time_table.reset_index(inplace=True)

    # Save the running time table
    output_file = file_path.replace('.csv', '_running_times.csv')
    time_table.to_csv(output_file, index=False)
    print(f"  - Running time table saved as: {output_file}")

# Detect and display all CSV files in the current directory
csv_files = glob.glob("*.csv")
print("Detected CSV files:", csv_files)

# Main processing for multiple CSV files
if not csv_files:
    print("No CSV files found in the current directory.")
else:
    for csv_file in csv_files:
        process_csv(csv_file)

# Test write permissions
try:
    with open("test_file.txt", "w") as f:
        f.write("Test output")
    print("Write permission test successful. 'test_file.txt' created.")
    os.remove("test_file.txt")  # Cleanup test file
except Exception as e:
    print("Write permission test failed:", e)
