import math
from matplotlib.ticker import ScalarFormatter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib_venn import venn2
import glob
import matplotlib.patches as patches

folder_path = "CSVs"
CF_Alpha_threshold_low = 0 #(0-1)
AS_Alpha_threshold = 0 #(200-1500)
CF_Beta_threshold_low = 0 #(0-1)
AS_Beta_threshold = 0 #(200-1500)

# Check if the folder path exists
if not os.path.isdir(folder_path):
    raise FileNotFoundError("The provided folder path does not exist.")

# Load and merge CSV files into df_0 and df_6 based on file names
def load_and_merge_csvs(folder_path, keyword):
    # Use glob to find all CSV files with the specified keyword
    file_pattern = os.path.join(folder_path, f'*{keyword}*.csv')
    file_paths = glob.glob(file_pattern)
    
    if not file_paths:
        raise FileNotFoundError(f"No CSV files found with '{keyword}' in their name.")
    
    # Read and concatenate all matching CSV files
    df_list = [pd.read_csv(file) for file in file_paths]
    combined_df = pd.concat(df_list, ignore_index=True)
    
    return combined_df

# Load data
df_0 = load_and_merge_csvs(folder_path, 'D0')
df_6 = load_and_merge_csvs(folder_path, 'D6')
df_0 = df_0[["cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore"]]
df_6 = df_6[["cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore"]]

df_0[['GeneSegment', 'Score']] = df_0['allVHitsWithScore'].str.extract(r'([^()]+)\(([\d.]+)\)')
df_6[['GeneSegment', 'Score']] = df_6['allVHitsWithScore'].str.extract(r'([^()]+)\(([\d.]+)\)')
df_0['Score'] = df_0['Score'].astype(float)
df_6['Score'] = df_6['Score'].astype(float)
df_0 = df_0.dropna(subset=['GeneSegment'])
df_6 = df_6.dropna(subset=['GeneSegment'])

df_0_filtered_Alpha = df_0[df_0['cloneFraction'] > CF_Alpha_threshold_low]
df_6_filtered_Alpha = df_6[df_6['cloneFraction'] > CF_Alpha_threshold_low]
df_0_filtered_Alpha = df_0_filtered_Alpha[df_0_filtered_Alpha['Score'] > AS_Alpha_threshold]
df_6_filtered_Alpha = df_6_filtered_Alpha[df_6_filtered_Alpha['Score'] > AS_Alpha_threshold]

df_0_filtered_Beta = df_0[df_0['cloneFraction'] > CF_Beta_threshold_low]
df_6_filtered_Beta = df_6[df_6['cloneFraction'] > CF_Beta_threshold_low]
df_0_filtered_Beta = df_0_filtered_Beta[df_0_filtered_Beta['Score'] > AS_Beta_threshold]
df_6_filtered_Beta = df_6_filtered_Beta[df_6_filtered_Beta['Score'] > AS_Beta_threshold]


total_TRAV_df_0 = df_0_filtered_Alpha[df_0_filtered_Alpha['GeneSegment'].str.contains('TRA', na=False)].shape[0]
total_TRBV_df_0 = df_0_filtered_Beta[df_0_filtered_Beta['GeneSegment'].str.contains('TRB', na=False)].shape[0]
# Filter and count TRAVs and TRBVs in df_6 (assuming df_6 corresponds to df_1)
total_TRAV_df_6 = df_6_filtered_Alpha[df_6_filtered_Alpha['GeneSegment'].str.contains('TRA', na=False)].shape[0]
total_TRBV_df_6 = df_6_filtered_Beta[df_6_filtered_Beta['GeneSegment'].str.contains('TRB', na=False)].shape[0]
# Print the results
print(f"Total TRAVs in df_0: {total_TRAV_df_0}")
print(f"Total TRBVs in df_0: {total_TRBV_df_0}")
print(f"Total TRAVs in df_6: {total_TRAV_df_6}")
print(f"Total TRBVs in df_6: {total_TRBV_df_6}")


plot_dir = "Plots"
os.makedirs(plot_dir, exist_ok=True)

# Find common GeneSegments in both datasets and remove anomalies
common_segments_Alpha = set(df_0_filtered_Alpha['GeneSegment'].unique()) & set(df_6_filtered_Alpha['GeneSegment'].unique())
common_segments_Alpha = [segment for segment in common_segments_Alpha if pd.notna(segment) and segment.strip()]
combined_df_all_Alpha = pd.concat([df_0_filtered_Alpha, df_6_filtered_Alpha])
frequency_counts_Alpha = combined_df_all_Alpha['GeneSegment'].value_counts()
common_segments_Alpha = sorted(common_segments_Alpha, key=lambda x: frequency_counts_Alpha.get(x, 0), reverse=True)

# Find common GeneSegments in both datasets and remove anomalies
common_segments_Beta = set(df_0_filtered_Beta['GeneSegment'].unique()) & set(df_6_filtered_Beta['GeneSegment'].unique())
common_segments_Beta = [segment for segment in common_segments_Beta if pd.notna(segment) and segment.strip()]
combined_df_all_Beta = pd.concat([df_0_filtered_Beta, df_6_filtered_Beta])
frequency_counts_Beta = combined_df_all_Beta['GeneSegment'].value_counts()
common_segments_Beta = sorted(common_segments_Beta, key=lambda x: frequency_counts_Beta.get(x, 0), reverse=True)



# Function to calculate Gini coefficient
def gini_coefficient(frequencies):
    sorted_freqs = np.sort(frequencies)
    n = len(frequencies)
    cumulative_freqs = np.cumsum(sorted_freqs)
    gini = (2 * np.sum((i + 1) * freq for i, freq in enumerate(sorted_freqs))) / (n * np.sum(sorted_freqs)) - (n + 1) / n
    return gini

# Function to calculate d50 coefficient
def lod10_d50_coefficient(frequencies):
    sorted_freqs = np.sort(frequencies)[::-1]
    cumulative_freqs = np.cumsum(sorted_freqs)
    total_reads = np.sum(sorted_freqs)
    d50 = np.searchsorted(cumulative_freqs, 0.5 * total_reads) + 1
    d50 = np.log(d50 + 1e-10)
    return d50

# Function to calculate Gini-Simpson index
def gini_simpson_index(frequencies):
    total = np.sum(frequencies)
    prob = frequencies / total
    index = 1 - np.sum(prob ** 2)
    return index

def calculate_metrics(df_0, df_6, name):
    results_0 = []
    results_6 = []
    if name == "Alpha":
        common_segments_filtered = [segment for segment in common_segments_Alpha if "TRA" in segment]
    else:
        common_segments_filtered = [segment for segment in common_segments_Beta if "TRB" in segment]
    
    for segment in common_segments_filtered:
        if name == "Alpha":
            freqs_0 = df_0_filtered_Alpha[df_0_filtered_Alpha['GeneSegment'] == segment]['cloneFraction'].values
            freqs_6 = df_6_filtered_Alpha[df_6_filtered_Alpha['GeneSegment'] == segment]['cloneFraction'].values
        else:
            freqs_0 = df_0_filtered_Beta[df_0_filtered_Beta['GeneSegment'] == segment]['cloneFraction'].values
            freqs_6 = df_6_filtered_Beta[df_6_filtered_Beta['GeneSegment'] == segment]['cloneFraction'].values
            
        gini_0 = gini_coefficient(freqs_0)
        gini_6 = gini_coefficient(freqs_6)

        d50_0 = lod10_d50_coefficient(freqs_0)
        d50_6 = lod10_d50_coefficient(freqs_6)

        simpson_0 = gini_simpson_index(freqs_0)
        simpson_6 = gini_simpson_index(freqs_6)

        results_0.append({'GeneSegment': segment, 'Gini Simp Coef': gini_0, 'log10 D50 Idx': d50_0, 'Gini Simp Idx': simpson_0, 'Timepoint': 'CTL'})
        results_6.append({'GeneSegment': segment, 'Gini Simp Coef': gini_6, 'log10 D50 Idx': d50_6, 'Gini Simp Idx': simpson_6, 'Timepoint': 'P3'})

    combined_results = pd.DataFrame(results_0 + results_6)
    return combined_results

def plot_metrics_all(combined_results_alpha, combined_results_beta, plot_dir):
    metrics = [['Gini Simp Coef', 'Gini Diversity'], ['log10 D50 Idx', 'D50 Diversity'], ['Gini Simp Idx', 'Gini Simpson']]
    
    # Set font properties
    plt.rcParams.update({
        'font.size': 20,              # Base font size
        'axes.titlesize': 26,         # Title font size
        'axes.labelsize': 24,         # Axis labels font size
        'xtick.labelsize': 24,        # X-tick labels font size
        'ytick.labelsize': 24,        # Y-tick labels font size
        'legend.fontsize': 20,        # Legend font size
        'legend.title_fontsize': 20   # Legend title font size
    })
    
    custom_palette = {
        'Alpha': '#87CEEB',   # Sky blue color in hexadecimal
        'Beta': '#90EE90',    # Light green color in hexadecimal
        'CTL': '#d3d8dd',     # Gray color in hexadecimal
        'P3': '#65a7e2'       # Blue color in hexadecimal
    }
    all_combined_results = []
    
    # Combine Alpha and Beta results for plotting
    combined_results_alpha['Type'] = 'Alpha'
    combined_results_beta['Type'] = 'Beta'
    combined_results_combined = pd.concat([combined_results_alpha, combined_results_beta])
    combined_results_combined.to_excel(f'{plot_dir}/Gene_Segment_Diversity_Metrics.xlsx', index=False)
    
    for metric, y_title in metrics:
        plt.figure(figsize=(6, 9))
        
        # Create the combined plot for Alpha and Beta
        ax = sns.boxplot(data=combined_results_combined, x='Type', y=metric, hue='Timepoint', palette=custom_palette, width=0.25, dodge=True, showfliers=False)
        
        plt.title(f'{y_title}')
        plt.ylabel(f'{metric}')
        plt.xlabel('')
        plt.legend(loc='lower left', fontsize='large', title_fontsize='17')

        plt.tight_layout()
        
        # Save the plot
        plt.savefig(f'{plot_dir}/{metric}.png', dpi=400)
        plt.close()

# Calculate metrics for Alpha and Beta
combined_results_Alpha = calculate_metrics(total_TRAV_df_0, total_TRAV_df_6, "Alpha")
combined_results_Beta = calculate_metrics(total_TRBV_df_0, total_TRBV_df_6, "Beta")

# Plot and save the results
plot_metrics_all(combined_results_Alpha, combined_results_Beta, plot_dir)

