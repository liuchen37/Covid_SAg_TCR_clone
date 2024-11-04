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
CF_Alpha_threshold_low = 1e-4 #(0-1)
AS_Alpha_threshold = 0 #(200-1500)
CF_Beta_threshold_low = 2e-5 #(0-1)
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

# Create box plots in batches of 10 Gene Segments per figure
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, ks_2samp, kruskal
from statsmodels.stats.weightstats import ttest_ind as ttest_ind_welch

def plot_GeneSegment_save(name):
    if name == "Alpha":
        common_segments_filtered = [segment for segment in common_segments_Alpha if "TRA" in segment]
    else:
        common_segments_filtered = [segment for segment in common_segments_Beta if "TRB" in segment]
    
    batch_size = int(np.ceil(len(common_segments_filtered)/3.0))
    num_batches = int(np.ceil(len(common_segments_filtered) / batch_size))
    
    plt.rcParams.update({
        'font.size': 14,               # Base font size
        'axes.titlesize': 16,          # Title font size
        'axes.labelsize': 16,          # Axis labels font size
        'xtick.labelsize': 14,         # X-tick labels font size
        'ytick.labelsize': 16,         # Y-tick labels font size
        'legend.fontsize': 16,         # Legend font size
        'legend.title_fontsize': 16    # Legend title font size
    })
    
    custom_palette = {
        'CTL': '#d3d8dd',  # Changed color for CTL
        'P3': '#65a7e2'   # Changed color for SIM
    }
    def get_max_ydata(ax, box_index):
        max_value = float('-inf')
        for line in ax.lines[10*box_index:10*(box_index+1)]:
            y_data = line.get_ydata()
            max_value = max(max_value, max(y_data))
        return max_value
    def print_boxplot_stats(data, label):
        """Prints the boxplot statistics for the given data."""
        # Ensure data is a numpy array
        data = np.array(data)
    
        # Calculate key statistics
        q1 = np.percentile(data, 25)  # 25th percentile
        q3 = np.percentile(data, 75)  # 75th percentile
        median = np.median(data)       # Median
        iqr = q3 - q1                 # Interquartile range
        lower_whisker = max(data.min(), q1 - 1.5 * iqr)  # Lower whisker
        upper_whisker = min(data.max(), q3 + 1.5 * iqr)  # Upper whisker
    
        # Print the statistics
        print(f"Statistics for {label}:")
        print(f"  Q1 (25th percentile): {q1}")
        print(f"  Median: {median}")
        print(f"  Q3 (75th percentile): {q3}")
        print(f"  IQR (Interquartile Range): {iqr}")
        print(f"  Lower Whisker: {lower_whisker}")
        print(f"  Upper Whisker: {upper_whisker}")
        print("")  # For readability
    stats_list = []
    for batch_num in range(num_batches):
        plt.figure(figsize=(15, 6))
    
        # Select a subset of gene segments for the current batch
        segments_batch = common_segments_filtered[batch_num * batch_size: (batch_num + 1) * batch_size]
    
        # Filter data for the selected segments
        if name == "Alpha":
            df_0_batch = df_0_filtered_Alpha[df_0_filtered_Alpha['GeneSegment'].isin(segments_batch)]
            df_6_batch = df_6_filtered_Alpha[df_6_filtered_Alpha['GeneSegment'].isin(segments_batch)]
        else:
            df_0_batch = df_0_filtered_Beta[df_0_filtered_Beta['GeneSegment'].isin(segments_batch)]
            df_6_batch = df_6_filtered_Beta[df_6_filtered_Beta['GeneSegment'].isin(segments_batch)]

        # Combine the data for plotting
        combined_df = pd.concat([
            df_0_batch.assign(Timepoint='CTL'),
            df_6_batch.assign(Timepoint='P3')
        ])
    
        # Create a box plot for the current batch
        ax = sns.boxplot(data=combined_df, x='GeneSegment', y='cloneFraction', hue='Timepoint', palette=custom_palette, width=0.35, 
                         showfliers=False, order=segments_batch)
            # Hide top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #plt.title(f'Clone Fraction for {name} Gene Segments')
        plt.xlabel('')  # Remove x-axis title
        plt.ylabel('Clonal Fraction')
        plt.xticks(rotation=65)
        
        y_ticks = ax.get_yticks()
        def calculate_upper_whisker(data):
            q3 = np.percentile(data, 75)  # 75th percentile
            iqr = np.percentile(data, 75) - np.percentile(data, 25)  # Interquartile range
            upper_whisker = q3 + 1.5 * iqr  # Upper whisker
            return min(np.max(data), upper_whisker)
        
        for i, segment in enumerate(segments_batch):
            if name == "Alpha":
                df_0_segment = df_0_filtered_Alpha[df_0_filtered_Alpha['GeneSegment'] == segment]['cloneFraction']
                df_6_segment = df_6_filtered_Alpha[df_6_filtered_Alpha['GeneSegment'] == segment]['cloneFraction']
            else:
                df_0_segment = df_0_filtered_Beta[df_0_filtered_Beta['GeneSegment'] == segment]['cloneFraction']
                df_6_segment = df_6_filtered_Beta[df_6_filtered_Beta['GeneSegment'] == segment]['cloneFraction']
                
            print_boxplot_stats(df_0_segment, segment +' df_0_segment')
            print_boxplot_stats(df_6_segment, segment +' df_1_segment')
            if len(df_0_segment) > 0 and len(df_6_segment) > 0:
                
                stat_mannwhitneyu, p_value_mannwhitneyu = mannwhitneyu(df_0_segment, df_6_segment, alternative='two-sided')
                
                q1_0 = np.percentile(df_0_segment, 25)
                median_0 = np.median(df_0_segment)
                q3_0 = np.percentile(df_0_segment, 75)
                iqr_0 = q3_0 - q1_0
                lower_whisker_0 = max(df_0_segment.min(), q1_0 - 1.5 * iqr_0)
                upper_whisker_0 = min(df_0_segment.max(), q3_0 + 1.5 * iqr_0)

                q1_6 = np.percentile(df_6_segment, 25)
                median_6 = np.median(df_6_segment)
                q3_6 = np.percentile(df_6_segment, 75)
                iqr_6 = q3_6 - q1_6
                lower_whisker_6 = max(df_6_segment.min(), q1_6 - 1.5 * iqr_6)
                upper_whisker_6 = min(df_6_segment.max(), q3_6 + 1.5 * iqr_6)
                
                stats_list.append({
                    'GeneSegment': segment,
                    'Q1_df_0': q1_0,
                    'Median_df_0': median_0,
                    'Q3_df_0': q3_0,
                    'Lower_Whisker_df_0': lower_whisker_0,
                    'Upper_Whisker_df_0': upper_whisker_0,
                    'Q1_df_6': q1_6,
                    'Median_df_6': median_6,
                    'Q3_df_6': q3_6,
                    'Lower_Whisker_df_6': lower_whisker_6,
                    'Upper_Whisker_df_6': upper_whisker_6,
                    'p_value_Mann-Whitney': p_value_mannwhitneyu,
                })

                #select between P-Values
                p_value = p_value_mannwhitneyu
                # Add stars based on p-value thresholds
                if p_value < 0.0001:
                    significance = '****'
                elif p_value < 0.001:
                    significance = '***'
                elif p_value < 0.01:
                    significance = '**'
                elif p_value < 0.05:
                    significance = '*'
                else:
                    significance = None

                y_min, y_max = ax.get_ylim()
                if significance:
                    upper_whisker = get_max_ydata(ax, i)

                    # Add significance stars at the calculated y_position
                    plt.text(i, upper_whisker, significance, ha='center', va='bottom', fontsize=13, color='black', fontweight='bold')
                    plt.text(i, upper_whisker+(y_max-y_min) / 20.0, f'{p_value:0.2e}', ha='center', va='bottom', fontsize=8, color='black')

        # Adjust legend position
        plt.legend(loc='upper right')
    
        plt.tight_layout(pad=2.0)
    
        # Save the figure
        plt.savefig(f'{plot_dir}/GS_CF_{name}_batch_{batch_num + 1}.png', dpi=400)
        plt.close()
    stats_df = pd.DataFrame(stats_list)
    stats_df.to_excel(f"{plot_dir}/GS_CF_{name}_allBatch.xlsx", index=False)
    

plot_GeneSegment_save("Alpha")
plot_GeneSegment_save("Beta")
