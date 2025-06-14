import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_age_histograms(obs_filtered, summary_filtered, fig=None, axes=None):
    bins = np.arange(0, 110, 10)
    if fig is None or axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sns.histplot(obs_filtered["age_int"].dropna(), bins=bins, ax=axes[0], discrete=False)
    axes[0].set_xlabel('Age (years, by decade)')
    axes[0].set_ylabel('Number of cells')
    axes[0].set_title('Cell Count per Age Bin (obs_filtered)')
    axes[0].set_xticks(bins)

    def parse_int_list(s):
        if isinstance(s, list):
            return s
        if isinstance(s, str):
            s = s.strip("[]")
            return [int(x) for x in s.split(",") if x.strip().isdigit()]
        return []

    age_list = summary_filtered["unique_age_ints_present"].apply(parse_int_list)
    all_ages = [item for sublist in age_list for item in sublist]

    sns.histplot(all_ages, bins=bins, ax=axes[1], discrete=False)
    axes[1].set_xlabel('Age (years, by decade)')
    axes[1].set_ylabel('Number of datasets')
    axes[1].set_title('Dataset Count per Age Bin (summary_filtered)')
    axes[1].set_xticks(bins)

    plt.tight_layout()
    return fig, axes
