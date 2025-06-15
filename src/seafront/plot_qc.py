import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd 


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

def plot_assay_counts_and_medians(obs_filtered, var="assay_simple", fig=None, axes=None, return_df=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    assay_counts = obs_filtered[var].value_counts()
    assay_medians = obs_filtered.groupby(var, observed=False)["raw_sum"].median()

    assay_df = pd.DataFrame({
        var: assay_counts.index,
        "n_cells": assay_counts.values,
        "median_raw_sum": assay_counts.index.map(assay_medians)
    }).set_index(var)
    assay_df["median_raw_sum"] = pd.to_numeric(assay_df["median_raw_sum"], errors="coerce")
    assay_df = assay_df.sort_index()

    if fig is None or axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(15, 0.5 * len(assay_df) + 2), sharey=True)

    order = assay_df.index

    sns.barplot(y=assay_df.index, x=assay_df["n_cells"], ax=axes[0], order=order)
    axes[0].set_ylabel("Assay")
    axes[0].set_xlabel("Number of cells")
    axes[0].set_title("Cell Count")
    axes[0].set_xscale('log')
    axes[0].tick_params(axis='y', labelrotation=0)

    sns.barplot(y=assay_df.index, x=assay_df["median_raw_sum"], ax=axes[1], order=order)
    axes[1].set_ylabel("")
    axes[1].set_xlabel("Median raw_sum")
    axes[1].set_title("Median raw_sum")
    axes[1].set_xscale('log')
    axes[1].tick_params(axis='y', labelleft=False)

    # Annotate values on bars
    for i, (count, y) in enumerate(zip(assay_df["n_cells"], range(len(assay_df)))):
        axes[0].text(count, y, f"{int(count):,}", va='center', ha='left', fontsize=9)

    for i, (median, y) in enumerate(zip(assay_df["median_raw_sum"], range(len(assay_df)))):
        if pd.notnull(median):
            if median > 0:
                display_value = f"{int(median):,}"
            else:
                display_value = "0"
            axes[1].text(median, y, display_value, va='center', ha='left', fontsize=9)

    plt.tight_layout()

    if return_df:
        return assay_df.reset_index()


def plot_ranked_median_counts_per_experiment(obs, fig=None, ax=None, return_df=False, cutoff=None):
    df = obs.copy()
    df["experiment"] = df["dataset_id"].astype(str) + "__" + df["assay_simple"].astype(str)
    median_counts = df.groupby("experiment", observed=True)["raw_sum"].median()
    median_counts = median_counts.sort_values(ascending=False).reset_index()

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    sns.lineplot(x=range(len(median_counts)), y=median_counts["raw_sum"], ax=ax)
    ax.set_xlabel("Experiment Rank")
    ax.set_ylabel("Median raw_sum per Experiment")
    ax.set_title("Rank Plot: Median Raw Counts per Dataset-Assay Combination")
    ax.set_yscale('log')

    if cutoff is not None:
        ax.axhline(y=cutoff, color="red", linestyle="--", label=f"cutoff={cutoff}")
        ax.legend()

    plt.tight_layout()
    if return_df:
        return median_counts
