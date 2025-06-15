import requests
import os
import pandas as pd
import cellxgene_census
import os
import random
import cellxgene_census



def get_cellxgene_paper_reference(experiment_id):
    url = f"https://api.cellxgene.cziscience.com/dp/v1/datasets/{experiment_id}"
    resp = requests.get(url)
    if not resp.ok:
        raise Exception(f"Failed to fetch metadata for {experiment_id}")

    data = resp.json()
    publication = None

    if "publication_doi" in data and data["publication_doi"]:
        doi = data["publication_doi"]
        title = data.get("publication_title", "")
        publication = f"DOI: {doi} {title}".strip()
    elif "links" in data and data["links"]:
        for link in data["links"]:
            if link.get("rel") == "publication":
                publication = link.get("href")
                break
    elif "source" in data:
        publication = data["source"]

    return publication or "No reference/publication found."


def load_census_obs():
    obs_file_path = "meta/census_obs.parquet"

    if os.path.exists(obs_file_path):
        obs = pd.read_parquet(obs_file_path)
    else:
        with cellxgene_census.open_soma(census_version="2025-01-30") as census:
            obs = census["census_data"]["homo_sapiens"].obs.read().concat().to_pandas()
        obs.to_parquet(obs_file_path, index=False)

    return obs

def filter_experiments_by_median_raw_sum(obs_filtered, threshold=3000):
    grouped = obs_filtered.groupby("experiment", observed=True)["raw_sum"].median()
    valid_experiments = grouped[grouped >= threshold].index
    obs_filtered = obs_filtered[obs_filtered["experiment"].isin(valid_experiments)].copy()
    summary_filtered = seafront.standardize.summarize_obs(obs_filtered)
    return obs_filtered, summary_filtered


def get_var_names_for_dataset(dataset_id):
    with cellxgene_census.open_soma(census_version="2025-01-30") as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_value_filter=f"dataset_id == '{dataset_id}'",
            column_names={"var": ["feature_name"]},
        )
        return adata.var["feature_name"].tolist()

def check_var_consistency_and_save(obs):
    dataset_ids = obs["dataset_id"].dropna().unique().tolist()
    if len(dataset_ids) < 3:
        raise ValueError("Need at least 3 datasets to sample from.")

    sampled_ids = random.sample(dataset_ids, 3)
    print(f"Checking gene consistency across: {sampled_ids}")

    reference_var = get_var_names_for_dataset(sampled_ids[0])
    all_match = True

    for ds_id in sampled_ids[1:]:
        current_var = get_var_names_for_dataset(ds_id)
        if current_var != reference_var:
            print(f"Dataset {ds_id} has different gene names than reference.")
            all_match = False
        else:
            print(f"Dataset {ds_id} matches reference.")

    if all_match:
        out_path = "meta/census_var.txt"
        with open(out_path, "w") as f:
            for gene in reference_var:
                f.write(f"{gene}\n")
        print(f"All datasets matched. Gene names written to: {out_path}")
    else:
        print("Gene name inconsistency detected. No file written.")
