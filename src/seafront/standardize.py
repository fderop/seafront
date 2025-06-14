import pandas as pd
import re

def summarize_obs(obs):
    exclude_cols = ["observation_joinid", "soma_joinid"]
    obs = obs[[c for c in obs.columns if c not in exclude_cols]]

    summary_dict = {}

    for col in obs.columns:
        if col == "dataset_id":
            continue
        elif col.startswith("n_") or col in {"raw_sum", "nnz", "raw_mean_nnz", "raw_variance_nnz"}:
            summary_dict[col] = ("median_" + col, "median")
        elif col in {"cell_type", "assay", "tissue", "disease", "self_reported_ethnicity",
                     "sex", "development_stage", "tissue_general", "suspension_type", "donor_id"}:
            summary_dict[col] = ("unique_" + col + "s_present", lambda x: ", ".join(sorted(set(map(str, x.dropna().unique())))))
        elif col.endswith("_ontology_term_id"):
            summary_dict[col] = ("unique_" + col + "s_present", lambda x: ", ".join(sorted(set(map(str, x.dropna().unique())))))
        else:
            summary_dict[col] = ("unique_" + col + "s_present", lambda x: ", ".join(sorted(set(map(str, x.dropna().unique())))))

    agg_dict = {k: v[1] for k, v in summary_dict.items()}
    rename_dict = {k: v[0] for k, v in summary_dict.items()}

    summary = obs.groupby("dataset_id", observed=False).agg(agg_dict).reset_index().rename(columns=rename_dict)
    return summary

def filter_age(stage_str, convert_map=None, drop_set=None):
    if convert_map is None:
        convert_map = {
            "child stage (1-4 yo)": 2,
            "juvenile stage (5-14 yo)": 10,
            "newborn stage (0-28 days)": 0,
            "fourth decade stage": 45,
            "fifth decade stage": 55,
            "sixth decade stage": 65,
            "seventh decade stage": 75,
            "eighth decade stage": 85,
            "ninth decade stage": 95,
            "third decade stage": 35,
        }
    if drop_set is None:
        drop_set = set([
            "adult stage", "late adult stage", "postnatal stage", "prime adult stage", "young adult stage",
            "pediatric stage", "infant stage", "middle aged stage", "organogenesis stage",
            "blastula stage", "embryonic stage", "unknown",
            "fifth LMP month stage", "fourth LMP month stage", "eighth LMP month stage", "ninth LMP month stage"
        ] + [f"Carnegie stage {str(i).zfill(2)}" for i in range(9, 24)])

    if stage_str in convert_map:
        return convert_map[stage_str]
    if stage_str in drop_set:
        return None
    if "LMP month stage" in stage_str:
        return None
    if stage_str.startswith("Carnegie stage"):
        return None
    match = re.match(r"(\d+)-year-old stage", stage_str)
    if match:
        return int(match.group(1))
    match = re.match(r"(\d+)-month-old stage", stage_str)
    if match:
        return 0
    if "week post-fertilization stage" in stage_str:
        return None
    return None

def filter_obs_with_age_int(obs):
    obs = obs.copy()
    obs["age_int"] = obs["development_stage"].map(lambda stage: filter_age(stage)).astype('Int64')
    obs = obs[obs["age_int"].notnull()].copy()
    return obs