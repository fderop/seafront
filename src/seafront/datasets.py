import os
import pandas as pd
import scanpy as sc
import anndata as ad
import tarfile, gzip, shutil
import warnings
from anndata import _core

# Suppress specific warnings about var_names not being unique
warnings.filterwarnings(
    "ignore",
    message="Variable names are not unique.*",
    category=UserWarning,
    module=_core.anndata.__name__
)


def download(name: str):
    if name == "ainciburu2023":
        return _download_ainciburu2023()
    else:
        raise ValueError(f"Unknown dataset: {name}")

def _download_ainciburu2023():
    data_dir = "raw/GSE180298"
    output_path = os.path.join(data_dir, "GSE180298_combined.h5ad")

    tar_file_path = os.path.join(data_dir, "GSE180298_RAW.tar")
    if os.path.exists(tar_file_path):
        with tarfile.open(tar_file_path, "r") as tar:
            tar.extractall(path=data_dir)

    gz_filenames = [
        "GSE180298_elderly_metadata.txt.gz",
        "GSE180298_mds_metadata.txt.gz",
        "GSE180298_young_metadata.txt.gz",
    ]

    for gz in gz_filenames:
        gz_path = os.path.join(data_dir, gz)
        out_path = gz_path[:-3]
        if not os.path.exists(out_path):
            with gzip.open(gz_path, 'rb') as f_in:
                with open(out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    metadata_files = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith("_metadata.txt")
    ])
    metadata_list = [
        pd.read_csv(os.path.join(data_dir, f), sep="\t", index_col=0)
        for f in metadata_files
    ]
    metadata = pd.concat(metadata_list, axis=0)
    print(f"Loaded metadata: {metadata.shape[0]} rows")

    h5_files = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith(".h5")
    ])
    adata_list = []
    for h5_fname in h5_files:
        full_path = os.path.join(data_dir, h5_fname)
        sample_id = h5_fname.replace("_filtered_feature_bc_matrix.h5", "")
        print(f"Reading {h5_fname}...")

        adata = sc.read_10x_h5(full_path)
        adata.var_names_make_unique()
        assert adata.var_names.is_unique

        adata.obs_names = [f"{bc}_{sample_id}" for bc in adata.obs_names]
        adata.obs["sample_id"] = sample_id
        adata_list.append(adata)

    adata_combined = ad.concat(adata_list, axis=0, label="batch", index_unique=None)

    patient_age_map = {
        "mds1": 71,
        "mds2": 54,
        "mds3": 83,
        "mds4": 83,
        "young1": 20,
        "young2": 23,
        "young3": 20,
        "young4": 19,
        "young5": 23,
        "elderly1": 61,
        "elderly2": 74,
        "elderly3": 72,
    }

    def extract_patient_id(index_str):
        return index_str.split("_")[-1].lower()

    adata_combined.obs["patient_id"] = adata_combined.obs.index.map(extract_patient_id)
    adata_combined.obs["patient_age"] = adata_combined.obs["patient_id"].map(patient_age_map)

    # Create a cleaned index for adata_combined: barcode_patientid
    adata_barcodes = adata_combined.obs_names.to_series()
    adata_clean_index = adata_barcodes.str.extract(r"^(?P<barcode>[^-]+)-1_(?P<patient>.+)$")
    adata_combined.obs["metadata_key"] = adata_clean_index["barcode"] + "_" + adata_clean_index["patient"]
    metadata_aligned = metadata.loc[adata_combined.obs["metadata_key"]].copy()
    adata_combined.obs["CellType"] = metadata_aligned["CellType"].values


    print(f"Combined AnnData shape: {adata_combined.shape}")
    adata_combined.write(output_path)
    print(f"Saved to {output_path}")
    return metadata, adata_combined
