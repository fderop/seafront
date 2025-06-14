import os
import anndata as ad
import cellxgene_census
import hashlib

def _file_checksum(path, hash_algo='sha256', blocksize=65536):
    hasher = hashlib.new(hash_algo)
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(blocksize), b''):
            hasher.update(block)
    return hasher.hexdigest()

def _update_checksums_file(dir_path, filename, checksum):
    checksum_file = os.path.join(dir_path, ".checksums")
    lines = []

    if os.path.exists(checksum_file):
        with open(checksum_file, "r") as f:
            lines = f.readlines()

    lines = [l for l in lines if not l.strip().endswith(f" {filename}")]
    lines.append(f"{checksum} {filename}\n")
    with open(checksum_file, "w") as f:
        f.writelines(lines)

def fetch_cellxgene(organism="Homo sapiens", dataset_id=None, raw_dir="raw_h5ad", checksum=None):
    """
    Downloads or loads an AnnData object from cellxgene_census. 
    Checks for an existing h5ad and validates checksum if provided.
    Automatically manages a .checksums file in the h5ad directory.
    """
    if not dataset_id:
        raise ValueError("dataset_id must be specified")
    os.makedirs(raw_dir, exist_ok=True)
    adata_fname = f"adata_{organism.replace(' ', '_').lower()}_{dataset_id}.h5ad"
    adata_path = os.path.join(raw_dir, adata_fname)

    if os.path.exists(adata_path):
        print(f"Found cached file: {adata_path}")
        if checksum:
            if _file_checksum(adata_path) != checksum:
                print("Checksum mismatch. Redownloading...")
                os.remove(adata_path)
                return fetch_and_save_data(organism, dataset_id, adata_path, raw_dir, adata_fname)
            print("Checksum OK. Loading AnnData...")
            return ad.read_h5ad(adata_path)
        print("Loading AnnData...")
        return ad.read_h5ad(adata_path)
    return fetch_and_save_data(organism, dataset_id, adata_path, raw_dir, adata_fname, checksum=checksum)

def fetch_and_save_data(organism, dataset_id, adata_path, raw_dir, adata_fname, checksum=None):
    print("Downloading AnnData from cellxgene_census...")
    with cellxgene_census.open_soma(census_version="2025-01-30") as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            obs_value_filter=f"dataset_id == '{dataset_id}'"
        )
    adata.write(adata_path)
    print(f"Saved to {adata_path}")

    sha = _file_checksum(adata_path)
    _update_checksums_file(raw_dir, adata_fname, sha)
    print(f"Checksum ({sha}) written to {os.path.join(raw_dir, '.checksums')}")

    if checksum:
        if sha != checksum:
            raise RuntimeError(f"Checksum failed after download. Got {sha}, expected {checksum}")
        print("Checksum OK.")

    return adata

