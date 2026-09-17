#!/usr/bin/env python3
import os
import torch
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import CAST


def read_tissues(vit_dir):
    """
    Read ViT features and coordinates for all tissues in vit_dir.

    Returns
    -------
    tissues : list
        List of tissue names
    adata_list : list
        List of AnnData objects with scaled features and coordinates
    """
    tissues = []
    adata_list = []
    scaler = StandardScaler()

    mf_files = sorted([f for f in os.listdir(vit_dir) if f.endswith("_vit_mf.csv")])
    if not mf_files:
        raise FileNotFoundError(f"No *_vit_mf.csv files found in {vit_dir}")

    for mf_file in mf_files:
        tissue_name = mf_file.replace("_vit_mf.csv", "")
        feature_file = os.path.join(vit_dir, mf_file)
        coord_file = os.path.join(vit_dir, f"{tissue_name}_vit_coord.csv")

        if not os.path.exists(coord_file):
            raise FileNotFoundError(f"Missing coordinate CSV for tissue {tissue_name}: {coord_file}")

        # Load and scale features
        features = pd.read_csv(feature_file, index_col=0)
        features = pd.DataFrame(scaler.fit_transform(features), index=features.index, columns=features.columns)

        # Load coordinates
        coords = pd.read_csv(coord_file, index_col=0)

        # Create AnnData
        adata = sc.AnnData(X=features.values, obs=coords)
        adata.obs['sample'] = tissue_name
        adata_list.append(adata)
        tissues.append(tissue_name)

        print(f"[INFO] Loaded tissue '{tissue_name}' with {adata.shape[0]} tiles")

    return tissues, adata_list


# -------------------------------
# Run CAST alignment
# -------------------------------
def run_cast(outdir):
    """
    Perform CAST alignment for all tissues.
    Assumes ViT features are in outdir/ViT_features.
    Saves aligned coordinates in outdir/CAST_output
    """
    vit_dir = os.path.join(outdir, "ViT_features")
    if not os.path.exists(vit_dir):
        raise FileNotFoundError(f"ViT features directory not found: {vit_dir}")

    tissues, adata_list = read_tissues(vit_dir)
    CAST_DIR = os.path.join(outdir, "CAST_output")
    os.makedirs(CAST_DIR, exist_ok=True)

    print("[INFO] Concatenating all tissues for CAST...")
    adata_all = adata_list[0].concatenate(*adata_list[1:])
    print(f"[INFO] Total tiles: {adata_all.shape[0]}")

    sample_list = np.unique(adata_all.obs['sample'])
    coords_raw, exps = CAST.extract_coords_exp(
        adata_all, batch_key='sample',
        cols=['x', 'y'],
        count_layer='.X',
        data_format='norm1e4'
    )

    coords_sub, exp_sub, sub_node_idxs = CAST.sub_data_extract(sample_list, coords_raw, exps, nodenum_t=20000)

    # Save preprocessed data
    torch.save(coords_raw, os.path.join(CAST_DIR, "coords_raw.pt"))
    torch.save(sub_node_idxs, os.path.join(CAST_DIR, "sub_node_idxs.pt"))
    torch.save(exp_sub, os.path.join(CAST_DIR, "exp_sub.pt"))
    torch.save(coords_sub, os.path.join(CAST_DIR, "coords_sub.pt"))

    # Reload embeddings for each tissue
    scaler = StandardScaler()
    emb_dict = {}
    for tissue in tissues:
        emb_file = os.path.join(vit_dir, f"{tissue}_vit_mf.csv")
        emb = pd.read_csv(emb_file, index_col=0)
        emb = scaler.fit_transform(emb)
        emb_dict[tissue] = emb[sub_node_idxs[tissue]]

    embed_dict = {t: torch.tensor(emb_dict[t]).float() for t in tissues}

    print("[INFO] Running CAST_STACK...")
    params_dist = CAST.reg_params(
        dataname=tissues[0],  # reference tissue
        diff_step=5,
        iterations=500,
        dist_penalty1=0,
        bleeding=500,
        d_list=[3, 2, 1, 0.5, 0.333],
        attention_params=[None, 3, 1, 0],
        dist_penalty2=[1],
        alpha_basis_bs=[500],
        meshsize=[8],
        iterations_bs=[400],
        attention_params_bs=[[None, 3, 1, 0]],
        mesh_weight=[None]
    )
    params_dist.alpha_basis = torch.Tensor([1/1000, 1/1000, 1/50, 5, 5]).reshape(5, 1).to(params_dist.device)

    coord_final = CAST.CAST_STACK(coords_raw, embed_dict, CAST_DIR, sample_list,
                                  params_dist, sub_node_idxs=sub_node_idxs, rescale=True)

    # Save aligned coordinates for each tissue
    for tissue in tissues:
        df_aligned = pd.DataFrame(coord_final[tissue], columns=["x_aligned", "y_aligned"])

        # Load original coordinate file to reuse its index
        coord_file = os.path.join(vit_dir, f"{tissue}_vit_coord.csv")
        coords = pd.read_csv(coord_file, index_col=0)
        df_aligned.index = coords.index  # keep same patch index names

        # Save with proper index
        df_aligned.to_csv(os.path.join(CAST_DIR, f"{tissue}_aligned.csv"))

    print(f"[INFO] CAST alignment completed. Aligned coordinates saved to {CAST_DIR}")


# -------------------------------
# CLI
# -------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="CAST alignment for multiple tissues using ViT features")
    parser.add_argument("--outdir", required=True, help="Base output directory for COAST_output")
    args = parser.parse_args()
    run_cast(args.outdir)
