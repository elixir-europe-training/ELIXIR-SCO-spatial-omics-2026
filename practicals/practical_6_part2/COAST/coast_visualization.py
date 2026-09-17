# coast_visualization.py
import glob
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def run_umap_clustering(outdir, n_clusters=9, spot_size=20, background=True):
    vit_feature_dir = os.path.join(outdir, "ViT_features")
    vis_outdir      = os.path.join(outdir, "visualization")
    os.makedirs(vis_outdir, exist_ok=True)

    # --- Load files ---
    vit_mf_files    = sorted(glob.glob(os.path.join(vit_feature_dir, "*_vit_mf.csv")))
    vit_coord_files = sorted(glob.glob(os.path.join(vit_feature_dir, "*_vit_coord.csv")))
    if not vit_mf_files or not vit_coord_files:
        raise FileNotFoundError("No _vit_mf.csv or _vit_coord.csv files found")

    feature_dfs, coord_dfs, tissue_labels = [], [], []
    for mf_file, coord_file in zip(vit_mf_files, vit_coord_files):
        df_feat  = pd.read_csv(mf_file,    index_col=0)
        df_coord = pd.read_csv(coord_file, index_col=0)
        if not background:
            mask     = ~df_feat.index.str.endswith("_grid")
            df_feat  = df_feat[mask]
            df_coord = df_coord[mask]
        feature_dfs.append(df_feat)
        coord_dfs.append(df_coord)
        tissue_labels.extend([os.path.splitext(os.path.basename(mf_file))[0]] * len(df_feat))

    df_features        = pd.concat(feature_dfs, axis=0)
    df_coords          = pd.concat(coord_dfs,   axis=0)
    df_coords["tissue"] = tissue_labels
    print(f"[INFO] Loaded {len(df_features)} patches from {len(vit_mf_files)} sections")

    # --- Scale + UMAP + KMeans ---
    features_scaled = StandardScaler().fit_transform(df_features.values)
    embedding       = umap.UMAP(random_state=42).fit_transform(features_scaled)
    df_features["UMAP1"] = embedding[:, 0]
    df_features["UMAP2"] = embedding[:, 1]
    df_features["cluster"] = KMeans(n_clusters=n_clusters, random_state=42).fit_predict(embedding)
    df_coords["cluster"]   = df_features["cluster"].values

    out_csv = os.path.join(vit_feature_dir, "umap_clusters_concatenated_sections.csv")
    df_features.to_csv(out_csv)
    print(f"[INFO] Saved: {out_csv}")

    # --- Plot 1: UMAP coloured by cluster ---
    fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
    sns.scatterplot(x="UMAP1", y="UMAP2", hue="cluster",
                    palette="tab10", data=df_features,
                    s=spot_size, linewidth=0, ax=ax, legend=False)
    from matplotlib.patches import Patch
    palette  = sns.color_palette("tab10", n_colors=n_clusters)
    handles  = [Patch(color=palette[c], label=str(c)) for c in range(n_clusters)]
    ax.legend(handles=handles, fontsize=6, markerscale=0.6, frameon=False,
              bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    out_umap = os.path.join(vis_outdir, "umap_clusters.png")
    plt.savefig(out_umap, dpi=200, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved: {out_umap}")

    # --- Plot 2: UMAP coloured by tissue ---
    section_colors = ["#F3CA40", "#5c8001"]
    _, ax = plt.subplots(figsize=(6, 5), dpi=200)
    for i, tissue in enumerate(df_coords["tissue"].unique()):
        mask = df_coords["tissue"] == tissue
        ax.scatter(df_features.loc[mask, "UMAP1"],
                   df_features.loc[mask, "UMAP2"],
                   c=section_colors[i % len(section_colors)],
                   s=spot_size, alpha=0.7, label=tissue)
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=6, markerscale=0.6, frameon=False, loc="upper right")
    plt.tight_layout()
    out_umap_section = os.path.join(vis_outdir, "umap_by_section.png")
    plt.savefig(out_umap_section, dpi=200, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved: {out_umap_section}")

    # --- Plot 3: spatial side-by-side ---
    unique_tissues = df_coords["tissue"].unique()
    _, axes = plt.subplots(1, len(unique_tissues),
                             figsize=(6 * len(unique_tissues), 6), dpi=200)
    if len(unique_tissues) == 1:
        axes = [axes]

    for ax, tissue in zip(axes, unique_tissues):
        df_sub = df_coords[df_coords["tissue"] == tissue].copy()
        df_sub["y_plot"] = -df_sub["y"]
        sns.scatterplot(x="x", y="y_plot", hue="cluster",
                        palette="tab10", data=df_sub,
                        s=spot_size, linewidth=0, ax=ax, legend=False)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")
        ax.set_title(f"{tissue} - spatial clusters", fontsize=9)

    # single legend on the last axis
    handles = [Patch(color=palette[c], label=str(c)) for c in range(n_clusters)]
    axes[-1].legend(handles=handles, fontsize=6, markerscale=0.6, frameon=False,
                    bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)

    plt.tight_layout()
    out_spatial = os.path.join(vis_outdir, "spatial_clusters.png")
    plt.savefig(out_spatial, dpi=200, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved: {out_spatial}")

    return df_features



def plot_final_aligned_scatter(outdir, spot_size, background=True):
    """
    Plot final aligned scatterplot of all *_aligned.csv files in CAST_output.

    Args:
        outdir (str): Base output directory (contains CAST_output/).
        spot_size (int): Size of scatterplot dots.
        background (bool): If True, include all patches (intissue + grid).
                           If False, exclude patches ending in '_grid'.
    """
    cast_outdir = os.path.join(outdir, "CAST_output")
    vis_outdir = os.path.join(outdir, "visualization")
    os.makedirs(vis_outdir, exist_ok=True)

    # Find all *_aligned.csv files
    aligned_files = sorted(glob.glob(os.path.join(cast_outdir, "*_aligned.csv")))
    if len(aligned_files) < 2:
        raise FileNotFoundError(f"Need at least 2 *_aligned.csv files in {cast_outdir}")

    plt.figure(figsize=(8, 8), dpi=200)

    # Pick a palette (yellow, green, blue, red, etc.)
    colors = ["#F3CA40", "#5c8001", "#1d3557", "#e63946", "#457b9d", "#ff7f0e"]

    for i, f in enumerate(aligned_files):
        df = pd.read_csv(f, index_col=0)

        if not background:
            df = df[~df.index.str.endswith("_grid")]

        if not {"x_aligned", "y_aligned"}.issubset(df.columns):
            print(f"[WARNING] Skipping {f} (missing x_aligned / y_aligned)")
            continue

        name = os.path.basename(f).replace("_aligned.csv", "")
        plt.scatter(
            df["x_aligned"], -df["y_aligned"],
            c=colors[i % len(colors)], s=spot_size, alpha=0.7,
            edgecolors="k", linewidth=0.01,
            label=name
        )

    plt.gca().set_aspect("equal", adjustable="box")
    plt.axis("off")
    plt.legend(markerscale=1.5, fontsize=8, loc="upper right")
    plt.tight_layout()

    suffix = "" if background else "_nobg"
    out_final = os.path.join(vis_outdir, f"aligned_tissues{suffix}.png")
    plt.savefig(out_final, dpi=200)
    plt.show()
    print(f"[INFO] Final aligned scatterplot saved in: {out_final}")