import os
import pandas as pd
import torch
from PIL import Image
from torchvision import transforms
from tqdm import tqdm
import scanpy as sc
import numpy as np

# -------------------------------
# Feature extraction from images
# -------------------------------
def extract_features(coord_path, img_path, output_name, outdir,
                     patch_size=224, step_size=224, pad="none", mode="moldata"):
    """
    Extract ViT features from tissue image patches.
    Supports:
      - moldata: tissue-guided patches + grid patches
      - nomoldata: only grid patches
    """
    VIT_DIR = os.path.join(outdir, "ViT_features")
    os.makedirs(VIT_DIR, exist_ok=True)

    # Load coordinates (supports .csv, .h5ad, .h5)
    if coord_path.endswith(".csv"):
        coords_df = pd.read_csv(coord_path, index_col=0)
    elif coord_path.endswith(".h5ad"):
        ad = sc.read_h5ad(coord_path)
        coords_df = ad.obs.copy()
    elif coord_path.endswith(".h5"):
        ad = sc.read_10x_h5(coord_path)
        coords_df = ad.obs.copy()
    else:
        raise ValueError(f"Unsupported coordinate file type: {coord_path}")

    # Ensure x and y exist
    if "x" not in coords_df.columns or "y" not in coords_df.columns:
        raise ValueError(f"{coord_path} missing 'x' and 'y' columns")

    # Safely cast coords: keep floats if there are NaNs, otherwise cast to int
    if coords_df[["x", "y"]].isna().any().any():
        coords_df[["x", "y"]] = coords_df[["x", "y"]].astype(float)
    else:
        coords_df[["x", "y"]] = coords_df[["x", "y"]].astype(int)

    print(f"[INFO] Loaded {len(coords_df)} coordinates from {coord_path}")

    # Load image
    img = Image.open(img_path).convert("RGB")
    print(f"[INFO] Loaded image {img_path} with size {img.size}")

    # Apply padding
    if pad != "none":
        new_w = img.width + patch_size if pad in ["left", "right"] else img.width
        new_h = img.height + patch_size if pad in ["top", "bottom"] else img.height
        new_img = Image.new("RGB", (new_w, new_h), (255, 255, 255))
        paste_x = patch_size if pad == "left" else 0
        paste_y = patch_size if pad == "top" else 0
        new_img.paste(img, (paste_x, paste_y))
        img = new_img
        print(f"[INFO] Applied {pad} padding, new image size {img.size}")

    # Initialize ViT model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = torch.hub.load('facebookresearch/dino:main', 'dino_vits16')
    model.eval().to(device)
    print(f"[INFO] ViT model loaded on {device}")

    transform = transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize([0.5]*3, [0.5]*3)
    ])

    # Prepare extraction
    img_width, img_height = img.size
    half_patch = patch_size // 2
    features, patch_labels, coords_list = [], [], []
    patches = []

    # -------------------------------
    # Patch extraction logic
    # -------------------------------
    if mode == "moldata":
        # Tissue-guided patches
        print("[INFO] Extracting tissue-guided patches...")
        for i, (x_center, y_center) in enumerate(zip(coords_df['x'], coords_df['y'])):
            x_start = max(0, x_center - half_patch)
            y_start = max(0, y_center - half_patch)
            x_end = min(img_width, x_start + patch_size)
            y_end = min(img_height, y_start + patch_size)
            if x_end - x_start >= patch_size and y_end - y_start >= patch_size:
                patch = img.crop((x_start, y_start, x_end, y_end))
                patches.append(patch)
                coords_list.append((x_center, y_center))
                patch_labels.append(f"{coords_df.index[i]}_intissue")
        print(f"[INFO] Total tissue patches: {len(patches)}")

        # Grid patches
        print("[INFO] Extracting grid patches...")
        for y in range(0, img_height - patch_size + 1, step_size):
            for x in range(0, img_width - patch_size + 1, step_size):
                patch = img.crop((x, y, x + patch_size, y + patch_size))
                patches.append(patch)
                coords_list.append((x + half_patch, y + half_patch))
                patch_labels.append(f"{len(patches)-1}_grid")
        print(f"[INFO] Total grid patches: {len(patches)}")

    elif mode == "nomoldata":
        # Grid patches only
        print("[INFO] Extracting grid patches...")
        for y in range(0, img_height - patch_size + 1, step_size):
            for x in range(0, img_width - patch_size + 1, step_size):
                patch = img.crop((x, y, x + patch_size, y + patch_size))
                patches.append(patch)
                coords_list.append((x + half_patch, y + half_patch))
                patch_labels.append(f"{len(patches)-1}_grid")
        print(f"[INFO] Total grid patches: {len(patches)}")

    # -------------------------------
    # Feature extraction
    # -------------------------------
    print(f"[INFO] Extracting features from {len(patches)} patches...")
    for patch in tqdm(patches, desc="Patches"):
        patch_tensor = transform(patch).unsqueeze(0).to(device)
        with torch.no_grad():
            feat = model(patch_tensor)
            if isinstance(feat, (tuple, list)):
                feat = feat[0]
            features.append(feat.flatten().cpu().numpy())
    
    # Stack features into a 2D array
    features_array = pd.DataFrame(
        data=np.vstack(features),           
        index=patch_labels,                
        columns=[f"feat_{i}" for i in range(len(features[0]))]
    )

    coords_array = pd.DataFrame(coords_list, columns=['x','y'], index=patch_labels)

    # Save CSVs with index
    features_array.to_csv(os.path.join(VIT_DIR, f"{output_name}_vit_mf.csv"), index=True)
    coords_array.to_csv(os.path.join(VIT_DIR, f"{output_name}_vit_coord.csv"), index=True)

    print(f"[INFO] Saved ViT features to {VIT_DIR}/{output_name}_vit_mf.csv")
    print(f"[INFO] Saved patch coordinates to {VIT_DIR}/{output_name}_vit_coord.csv")

    return features_array, coords_array

