# Post-course change log

Below is a list of changes that have been made to this branch after the course:

## 2026-11-06

### Added

- Notebook includes an option to subset the data for practical 1 to reduce Cellpose run time

## 2026-11-03

### Fixed

- All images are now visible in html-files for practical exercises

## 2026-11-02

### Fixed

- Xenium boundary cropping in p3 has been fixed to align with sainsc crop
- An option to recompute steps in p9 has been added

## 2026-10-30

### Added 

- This changelog has been added

### Fixed

- A smaller version of the Xenium dataset has been uploaded to <https://zenodo.org/records/23766246>.
    - This link replaces the original link in `docker/practical_data.tsv` 
- Practical 0 now creates a symlink to avoid downloading Practical1_cropped.zarr if the cropped dataset has already been created
- `docker/download_labs.sh` can now download from branches with `/` in the branch name

## 2026-10-29

### Added

- Functionality in `docker/download_labs.sh` to delete data from individual practicals when resetting (`reset --data`)

### Fixed

- README.md has been updated to reflect course status
- The data file for practial 6 part 1 dataset is no longer downloaded
    - The paths in the practical 6 part 1 notebook now link to the identical files in the part 2 download
    - The link for part 1 in `docker/practical_data.tsv` is now the same as part 2

