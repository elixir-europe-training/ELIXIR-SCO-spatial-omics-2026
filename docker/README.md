# Instructions

## How to create Pixi environments

Each practical has its own Pixi environment defined in the `practicals/practical_X/env/` directory (replace `X` with the practical number). This right now contains a basic environment that enables `Jupyterlab` or `Rstudio` to run, but you can customize it by adding any additional packages you need for the practical. Feel free to modify the `pixi.toml` file to remove conflicting packages but highlight them during the pull request review so we can discuss them.

To create the environment for a specific practical, navigate to that directory and run:

```bash
cd practicals/practical_X/env/
pixi add --no-install <package-name-a> <package-name-b> ...
```
Replace `<package-name-a>`, `<package-name-b>`, etc. with the actual package names you want to include in the environment. After adding the necessary packages

## How to create your containers

Go to the main repo directory (not this `docker` directory) and run the following command to build the Jupyter container for a specific practical (replace `X` with the practical number):

```bash
docker build --platform linux/amd64 --build-arg PRAC=X -t <YOUR-DOCKERHUB-ID>/elixir-prac-0:<version> -f docker/jupyter/Dockerfile .

# For practical 7 : 
# docker build --platform linux/amd64 --build-arg PRAC=7 -t <YOUR-DOCKERHUB-ID>/elixir-prac-7:<version> -f docker/rstudio/Dockerfile .
```

Replace `<YOUR-DOCKERHUB-ID>` with your actual Docker Hub username and `<version>` with the version tag you want to use (e.g., `v1.0`) and `X` with the practical number (e.g., `0` for Practical 0). This command will build the Docker image for the specified practical and tag it appropriately.

# Running containers

The final containers will be run on SciLifeLab Serve. To test on Serve, create a Custom App linking to your image.

The password for the jupyter labs is `spatial`.

# Practical Data

If your practical requires large datasets to be downloaded, we recommend hosting them on `Zenodo` and letting us know the DOI and your practical number in the pull request description. We will consolidate all and add them to the download script. In case you need it beforehand, you can simply open an issue asking for this to be added and we will do it as soon as possible.

# Download Script

`download_labs.sh` lets participants fetch practicals from GitHub into the container's `work/` directory. It can be run interactively or non-interactively via command-line arguments.

## Usage

```
bash download_labs.sh [OPTIONS] [COMMAND]
```

### Options

All options are optional and fall back to the defaults shown below.

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| `-r` | `--repo` | `elixir-europe-training/ELIXIR-SCO-spatial-omics-2026` | GitHub repository (`OWNER/NAME`) |
| `-b` | `--branch` | `main` | Git branch to download from |
| `-p` | `--parent` | `practicals` | Folder inside the repo that contains the practicals |
| `-n` | `--n-practicals` | `10` | Total number of practicals (used by the status view and reset-all) |
| `-d` | `--dest` | `work/` | Local destination directory |

Options must be placed **before** the command:

```bash
bash download_labs.sh --dest ~/mywork --branch dev all
```

### Commands

| Command | Description |
|---------|-------------|
| *(none)* | Open the interactive menu |
| `all` | Download all practicals |
| `0 3 7` | Download specific practicals by number |
| `reset 2` | Delete practical 2 |
| `reset all` | Delete all practicals |

### Examples

```bash
# Interactive menu (default)
bash download_labs.sh

# Download everything to a custom folder
bash download_labs.sh --dest ~/practicals all

# Download practicals 0 and 3 from a feature branch
bash download_labs.sh -b dev 0 3

# Reset practical 5
bash download_labs.sh reset 5
```

# GitHub Actions

GitHub Actions are set up to automatically build and push the Docker images to GitHub container registry (GHCR) whenever changes are made to the `docker/` directory or the `practicals/*/env/` directories. This ensures that the latest versions of the images are always available for use. The workflow is defined in the `.github/workflows/build_push.yml` file, which specifies the conditions under which the images should be built and pushed.
