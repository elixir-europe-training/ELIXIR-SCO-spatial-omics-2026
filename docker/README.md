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
# docker build --platform linux/amd64 --build-arg PRAC=6 -t <YOUR-DOCKERHUB-ID>/elixir-prac-6:<version> -f docker/rstudio/Dockerfile 
```

Replace `<YOUR-DOCKERHUB-ID>` with your actual Docker Hub username and `<version>` with the version tag you want to use (e.g., `v1.0`) and `X` with the practical number (e.g., `0` for Practical 0). This command will build the Docker image for the specified practical and tag it appropriately.

# Practical Data

If your practical requires large datasets to be downloaded, we recommend hosting them on `Zenodo` and letting us know the DOI and your practical number in the pull request description. We will consolidate all and add them to the download script. In case you need it beforehand, you can simply open an issue asking for this to be added and we will do it as soon as possible.

# GitHub Actions

GitHub Actions are set up to automatically build and push the Docker images to GitHub container registry (GHCR) whenever changes are made to the `docker/` directory or the `practicals/*/env/` directories. This ensures that the latest versions of the images are always available for use. The workflow is defined in the `.github/workflows/build_push.yml` file, which specifies the conditions under which the images should be built and pushed.
