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
```

Replace `<YOUR-DOCKERHUB-ID>` with your actual Docker Hub username and `<version>` with the version tag you want to use (e.g., `v1.0`) and `X` with the practical number (e.g., `0` for Practical 0). This command will build the Docker image for the specified practical and tag it appropriately.

# GitHub Actions

**TO DO**: We will set up a GitHub Actions workflow in the `.github/workflows/` directory to automate the building and pushing of the Docker images to Docker Hub whenever changes are made to the `practicals/practical_X/env/` directories or the `docker/jupyter/Dockerfile`. The workflow will be triggered on push events to the main branch and will build the Docker images for all practicals, tagging them with the appropriate version and pushing them to Docker Hub under your account.
