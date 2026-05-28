FROM ghcr.io/prefix-dev/pixi:latest AS build
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG PRAC
ENV PIXI_ENV_PATH=practicals/practical_${PRAC}/env

WORKDIR /app
# Install git and curl to fetch the environment files from the repo, then clean up apt cache to reduce image size
RUN apt-get update && apt-get install -y git curl && rm -rf /var/lib/apt/lists/*

COPY ${PIXI_ENV_PATH}/pixi.toml /app/
COPY ${PIXI_ENV_PATH}/pixi.lock /app/

# Install the selected environment (locked for reproducibility)
RUN pixi install --locked

# Create an activation script that sets up the env vars
RUN pixi shell-hook -s bash > /shell-hook.sh && \
    echo 'exec "$@"' >> /shell-hook.sh

FROM ubuntu:22.04 AS production
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG NB_USER=nbis
ARG NB_UID=1000
ARG NB_GID=1000

ENV NB_USER=${NB_USER} \
    NB_UID=${NB_UID} \
    HOME=/home/${NB_USER} \
    DEBIAN_FRONTEND=noninteractive \
    SHELL=/bin/bash \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    LANGUAGE=en_US:en

# Install minimal runtime deps
RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        tini \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user with UID 1000 (required for SciLifeLab Serve)
RUN userdel -r ubuntu 2>/dev/null || true && \
    groupadd --gid 1000 ${NB_USER} && \
    useradd --uid 1000 --gid 1000 --create-home --shell /bin/bash ${NB_USER}

# Copy pixi-installed environment + activation hook
COPY --from=build /app /app
COPY --from=build /shell-hook.sh /shell-hook.sh
ENV PATH="/app/.pixi/envs/default/bin:${PATH}"

# Start script in WORKDIR (Serve requirement)
RUN cat > /app/start-script.sh <<'EOF'
#!/bin/bash
set -e

HOME=/home/nbis
cd "$HOME"
source /shell-hook.sh && \
    exec jupyter lab \
        --no-browser \
        --ip=0.0.0.0 \
        --port=8888 \
        --PasswordIdentityProvider.hashed_password='argon2:$argon2id$v=19$m=10240,t=10,p=8$lELOJTD7tb0yBhYrLjRMEg$NKicFPxobrYUdyU7EdTVjeZhJsv85U6uzJJLHa8MfhE' \
        --IdentityProvider.token='' \
        --MappingKernelManager.kernel_info_timeout=300
EOF
RUN mkdir -p ${HOME}/.local/share/jupyter \
             ${HOME}/.jupyter \
             ${HOME}/.cache \
             ${HOME}/work
# Make sure user owns its home and app dirs
RUN chmod +x /app/start-script.sh /shell-hook.sh && \
    chown -R ${NB_UID}:${NB_GID} ${HOME} /app

# Copy helper script to download labs
COPY docker/download_labs.sh /home/${NB_USER}/download-labs.sh
RUN chmod +x /home/${NB_USER}/download-labs.sh && \
    chown ${NB_UID}:${NB_GID} /home/${NB_USER}/download-labs.sh

# Make /root/.pixi/bin accessible to the non-root user for pixi shell-hook to work (since pixi installs in /root by default) and it's in the PATH, this is a simple workaround
RUN mkdir -p /root/.pixi/bin && chmod -R 755 /root && chown -R ${NB_UID}:${NB_GID} /root

# Switch to non-root user
USER ${NB_UID}:${NB_GID}
RUN /shell-hook.sh python -m ipykernel install --user --name=ELIXIR_SP --display-name "ELIXIR Spatial Omics Kernel"
WORKDIR /app

EXPOSE 8888

# tini handles signals; shell-hook activates pixi env then runs the start script
ENTRYPOINT ["tini", "-g", "--", "./start-script.sh"]