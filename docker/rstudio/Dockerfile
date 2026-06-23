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
RUN pixi install --locked --run-post-link-scripts


ARG NB_UID=1000
ARG NB_GID=1000

RUN chown -R ${NB_UID}:${NB_GID} /app

# Create an activation script that sets up the env vars
RUN pixi shell-hook -s bash > /shell-hook.sh && \
    echo 'exec "$@"' >> /shell-hook.sh

FROM ubuntu:latest AS production
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG NB_USER=rstudio
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
        curl \
        libglx-mesa0 \
        psmisc \
        gdebi-core \
        libssl-dev\
        libclang-dev \
        lsb-release \
        sudo \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user with UID 1000 (required for SciLifeLab Serve)
RUN userdel -r ubuntu 2>/dev/null || true && \
    groupadd --gid 1000 ${NB_USER} && \
    useradd --system --uid 1000 --gid 1000 --create-home --shell /bin/bash ${NB_USER}

# Copy pixi-installed environment + activation hook
COPY --from=build /app /app
COPY --from=build /shell-hook.sh /shell-hook.sh
ENV PATH="/app/.pixi/envs/default/bin:${PATH}"

RUN ARCH="$(dpkg --print-architecture)" && \
    echo "Detected architecture: $ARCH" && \
    if [ "$ARCH" = "amd64" ]; then \
        RSTUDIO_SERVER_DEB="https://s3.amazonaws.com/rstudio-ide-build/server/jammy/amd64/rstudio-server-2026.05.0-daily-128-amd64.deb"; \
    elif [ "$ARCH" = "arm64" ]; then \
        RSTUDIO_SERVER_DEB="https://s3.amazonaws.com/rstudio-ide-build/server/jammy/arm64/rstudio-server-2026.05.0-daily-128-arm64.deb"; \
    else \
        echo "Unsupported architecture: $ARCH"; \
        exit 1; \
    fi && \
    wget $RSTUDIO_SERVER_DEB -O /tmp/rstudio-server.deb && \
    gdebi -n /tmp/rstudio-server.deb && \
    rm /tmp/rstudio-server.deb
# Set a password for the RStudio user created above.
RUN echo "${NB_USER}:rstudio" | chpasswd && \
    usermod -aG sudo ${NB_USER}

# Set the R working directory for the RStudio
RUN echo 'session-default-working-dir=/work/' >> /etc/rstudio/rsession.conf && \
    echo 'session-default-new-project-dir=/work/' >> /etc/rstudio/rsession.conf

RUN mkdir -p /home/rstudio/.rserver/db /home/rstudio/.rserver/datadir && \
    chown -R ${NB_UID}:${NB_GID} /home/rstudio

RUN cat > /etc/rstudio/rserver.conf <<'EOF'
server-user=rstudio
auth-none=1
server-data-dir=/home/rstudio/.rserver/datadir
auth-timeout-minutes=0
auth-stay-signed-in-days=30
EOF

RUN cat > /etc/rstudio/database.conf <<'EOF'
directory=/home/rstudio/.rserver/db
auto-create=1
EOF

# Set the download method to wget
ENV R_LIBS_USER="/${NB_USER}/Rlibs"
# Apptainer compatibility layer

# Fix permissions for the R libraries directory
RUN chown -R ${NB_USER}:${NB_GID} /etc/rstudio && \
    mkdir -p ${R_LIBS_USER} && \
    chown -R ${NB_USER}:${NB_GID} ${R_LIBS_USER} && \
    chown -R ${NB_USER}:${NB_GID} /var/lib/rstudio-server
EXPOSE 8787
# Make /root/.pixi/bin accessible to the non-root user for pixi shell-hook to work (since pixi installs in /root by default) and it's in the PATH, this is a simple workaround
RUN mkdir -p /root/.pixi/bin && chmod -R 755 /root && chown -R ${NB_UID}:${NB_GID} /root

# Start script in WORKDIR (Serve requirement)
RUN cat > /app/start-script.sh <<'EOF'
#!/bin/bash
set -e

export USER=rstudio
export HOME=/home/rstudio
export USERNAME=rstudio

cd "$HOME"
source /shell-hook.sh
exec rstudio-server start
EOF

RUN mkdir -p ${HOME}/.cache \
             ${HOME}/work

# Make sure user owns its home and app dirs
RUN chmod +x /app/start-script.sh /shell-hook.sh && \
    chown -R ${NB_UID}:${NB_GID} ${HOME}

# Copy helper script to download labs
COPY docker/download_labs.sh /home/${NB_USER}/download-labs.sh
RUN chmod +x /home/${NB_USER}/download-labs.sh && \
    chown ${NB_UID}:${NB_GID} /home/${NB_USER}/download-labs.sh

# Switch to non-root user
USER ${NB_UID}:${NB_GID}
WORKDIR /app


# tini handles signals; shell-hook activates pixi env then runs the start script
ENTRYPOINT ["tini", "-g", "--", "./start-script.sh"]
