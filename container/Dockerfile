# OpenFOAM-14 + OpenMPI + PARPACK + ofLSA solvers, ready to run.
#
# Built automatically by .github/workflows/container.yml and pushed to GHCR.
# On a cluster (no root) run it with Apptainer:
#   apptainer pull oflsa.sif docker://ghcr.io/yurifeld/oflsa:latest
#   apptainer exec oflsa.sif oflsa-run blockMesh
#   mpirun -np 8 apptainer exec oflsa.sif oflsa-run StabilitySolver -parallel

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates wget gnupg \
 && wget -qO /etc/apt/trusted.gpg.d/openfoam.asc https://dl.openfoam.org/gpg.key \
 && echo "deb http://dl.openfoam.org/ubuntu noble main" > /etc/apt/sources.list.d/openfoam.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
        openfoam14 \
        libarpack2-dev libparpack2-dev gfortran \
 && rm -rf /var/lib/apt/lists/*

COPY . /opt/ofLSA

# Build into a fixed location inside the image (not $HOME, which the cluster mounts over)
RUN bash -c '\
    . /opt/openfoam14/etc/bashrc || true; \
    export FOAM_USER_LIBBIN=/opt/ofLSA/platforms/lib FOAM_USER_APPBIN=/opt/ofLSA/platforms/bin; \
    mkdir -p $FOAM_USER_LIBBIN $FOAM_USER_APPBIN; \
    . /opt/ofLSA/etc/bashrc; \
    cd /opt/ofLSA && ./Allwmake -j' \
 && install -m 755 /opt/ofLSA/container/oflsa-run /usr/local/bin/oflsa-run

ENTRYPOINT ["/usr/local/bin/oflsa-run"]
CMD ["bash"]
