#!/bin/bash
# One-step setup on Ubuntu (needs sudo): installs what is missing — OpenFOAM-14 from
# openfoam.org, OpenMPI, gfortran, ARPACK/PARPACK — and builds the ofLSA solvers.
# Already installed packages are left alone; the script can be re-run safely.
#
#   scripts/install_ubuntu.sh            # install what is missing, then build
#   scripts/install_ubuntu.sh --check    # only report what is installed / missing

set -u

CHECK_ONLY=no
[ "${1:-}" = "--check" ] && CHECK_ONLY=yes

OFLSA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FOAM_BASHRC=/opt/openfoam14/etc/bashrc
PACKAGES="openfoam14 openmpi-bin libopenmpi-dev gfortran libarpack2-dev libparpack2-dev"

die() { echo "Error: $*" >&2; exit 1; }

# ---- system -----------------------------------------------------------------
[ -r /etc/os-release ] || die "cannot identify the OS (no /etc/os-release)"
. /etc/os-release
[ "${ID:-}" = ubuntu ] || die "this script is for Ubuntu (found: ${PRETTY_NAME:-unknown}).
On other systems see README section 'Server without sudo'."
CODENAME="${VERSION_CODENAME:-${UBUNTU_CODENAME:-}}"
[ "$(dpkg --print-architecture)" = amd64 ] || echo "Warning: openfoam.org packages are built for amd64."
echo "System: $PRETTY_NAME ($CODENAME)"

# ---- what is missing ---------------------------------------------------------
installed() { dpkg-query -W -f='${Status}' "$1" 2>/dev/null | grep -q "install ok installed"; }

MISSING=""
for p in $PACKAGES; do
    if installed "$p"; then echo "  [ok]      $p"
    else                    echo "  [missing] $p"; MISSING="$MISSING $p"
    fi
done

if [ "$CHECK_ONLY" = yes ]; then
    [ -z "$MISSING" ] && echo "Everything is installed." || echo "Missing:$MISSING"
    exit 0
fi

# ---- install ------------------------------------------------------------------
if [ -n "$MISSING" ]; then
    command -v sudo >/dev/null || die "sudo not available. See README section 'Server without sudo'."
    echo "Installing:$MISSING  (sudo will ask for your password)"

    case " $MISSING " in
    *" openfoam14 "*)
        if ! grep -rqs "dl.openfoam.org/ubuntu" /etc/apt/sources.list /etc/apt/sources.list.d/; then
            curl -fsI "http://dl.openfoam.org/ubuntu/dists/$CODENAME/Release" >/dev/null 2>&1 \
                || die "openfoam.org has no packages for Ubuntu '$CODENAME'. Supported: see https://openfoam.org/download/14-ubuntu/"
            echo "Adding the openfoam.org package repository"
            sudo apt-get update -qq
            sudo apt-get install -y -qq ca-certificates curl >/dev/null || die "apt-get failed"
            curl -fsSL https://dl.openfoam.org/gpg.key | sudo tee /etc/apt/trusted.gpg.d/openfoam.asc >/dev/null \
                || die "could not download the openfoam.org key"
            echo "deb http://dl.openfoam.org/ubuntu $CODENAME main" \
                | sudo tee /etc/apt/sources.list.d/openfoam.list >/dev/null
        fi
        ;;
    esac

    sudo apt-get update -qq || die "apt-get update failed"
    # shellcheck disable=SC2086
    sudo apt-get install -y $MISSING || die "apt-get install failed"
fi

[ -f "$FOAM_BASHRC" ] || die "$FOAM_BASHRC not found after installation"

# ---- build --------------------------------------------------------------------
# Clean environment: sourcing OpenFOAM twice in one shell breaks it, and on some
# desktops its ParaView check hangs (same workaround as scripts/oflsa-shell).
echo "Building ofLSA ..."
env -i HOME="$HOME" USER="${USER:-}" PATH=/usr/local/bin:/usr/bin:/bin \
    bash -c "
        source '$FOAM_BASHRC' ParaView_TYPE=none < /dev/null > /dev/null 2>&1
        cd '$OFLSA_DIR' && ./Allwmake -j
    " > "$OFLSA_DIR/log.Allwmake" 2>&1
STATUS=$?
tail -4 "$OFLSA_DIR/log.Allwmake"
[ $STATUS -eq 0 ] || die "build failed, see $OFLSA_DIR/log.Allwmake"

cat <<EOF

ofLSA is installed. To run the example:

    $OFLSA_DIR/scripts/oflsa-shell
    cd tutorials/lidDrivenCavity3D
    ./Allrun -np 4
EOF
