#!/usr/bin/env bash
###############################################################################
# Download practicals from the ELIXIR Spatial Omics 2026 course
#
# Usage:
#   bash download_practicals.sh              # interactive menu
#   bash download_practicals.sh all          # download all practicals
#   bash download_practicals.sh 0 3 7        # download practical_0, _3, _7
#   bash download_practicals.sh reset 2      # remove practical_2
#   bash download_practicals.sh reset all    # remove all practicals
###############################################################################

set -euo pipefail

# ---------- Config ----------
REPO="addityea/ELIXIR-SCO-spatial-omics-2026"
BRANCH="main"
PARENT="practicals"
N_PRACTICALS=10
DEST="work/"
# Define Zenodo DOI links to download and unzip in data/ directory for certain practicals
declare -A practical_dois=(
    [practical_0]="https://zenodo.org/records/17641420/files/Practical0.zip?download=1"
)

REPO_NAME="${REPO##*/}"
TARBALL_URL="https://github.com/${REPO}/archive/refs/heads/${BRANCH}.tar.gz"

# ---------- Colors ----------
if [ -t 1 ]; then
    BOLD="\033[1m"; GREEN="\033[32m"; BLUE="\033[34m"
    YELLOW="\033[33m"; RED="\033[31m"; RESET="\033[0m"
else
    BOLD=""; GREEN=""; BLUE=""; YELLOW=""; RED=""; RESET=""
fi

log()  { echo -e "${BLUE}ℹ️  $*${RESET}"; }
ok()   { echo -e "${GREEN}✅ $*${RESET}"; }
warn() { echo -e "${YELLOW}⚠️  $*${RESET}"; }
err()  { echo -e "${RED}❌ $*${RESET}"; }

# ---------- Functions ----------
download_one() {
    local n="$1"
    local name="practical_${n}"
    local target="${DEST}/${name}"

    if [ -d "$target" ]; then
        warn "${name} already exists — removing old copy first..."
        rm -rf "$target"
    fi

    mkdir -p "$DEST"
    log "Downloading ${BOLD}${name}${RESET}..."

    if curl -sL "$TARBALL_URL" \
        | tar -xz -C "$DEST" --strip-components=2 \
            "${REPO_NAME}-${BRANCH}/${PARENT}/${name}" 2>/dev/null; then
        local n_files n_nb
        n_files=$(find "$target" -type f | wc -l | tr -d ' ')
        n_nb=$(find "$target" -name "*.ipynb" | wc -l | tr -d ' ')
        ok "${name} → ${target} (${n_files} files, ${n_nb} notebooks)"
        if [[ -n "${practical_dois[$name]:-}" ]]; then
            log "Downloading data for ${name} from Zenodo DOI: ${practical_dois[$name]}"
            curl --cookie zenodo-cookies.txt -L "${practical_dois[$name]}" -o "work/${name}.zip"
            unzip -o "work/${name}.zip" -d "$target/data/"
            rm "work/${name}.zip"
            # Clean up any __MACOSX folders that may be included in the zip
            find "$target/data/" -name "__MACOSX" -type d -exec rm -rf {} + >/dev/null 2>&1 || true
            ok "Data for ${name} downloaded and extracted to ${target}/data/"
        fi
    else
        err "Failed to download ${name}. Does it exist in the repo?"
        return 1
    fi
}

download_all() {
    mkdir -p "$DEST"
    log "Downloading ${BOLD}ALL${RESET} practicals..."
    if curl -sL "$TARBALL_URL" \
        | tar -xz -C "$DEST" --strip-components=2 \
            "${REPO_NAME}-${BRANCH}/${PARENT}"; then
        ok "All practicals downloaded → ${DEST}"

    else
        err "Download failed."
        return 1
    fi
}

reset_one() {
    local n="$1"
    local target="${DEST}/practical_${n}"
    if [ -d "$target" ]; then
        rm -rf "$target"
        ok "Removed practical_${n}"
    else
        warn "practical_${n} not present."
    fi
}

reset_all() {
    local removed=0
    for i in $(seq 0 $((N_PRACTICALS - 1))); do
        local target="${DEST}/practical_${i}"
        if [ -d "$target" ]; then
            rm -rf "$target"
            removed=$((removed + 1))
        fi
    done
    ok "Removed ${removed} practicals."
}

show_status() {
    echo
    echo -e "${BOLD}Current status in ${DEST}:${RESET}"
    for i in $(seq 0 $((N_PRACTICALS - 1))); do
        local target="${DEST}/practical_${i}"
        if [ -d "$target" ]; then
            echo -e "  ${GREEN}✅ practical_${i}${RESET}"
        else
            echo -e "  ${YELLOW}⬜ practical_${i}${RESET}  (not downloaded)"
        fi
    done
    echo
}

interactive_menu() {
    while true; do
        show_status
        echo -e "${BOLD}Options:${RESET}"
        echo "  [0-9]   Download a specific practical (e.g. type '3')"
        echo "  a       Download ALL practicals"
        echo "  r N     Reset (delete) practical N (e.g. 'r 2')"
        echo "  R       Reset ALL practicals"
        echo "  q       Quit"
        echo
        read -rp "Choice: " choice

        case "$choice" in
            [0-9])         download_one "$choice" ;;
            a|A)           download_all ;;
            r\ [0-9])      reset_one "${choice#r }" ;;
            R)             reset_all ;;
            q|Q|exit|quit) ok "Bye!"; exit 0 ;;
            *)             warn "Unknown choice: $choice" ;;
        esac
    done
}

# ---------- Main ----------
echo -e "${BOLD}🧬 ELIXIR Spatial Omics 2026 — Practicals downloader${RESET}"
echo -e "Repo: ${BLUE}https://github.com/${REPO}/tree/${BRANCH}/${PARENT}${RESET}"
echo -e "Destination: ${BLUE}${DEST}${RESET}"
echo

# No args → interactive menu
if [ $# -eq 0 ]; then
    interactive_menu
    exit 0
fi

# Handle reset subcommand
if [ "$1" = "reset" ]; then
    shift
    if [ $# -eq 0 ] || [ "$1" = "all" ]; then
        reset_all
    else
        for n in "$@"; do
            reset_one "$n"
        done
    fi
    exit 0
fi

# Handle "all"
if [ "$1" = "all" ]; then
    download_all
    exit 0
fi

# Handle list of numbers
for n in "$@"; do
    if [[ "$n" =~ ^[0-9]+$ ]]; then
        download_one "$n"
    else
        warn "Skipping invalid argument: $n"
    fi
done

ok "Done! 👉 Refresh the file browser to see your files."
