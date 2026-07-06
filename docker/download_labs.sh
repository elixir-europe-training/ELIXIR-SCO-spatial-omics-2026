#!/usr/bin/env bash
###############################################################################
# Download practicals from the ELIXIR Spatial Omics 2026 course
#
# Usage:
#   bash download_labs.sh [OPTIONS] [COMMAND]
#
# Options (all optional, shown with defaults):
#   -r, --repo          OWNER/NAME  GitHub repository  (elixir-europe-training/ELIXIR-SCO-spatial-omics-2026)
#   -b, --branch        BRANCH      Git branch          (main)
#   -p, --parent        DIR         Folder in repo      (practicals)
#   -n, --n-practicals  N           Number of practicals (10)
#   -d, --dest          PATH        Local destination   (work/)
#
# Commands:
#   (none)              Interactive menu
#   all                 Download all practicals
#   0 3 7               Download practical_0, _3, _7
#   reset 2             Remove practical_2
#   reset all           Remove all practicals
###############################################################################

set -euo pipefail

# ---------- Config ----------
REPO="elixir-europe-training/ELIXIR-SCO-spatial-omics-2026"
BRANCH="main"
PARENT="practicals"
N_PRACTICALS=10
DEST="work/"
# Zenodo archives to download and unzip into the central <DEST>/data/ directory.
# Each dataset is stored under its ORIGINAL archive name (e.g. data/Practical0),
# and reuse is keyed on that name: point several practicals at the SAME URL and
# the (multi-GB) archive is downloaded/extracted only once and shared.
declare -A practical_dois=(
    [practical_0]="https://zenodo.org/records/17641420/files/Practical0.zip?download=1"
)

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
# The repo archive is large (~230 MB) because data/images are committed in git.
# Download it ONCE per run and extract every requested practical from the cached
# copy, instead of re-streaming the whole tarball for each practical.
TARBALL_FILE=""
cleanup() { [ -n "$TARBALL_FILE" ] && rm -f "$TARBALL_FILE"; }
trap cleanup EXIT

fetch_tarball() {
    [ -n "$TARBALL_FILE" ] && [ -f "$TARBALL_FILE" ] && return 0
    TARBALL_FILE="$(mktemp "${TMPDIR:-/tmp}/elixir-sco.XXXXXX")"
    log "Fetching repository archive..."
    if ! curl -fL --progress-bar "$TARBALL_URL" -o "$TARBALL_FILE"; then
        err "Failed to download archive from ${TARBALL_URL}"
        rm -f "$TARBALL_FILE"; TARBALL_FILE=""
        return 1
    fi
}

download_data() {
    # Download + extract Zenodo data for a practical, if one is registered.
    #
    # Datasets live in a single central ${DEST}/data dir and keep their ORIGINAL
    # archive name (e.g. work/data/Practical0), not a per-practical name. The
    # existence check is keyed on that dataset name, so if several practicals
    # point at the same Zenodo archive the (multi-GB) download happens only once
    # and the rest reuse the already-extracted folder.
    local name="$1"
    [[ -n "${practical_dois[$name]:-}" ]] || return 0
    local url="${practical_dois[$name]}"

    # Dataset name = the archive's own filename stem (Practical0.zip -> Practical0).
    local base="${url##*/}"; base="${base%%\?*}"   # strip path + query string
    local data_name="${base%.*}"                    # strip .zip
    local data_dir="${DEST}/data"
    local dest="${data_dir}/${data_name}"

    if [ -d "$dest" ]; then
        ok "Dataset '${data_name}' already present in ${data_dir}/ — reusing it (no download for ${name})."
        return 0
    fi

    mkdir -p "$data_dir"
    local zip stage
    zip="$(mktemp "${TMPDIR:-/tmp}/${data_name}.XXXXXX")"
    stage="$(mktemp -d "${data_dir}/.stage.XXXXXX")"   # same FS as dest -> atomic move
    log "Downloading dataset ${BOLD}${data_name}${RESET} for ${name} from Zenodo: ${url}"
    if ! curl -fL --progress-bar "$url" -o "$zip"; then
        err "Failed to download dataset '${data_name}'."
        rm -f "$zip"; rm -rf "$stage"; return 1
    fi
    if ! unzip -oq "$zip" -d "$stage"; then
        err "Failed to unzip dataset '${data_name}'."
        rm -f "$zip"; rm -rf "$stage"; return 1
    fi
    rm -f "$zip"
    # Clean up any __MACOSX folders that may be included in the zip
    find "$stage" -name "__MACOSX" -type d -exec rm -rf {} + 2>/dev/null || true

    # Publish under the original name. Only stage -> dest at the end, so an
    # interrupted download never leaves a half-populated dataset folder behind.
    local tops only
    tops=$(find "$stage" -mindepth 1 -maxdepth 1 | wc -l | tr -d ' ')
    only=$(find "$stage" -mindepth 1 -maxdepth 1)
    if [ "$tops" -eq 1 ] && [ -d "$only" ]; then
        mv "$only" "$dest"          # zip already had a single top-level folder
        rm -rf "$stage"
    else
        mv "$stage" "$dest"         # loose files: keep them under data_name/
    fi
    ok "Dataset '${data_name}' extracted to ${dest}/"
}

download_one() {
    local n="$1"
    local name="practical_${n}"
    local target="${DEST}/${name}"

    fetch_tarball || return 1

    if [ -d "$target" ]; then
        warn "${name} already exists — removing old copy first..."
        rm -rf "$target"
    fi

    mkdir -p "$DEST"
    log "Extracting ${BOLD}${name}${RESET}..."

    if tar -xz -f "$TARBALL_FILE" -C "$DEST" --strip-components=2 \
            "${REPO_NAME}-${BRANCH}/${PARENT}/${name}"; then
        local n_files n_nb
        n_files=$(find "$target" -type f | wc -l | tr -d ' ')
        n_nb=$(find "$target" -name "*.ipynb" | wc -l | tr -d ' ')
        ok "${name} → ${target} (${n_files} files, ${n_nb} notebooks)"
        download_data "$name"
    else
        err "Failed to extract ${name}. Does it exist in the repo?"
        return 1
    fi
}

download_all() {
    fetch_tarball || return 1
    mkdir -p "$DEST"
    log "Extracting ${BOLD}ALL${RESET} practicals..."
    if tar -xz -f "$TARBALL_FILE" -C "$DEST" --strip-components=2 \
            "${REPO_NAME}-${BRANCH}/${PARENT}"; then
        ok "All practicals downloaded → ${DEST}"
    else
        err "Extraction failed."
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
            a|A)           download_all ;;
            r\ *)          reset_one "${choice#r }" ;;
            R)             reset_all ;;
            q|Q|exit|quit) ok "Bye!"; exit 0 ;;
            '')            : ;;
            *)
                if [[ "$choice" =~ ^[0-9]+$ ]]; then
                    download_one "$choice"
                else
                    warn "Unknown choice: $choice"
                fi
                ;;
        esac
    done
}

# ---------- Parse options ----------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--repo)           REPO="$2";          shift 2 ;;
        -b|--branch)         BRANCH="$2";        shift 2 ;;
        -p|--parent)         PARENT="$2";        shift 2 ;;
        -n|--n-practicals)   N_PRACTICALS="$2";  shift 2 ;;
        -d|--dest)           DEST="$2";          shift 2 ;;
        --) shift; break ;;
        -*) err "Unknown option: $1"; exit 1 ;;
        *) break ;;
    esac
done

DEST="${DEST%/}"   # strip trailing slash to avoid work//practical_0
REPO_NAME="${REPO##*/}"
TARBALL_URL="https://github.com/${REPO}/archive/refs/heads/${BRANCH}.tar.gz"

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
