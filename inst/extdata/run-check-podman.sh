#!/usr/bin/env bash
# Run `devtools::check()` locally inside the mosuite-minimal container image
# via podman (e.g. on macOS with Podman Desktop / podman machine).
# Note: you may need to disable the VPN in order to pull the image.
#
# Usage:
#   ./inst/extdata/run-check-podman.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

IMAGE="docker.io/nciccbr/mosuite-minimal:latest"

if ! command -v podman >/dev/null 2>&1; then
  echo "Error: podman is not installed or not on \$PATH. Install it from https://podman.io/" >&2
  exit 1
fi

if ! podman info >/dev/null 2>&1; then
  echo "Podman machine is not running. Starting it..." >&2
  podman machine set --memory 8192 && podman machine start
fi

echo "Using image: $IMAGE"
podman pull "$IMAGE"

podman run --rm \
  -v "$REPO_ROOT":/workspace:Z \
  -w /workspace \
  "$IMAGE" \
  R -e 'devtools::check()'
