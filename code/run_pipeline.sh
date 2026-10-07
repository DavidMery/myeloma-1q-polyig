#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAPSULE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$CAPSULE_ROOT"

exec Rscript "$SCRIPT_DIR/run_pipeline.R"