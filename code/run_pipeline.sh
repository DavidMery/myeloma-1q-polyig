#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")/.."

Rscript code/run_pipeline.R