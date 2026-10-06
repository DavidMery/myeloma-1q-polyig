#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")"

Rscript scripts/master-pipeline/run-synthetic-pipeline.R