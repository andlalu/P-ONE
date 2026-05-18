#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "usage: bash scripts/sync_run_to_s3.sh <local-run-root> <s3-prefix>" >&2
  exit 2
fi

LOCAL_RUN_ROOT="$1"
S3_PREFIX="$2"

aws s3 sync "$LOCAL_RUN_ROOT" "$S3_PREFIX"
