#!/usr/bin/env bash
set -euo pipefail

RECORD_ID="19248700"
FILE="Datasets.tar.gz"
BASE_URL="https://zenodo.org/records/${RECORD_ID}/files"

echo "Record ID: ${RECORD_ID}"
echo "Downloading: ${FILE}"
wget -c "${BASE_URL}/${FILE}"

echo "Downloading: ${FILE}.sha256"
wget -c "${BASE_URL}/${FILE}.sha256"

echo "Verifying checksum..."
sha256sum -c "${FILE}.sha256"

echo "Extracting archive..."
tar -xzf "${FILE}"

echo "Removing downloaded archive: ${FILE}"
rm -f "${FILE}"
echo "Removing downloaded archive: ${FILE}.sha256"
rm -f "${FILE}.sha256"

echo "Done."
