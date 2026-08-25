#!/usr/bin/env bash
# Download every file of an ALMA EC (external correlation) delivery listing.
# Default: the EHT 2021 3C273 delivery (project 2019.1.01218.V), ~1.2 TB.
#
# Usage:
#   scripts/download_alma_ec_delivery.sh [--sizes] [URL] [DEST_DIR]
#
#   --sizes   Only print the per-file size manifest and grand total, then exit.
#
# Resumable and re-runnable: interrupted transfers continue where they left
# off (curl -C -) and files already complete (local size == remote size) are
# skipped. Set CHECK_MD5=1 to verify against the delivery's md5sum.txt at the
# end (reads every downloaded byte, so it takes a while on the full ~1.2 TB).
set -euo pipefail

SIZES_ONLY=0
if [ "${1:-}" = "--sizes" ]; then
    SIZES_ONLY=1
    shift
fi
URL="${1:-https://almascience.eso.org/almadata/ec/eht-2021/2019.1.01218.V/}"
DEST="${2:-$(basename "${URL%/}")}"

command -v curl >/dev/null || { echo "error: curl is required" >&2; exit 1; }

echo "listing $URL"
files=$(curl -fsS "$URL" | grep -oE 'href="[^"?/][^"]*"' | sed 's/href="//;s/"//' | sort -u)
[ -n "$files" ] || { echo "error: no files found in the listing" >&2; exit 1; }

manifest=$(mktemp)
trap 'rm -f "$manifest"' EXIT
total=0
count=0
for f in $files; do
    size=$(curl -fsSI --max-time 30 "$URL$f" | tr -d '\r' \
           | awk 'tolower($1)=="content-length:"{print $2}')
    size=${size:-0}
    printf '%s %s\n' "$size" "$f" >> "$manifest"
    total=$((total + size))
    count=$((count + 1))
done
awk -v t="$total" -v n="$count" \
    'BEGIN { printf "%d files, %.1f GB total (%.2f TB)\n", n, t/1e9, t/1e12 }'
if [ "$SIZES_ONLY" = 1 ]; then
    sort -rn "$manifest" | awk '{printf "%10.2f GB  %s\n", $1/1e9, $2}'
    exit 0
fi

mkdir -p "$DEST"
avail=$(df -k "$DEST" | awk 'NR==2 {print $4 * 1024}')
if [ "$avail" -lt "$total" ]; then
    awk -v a="$avail" -v t="$total" 'BEGIN {
        printf "WARNING: %.1f GB free in destination, %.1f GB needed\n", a/1e9, t/1e9 }'
fi

while read -r size f; do
    local_size=$(stat -f %z "$DEST/$f" 2>/dev/null \
                 || stat -c %s "$DEST/$f" 2>/dev/null || echo 0)
    if [ "$local_size" = "$size" ] && [ "$size" != 0 ]; then
        echo "complete   $f"
        continue
    fi
    awk -v s="$size" -v f="$f" 'BEGIN {printf "fetching   %s (%.1f GB)\n", f, s/1e9}'
    curl -fS --retry 5 --retry-delay 15 -C - -o "$DEST/$f" "$URL$f"
done < "$manifest"

md5file=$(find "$DEST" -maxdepth 1 -name "*md5sum.txt" | head -1)
if [ -n "$md5file" ] && [ "${CHECK_MD5:-0}" = 1 ]; then
    if command -v md5sum >/dev/null; then
        echo "verifying checksums (reads every byte)"
        (cd "$DEST" && md5sum -c "$(basename "$md5file")")
    else
        echo "note: md5sum not available (macOS: brew install coreutils); skipping verification"
    fi
fi
echo "delivery complete in $DEST"
