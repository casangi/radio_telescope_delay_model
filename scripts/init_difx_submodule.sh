#!/usr/bin/env bash
# Initialize the pinned DiFX submodule (extern/difx) the way this package
# expects it: shallow (single commit), blob-filtered and sparse, so only
# applications/difxcalc11 is materialized (~60 MB) instead of the full
# ~200 MB DiFX worktree.
#
# A plain `git submodule update --init extern/difx` also works if you do not
# mind downloading the whole tree; this script exists so no machine has to.
# Idempotent: safe to re-run, and after a deliberate SHA bump it moves the
# checkout to the newly recorded commit.
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "$0")/.." && pwd)
path=extern/difx
url=$(git -C "$repo_root" config -f .gitmodules "submodule.${path}.url")
sha=$(git -C "$repo_root" ls-files -s -- "$path" | awk '{print $2}')
if [ -z "$sha" ]; then
    echo "error: no gitlink for ${path} in the index" >&2
    exit 1
fi

dir="$repo_root/$path"
mkdir -p "$dir"
# Keep Dropbox away from the submodule's git internals (macOS Dropbox honors
# this xattr; a no-op elsewhere). Machines syncing the repo through Dropbox
# get the submodule by running this script, not through sync.
xattr -w com.dropbox.ignored 1 "$dir" 2>/dev/null || true

if [ ! -e "$dir/.git" ]; then
    git init -q "$dir"
    git -C "$dir" remote add origin "$url"
fi
git -C "$dir" config remote.origin.promisor true
git -C "$dir" config remote.origin.partialclonefilter blob:none
git -C "$dir" sparse-checkout set --cone applications/difxcalc11
git -C "$dir" fetch --depth 1 --filter=blob:none origin "$sha"
git -C "$dir" checkout -q --detach "$sha"
git -C "$repo_root" submodule init -- "$path" >/dev/null 2>&1 || true

echo "extern/difx ready at $(git -C "$dir" log --oneline -1)"
