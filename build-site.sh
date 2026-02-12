#!/usr/bin/env bash
#
# build-site.sh -- Assemble and render the GitHub Pages site.
#
# Usage:
#   ./build-site.sh              # Build the site locally
#   ./build-site.sh --publish    # Build and publish to gh-pages branch
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
SITE_DIR="$REPO_ROOT/site"

echo "==> Syncing content into site/ ..."

# Analysis notebooks
mkdir -p "$SITE_DIR/analysis"
rsync -a --delete \
  --include='*.qmd' \
  --include='_notebook_setup.R' \
  --exclude='*' \
  "$REPO_ROOT/analysis/" "$SITE_DIR/analysis/"

# Docs (excluding agent_work)
mkdir -p "$SITE_DIR/docs/diagrams"
rsync -a --delete --exclude='agent_work' \
  "$REPO_ROOT/docs/" "$SITE_DIR/docs/"

# README
cp "$REPO_ROOT/README.md" "$SITE_DIR/README.md"

# Snakemake report (if present)
if ls "$REPO_ROOT"/report*.html 2>/dev/null 1>&2; then
  cp "$REPO_ROOT"/report*.html "$SITE_DIR/"
  echo "    Copied Snakemake report."
else
  echo "    No Snakemake report found (skipping)."
fi

# Create symlinks for data dependencies so notebooks can render
for dir in R config results resources manuscript; do
  target="$REPO_ROOT/$dir"
  link="$SITE_DIR/$dir"
  if [ -d "$target" ]; then
    ln -sfn "$target" "$link"
  fi
done

# Symlink styles.css and bibliography for notebooks that reference them
[ -f "$REPO_ROOT/styles.css" ] && ln -sfn "$REPO_ROOT/styles.css" "$SITE_DIR/styles.css"

# Sentinel file so here::here() resolves inside site/
touch "$SITE_DIR/.here"

echo "==> Rendering site ..."
cd "$SITE_DIR"
quarto render

echo "==> Site built at $SITE_DIR/_site/"

if [ "${1:-}" = "--publish" ]; then
  echo "==> Publishing to gh-pages branch ..."
  quarto publish gh-pages --no-browser
fi
