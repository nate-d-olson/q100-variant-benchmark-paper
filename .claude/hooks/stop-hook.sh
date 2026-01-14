#!/usr/bin/env bash
#
# Stop hook for Claude Code - runs after each AI response
# This project-local hook overrides the global stop-hook.sh
#
# Purpose: Check if codebase improvements task is complete,
#          and only request if not done

set -euo pipefail

# Check if the improvements documentation exists
if [ -f "CODEBASE_IMPROVEMENTS.md" ]; then
    # Check if the implementation branch exists
    if git rev-parse --verify feature/implement-improvements >/dev/null 2>&1; then
        # Task is complete - no feedback needed
        exit 0
    fi
fi

# Task not complete - request it
echo "Identify 5 features or whys this codebase can be improved. Create a markdown file with the suggested features then create a new branch to implement the features."
exit 1
