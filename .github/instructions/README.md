# Copilot Instructions Directory

This directory contains path-specific instructions for GitHub Copilot Coding Agent.

## File Organization

- **snakemake.instructions.md** - Rules and templates for Snakemake workflow development
  - Applies to: `workflow/**/*.smk`, `workflow/Snakefile`
  - Contains: Mandatory rule directives, conda environment standards, rule templates

- **python-scripts.instructions.md** - Guidelines for Python scripts in Snakemake rules
  - Applies to: `workflow/scripts/**/*.py`
  - Contains: Script templates, resource management, Truvari patterns, error handling

- **r-analysis.instructions.md** - R and Quarto analysis code guidelines
  - Applies to: `R/**/*.R`, `analysis/**/*.{R,qmd}`
  - Contains: Code style, caching system, data loading, factor levels, Quarto patterns

- **development.instructions.md** - General development workflow and tooling
  - Applies to: `**` (all files)
  - Contains: Environment setup, git workflow, tool documentation, commit conventions

## How These Files Work

Each `.instructions.md` file uses YAML frontmatter to specify which files it applies to:

```yaml
---
applyTo: "path/pattern/**/*.ext"
---
```

When GitHub Copilot Coding Agent works on a file, it will automatically load and follow the relevant instructions based on the path patterns.

## Related Files

- **../.github/copilot-instructions.md** - Repository-wide instructions (always loaded)
- **../.github/AGENTS.md** - Quick reference guide for AI agents (always loaded)

## Best Practices

1. Keep instructions focused and actionable
2. Use specific examples and templates
3. Document project-specific conventions and patterns
4. Update instructions when conventions change
5. Test instructions by having Copilot work on relevant files

## Adding New Instructions

To add new path-specific instructions:

1. Create a new file: `<topic>.instructions.md`
2. Add YAML frontmatter with `applyTo` pattern
3. Write clear, specific guidance for that code area
4. Update this README to document the new file

## References

- [GitHub Copilot Custom Instructions Documentation](https://docs.github.com/en/copilot/how-tos/configure-custom-instructions/add-repository-instructions)
- [AGENTS.md Specification](https://agents.md/)
