# Copilot Instructions Setup - Verification Guide

This document helps verify that the Copilot instructions are properly configured.

## Files Created/Modified

### Core Instruction Files

1. **`.github/copilot-instructions.md`** (existing, verified)
   - Repository-wide instructions for all files
   - Always loaded by Copilot
   - Contains project overview, key commands, data flow, debugging patterns

2. **`.github/AGENTS.md`** (renamed from AGENT.md)
   - Quick reference for AI agents
   - Follows agents.md specification (https://agents.md/)
   - Contains prescriptive steps for safe, minimal changes

### Path-Specific Instructions

Located in `.github/instructions/`:

1. **`snakemake.instructions.md`** (existing, verified)
   - Applies to: `workflow/**/*.smk`, `workflow/Snakefile`
   - Mandatory rule directives, templates, best practices

2. **`python-scripts.instructions.md`** (new)
   - Applies to: `workflow/scripts/**/*.py`
   - Script templates, Truvari patterns, error handling

3. **`r-analysis.instructions.md`** (new)
   - Applies to: `R/**/*.R`, `analysis/**/*.{R,qmd}`
   - Code style, caching system, data loading patterns

4. **`development.instructions.md`** (existing, verified)
   - Applies to: `**` (all files)
   - Environment setup, git workflow, tool documentation

5. **`README.md`** (new)
   - Documents the instruction file organization
   - Explains how path-specific instructions work

### Automation Workflows

1. **`.github/workflows/copilot-setup-steps.yml`** (new)
   - Automates environment setup for Copilot Coding Agent
   - Installs dependencies, validates workflow, runs linting
   - Follows GitHub's recommended best practices

## Verification Checklist

### Structure Verification

- [x] `.github/copilot-instructions.md` exists and is comprehensive
- [x] `.github/AGENTS.md` exists (not AGENT.md)
- [x] All `.instructions.md` files have proper YAML frontmatter with `applyTo`
- [x] Instructions directory has a README.md
- [x] copilot-setup-steps.yml workflow exists

### Content Verification

Check that each instruction file contains:

**snakemake.instructions.md:**
- [ ] Required directives (log, ensure(), threads, resources, conda)
- [ ] Download rules pattern
- [ ] Conda environment standards
- [ ] Complete rule template

**python-scripts.instructions.md:**
- [ ] Resource management (with statements)
- [ ] Error handling (raise exceptions)
- [ ] Truvari patterns with .name extraction
- [ ] Common patterns for CSV, Parquet, BED files

**r-analysis.instructions.md:**
- [ ] Code style (air formatter, lintr)
- [ ] Caching system documentation
- [ ] Factor levels and schema patterns
- [ ] Quarto rendering instructions

**development.instructions.md:**
- [ ] Environment setup commands
- [ ] Git workflow and commit conventions
- [ ] Tool documentation links
- [ ] Testing and validation commands

### Functional Verification

To test that instructions are working:

1. **Test Snakemake instruction application:**
   - Open a `.smk` file
   - Ask Copilot to create a new rule
   - Verify it includes all required directives

2. **Test Python instruction application:**
   - Open a Python script in `workflow/scripts/`
   - Ask Copilot to add error handling
   - Verify it uses raise instead of returning empty values

3. **Test R instruction application:**
   - Open an R file or Quarto document
   - Ask Copilot about caching
   - Verify it references the cache.R patterns

4. **Test setup workflow:**
   - Trigger the copilot-setup-steps workflow
   - Verify it completes successfully
   - Check that environment is properly configured

## Best Practices Followed

✅ **Repository-wide instructions** in `.github/copilot-instructions.md`
✅ **Path-specific instructions** in `.github/instructions/*.instructions.md`
✅ **AGENTS.md** for AI agent targeting (agents.md spec)
✅ **Setup workflow** for environment automation
✅ **YAML frontmatter** with `applyTo` patterns
✅ **Comprehensive documentation** with examples and templates
✅ **Clear organization** with README in instructions directory

## References

- [GitHub Copilot Custom Instructions](https://docs.github.com/en/copilot/how-tos/configure-custom-instructions/add-repository-instructions)
- [AGENTS.md Specification](https://agents.md/)
- [Copilot Setup Steps Blog](https://github.blog/changelog/2025-07-30-copilot-coding-agent-custom-setup-steps-are-more-reliable-and-easier-to-debug/)
- [Onboarding AI Peer Programmer](https://github.blog/ai-and-ml/github-copilot/onboarding-your-ai-peer-programmer-setting-up-github-copilot-coding-agent-for-success/)

## Maintenance

To keep instructions effective:

1. **Update when conventions change** - Keep instructions synchronized with actual practices
2. **Add examples** - Real examples are more helpful than abstract rules
3. **Test regularly** - Verify Copilot follows instructions correctly
4. **Remove outdated info** - Delete or update deprecated patterns
5. **Link to docs** - Reference authoritative documentation
6. **Be specific** - Concrete patterns work better than vague guidance

## Next Steps

1. Test the instructions by having Copilot work on various files
2. Gather feedback from team members on instruction effectiveness
3. Iterate and improve based on real usage
4. Consider adding more specialized instruction files if needed (e.g., for tests/)
5. Keep instructions in sync with codebase changes
