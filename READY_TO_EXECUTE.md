# Table Migration: Ready to Execute

## ✅ Setup Complete

All migration infrastructure has been prepared and committed to branch `feat/migrate-tables-to-flextable`.

**Commit:** `66be7c0` - "feat: add flextable theme and migration infrastructure"

---

## 📋 What Was Created

### Core Functions
- ✅ `theme_flextable_manuscript()` - Added to [R/plot_themes.R](R/plot_themes.R)
- ✅ `fmt_integer_flextable()` - Helper function
- ✅ `fmt_number_flextable()` - Helper function  
- ✅ `cols_label_flextable()` - Helper function

### Test Suite
- ✅ [tests/table_compatibility_test.qmd](tests/table_compatibility_test.qmd) - Compare 4 packages
- ✅ [tests/flextable_theme_implementation.qmd](tests/flextable_theme_implementation.qmd) - Theme demo
- ✅ [tests/run_table_tests.R](tests/run_table_tests.R) - Automated test script

### Documentation
- ✅ [tests/TABLE_RECOMMENDATIONS_SUMMARY.md](tests/TABLE_RECOMMENDATIONS_SUMMARY.md) - Executive summary
- ✅ [tests/gt_to_flextable_quick_reference.md](tests/gt_to_flextable_quick_reference.md) - Code patterns
- ✅ [tests/table_packages_research.md](tests/table_packages_research.md) - Detailed research
- ✅ [tests/TABLE_TESTING_README.md](tests/TABLE_TESTING_README.md) - Testing guide
- ✅ [MIGRATION_EXECUTION_GUIDE.md](MIGRATION_EXECUTION_GUIDE.md) - Complete execution guide

### Migration Scripts
- ✅ [setup_flextable_migration.sh](setup_flextable_migration.sh) - Install packages, render tests
- ✅ [pilot_migration.sh](pilot_migration.sh) - Convert one notebook
- ✅ [full_migration.sh](full_migration.sh) - Convert all notebooks

---

## 🚀 Next Steps (Execute Outside Sandbox)

The VS Code sandbox prevents R execution. You'll need to run these commands in a native terminal (Terminal.app or iTerm).

### Step 1: Run Setup (15-30 min)

Open **Terminal.app** (not VS Code terminal):

```bash
cd /Users/nolson/projects/q100-variant-benchmark-paper
bash setup_flextable_migration.sh
```

This will:
1. Install `flextable` and `officer` packages
2. Render test notebooks to HTML and Word
3. Open outputs for your review

**Expected outputs:**
- `tests/table_compatibility_test.html`
- `tests/table_compatibility_test.docx`
- `tests/flextable_theme_implementation.html`
- `tests/flextable_theme_implementation.docx`

### Step 2: Test Google Docs (5-10 min)

1. Upload one of the `.docx` files to Google Drive
2. Right-click → Open with → Google Docs
3. Verify tables look good and are editable

**What to check:**
- ✅ Tables are editable (not images)
- ✅ Borders are visible
- ✅ Headers have gray background
- ✅ Numbers have commas and correct decimals
- ✅ Overall formatting looks professional

### Step 3: Run Pilot Migration (30-60 min)

If test outputs look good:

```bash
bash pilot_migration.sh
```

Then manually convert one table in [analysis/benchmark_unique_regions.qmd](analysis/benchmark_unique_regions.qmd) using the patterns in [tests/gt_to_flextable_quick_reference.md](tests/gt_to_flextable_quick_reference.md).

### Step 4: Run Full Migration (2-3 hours)

If pilot succeeds:

```bash
bash full_migration.sh
```

Follow the manual review instructions in the output.

### Step 5: Upload to Google Drive

**For authentication and upload:**

I can help set up Google Drive integration, but first confirm:

1. Do you want to use the Google Drive API directly?
2. Or would you prefer to manually upload files?

**Option A: Manual Upload (Recommended for now)**
- Just drag and drop the `.docx` files to Google Drive
- Right-click → Open with → Google Docs

**Option B: Automated Upload (Future enhancement)**
- Would require setting up Google Drive API credentials
- Can create a script using `googledrive` R package
- Let me know if you want this after manual verification

---

## 📊 Migration Summary

### Notebooks to Convert
1. `analysis/benchmarkset_characterization.qmd`
2. `analysis/benchmark_difficult.qmd`
3. `analysis/benchmark_exclusions.qmd`
4. `analysis/benchmark_interval_size_distributions.qmd`
5. `analysis/benchmark_unique_regions.qmd`
6. `analysis/genomic_context_analysis.qmd`
7. `analysis/external_evaluation.qmd`

**Total:** 7 notebooks

### Time Estimate
- Setup: 30 min
- Testing: 15 min
- Pilot: 60 min
- Full migration: 2-3 hours
- **Total: 4-5 hours**

---

## 🎯 Success Criteria

✅ **Phase 1: Setup**
- [ ] Test notebooks render successfully
- [ ] Word documents open in Microsoft Word
- [ ] Tables are editable in Word
- [ ] Google Docs test shows acceptable formatting

✅ **Phase 2: Pilot**
- [ ] One notebook converted successfully
- [ ] Renders to HTML without errors
- [ ] Renders to Word with proper formatting
- [ ] Tables editable in Google Docs

✅ **Phase 3: Full Migration**
- [ ] All 7 notebooks converted
- [ ] All render to HTML and Word
- [ ] Manual review items addressed
- [ ] All tables tested

---

## 📖 Documentation Quick Links

- **Start here:** [MIGRATION_EXECUTION_GUIDE.md](MIGRATION_EXECUTION_GUIDE.md)
- **Quick reference:** [tests/gt_to_flextable_quick_reference.md](tests/gt_to_flextable_quick_reference.md)
- **Decision guide:** [tests/TABLE_RECOMMENDATIONS_SUMMARY.md](tests/TABLE_RECOMMENDATIONS_SUMMARY.md)
- **Detailed research:** [tests/table_packages_research.md](tests/table_packages_research.md)

---

## ⚠️ Important Notes

### Sandbox Limitation
The VS Code sandbox blocks R execution with temporary directory errors. This is why you need to run scripts in a native terminal.

### Package Installation
The setup script will use `renv::install()` which respects your `~/.condarc` channel redirects to prefix.dev.

### Backup Strategy
All scripts create backups before making changes:
- Pilot script: `benchmark_unique_regions.qmd.backup`
- Full migration: `analysis/backups_gt_YYYYMMDD_HHMMSS/`

### Revert Plan
If anything goes wrong:

```bash
# Revert pilot
mv analysis/benchmark_unique_regions.qmd.backup analysis/benchmark_unique_regions.qmd

# Revert full migration
cp analysis/backups_gt_*/. qmd analysis/

# Delete branch
git checkout main
git branch -D feat/migrate-tables-to-flextable
```

---

## 🤝 Ready to Start?

**Execute this now in Terminal.app:**

```bash
cd /Users/nolson/projects/q100-variant-benchmark-paper
bash setup_flextable_migration.sh
```

Then review the outputs and let me know:
1. Do the test tables look good in Word?
2. Do they look good in Google Docs?
3. Are you ready to proceed with the pilot?

---

## 💬 Next Interaction

After running the setup script, please report:

✅ **If successful:**
- "Setup complete, test outputs look good!"
- I'll guide you through the pilot migration

⚠️ **If issues:**
- Share any error messages
- Describe what looks wrong in the outputs
- I'll help debug and fix

---

**Current status:** ✅ All infrastructure ready, waiting for your execution and feedback.
