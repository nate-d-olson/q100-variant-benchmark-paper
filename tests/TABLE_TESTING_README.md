# Table Compatibility Testing

This directory contains test notebooks for evaluating R table packages for Word and Google Docs compatibility.

## Problem Statement

The `gt` package creates beautiful tables for HTML output but loses formatting when:

- Exported to Word (.docx) via Quarto/RMarkdown
- Copied and pasted into Google Docs
- Shared with collaborators who use Office applications

This is a significant issue for manuscript preparation and collaborative research.

## Test Files

### 1. `table_compatibility_test.qmd`

**Purpose**: Comprehensive comparison of table packages (flextable, kableExtra, huxtable, gt)

**Contents**:

- Sample data mimicking project tables
- Basic and styled examples for each package
- Grouped tables
- Comparison matrix
- Testing instructions

**Renders to**:

- `table_compatibility_test.html` - View in browser
- `table_compatibility_test.docx` - Test Word compatibility

### 2. `flextable_theme_implementation.qmd`

**Purpose**: Ready-to-use flextable theme functions replicating `theme_gt_manuscript()`

**Contents**:

- `theme_flextable_manuscript()` function definition
- Helper functions matching `gt` API
- Usage examples
- Migration guide
- Direct Word document creation examples

**Use**: Copy functions to `R/flextable_themes.R` for production use

### 3. `table_packages_research.md`

**Purpose**: Complete research documentation and recommendations

**Contents**:

- Detailed package comparisons
- Pros/cons for each approach
- Implementation guide
- Migration checklist
- Links to documentation and resources

## Running the Tests

### Prerequisites

```r
# Install required packages (in R with renv activated)
renv::install(c("flextable", "kableExtra", "huxtable", "officer"))

# If using system Quarto:
# brew install quarto  # macOS
# Or download from: https://quarto.org/docs/get-started/
```

### Render Test Notebooks

```r
# From R console
quarto::quarto_render("tests/table_compatibility_test.qmd", output_format = "all")
quarto::quarto_render("tests/flextable_theme_implementation.qmd")

# From terminal
quarto render tests/table_compatibility_test.qmd
quarto render tests/flextable_theme_implementation.qmd
```

## Testing Workflow

### 1. Visual Inspection (HTML)

```bash
# Open HTML output in browser
open tests/table_compatibility_test.html
open tests/flextable_theme_implementation.html
```

**Check**:

- ✓ Table formatting (borders, colors, fonts)
- ✓ Number formatting (commas, decimals)
- ✓ Column alignment
- ✓ Grouped table structure
- ✓ Overall appearance

### 2. Word Export Testing

```bash
# Open Word document
open tests/table_compatibility_test.docx
```

**Checklist**:

- [ ] **Borders**: Top/bottom 2px black borders preserved?
- [ ] **Header formatting**: Bold text, gray background (#F5F5F5)?
- [ ] **Fonts**: Roboto (or fallback Arial) at correct sizes?
- [ ] **Row striping**: Alternating row colors visible?
- [ ] **Number formatting**: Commas in large numbers? Correct decimals?
- [ ] **Column alignment**: Headers centered, numeric columns right-aligned?
- [ ] **Grouping**: Merged cells for grouped columns work correctly?

**Known issues to watch for**:

- kableExtra: Most Bootstrap styling disappears
- gt: Tables may export as images or lose all formatting
- huxtable: Occasional border inconsistencies

### 3. Google Docs Testing

**Method 1: Upload .docx to Google Drive**

1. Upload `table_compatibility_test.docx` to Google Drive
2. Right-click → Open with → Google Docs
3. Check formatting preservation

**Method 2: Copy from Word, Paste to Google Docs**

1. Open `table_compatibility_test.docx` in Microsoft Word
2. Select table → Copy (Cmd+C)
3. Open new Google Doc → Paste (Cmd+V)
4. Check formatting

**Method 3: Copy from HTML, Paste to Google Docs**

1. Open `table_compatibility_test.html` in browser
2. Select table → Copy
3. Paste into Google Doc
4. Check formatting

**Compare results**: Which method preserves formatting best?

### 4. Editing Testing

In Word or Google Docs:

1. Try editing table content
2. Try adding/removing rows
3. Try changing column widths
4. Verify no corruption or strange behavior

### 5. Performance Testing

For large tables (100+ rows):

```r
# Generate large test dataset
large_data <- tibble(
  id = 1:500,
  value1 = rnorm(500),
  value2 = rpois(500, 10),
  category = sample(letters[1:5], 500, replace = TRUE)
)

# Test rendering time
system.time({
  large_data %>%
    flextable() %>%
    theme_flextable_manuscript()
})
```

## Expected Results

### flextable (Best for Word/Google Docs)

**HTML output**: ⭐⭐⭐⭐ Good (clean, professional)  
**Word export**: ⭐⭐⭐⭐⭐ Excellent (formatting preserved)  
**Google Docs**: ⭐⭐⭐⭐ Good (most formatting preserved)  
**Editability**: ⭐⭐⭐⭐⭐ Excellent (native Word table)

### kableExtra

**HTML output**: ⭐⭐⭐⭐⭐ Excellent (Bootstrap styling)  
**Word export**: ⭐⭐ Poor (most styling lost)  
**Google Docs**: ⭐⭐ Poor (basic text only)  
**Editability**: ⭐⭐⭐ Fair (loses styling when edited)

### huxtable

**HTML output**: ⭐⭐⭐ Fair (functional but plain)  
**Word export**: ⭐⭐⭐⭐ Good (mostly preserved)  
**Google Docs**: ⭐⭐⭐ Fair (some formatting lost)  
**Editability**: ⭐⭐⭐⭐ Good (editable table)

### gt (Current Package)

**HTML output**: ⭐⭐⭐⭐⭐ Excellent (beautiful design)  
**Word export**: ⭐ Very Poor (formatting lost or image export)  
**Google Docs**: ⭐ Very Poor (loses all styling)  
**Editability**: ⭐ Very Poor (image or unformatted text)

## Decision Criteria

Choose **flextable** if:

- ✓ Tables need Word export
- ✓ Collaborators use Google Docs
- ✓ Cross-format consistency required
- ✓ Programmatic Word document generation needed

Keep **gt** if:

- ✓ Output is HTML-only (website, supplementary materials)
- ✓ Maximum aesthetic control needed for web
- ✓ Interactive features required

Use **hybrid approach** if:

- ✓ Different tables have different requirements
- ✓ Want to maintain existing `gt` code
- ✓ Need both HTML beauty and Word compatibility

## Troubleshooting

### "Package not found" error

```r
# Ensure packages are installed
renv::install("flextable")
renv::restore()  # Restore from lockfile
```

### Fonts not rendering

```r
# flextable uses system fonts
# If Roboto not available, install it or use fallback:
theme_flextable_manuscript(ft, font_family = "Arial")
```

### Quarto rendering fails

```r
# Check Quarto version
system("quarto --version")  # Should be 1.3+

# Update Quarto if needed:
# https://quarto.org/docs/download/
```

### Word document won't open

- Check that `officer` package is installed
- Ensure write permissions to output directory
- Try explicitly specifying output path:

```r
quarto::quarto_render(
  "tests/table_compatibility_test.qmd",
  output_format = "docx",
  output_file = "test_output.docx"
)
```

### Tables look different in Word vs HTML

This is expected! Different rendering engines handle styling differently. Key is that:

- Core formatting (borders, alignment) is preserved
- Numbers display correctly
- Table is editable in Word

## Next Steps

After testing:

1. **Review results** - Which package best meets your needs?
2. **Read `table_packages_research.md`** - Understand tradeoffs
3. **Copy theme function** - From `flextable_theme_implementation.qmd` to `R/flextable_themes.R`
4. **Test with real data** - Use actual project tables
5. **Pilot migration** - Update one analysis notebook
6. **Full migration** - Update all tables if pilot successful

## Additional Resources

- **flextable book**: <https://ardata-fr.github.io/flextable-book/>
- **officer documentation**: <https://ardata-fr.github.io/officeverse/>
- **Quarto tables guide**: <https://quarto.org/docs/authoring/tables.html>
- **R Markdown Word output**: <https://bookdown.org/yihui/rmarkdown/word-document.html>

## Questions?

Common questions and answers:

**Q: Do I need to refactor all existing tables?**  
A: Only if you need Word/Google Docs export. HTML-only files can keep `gt`.

**Q: Can I use both `gt` and `flextable` in same project?**  
A: Yes! Use conditional rendering or separate files.

**Q: Will this work with Quarto websites?**  
A: Yes, flextable renders fine to HTML. Not as polished as `gt`, but functional.

**Q: What about PDF output?**  
A: Both `gt` and `flextable` work with PDF. Test your specific tables.

**Q: Can I customize the theme further?**  
A: Yes! See `flextable_theme_implementation.qmd` for examples of theme customization.
