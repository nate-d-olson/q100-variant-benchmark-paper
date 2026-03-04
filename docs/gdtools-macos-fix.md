# gdtools Compilation Troubleshooting on macOS

## Problem

When installing flextable or its dependencies, you may encounter this error:

```
library not found for -lintl
arm64-apple-darwin20.0: error: linker command failed with exit code 1
```

## Cause

The `gdtools` package (a dependency of flextable) requires:
- **gettext** library (provides libintl for internationalization)
- Proper compiler flags to find Homebrew-installed libraries

On Apple Silicon Macs, Homebrew installs to `/opt/homebrew/`, and R's compiler may not find libraries there by default.

## Solutions

### Solution 1: Install gettext (Recommended)

**Using Homebrew (automatic in setup script):**

```bash
brew install gettext
```

The `.Renviron` file in the project root configures compiler flags automatically.

### Solution 2: Manual Compiler Configuration

If Solution 1 doesn't work, configure R's compiler flags manually:

**For Apple Silicon (M1/M2/M3):**

Create or edit `~/.R/Makevars`:

```makefile
# Homebrew paths for Apple Silicon
CPPFLAGS=-I/opt/homebrew/include -I/opt/homebrew/opt/gettext/include
LDFLAGS=-L/opt/homebrew/lib -L/opt/homebrew/opt/gettext/lib
```

**For Intel Macs:**

```makefile
# Homebrew paths for Intel Macs
CPPFLAGS=-I/usr/local/include -I/usr/local/opt/gettext/include
LDFLAGS=-L/usr/local/lib -L/usr/local/opt/gettext/lib
```

### Solution 3: Use Pre-compiled Binaries

Instead of compiling from source, use precompiled binaries:

```r
# In R console
options(pkgType = "both")  # Prefer binaries when available
renv::install("gdtools")
renv::install("flextable")
```

### Solution 4: System Package Manager

**Using MacPorts (if you don't use Homebrew):**

```bash
sudo port install gettext
```

Then update compiler flags in `~/.R/Makevars`:

```makefile
CPPFLAGS=-I/opt/local/include
LDFLAGS=-L/opt/local/lib
```

## Verification

After applying a solution, test the installation:

```r
# In R console
Sys.which("pkg-config")  # Should find pkg-config
system("pkg-config --libs libintl")  # Should show library flags

# Try installing gdtools
install.packages("gdtools")

# If successful, install flextable
renv::install("flextable")
```

## Additional Dependencies

gdtools may also require:

- **cairo** - Graphics library
- **freetype** - Font rendering
- **fontconfig** - Font configuration

Install all at once:

```bash
brew install gettext cairo freetype fontconfig
```

## Project Configuration

This project includes a `.Renviron` file with compiler flags for Apple Silicon. If you're on an Intel Mac, edit `.Renviron` and uncomment the Intel-specific lines:

```r
# .Renviron
# For Intel Macs - uncomment these lines:
CPPFLAGS=-I/usr/local/include -I/usr/local/opt/gettext/include
LDFLAGS=-L/usr/local/lib -L/usr/local/opt/gettext/lib
```

## Still Having Issues?

1. **Check Homebrew installation:**
   ```bash
   brew doctor
   brew update
   ```

2. **Verify gettext installation:**
   ```bash
   brew list gettext
   ls -la /opt/homebrew/opt/gettext/lib/libintl.*
   ```

3. **Check R compiler configuration:**
   ```r
   # In R console
   Sys.getenv("CPPFLAGS")
   Sys.getenv("LDFLAGS")
   ```

4. **Try reinstalling R build tools:**
   ```bash
   xcode-select --install
   ```

5. **Use binary packages only:**
   ```r
   # In R console
   options(pkgType = "mac.binary")
   renv::install("flextable")
   ```

## Alternative: Skip gdtools

If you absolutely cannot compile gdtools, you can use flextable without it (some features may be limited):

```r
# Install without gdtools-dependent features
Sys.setenv(R_COMPILE_PKGS = "0")
renv::install("flextable")
```

Note: This may limit font rendering capabilities.

## Reference Links

- [gdtools GitHub Issues](https://github.com/davidgohel/gdtools/issues)
- [Homebrew gettext](https://formulae.brew.sh/formula/gettext)
- [RStudio Community: macOS compilation issues](https://community.rstudio.com/tag/macos)
- [R macOS Developer Tools](https://mac.r-project.org/tools/)
