# S3 Object Checker

Script to verify and fix S3 object permissions and storage tiers for URLs in `config/config.yaml`.

## Features

- **Object Existence**: Verifies that objects exist in S3
- **ACL Checking**: Ensures objects have `public-read` ACL
- **Storage Tier**: Checks that objects are in a storage class that allows immediate download
- **Automatic Fixes**: Can automatically set ACLs and initiate Glacier restores
- **Dry Run Mode**: Preview changes before making them

## Installation

Install required dependencies:

```bash
pip install -r scripts/requirements-s3-check.txt
```

## AWS Credentials

The script uses boto3 and requires AWS credentials. Configure them using:

```bash
# Option 1: AWS CLI
aws configure

# Option 2: Environment variables
export AWS_ACCESS_KEY_ID="your-access-key"
export AWS_SECRET_ACCESS_KEY="your-secret-key"
export AWS_DEFAULT_REGION="us-east-1"

# Option 3: Use a specific AWS profile
python3 scripts/check_s3_objects.py config/config.yaml --profile myprofile
```

## Usage

### Check-only mode (no changes)

```bash
# Check all URLs in config
python3 scripts/check_s3_objects.py config/config.yaml

# Verbose output
python3 scripts/check_s3_objects.py config/config.yaml --verbose
```

### Fix mode (makes changes)

```bash
# See what would be fixed (dry run)
python3 scripts/check_s3_objects.py config/config.yaml --fix --dry-run

# Actually fix issues
python3 scripts/check_s3_objects.py config/config.yaml --fix
```

### Check a specific URL

```bash
# Check single URL
python3 scripts/check_s3_objects.py --url s3://giab-data/path/to/file.bed --fix
```

## What the Script Does

### 1. Object Existence
- Uses `HeadObject` to verify the object exists
- Reports missing objects

### 2. ACL Permissions
- Checks if object has `public-read` ACL
- If not, sets ACL to `public-read` (when `--fix` is used)
- Looks for AllUsers group with READ permission

### 3. Storage Class
- Checks the object's storage class
- For Glacier storage classes (GLACIER, DEEP_ARCHIVE, GLACIER_IR):
  - Checks if restore is in progress or complete
  - Initiates restore request if needed (when `--fix` is used)
  - Uses Expedited tier for GLACIER/GLACIER_IR
  - Uses Standard tier for DEEP_ARCHIVE
  - Keeps restored objects available for 30 days

## Storage Classes

The script handles these storage classes:

- **STANDARD** - Immediate download (OK)
- **STANDARD_IA** - Immediate download (OK)
- **ONEZONE_IA** - Immediate download (OK)
- **INTELLIGENT_TIERING** - Immediate download (OK)
- **GLACIER** - Needs restore (1-5 minutes for Expedited)
- **GLACIER_IR** - Needs restore (minutes)
- **DEEP_ARCHIVE** - Needs restore (12 hours)

## Example Output

```
Parsing config: config/config.yaml
Found 125 S3 URLs

*** FIX MODE - Issues will be corrected ***

https://giab-data.s3.amazonaws.com/defrabb_runs/.../file.bed
  ✓ Object exists
  ✗ ACL is NOT public-read
  ✓ Set ACL to public-read
  ✓ Storage class STANDARD allows immediate download

================================================================================
SUMMARY
================================================================================
Total URLs checked: 125
Objects exist: 125/125
Public-readable: 125/125
Need restore: 0
Errors: 0
Objects with actions taken: 15
```

## Exit Codes

- `0` - Success (or all issues fixed)
- `1` - Issues found (when not using `--fix`)

## Notes

- S3 restore requests for Glacier storage may take time to complete
- Changing ACLs requires appropriate IAM permissions
- The script only processes S3 URLs (filters out FTP and other URLs)
- Supported S3 URL formats:
  - `s3://bucket/key`
  - `https://bucket.s3.amazonaws.com/key`
  - `https://s3.amazonaws.com/bucket/key`
  - `https://s3-region.amazonaws.com/bucket/key`
