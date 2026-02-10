#!/usr/bin/env python3
"""
Check and fix S3 object permissions and storage tiers for URLs in config.yaml.

This script:
1. Parses config.yaml to find all S3 URLs
2. Verifies objects exist in S3
3. Checks object ACLs (should be public-read)
4. Checks storage tiers (should allow immediate download)
5. Updates ACL and storage tier as needed
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from urllib.parse import urlparse

import boto3
import yaml
from botocore.exceptions import ClientError, NoCredentialsError


def parse_s3_url(url: str) -> Tuple[str, str]:
    """
    Parse S3 URL to extract bucket and key.

    Supports formats:
    - https://bucket.s3.amazonaws.com/key
    - https://s3.amazonaws.com/bucket/key
    - https://s3-region.amazonaws.com/bucket/key
    - s3://bucket/key

    Returns:
        Tuple of (bucket_name, object_key)
    """
    # Handle s3:// protocol
    if url.startswith("s3://"):
        parsed = urlparse(url)
        return parsed.netloc, parsed.path.lstrip("/")

    # Handle https URLs
    parsed = urlparse(url)
    hostname = parsed.hostname
    path = parsed.path.lstrip("/")

    # Format: bucket.s3.amazonaws.com or bucket.s3-region.amazonaws.com
    if hostname and ".s3" in hostname and hostname.endswith(".amazonaws.com"):
        bucket = hostname.split(".s3")[0]
        key = path
        return bucket, key

    # Format: s3.amazonaws.com/bucket/key or s3-region.amazonaws.com/bucket/key
    if hostname and hostname.startswith("s3") and hostname.endswith(".amazonaws.com"):
        parts = path.split("/", 1)
        if len(parts) == 2:
            return parts[0], parts[1]
        elif len(parts) == 1:
            return parts[0], ""

    raise ValueError(f"Cannot parse S3 URL: {url}")


def is_s3_url(url: str) -> bool:
    """Check if URL is an S3 URL."""
    return (
        url.startswith("s3://")
        or ".s3.amazonaws.com" in url
        or ".s3-" in url
        and ".amazonaws.com" in url
        or url.startswith("https://s3.amazonaws.com/")
        or url.startswith("https://s3-")
    )


def extract_urls_from_dict(data: dict, urls: Set[str]) -> None:
    """Recursively extract URLs from nested dictionary."""
    if isinstance(data, dict):
        for key, value in data.items():
            if key == "url" and isinstance(value, str):
                urls.add(value)
            elif isinstance(value, (dict, list)):
                extract_urls_from_dict(value, urls)
    elif isinstance(data, list):
        for item in data:
            extract_urls_from_dict(item, urls)


def get_s3_urls_from_config(config_path: Path) -> List[str]:
    """Extract all S3 URLs from config YAML file."""
    with open(config_path) as f:
        config = yaml.safe_load(f)

    all_urls = set()
    extract_urls_from_dict(config, all_urls)

    # Filter to only S3 URLs
    s3_urls = [url for url in all_urls if is_s3_url(url)]
    return sorted(s3_urls)


def check_object_exists(s3_client, bucket: str, key: str) -> bool:
    """Check if S3 object exists."""
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        if e.response["Error"]["Code"] == "404":
            return False
        raise


def get_object_info(s3_client, bucket: str, key: str) -> Dict:
    """Get object metadata including storage class."""
    try:
        response = s3_client.head_object(Bucket=bucket, Key=key)
        return {
            "storage_class": response.get("StorageClass", "STANDARD"),
            "restore": response.get("Restore"),
            "content_length": response.get("ContentLength", 0),
        }
    except ClientError as e:
        raise Exception(f"Failed to get object info: {e}")


def get_object_acl(s3_client, bucket: str, key: str) -> Dict:
    """Get object ACL."""
    try:
        response = s3_client.get_object_acl(Bucket=bucket, Key=key)
        return response
    except ClientError as e:
        raise Exception(f"Failed to get object ACL: {e}")


def is_publicly_readable(acl: Dict) -> bool:
    """Check if object has public-read ACL."""
    for grant in acl.get("Grants", []):
        grantee = grant.get("Grantee", {})
        permission = grant.get("Permission")

        # Check for AllUsers group with READ permission
        if (
            grantee.get("Type") == "Group"
            and grantee.get("URI") == "http://acs.amazonaws.com/groups/global/AllUsers"
            and permission == "READ"
        ):
            return True

    return False


def set_public_read_acl(
    s3_client, bucket: str, key: str, dry_run: bool = False
) -> None:
    """Set object ACL to public-read."""
    if dry_run:
        print("  [DRY RUN] Would set ACL to public-read")
        return

    try:
        s3_client.put_object_acl(Bucket=bucket, Key=key, ACL="public-read")
        print("  ✓ Set ACL to public-read")
    except ClientError as e:
        raise Exception(f"Failed to set ACL: {e}")


def needs_restore(storage_class: str, restore_status: str) -> bool:
    """Check if object needs to be restored from Glacier."""
    glacier_classes = ["GLACIER", "DEEP_ARCHIVE", "GLACIER_IR"]

    if storage_class not in glacier_classes:
        return False

    # If already restored or in progress
    if restore_status:
        if 'ongoing-request="true"' in restore_status:
            return False  # Restoration in progress
        elif 'ongoing-request="false"' in restore_status:
            return False  # Already restored

    return True


def restore_object(
    s3_client, bucket: str, key: str, storage_class: str, dry_run: bool = False
) -> None:
    """Restore object from Glacier or change storage class."""
    if dry_run:
        print(f"  [DRY RUN] Would restore from {storage_class}")
        return

    glacier_classes = ["GLACIER", "DEEP_ARCHIVE", "GLACIER_IR"]

    if storage_class in glacier_classes:
        # Initiate restore request
        try:
            restore_config = {
                "Days": 30,  # Keep restored copy for 30 days
                "GlacierJobParameters": {
                    "Tier": (
                        "Expedited" if storage_class != "DEEP_ARCHIVE" else "Standard"
                    )
                },
            }

            # Deep Archive doesn't support Expedited
            if storage_class == "DEEP_ARCHIVE":
                restore_config["GlacierJobParameters"]["Tier"] = "Standard"

            s3_client.restore_object(
                Bucket=bucket, Key=key, RestoreRequest=restore_config
            )
            print(f"  ✓ Initiated restore from {storage_class}")
        except ClientError as e:
            if e.response["Error"]["Code"] == "RestoreAlreadyInProgress":
                print("  ℹ Restore already in progress")
            else:
                raise Exception(f"Failed to restore object: {e}")
    else:
        # For non-Glacier classes that aren't STANDARD, consider copying to STANDARD
        print(f"  ℹ Object in {storage_class} - can be downloaded immediately")


def check_and_fix_object(
    s3_client, url: str, dry_run: bool = False, verbose: bool = False
) -> Dict:
    """
    Check and fix a single S3 object.

    Returns:
        Dictionary with status information
    """
    result = {
        "url": url,
        "exists": False,
        "public": False,
        "storage_class": None,
        "needs_restore": False,
        "actions_taken": [],
    }

    try:
        bucket, key = parse_s3_url(url)
    except ValueError as e:
        result["error"] = str(e)
        return result

    if verbose:
        print(f"\nChecking: {url}")
        print(f"  Bucket: {bucket}")
        print(f"  Key: {key}")
    else:
        print(f"\n{url}")

    # Check if object exists
    try:
        exists = check_object_exists(s3_client, bucket, key)
        result["exists"] = exists

        if not exists:
            print("  ✗ Object does not exist")
            result["error"] = "Object not found"
            return result

        if verbose:
            print("  ✓ Object exists")

    except Exception as e:
        print(f"  ✗ Error checking existence: {e}")
        result["error"] = str(e)
        return result

    # Get object info
    try:
        info = get_object_info(s3_client, bucket, key)
        result["storage_class"] = info["storage_class"]

        if verbose:
            print(f"  Storage class: {info['storage_class']}")
            print(f"  Size: {info['content_length']:,} bytes")

    except Exception as e:
        print(f"  ✗ Error getting object info: {e}")
        result["error"] = str(e)
        return result

    # Check ACL
    try:
        acl = get_object_acl(s3_client, bucket, key)
        is_public = is_publicly_readable(acl)
        result["public"] = is_public

        if is_public:
            if verbose:
                print("  ✓ ACL is public-read")
        else:
            print("  ✗ ACL is NOT public-read")
            set_public_read_acl(s3_client, bucket, key, dry_run)
            result["actions_taken"].append("set_acl")

    except Exception as e:
        print(f"  ⚠ Error checking/setting ACL: {e}")
        result["acl_error"] = str(e)

    # Check storage class and restore if needed
    if needs_restore(info["storage_class"], info.get("restore")):
        print(f"  ✗ Object in {info['storage_class']} - needs restore")
        result["needs_restore"] = True
        try:
            restore_object(s3_client, bucket, key, info["storage_class"], dry_run)
            result["actions_taken"].append("restore")
        except Exception as e:
            print(f"  ⚠ Error restoring object: {e}")
            result["restore_error"] = str(e)
    elif info["storage_class"] in ["GLACIER", "DEEP_ARCHIVE", "GLACIER_IR"]:
        if info.get("restore"):
            print(
                f"  ℹ Object in {info['storage_class']} but restore is in progress/complete"
            )
    else:
        if verbose:
            print(
                f"  ✓ Storage class {info['storage_class']} allows immediate download"
            )

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Check and fix S3 object permissions and storage tiers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check all objects (read-only)
  %(prog)s config/config.yaml

  # Check and fix objects
  %(prog)s config/config.yaml --fix

  # Show what would be fixed without making changes
  %(prog)s config/config.yaml --fix --dry-run

  # Check specific URL
  %(prog)s --url s3://bucket/path/to/object --fix
        """,
    )

    parser.add_argument("config", type=Path, nargs="?", help="Path to config.yaml file")
    parser.add_argument(
        "--url", help="Check a specific S3 URL instead of parsing config"
    )
    parser.add_argument(
        "--fix", action="store_true", help="Fix issues (set ACL and restore objects)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--profile", help="AWS profile to use")

    args = parser.parse_args()

    if not args.url and not args.config:
        parser.error("Either config file or --url must be provided")

    # Initialize S3 client
    try:
        session_kwargs = {}
        if args.profile:
            session_kwargs["profile_name"] = args.profile

        session = boto3.Session(**session_kwargs)
        s3_client = session.client("s3")
    except NoCredentialsError:
        print("Error: AWS credentials not found.", file=sys.stderr)
        print(
            "Please configure credentials using 'aws configure' or set environment variables.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Get URLs to check
    if args.url:
        urls = [args.url]
    else:
        if not args.config.exists():
            print(f"Error: Config file not found: {args.config}", file=sys.stderr)
            sys.exit(1)

        print(f"Parsing config: {args.config}")
        urls = get_s3_urls_from_config(args.config)
        print(f"Found {len(urls)} S3 URLs")

    if args.dry_run:
        print("\n*** DRY RUN MODE - No changes will be made ***\n")
    elif args.fix:
        print("\n*** FIX MODE - Issues will be corrected ***\n")
    else:
        print("\n*** CHECK MODE - Read-only ***\n")

    # Check each URL
    results = []
    for url in urls:
        result = check_and_fix_object(
            s3_client, url, dry_run=args.dry_run or not args.fix, verbose=args.verbose
        )
        results.append(result)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    total = len(results)
    exists = sum(1 for r in results if r["exists"])
    public = sum(1 for r in results if r.get("public"))
    needs_restore = sum(1 for r in results if r.get("needs_restore"))
    errors = sum(1 for r in results if "error" in r)
    actions = sum(1 for r in results if r.get("actions_taken"))

    print(f"Total URLs checked: {total}")
    print(f"Objects exist: {exists}/{total}")
    print(f"Public-readable: {public}/{exists if exists > 0 else total}")
    print(f"Need restore: {needs_restore}")
    print(f"Errors: {errors}")
    print(f"Objects with actions taken: {actions}")

    # List objects that need attention
    issues = [
        r
        for r in results
        if not r.get("exists")
        or not r.get("public")
        or r.get("needs_restore")
        or "error" in r
    ]

    if issues:
        print(f"\n{len(issues)} object(s) need attention:")
        for r in issues:
            print(f"  - {r['url']}")
            if not r.get("exists"):
                print("    Issue: Object does not exist")
            if r.get("exists") and not r.get("public"):
                print("    Issue: Not publicly readable")
            if r.get("needs_restore"):
                print(f"    Issue: Needs restore from {r.get('storage_class')}")
            if "error" in r:
                print(f"    Error: {r['error']}")

    # Exit with error code if there are issues and we're not fixing them
    if issues and not args.fix:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
