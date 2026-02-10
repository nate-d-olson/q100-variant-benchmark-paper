import yaml

print("Script starting...", flush=True)
import hashlib
import urllib.request
import os
import sys

CONFIG_FILE = "config/config.yaml"
OUTPUT_FILE = "checksum_verification.txt"


def calculate_checksum(file_path, algorithm):
    hash_func = getattr(hashlib, algorithm)()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_func.update(chunk)
    return hash_func.hexdigest()


def verify_url(url, expected_checksum, algorithm):
    print(f"Checking {url}...", flush=True)
    try:
        temp_filename = "temp_download_file"
        # Set a user agent to avoid 403 errors from some servers
        opener = urllib.request.build_opener()
        opener.addheaders = [("User-agent", "Mozilla/5.0")]
        urllib.request.install_opener(opener)

        urllib.request.urlretrieve(url, temp_filename)
        calculated_checksum = calculate_checksum(temp_filename, algorithm)
        if os.path.exists(temp_filename):
            os.remove(temp_filename)

        if calculated_checksum == expected_checksum:
            return True, f"PASS:    ({algorithm}: {calculated_checksum}) {url}"
        else:
            return (
                False,
                f"FAIL:   (Expected {algorithm}: {expected_checksum}, Got: {calculated_checksum}) {url}",
            )
    except Exception as e:
        if os.path.exists("temp_download_file"):
            os.remove("temp_download_file")
        return False, f"ERROR:  ({str(e)})  {url}"


def traverse_and_verify(data, results):
    if isinstance(data, dict):
        if "url" in data and ("sha256" in data or "md5" in data):
            url = data["url"]
            if "sha256" in data:
                algorithm = "sha256"
                expected = data["sha256"]
            else:
                algorithm = "md5"
                expected = data["md5"]

            success, message = verify_url(url, expected, algorithm)
            results.append(message)
            print(message)
            # Continue traversing in case there are nested objects, though unlikely in this specific schema
            # But wait, looking at the yaml, 'stratifications' is inside the object that has a url (for references).
            # The structure is:
            # key:
            #   url: ...
            #   md5: ...
            #   stratifications: ...
            # So I should continue to traverse children even if I found a URL here.

        for key, value in data.items():
            # Avoid re-checking the 'url' string itself, but check other keys
            if key != "url" and key != "sha256" and key != "md5":
                traverse_and_verify(value, results)

    elif isinstance(data, list):
        for item in data:
            traverse_and_verify(item, results)


def main():
    if not os.path.exists(CONFIG_FILE):
        print(f"Config file not found: {CONFIG_FILE}")
        sys.exit(1)

    with open(CONFIG_FILE, "r") as f:
        config = yaml.safe_load(f)

    results = []
    traverse_and_verify(config, results)

    with open(OUTPUT_FILE, "w") as f:
        f.write("\n".join(results) + "\n")

    print(f"Verification complete. Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
