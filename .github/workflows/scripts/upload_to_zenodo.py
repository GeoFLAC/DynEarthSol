import os
import sys
import json
import glob
import re
import requests

TOKEN = os.environ.get("ZENODO_TOKEN")
if not TOKEN:
    print("Error: ZENODO_TOKEN is missing in environment variables.")
    sys.exit(1)

TAG_NAME = os.environ.get("GITHUB_REF_NAME", "unknown-version")

CFF_FILE = "CITATION.cff"
ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"
HEADERS = {"Content-Type": "application/json"}
PARAMS = {"access_token": TOKEN}
TIMEOUT = 60         # seconds, for the JSON calls
UPLOAD_TIMEOUT = 600 # seconds, for the file upload

def extract_zenodo_id_from_cff(filepath):
    """
    Reads CITATION.cff and extracts the trailing integer of the Zenodo DOI.
    Matches formats like: doi: 10.5281/zenodo.20293558 or doi: "10.5281/zenodo.20293558"
    """
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found in the repository root.")
        sys.exit(1)
        
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
        
    # Regex to capture the digits after '10.5281/zenodo.'
    match = re.search(r"doi:\s*\"?10\.5281/zenodo\.(\d+)\"?", content)
    if match:
        extracted_id = match.group(1)
        print(f"-> Extracted Zenodo Record ID from {filepath}: {extracted_id}")
        return extracted_id
    else:
        print(f"Error: Could not find a valid Zenodo DOI pattern in {filepath}.")
        print("Expected pattern: doi: 10.5281/zenodo.XXXXX")
        sys.exit(1)

def main():
    print(f"--- Starting Zenodo Upload Process for Release: {TAG_NAME} ---")

    # 1. Dynamically extract the Target ID from CITATION.cff
    TARGET_ID = extract_zenodo_id_from_cff(CFF_FILE)

    # 2. Collect what to upload. This precedes the delete in step 5, which
    # would otherwise clear the draft and leave it with no files at all.
    upload_dir = "output_files"
    files_to_upload = sorted(
        f for f in glob.glob(f"{upload_dir}/*") if os.path.isfile(f))
    if not files_to_upload:
        print(f"Error: no file to upload in '{upload_dir}/'.")
        sys.exit(1)
    print("-> Staged: " + ", ".join(os.path.basename(f) for f in files_to_upload))

    # 3. Fetch the metadata of this ID from Zenodo
    res = requests.get(f"{ZENODO_API_URL}/{TARGET_ID}", params=PARAMS,
                       timeout=TIMEOUT)
    res.raise_for_status()
    draft_data = res.json()
    
    # 4. Double-check the status to enforce workflow integrity
    is_submitted = draft_data.get("submitted", False)
    links = draft_data.get("links", {})

    if is_submitted:
        print(f"Workflow Blocked: The DOI {TARGET_ID} in {CFF_FILE} is already PUBLISHED!")
        print("This means you forgot to update the CITATION.cff file with the new reserved DOI before merging.")
        if "latest_draft" in links:
            # Inform the user what the actual correct draft ID should be
            actual_draft_id = links["latest_draft"].split("/")[-1]
            print(f"Hint: A new draft exists on Zenodo with ID: {actual_draft_id}")
            print(f"Please update the 'doi' field in {CFF_FILE} to point to this new ID, merge it, and re-run.")
        sys.exit(1)

    # If not submitted, it's a valid open draft
    draft_id = draft_data["id"]
    bucket_url = draft_data["links"]["bucket"]

    print(f"Successfully verified target as an active Draft (ID: {draft_id})")

    # 5. Clear old files that Zenodo automatically inherited from the previous version
    existing_files = draft_data.get("files", [])
    for f in existing_files:
        requests.delete(f"{ZENODO_API_URL}/{draft_id}/files/{f['id']}",
                        params=PARAMS, timeout=TIMEOUT).raise_for_status()
    if existing_files:
        print("-> Cleared inherited old files from the draft.")

    # 6. Upload the files staged in step 2
    for filepath in files_to_upload:
        filename = os.path.basename(filepath)
        with open(filepath, "rb") as fp:
            r = requests.put(f"{bucket_url}/{filename}", data=fp, params=PARAMS,
                             timeout=UPLOAD_TIMEOUT)
            r.raise_for_status()
            print(f"-> Uploaded: {filename}")

    # 7. Read .zenodo.json (if exists) and forcefully overwrite the version with GitHub Tag
    metadata_payload = {}
    if os.path.exists(".zenodo.json"):
        with open(".zenodo.json", "r", encoding="utf-8") as f:
            metadata_payload = json.load(f)
        print("-> Successfully loaded base metadata from .zenodo.json")
    else:
        print("-> Warning: .zenodo.json not found, using default metadata.")
        metadata_payload = {
            "title": "Automated Release",
            "upload_type": "software",
            "creators": [{"name": "Automated Uploader"}]
        }

    # Force synchronizing version name with GitHub Release Tag
    metadata_payload["version"] = TAG_NAME

    # Point the deposit at the tag it was built from, the way Zenodo's own
    # GitHub integration does.
    repository = os.environ.get("GITHUB_REPOSITORY")
    if repository:
        tree_url = f"https://github.com/{repository}/tree/{TAG_NAME}"
        metadata_payload.setdefault("related_identifiers", []).append({
            "relation": "isSupplementTo",
            "identifier": tree_url,
            "resource_type": "software",
        })
        print(f"-> Linked as a supplement to {tree_url}")
    else:
        print("WARNING: GITHUB_REPOSITORY is unset, so the deposit will not "
              "name the tag it came from.")

    meta_res = requests.put(
        f"{ZENODO_API_URL}/{draft_id}", 
        json={"metadata": metadata_payload}, 
        params=PARAMS, 
        headers=HEADERS,
        timeout=TIMEOUT
    )
    meta_res.raise_for_status()
    
    print(f"\nSUCCESS! Files have been uploaded to the manual Draft, and version updated to {TAG_NAME}.")
    print(f"🔗 Please review and publish manually at: https://zenodo.org/deposit/{draft_id}")

if __name__ == "__main__":
    main()