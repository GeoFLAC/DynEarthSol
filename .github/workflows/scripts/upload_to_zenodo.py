import os
import sys
import json
import glob
import requests

TOKEN = os.environ.get("ZENODO_TOKEN")
if not TOKEN:
    print("Error: ZENODO_TOKEN is missing in environment variables.")
    sys.exit(1)

TAG_NAME = os.environ.get("GITHUB_REF_NAME", "unknown-version")

# [TODO]: Replace this with the Version Record ID of the latest published release in this concept series.
LAST_DEPOSITION_ID = "20293558" 

ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"
HEADERS = {"Content-Type": "application/json"}
PARAMS = {"access_token": TOKEN}

def main():
    print(f"--- Starting Zenodo Upload Process for Release: {TAG_NAME} ---")

    # 1. Check if a manual draft has already been created on Zenodo
    res = requests.get(f"{ZENODO_API_URL}/{LAST_DEPOSITION_ID}", params=PARAMS)
    res.raise_for_status()
    links = res.json().get("links", {})

    if "latest_draft" not in links:
        print("Error: No unpublished Draft found on Zenodo!")
        print("Failsafe triggered: Please ensure you have manually clicked 'New version' "
              "on the Zenodo website BEFORE publishing this GitHub Release.")
        sys.exit(1)

    draft_url = links["latest_draft"]
    draft_res = requests.get(draft_url, params=PARAMS)
    draft_res.raise_for_status()
    draft_data = draft_res.json()
    
    draft_id = draft_data["id"]
    bucket_url = draft_data["links"]["bucket"]

    print(f"Successfully found the manually opened Draft (ID: {draft_id})")

    # 2. Clear old files that Zenodo automatically inherited from the previous version
    existing_files = draft_data.get("files", [])
    for f in existing_files:
        requests.delete(f"{ZENODO_API_URL}/{draft_id}/files/{f['id']}", params=PARAMS).raise_for_status()
    if existing_files:
        print("-> Cleared inherited old files from the draft.")

    # 3. Upload all files from the output_files directory
    upload_dir = "output_files"
    files_to_upload = glob.glob(f"{upload_dir}/*")
    
    if not files_to_upload:
        print(f"Warning: No files found in '{upload_dir}/' to upload.")
    
    for filepath in files_to_upload:
        filename = os.path.basename(filepath)
        with open(filepath, "rb") as fp:
            r = requests.put(f"{bucket_url}/{filename}", data=fp, params=PARAMS)
            r.raise_for_status()
            print(f"-> Uploaded: {filename}")

    # 4. Read .zenodo.json (if exists) and forcefully overwrite the version with GitHub Tag
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

    metadata_payload["version"] = TAG_NAME
    
    meta_res = requests.put(
        draft_url, 
        json={"metadata": metadata_payload}, 
        params=PARAMS, 
        headers=HEADERS
    )
    meta_res.raise_for_status()
    
    print(f"\nSUCCESS! Files have been uploaded to the manual Draft, and version updated to {TAG_NAME}.")
    print(f"🔗 Please review and publish manually at: https://zenodo.org/deposit/{draft_id}")

if __name__ == "__main__":
    main()