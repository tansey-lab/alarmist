#!/usr/bin/env python3
"""
Sync a pipeline to Cirro from CI.

This script uses a refresh token (from GitHub secrets) to authenticate
without interactive login and sync the pipeline configuration to Cirro.

Required environment variables:
    CIRRO_REFRESH_TOKEN: The refresh token from an authenticated session
    CIRRO_BASE_URL: The Cirro instance URL (e.g., breakthroughcancer.cirro.bio)

Optional environment variables:
    CIRRO_PROCESS_ID: The process ID to sync (if known)
    CIRRO_PROCESS_NAME: The process name to sync (alternative to ID)

To get the refresh token for CI:
    1. Run `cirro configure` locally and save login
    2. Run `python scripts/check_token.py` to verify auth works
    3. Extract the refresh_token from ~/.cirro/*.token.dat
       (requires keychain access on macOS)
    4. Store it as a GitHub secret named CIRRO_REFRESH_TOKEN
"""

import os
import sys

import boto3
import requests
from botocore.exceptions import ClientError
from cirro import CirroApi
from cirro.auth.access_token import AccessTokenAuth


def get_cirro_system_info(base_url: str) -> dict:
    """Fetch Cirro system info to get Cognito configuration."""
    rest_endpoint = f"https://{base_url}/api"
    resp = requests.get(f"{rest_endpoint}/info/system")
    resp.raise_for_status()
    return resp.json()


def refresh_access_token(refresh_token: str, client_id: str, region: str) -> str:
    """Use the refresh token to get a new access token from Cognito."""
    cognito = boto3.client("cognito-idp", region_name=region)

    try:
        resp = cognito.initiate_auth(
            ClientId=client_id,
            AuthFlow="REFRESH_TOKEN_AUTH",
            AuthParameters={"REFRESH_TOKEN": refresh_token},
        )
        return resp["AuthenticationResult"]["AccessToken"]
    except ClientError as e:
        error_code = e.response.get("Error", {}).get("Code", "")
        if error_code == "NotAuthorizedException":
            print("ERROR: Refresh token is invalid or expired.", file=sys.stderr)
            print(
                "Re-authenticate and update the CIRRO_REFRESH_TOKEN secret.",
                file=sys.stderr,
            )
        raise


def find_process_by_name(cirro: CirroApi, name: str):
    """Find a process by its display name."""
    processes = cirro.processes.list()
    for p in processes:
        if p.name == name:
            return p
    return None


def main():
    # Required env vars
    refresh_token = os.environ.get("CIRRO_REFRESH_TOKEN")
    base_url = os.environ.get("CIRRO_BASE_URL")

    if not refresh_token:
        print(
            "ERROR: CIRRO_REFRESH_TOKEN environment variable is required",
            file=sys.stderr,
        )
        sys.exit(1)

    if not base_url:
        print("ERROR: CIRRO_BASE_URL environment variable is required", file=sys.stderr)
        sys.exit(1)

    # Optional: process identification
    process_id = os.environ.get("CIRRO_PROCESS_ID")
    process_name = os.environ.get("CIRRO_PROCESS_NAME")

    # Get Cirro system info (includes Cognito config)
    print(f"Connecting to Cirro at {base_url}...")
    system_info = get_cirro_system_info(base_url)

    client_id = system_info["auth"]["sdkAppId"]
    region = system_info["region"]

    # Get access token using refresh token
    print("Refreshing access token...")
    access_token = refresh_access_token(refresh_token, client_id, region)

    # Create authenticated Cirro client
    auth = AccessTokenAuth(token=access_token)
    cirro = CirroApi(auth_info=auth, base_url=base_url)

    print(f"Authenticated as: {auth.get_current_user()}")

    # Find the process to sync
    if process_id:
        print(f"Using process ID: {process_id}")
    elif process_name:
        print(f"Looking up process by name: {process_name}")
        process = find_process_by_name(cirro, process_name)
        if not process:
            print(f"ERROR: Process '{process_name}' not found", file=sys.stderr)
            print("\nAvailable processes:", file=sys.stderr)
            for p in cirro.processes.list():
                print(f"  - {p.name} (ID: {p.id})", file=sys.stderr)
            sys.exit(1)
        process_id = process.id
        print(f"Found process ID: {process_id}")
    else:
        print(
            "ERROR: Either CIRRO_PROCESS_ID or CIRRO_PROCESS_NAME must be set",
            file=sys.stderr,
        )
        sys.exit(1)

    # Sync the pipeline
    print(f"Syncing pipeline {process_id}...")
    result = cirro.processes.sync_custom_process(process_id)

    print("Sync completed successfully!")
    print(f"  Repository: {result.repository}")
    print(f"  Branch: {result.branch}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
