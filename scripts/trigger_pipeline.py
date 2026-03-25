#!/usr/bin/env python3
"""
Trigger and monitor Cirro pipeline runs.

This script supports two subcommands:
    trigger: Start a new pipeline run and output the dataset ID/URL
    wait: Wait for an existing pipeline run to complete

Authentication modes:
1. CI mode: Uses CIRRO_REFRESH_TOKEN (from GitHub secrets) to authenticate
2. Local mode: Uses the SDK's default auth from ~/.cirro/config.ini

Required environment variables:
    CIRRO_BASE_URL: The Cirro instance URL (e.g., breakthroughcancer.cirro.bio)
    CIRRO_PROJECT_ID: The project ID containing the dataset

For 'trigger' command:
    CIRRO_PROCESS_ID: The process ID to run
    CIRRO_DATASET_ID: The dataset ID to use as input

For 'wait' command:
    CIRRO_OUTPUT_DATASET_ID: The dataset ID to monitor (from trigger output)

Optional environment variables:
    CIRRO_REFRESH_TOKEN: The refresh token (if not set, uses local SDK auth)
    CIRRO_RUN_NAME: Custom name for the pipeline run (trigger only)
    CIRRO_POLL_INTERVAL: Poll interval in seconds (default: 60)
    CIRRO_TIMEOUT: Timeout in seconds (default: 7200 = 2 hours)
"""

import argparse
import os
import sys
import time
from datetime import datetime

import boto3
import requests
from botocore.exceptions import ClientError
from cirro import CirroApi
from cirro.auth.access_token import AccessTokenAuth
from cirro_api_client.v1.errors import UnexpectedStatus
from cirro_api_client.v1.models import RunAnalysisRequest, RunAnalysisRequestParams
from cirro_api_client.v1.models.status import Status


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


def get_dashboard_url(base_url: str, project_id: str, dataset_id: str) -> str:
    """Generate the Cirro dashboard URL for a dataset."""
    return f"https://{base_url}/project/{project_id}/dataset/{dataset_id}"


def get_cirro_client(
    base_url: str, refresh_token: str | None
) -> tuple[CirroApi, str | None, str | None]:
    """
    Create a Cirro API client.

    Returns:
        tuple of (cirro_client, client_id, region)
        client_id and region are only set if using refresh token auth
    """
    client_id = None
    region = None

    if refresh_token:
        # CI mode: use refresh token authentication
        system_info = get_cirro_system_info(base_url)
        client_id = system_info["auth"]["sdkAppId"]
        region = system_info["region"]

        print("Refreshing access token...")
        access_token = refresh_access_token(refresh_token, client_id, region)

        auth = AccessTokenAuth(token=access_token)
        cirro = CirroApi(auth_info=auth, base_url=base_url)
        print(f"Authenticated as: {auth.get_current_user()}")
    else:
        # Local mode: use SDK default auth from ~/.cirro/config.ini
        print("Using local SDK authentication...")
        cirro = CirroApi(base_url=base_url)

    return cirro, client_id, region


def wait_for_completion(
    cirro: CirroApi,
    project_id: str,
    dataset_id: str,
    poll_interval: int = 60,
    timeout: int = 7200,
    refresh_token: str | None = None,
    client_id: str | None = None,
    region: str | None = None,
    base_url: str | None = None,
) -> tuple[Status, str]:
    """
    Poll the dataset status until it reaches a terminal state.

    If refresh_token, client_id, region, and base_url are provided,
    the function will refresh the access token when it expires (401 error).

    Returns:
        tuple of (final_status, status_message)
    """
    terminal_states = {Status.COMPLETED, Status.FAILED, Status.DELETED, Status.ARCHIVED}
    start_time = time.time()

    while True:
        elapsed = time.time() - start_time
        if elapsed > timeout:
            return Status.UNKNOWN, f"Timeout after {timeout} seconds"

        try:
            dataset = cirro.datasets.get(project_id, dataset_id)
        except UnexpectedStatus as e:
            if (
                e.status_code == 401
                and refresh_token
                and client_id
                and region
                and base_url
            ):
                print("  Access token expired, refreshing...")
                access_token = refresh_access_token(refresh_token, client_id, region)
                auth = AccessTokenAuth(token=access_token)
                cirro = CirroApi(auth_info=auth, base_url=base_url)
                print("  Token refreshed, retrying...")
                continue
            raise

        if not dataset:
            return Status.UNKNOWN, f"Dataset {dataset_id} not found"

        status = dataset.status
        print(f"  Status: {status.value} (elapsed: {int(elapsed)}s)")

        if status in terminal_states:
            if status == Status.COMPLETED:
                return status, "Pipeline completed successfully"
            elif status == Status.FAILED:
                return status, "Pipeline failed"
            else:
                return status, f"Pipeline ended with status: {status.value}"

        time.sleep(poll_interval)


def set_github_output(name: str, value: str) -> None:
    """Set a GitHub Actions output variable."""
    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a") as f:
            f.write(f"{name}={value}\n")


def cmd_trigger(args: argparse.Namespace) -> int:
    """Trigger a new pipeline run."""
    refresh_token = os.environ.get("CIRRO_REFRESH_TOKEN")
    base_url = os.environ.get("CIRRO_BASE_URL", "breakthroughcancer.cirro.bio")
    process_id = os.environ.get("CIRRO_PROCESS_ID", "quinnj2-alarmist")
    dataset_id = os.environ.get("CIRRO_DATASET_ID")
    project_id = os.environ.get("CIRRO_PROJECT_ID")
    run_name = os.environ.get(
        "CIRRO_RUN_NAME", f"CI-triggered-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    )

    if not dataset_id:
        print("ERROR: CIRRO_DATASET_ID is required", file=sys.stderr)
        return 1

    if not project_id:
        print("ERROR: CIRRO_PROJECT_ID is required", file=sys.stderr)
        return 1

    print(f"Connecting to Cirro at {base_url}...")
    cirro, _, _ = get_cirro_client(base_url, refresh_token)

    # Verify dataset exists
    print(f"Fetching dataset info for: {dataset_id}")
    dataset = cirro.datasets.get(project_id, dataset_id)
    if not dataset:
        print(
            f"ERROR: Dataset {dataset_id} not found in project {project_id}",
            file=sys.stderr,
        )
        return 1
    print(f"Found dataset: {dataset.name}")

    # Trigger the pipeline
    print(f"Triggering pipeline {process_id}...")
    print(f"  Dataset: {dataset_id}")
    print(f"  Project: {project_id}")
    print(f"  Run name: {run_name}")

    # Build params
    params = RunAnalysisRequestParams()

    request = RunAnalysisRequest(
        name=run_name,
        description="Automated CI run triggered by GitHub Actions",
        process_id=process_id,
        source_dataset_ids=[dataset_id],
        params=params,
        notification_emails=[],
    )
    result = cirro.execution.run_analysis(project_id, request)

    print("Pipeline triggered successfully!")
    print(f"  Dataset ID: {result.id}")
    if result.message:
        print(f"  Message: {result.message}")

    # Generate dashboard URL for the new dataset (pipeline output)
    output_dataset_id = result.id
    dashboard_url = get_dashboard_url(base_url, project_id, output_dataset_id)
    print(f"  Dashboard URL: {dashboard_url}")

    # Set GitHub Actions outputs
    set_github_output("dataset_id", output_dataset_id)
    set_github_output("dashboard_url", dashboard_url)
    set_github_output("run_name", run_name)
    set_github_output("status", "TRIGGERED")
    set_github_output("status_message", "Pipeline triggered successfully")

    return 0


def cmd_wait(args: argparse.Namespace) -> int:
    """Wait for a pipeline run to complete."""
    refresh_token = os.environ.get("CIRRO_REFRESH_TOKEN")
    base_url = os.environ.get("CIRRO_BASE_URL", "breakthroughcancer.cirro.bio")
    project_id = os.environ.get("CIRRO_PROJECT_ID")
    output_dataset_id = os.environ.get("CIRRO_OUTPUT_DATASET_ID")
    poll_interval = int(os.environ.get("CIRRO_POLL_INTERVAL", "60"))
    timeout = int(os.environ.get("CIRRO_TIMEOUT", "7200"))

    if not output_dataset_id:
        print(
            "ERROR: CIRRO_OUTPUT_DATASET_ID is required for wait command",
            file=sys.stderr,
        )
        return 1

    if not project_id:
        print("ERROR: CIRRO_PROJECT_ID is required", file=sys.stderr)
        return 1

    print(f"Connecting to Cirro at {base_url}...")
    cirro, client_id, region = get_cirro_client(base_url, refresh_token)

    dashboard_url = get_dashboard_url(base_url, project_id, output_dataset_id)

    print(f"Waiting for dataset {output_dataset_id} to complete...")
    print(f"  Dashboard URL: {dashboard_url}")
    print(f"  Poll interval: {poll_interval}s, Timeout: {timeout}s")

    final_status, status_message = wait_for_completion(
        cirro,
        project_id,
        output_dataset_id,
        poll_interval,
        timeout,
        refresh_token=refresh_token,
        client_id=client_id,
        region=region,
        base_url=base_url,
    )

    print(f"\n{status_message}")
    print(f"  Final status: {final_status.value}")
    print(f"  Dashboard URL: {dashboard_url}")

    # Set final status outputs
    set_github_output("status", final_status.value)
    set_github_output("status_message", status_message)
    set_github_output("dashboard_url", dashboard_url)

    if final_status == Status.COMPLETED:
        return 0
    else:
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Trigger and monitor Cirro pipeline runs"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # trigger subcommand
    trigger_parser = subparsers.add_parser("trigger", help="Trigger a new pipeline run")
    trigger_parser.set_defaults(func=cmd_trigger)

    # wait subcommand
    wait_parser = subparsers.add_parser(
        "wait", help="Wait for a pipeline run to complete"
    )
    wait_parser.set_defaults(func=cmd_wait)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
