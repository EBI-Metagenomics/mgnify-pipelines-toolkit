#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Custom Exceptions
class WebinCredentialsError(Exception):
    """Raised when Webin credentials are missing or invalid."""

    pass


class ManifestValidationError(Exception):
    """Raised when manifest file validation fails."""

    pass


class WebinCLINotFoundError(Exception):
    """Raised when webin-cli executable cannot be found."""

    pass


class WebinCLIExecutionError(Exception):
    """Raised when webin-cli execution fails after all retries."""

    pass


class DownloadError(Exception):
    """Raised when downloading webin-cli fails."""

    pass


ENA_WEBIN = "ENA_WEBIN"
ENA_WEBIN_PASSWORD = "ENA_WEBIN_PASSWORD"
REPORT_FILE = "webin-cli.report"
WEBIN_SUBMISSION_RESULT_FILES = ["analysis.xml", "receipt.xml", "submission.xml", "webin-submission.xml"]
RETRIES = 3
RETRY_DELAY = 5
INSDC_CENTRE_PREFIXES = "EDS"
ENA_ASSEMBLY_ACCESSION_REGEX = f"([{INSDC_CENTRE_PREFIXES}]RZ[0-9]{{6,}})"


logger = logging.getLogger(__name__)


def configure_logging(debug: bool = False) -> None:
    """Configure process-wide logging level and format."""
    logging.basicConfig(level=logging.DEBUG if debug else logging.INFO, format="[%(asctime)s] - %(levelname)s - %(message)s", force=True)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the Webin-CLI submission tool.

    Returns:
        argparse.Namespace: Parsed arguments object containing attributes:
            - manifest (str): Path to a manifest file or a directory of manifest files
            - context (str): Submission context: genome, transcriptome, etc.
            - mode (str): 'submit' or 'validate'
            - test (bool): Whether to use Webin test server
            - workdir (Optional[str]): Working directory for temporary output
            - download_webin_cli (bool): Whether to download Webin-CLI jar
            - download_webin_cli_directory (str): Path to store downloaded Webin-CLI
            - download_webin_cli_version (Optional[str]): Specific Webin-CLI version
            - webin_cli_jar (Optional[str]): Path to pre-downloaded jar file
            - retries (int): Number of retry attempts
            - retry_delay (int): Initial retry delay in seconds
            - java_heap_size_initial (Optional[int]): Java initial heap size in GB
            - java_heap_size_max (Optional[int]): Java maximum heap size in GB
            - debug (bool): Enable verbose debug logging
    """

    def positive_int(value: str) -> int:
        """Argparse type: integer >= 1."""
        parsed_value = int(value)
        if parsed_value < 1:
            raise argparse.ArgumentTypeError("Value must be >= 1")
        return parsed_value

    def non_negative_int(value: str) -> int:
        """Argparse type: integer >= 0."""
        parsed_value = int(value)
        if parsed_value < 0:
            raise argparse.ArgumentTypeError("Value must be >= 0")
        return parsed_value

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--manifest",
        required=True,
        type=str,
        help="Path to a single manifest file or a directory containing manifest files",
    )
    parser.add_argument(
        "-o",
        "--output-accessions",
        required=False,
        default="ena_accessions.tsv",
        type=str,
        help="File to write assigned accessions to (TSV, default: ena_accessions.tsv)",
    )
    parser.add_argument(
        "-c",
        "--context",
        required=True,
        type=str,
        help="Submission type: genome, transcriptome, sequence, polysample, reads, taxrefset",
        choices=["genome", "transcriptome", "sequence", "polysample", "reads", "taxrefset"],
    )
    parser.add_argument("--mode", required=True, type=str, help="submit or validate", choices=["submit", "validate"])
    parser.add_argument("--test", required=False, action="store_true", help="Specify to use test server instead of live")
    parser.add_argument("--workdir", required=False, help="Path to working directory")
    parser.add_argument("--download-webin-cli", required=False, action="store_true", help="Specify if you do not have ena-webin-cli installed")
    parser.add_argument("--download-webin-cli-directory", required=False, default=".", type=str, help="Path to save webin-cli into")
    parser.add_argument("--download-webin-cli-version", required=False, type=str, help="Version of ena-webin-cli to download, default: latest")
    parser.add_argument("--webin-cli-jar", required=False, type=str, help="Path to pre-downloaded webin-cli.jar file to execute")
    parser.add_argument(
        "--retries",
        required=False,
        type=positive_int,
        default=RETRIES,
        help=f"Number of retry attempts (must be >= 1, default: {RETRIES})",
    )
    parser.add_argument(
        "--retry-delay",
        required=False,
        type=non_negative_int,
        default=RETRY_DELAY,
        help=f"Initial retry delay in seconds (must be >= 0, default: {RETRY_DELAY})",
    )
    parser.add_argument(
        "--java-heap-size-initial",
        required=False,
        type=positive_int,
        default=None,
        help="Java initial heap size in GB (-Xms); only added when explicitly provided",
    )
    parser.add_argument(
        "--java-heap-size-max",
        required=False,
        type=positive_int,
        default=None,
        help="Java maximum heap size in GB (-Xmx); only added when explicitly provided",
    )
    parser.add_argument("--debug", required=False, action="store_true", help="Enable debug logging")
    return parser.parse_args()


def get_latest_webin_cli_version() -> str:
    """
    Fetch the latest release version tag of ena-webin-cli from GitHub API.

    Returns:
        str: Latest version string (e.g. '9.0.1')

    Raises:
        DownloadError: If fetching the version fails.
    """
    api_url = "https://api.github.com/repos/enasequence/webin-cli/releases/latest"
    try:
        with urllib.request.urlopen(api_url, timeout=10) as response:
            data = json.load(response)
            return data["tag_name"]
    except (urllib.error.URLError, KeyError, json.JSONDecodeError) as e:
        logger.error(f"Failed to fetch latest webin-cli version: {e}")
        raise DownloadError(f"Failed to fetch latest webin-cli version: {e}")


def download_webin_cli(version: Optional[str] = None, dest_dir: str = "") -> None:
    """
    Download the Webin-CLI `.jar` file.

    Args:
        version (str, optional): Version to download. If None, fetch latest.
        dest_dir (str): Directory to save the jar file.

    Raises:
        DownloadError: If download fails.
    """
    if version is None:
        version = get_latest_webin_cli_version()

    jar_url = f"https://github.com/enasequence/webin-cli/releases/download/{version}/webin-cli-{version}.jar"
    dest_path = Path(dest_dir) / "webin-cli.jar"

    try:
        logger.info(f"Downloading webin-cli.jar version {version} to {dest_path} ...")
        urllib.request.urlretrieve(jar_url, str(dest_path))
        logger.info(f"Successfully downloaded webin-cli.jar to {dest_path}")
    except (urllib.error.URLError, OSError) as e:
        logger.error(f"Failed to download webin-cli.jar: {e}")
        raise DownloadError(f"Failed to download webin-cli.jar: {e}")


def ensure_webin_credentials_exist() -> None:
    """
    Verify that required Webin credentials are present in environment variables.

    Raises:
        WebinCredentialsError: If credentials are missing.
    """
    logger.debug("Checking required Webin credential environment variables")
    if ENA_WEBIN not in os.environ:
        raise WebinCredentialsError(f"The variable {ENA_WEBIN} is missing from the environment.")
    if ENA_WEBIN_PASSWORD not in os.environ:
        raise WebinCredentialsError(f"The variable {ENA_WEBIN_PASSWORD} is missing from the environment.")


def get_webin_credentials() -> str:
    """
    Retrieve Webin credentials from environment variables.

    Returns:
        str: Webin username.

    Raises:
        WebinCredentialsError: If credentials are missing.
    """
    ensure_webin_credentials_exist()
    webin = os.environ.get(ENA_WEBIN)
    password = os.environ.get(ENA_WEBIN_PASSWORD)
    if not password:
        raise WebinCredentialsError(f"The variable {ENA_WEBIN_PASSWORD} is missing password.")
    logger.debug(f"Webin credentials retrieved from environment: {webin}")
    return webin


def parse_manifest(manifest_path: Path) -> dict[str, str]:
    """
    Parse manifest file to extract ASSEMBLYNAME, FASTA file path, and other fields into dict.
    Example:
    STUDY	PRJ1
    SAMPLE	SAMEA7687881
    RUN_REF	ERR4918394
    ASSEMBLYNAME	ERR4918394_d41d8cd98f00
    ASSEMBLY_TYPE	primary metagenome
    COVERAGE	20.0
    PROGRAM	metaSPADES v3.12.1
    PLATFORM	DNBSEQ-G400
    FASTA	tests/fixtures/ERR4918394.fasta.gz
    TPA	true
    Args:
        manifest_path (Path): Path to the manifest file.
    Returns:
        Dict[str, str]:
            - keys: field names from the manifest
            - values: corresponding values from the manifest
    """
    manifest_dict = {}
    with open(manifest_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            parts = line.split(None, 1)  # split on any whitespace
            if len(parts) != 2:
                continue  # skip malformed lines
            field, value = parts
            manifest_dict[field] = value
    return manifest_dict


def check_manifest(manifest_file: str, workdir: Optional[str]) -> Tuple[str, str]:
    """
    Validate a manifest file and create an update one with absolute paths if needed.
    Args:
        manifest_file (str): Path to the original manifest.
        workdir (str or None): Optional working directory to prefix FASTA paths.
    Returns:
        Tuple[str, str]:
            - assembly_name (str): Parsed value of ASSEMBLYNAME.
            - manifest_for_submission (str): Path to the (possibly updated) manifest file.
    Raises:
        ManifestValidationError: If manifest validation fails or FASTA file not found.
    """
    manifest_path = Path(manifest_file)
    logger.debug(f"Checking manifest file {manifest_path} with workdir={workdir}")
    if not manifest_path.exists():
        logger.error(f"Manifest file does not exist: {manifest_path}")
        raise ManifestValidationError(f"Manifest file does not exist: {manifest_path}")
    manifest_dict = parse_manifest(manifest_path)
    assembly_name = manifest_dict.get("ASSEMBLYNAME")
    fasta_path = Path(manifest_dict.get("FASTA")) if manifest_dict.get("FASTA") else None
    if not assembly_name:
        logger.error(f"ASSEMBLYNAME field not found in {manifest_path}")
        raise ManifestValidationError(f"ASSEMBLYNAME field not found in {manifest_path}")
    if not fasta_path:
        logger.error(f"FASTA field not found in {manifest_path}")
        raise ManifestValidationError(f"FASTA field not found in {manifest_path}")
    # Validate FASTA file exists (if absolute path or workdir provided)
    if fasta_path.is_absolute():
        logger.debug(f"FASTA path is absolute: {fasta_path}")
        if not fasta_path.exists():
            logger.error(f"FASTA file not found: {fasta_path}")
            raise ManifestValidationError(f"FASTA file not found: {fasta_path}")
    elif workdir:
        # Check if FASTA exists relative to workdir
        full_fasta_path = Path(workdir) / fasta_path
        logger.debug(f"Resolved FASTA path using workdir: {full_fasta_path}")
        if not full_fasta_path.exists():
            logger.error(f"FASTA file not found: {full_fasta_path}")
            raise ManifestValidationError(f"FASTA file not found: {full_fasta_path}")

        manifest_dict["FASTA"] = str(full_fasta_path)
        updated_manifest_path = manifest_path.parent / f"updated_{manifest_path.name}"
        with open(updated_manifest_path, "w") as f:
            for key, value in manifest_dict.items():
                f.write(f"{key}\t{value}\n")
        logger.info(f"New manifest {updated_manifest_path} was created with absolute paths")
        return assembly_name, str(updated_manifest_path)
    logger.debug(f"Manifest check complete without update: {manifest_file}")
    return assembly_name, manifest_file


def get_webin_cli_command(
    manifest: str,
    context: str,
    webin: str,
    mode: str,
    test: bool,
    output_dir: str,
    jar: Optional[str] = None,
    java_heap_size_initial: Optional[int] = None,
    java_heap_size_max: Optional[int] = None,
) -> List[str]:
    """
    Build the webin-cli command list based on execution method.

    Args:
        manifest (str): Path to manifest file.
        context (str): Submission context (genome, assembly, etc.).
        webin (str): ENA Webin username.
        password (str): ENA Webin password.
        mode (str): Execution mode ('submit' or 'validate').
        test (bool): Whether to use test server.
        output_dir (str): Directory where webin-cli writes output files.
        jar (Optional[str]): Path to webin-cli jar file, if using jar execution.
        java_heap_size_initial (Optional[int]): Java initial heap size in GB (-Xms).
        java_heap_size_max (Optional[int]): Java maximum heap size in GB (-Xmx).

    Returns:
        List[str]: Complete command list ready for execution.

    Raises:
        WebinCLINotFoundError: If neither jar is provided nor ena-webin-cli is installed.
    """
    if jar:
        # Use jar file execution
        logger.info(f"Using java execution with jar: {jar}")
        logback_config = Path(__file__).parent / "webincli_logback.xml"
        cmd = ["java"]
        if logback_config.exists():
            cmd.append(f"-Dlogback.configurationFile={logback_config}")
        else:
            logger.warning(f"Logback config not found at {logback_config}, continuing without it")
        if java_heap_size_initial is not None:
            cmd.append(f"-Xms{java_heap_size_initial}g")
        if java_heap_size_max is not None:
            cmd.append(f"-Xmx{java_heap_size_max}g")

        cmd += ["-jar", jar]
    else:
        # Use mamba/conda installation
        webin_cli_path = shutil.which("ena-webin-cli")
        if webin_cli_path:
            logger.info(f"Using ena-webin-cli from: {webin_cli_path}")
            cmd = ["ena-webin-cli"]
        else:
            logger.error("ena-webin-cli was not found in PATH. Install it with mamba/conda or use --download-webin-cli / --webin-cli-jar options.")
            raise WebinCLINotFoundError(
                "ena-webin-cli was not found in PATH. Install it with mamba/conda or use --download-webin-cli / --webin-cli-jar options."
            )

    # Add webin-cli arguments
    cmd += [
        f"-context={context}",
        f"-manifest={manifest}",
        f"-userName={webin}",
        f"-passwordEnv={ENA_WEBIN_PASSWORD}",
        f"-outputDir={output_dir}",
        f"-{mode}",
    ]

    if test:
        cmd.append("-test")

    return cmd


def prepare_output_dir(manifest: str, assembly_name: str) -> str:
    """Prepare webin-cli output directory named after assembly name."""
    output_dir = Path(manifest).parent / assembly_name
    if output_dir.exists():
        if output_dir.is_dir():
            logger.warning(f"Output directory already exists and will be reused: {output_dir}")
        else:
            raise ManifestValidationError(f"Output path exists and is not a directory: {output_dir}")
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created output directory for webin-cli: {output_dir}")
    return str(output_dir)


def run_webin_cli(
    manifest: str,
    context: str,
    webin: str,
    mode: str,
    test: bool,
    output_dir: str,
    jar: Optional[str] = None,
    retries: int = RETRIES,
    retry_delay: int = RETRY_DELAY,
    java_heap_size_initial: Optional[int] = None,
    java_heap_size_max: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """
    Execute webin-cli with retry logic and error handling.

    Args:
        manifest (str): Path to manifest file.
        context (str): Submission context (genome, assembly, etc.).
        webin (str): ENA Webin username.
        password (str): ENA Webin password.
        mode (str): Execution mode ('submit' or 'validate').
        test (bool): Whether to use test server.
        output_dir (str): Directory where webin-cli writes output files.
        jar (Optional[str]): Path to webin-cli jar file, if using jar execution.
        retries (int): Number of retry attempts.
        retry_delay (int): Initial retry delay in seconds.
        java_heap_size_initial (Optional[int]): Java initial heap size in GB (-Xms).
        java_heap_size_max (Optional[int]): Java maximum heap size in GB (-Xmx).

    Returns:
        subprocess.CompletedProcess: Result of successful webin-cli execution.

    Raises:
        WebinCredentialsError: If credentials are invalid.
        WebinCLIExecutionError: If all retry attempts fail.
    """

    def get_combined_log(run_result: subprocess.CompletedProcess) -> str:
        """Combine stdout and stderr from the subprocess result for logging."""
        output_parts = []
        if run_result.stdout and run_result.stdout.strip():
            output_parts.append(run_result.stdout.strip())
        if run_result.stderr and run_result.stderr.strip():
            output_parts.append(run_result.stderr.strip())
        if not output_parts:
            return "webin-cli produced no stdout/stderr output"
        return "\n".join(output_parts)

    def redact_command(command: List[str]) -> str:
        """Return a shell-quoted command string with credential values replaced by ***."""
        redacted = []
        for token in command:
            if token.startswith("-userName=") or token.startswith("-passwordEnv="):
                key, _ = token.split("=", 1)
                redacted.append(f"{key}=***")
            else:
                redacted.append(token)
        return " ".join(redacted)

    cmd = get_webin_cli_command(manifest, context, webin, mode, test, output_dir, jar, java_heap_size_initial, java_heap_size_max)
    logger.info(f"Command: {redact_command(cmd)}")

    for attempt in range(1, retries + 1):
        logger.info(f"🚀 Running ena-webin-cli (attempt {attempt}/{retries})...")

        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            env=os.environ.copy(),
        )
        combined_log = get_combined_log(result)
        logger.debug(
            f"webin-cli process finished with return code {result.returncode}, stdout_len={len(result.stdout or '')}, stderr_len={len(result.stderr or '')}"
        )
        logger.debug("ena-webin-cli output:\n" + combined_log)

        if result.returncode == 0:
            logger.info("✅ Command completed successfully")
            return result  # success

        logger.warning(f"❌ ena-webin-cli exited with non-zero code {result.returncode}")

        if "Invalid submission account user name or password." in combined_log:
            logger.error("💥 Invalid credentials for Webin account - not retrying")
            raise WebinCredentialsError("Invalid credentials for Webin account")

        if "The object being added already exists in the submission account with accession:" in combined_log:
            logger.info("📦 Submission already exists")
            logger.info("✅ Command completed successfully")
            return result

        logger.error("ena-webin-cli output:\n" + combined_log)

        if attempt < retries:
            sleep_time = retry_delay * (2 ** (attempt - 1))  # exponential backoff
            logger.warning(f"🔁 Retrying in {sleep_time} seconds...")
            logger.debug(f"Retry backoff details: attempt={attempt}, max_retries={retries}, sleep_time={sleep_time}")
            time.sleep(sleep_time)
        else:
            logger.error("💥 All retries failed.")
            raise WebinCLIExecutionError(f"Webin-CLI execution failed after {retries} attempts. See logs for details.")


def check_submission_status_test(report_text: str) -> Tuple[bool, bool]:
    """
    Check the test Webin-CLI submission report and determine submission status.

    On ENA's test server, when an object is resubmitted (already exists), webin-cli
    returns a minimal report containing only "This was a TEST submission(s)." with
    no additional details. For first-time submissions, the report contains this
    message plus additional submission details.

    Args:
        report_text (str): Contents of the webin-cli.report file.

    Returns:
        Tuple[bool, bool]: (success, is_resubmission)
            - success (bool): True if submission/validation succeeded, False otherwise.
            - is_resubmission (bool): True if object was previously submitted to test
              server (duplicate submission), False if this is first-time submission.

    Examples:
        First-time successful submission
        check_submission_status_test("This was a TEST submission(s).\\nsubmission has been completed successfully")
        (True, False)

        Resubmission (object already exists)
        check_submission_status_test("This was a TEST submission(s).")
        (True, True)

        Failed submission
        check_submission_status_test("Error: validation failed")
        (False, False)
    """
    logger.debug("Evaluating TEST server submission report content")
    if "This was a TEST submission(s)." not in report_text:
        logger.info("Submission failed on TEST server")
        return False, False

    # Check for resubmission: minimal report with only success message
    # For resubmissions, report contains only "This was a TEST submission(s)." without details
    is_minimal_report = report_text.strip() == "This was a TEST submission(s)." or (
        report_text.count("\n") <= 2 and "submission has been completed successfully" not in report_text
    )

    if is_minimal_report:
        logger.info("Submitted object already exists on TEST server")
        return True, True
    elif "submission has been completed successfully" in report_text:
        logger.info("Successfully submitted object for the first time on TEST server")
        return True, False
    else:
        logger.info("Submission failed on TEST server")
        return False, False


def check_submission_status_live(report_text: str) -> Tuple[bool, bool]:
    """
    Check the live Webin-CLI submission report and determine submission status.

    On ENA's live server, the webin-cli report explicitly indicates when an object
    being submitted already exists in the account. Unlike the test server, the live
    server provides clear messaging in both first-time and duplicate submissions.

    Args:
        report_text (str): Contents of the webin-cli.report file.

    Returns:
        Tuple[bool, bool]: (success, is_resubmission)
            - success (bool): True if submission succeeded (either new or duplicate),
              False if submission failed.
            - is_resubmission (bool): True if object was previously submitted to live
              server (duplicate), False if this is first-time submission.

    Examples:
        First-time successful submission
        check_submission_status_live("submission has been completed successfully")
        (True, False)

        Resubmission (object already exists)
        check_submission_status_live("object being added already exists in the submission account with accession ERZ123456")
        (True, True)

        Failed submission
        check_submission_status_live("Error: validation failed")
        (False, False)
    """
    is_resubmission = False
    logger.debug("Evaluating MAIN server submission report content")
    if "submission has been completed successfully" in report_text:
        logger.info("Successfully submitted object on MAIN server")
        return True, is_resubmission
    elif "object being added already exists in the submission account with accession" in report_text:
        logger.info("Submitted object already exists on MAIN server")
        is_resubmission = True
        return True, is_resubmission
    else:
        logger.info("Submission failed on MAIN server")
        return False, is_resubmission


def check_report(fasta_location: str, assembly_name: str, mode: str, test: bool) -> Tuple[bool, bool]:
    """
    Parse and validate the webin-cli report file to determine operation outcome.

    This function reads the webin-cli.report file generated after a submission or
    validation attempt, extracts any assigned accessions, and determines whether
    the operation succeeded and whether this was a resubmission of existing data.

    Args:
        fasta_location (str): Directory path where webin-cli.report is located
            (typically the directory containing the FASTA file).
        assembly_name (str): Assembly identifier
        mode (str): Operation mode - either 'submit' or 'validate'.
        test (bool): Whether operation was run against test server (True) or
            live server (False).

    Returns:
        Tuple[bool, bool]: (success, is_resubmission)
            - success (bool): True if the operation (submit/validate) succeeded,
              False otherwise.
            - is_resubmission (bool): True if this was a duplicate submission
              (object already existed), False if first-time submission. Always
              False for validation mode.

    Note:
        For submit mode, the behavior differs between test and live servers.
        See check_submission_status_test() and check_submission_status_live()
        for detailed explanations of server-specific behaviors.
    """
    report_path = Path(fasta_location) / REPORT_FILE
    logger.info(f"Checking webin report {report_path}")

    if not report_path.exists():
        logger.warning(f"⚠️ Report file not found: {report_path}")
        return False, False

    report_text = report_path.read_text()
    logger.debug(f"Report content:\n{report_text}")
    assembly_accession_match_list = re.findall(ENA_ASSEMBLY_ACCESSION_REGEX, report_text)
    if assembly_accession_match_list:
        logger.info(f"Assigned accessions: {','.join(assembly_accession_match_list)}")

    if mode == "submit":
        if test:
            return check_submission_status_test(report_text)
        else:
            return check_submission_status_live(report_text)

    elif mode == "validate":
        if "Submission(s) validated successfully." in report_text:
            logger.info("Submission validation succeeded")
            return True, False
        else:
            logger.info("Submission validation failed")
            return False, False


def check_result(result_location: str, context: str, assembly_name: str, mode: str, test: bool) -> bool:
    """
    Verify that webin-cli execution completed successfully and produced expected outputs.

    For first-time submissions (not resubmissions), this function validates that
    webin-cli created the expected directory structure and output files.

    Args:
        result_location (str): Directory where webin-cli was executed (contains report).
        context (str): Submission context (genome, assembly, etc.).
        assembly_name (str): Name of the assembly from manifest ASSEMBLYNAME field.
        mode (str): Operation mode - 'submit' or 'validate'.
        test (bool): Whether operation was against test server.

    Returns:
        bool: True if operation succeeded and all expected outputs exist,
              False if operation failed or outputs are missing.

    Note:
        For resubmissions (when object already exists), we only check the report
        status and don't verify output files, since webin-cli may not regenerate them.
    """
    success, is_resubmission = check_report(result_location, assembly_name, mode, test)
    logger.debug(f"Result status after report parsing: success={success}, is_resubmission={is_resubmission}, mode={mode}, context={context}")

    if not success:
        logger.info("Command failed. Check logs")
        return False

    # For validation mode or resubmissions, report status is sufficient
    if mode == "validate" or is_resubmission:
        return True

    # For first-time submissions, verify output directory structure and files exist
    result_path = Path(result_location)
    assembly_dir = result_path / context / assembly_name

    if not assembly_dir.exists():
        logger.info(f"No result directory found for {assembly_name} at {assembly_dir}")
        return False

    submission_dir = assembly_dir / mode
    if not submission_dir.exists():
        logger.info(f"No submission result directory found at {submission_dir}")
        return False

    # Verify all expected result files are present
    missing_files = []
    for filename in WEBIN_SUBMISSION_RESULT_FILES:
        file_path = submission_dir / filename
        if not file_path.exists():
            missing_files.append(filename)

    if missing_files:
        logger.info(f"Missing required files in {submission_dir}: {', '.join(missing_files)}")
        return False

    logger.info(f"Result directory {submission_dir} contains all expected files")
    return True


def resolve_manifests(manifest_input: str) -> List[str]:
    """
    Resolve a manifest argument to a list of manifest file paths.

    Accepts either a single manifest file or a directory containing manifest files
    (matched by the ``*.manifest`` glob pattern).

    Args:
        manifest_input (str): Path to a manifest file or a directory of manifests.

    Returns:
        List[str]: Sorted list of manifest file paths.

    Raises:
        ManifestValidationError: If the input is neither a file nor a directory,
            or if a directory contains no ``*.manifest`` files.
    """
    path = Path(manifest_input)
    logger.debug(f"Resolving manifest input: {manifest_input}")
    if path.is_file():
        return [str(path)]
    elif path.is_dir():
        manifests = sorted(path.glob("*.manifest"))
        if not manifests:
            raise ManifestValidationError(f"No .manifest files found in directory: {path}")
        return [str(m) for m in manifests]
    else:
        raise ManifestValidationError(f"Invalid manifest input '{manifest_input}': must be a single manifest file or a directory of manifest files.")


def write_assigned_accessions(accessions: Dict[str, str], output_accessions: str) -> None:
    """
    Write assigned assembly accessions for all successful submissions to a TSV file.

    Args:
        accessions (Dict[str, str]): Mapping of assembly name to accession.
        output_accessions (str): Path to output TSV file (overwritten each run).
    """
    logger.debug(f"Preparing to write accessions table to {output_accessions}")
    if not accessions:
        logger.warning("No accessions to write.")
        return
    output_file = Path(output_accessions)
    with open(output_file, "w", newline="") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t")
        writer.writerow(["alias", "accession"])
        for assembly_name, accession in accessions.items():
            writer.writerow([assembly_name, accession])
    logger.info(f"Assigned accessions written to {output_file}")


def main() -> int:
    """
    Main entry point for the webin-cli handler.

    Returns:
        int: Exit code (0 for success, non-zero for failure).
    """
    try:
        args = parse_arguments()
        configure_logging(debug=args.debug)
        logger.debug(
            f"Parsed arguments: manifest={args.manifest}, context={args.context}, mode={args.mode}, test={args.test}, "
            f"workdir={args.workdir}, download_webin_cli={args.download_webin_cli}, "
            f"download_webin_cli_directory={args.download_webin_cli_directory}, download_webin_cli_version={args.download_webin_cli_version}, "
            f"webin_cli_jar={args.webin_cli_jar}, retries={args.retries}, retry_delay={args.retry_delay}, "
            f"java_heap_size_initial={args.java_heap_size_initial}, java_heap_size_max={args.java_heap_size_max}, debug={args.debug}"
        )

        if args.download_webin_cli:
            download_webin_cli(version=args.download_webin_cli_version, dest_dir=args.download_webin_cli_directory)

        webin = get_webin_credentials()

        # pick execution method for webin-cli
        execute_jar = None
        if args.download_webin_cli:
            execute_jar = str(Path(args.download_webin_cli_directory) / "webin-cli.jar")
        elif args.webin_cli_jar:
            execute_jar = args.webin_cli_jar

        manifests = resolve_manifests(args.manifest)
        logger.info(f"Resolved {len(manifests)} manifest(s) for processing")

        failed_manifests = []
        accessions_to_write: Dict[str, str] = {}
        for manifest in manifests:
            logger.debug(f"Starting submission loop for manifest: {manifest}")
            assembly_name, manifest_for_submission = check_manifest(manifest, args.workdir)
            output_dir = prepare_output_dir(manifest_for_submission, assembly_name)
            run_webin_cli(
                manifest=manifest_for_submission,
                context=args.context,
                webin=webin,
                mode=args.mode,
                test=args.test,
                output_dir=output_dir,
                jar=execute_jar,
                retries=args.retries,
                retry_delay=args.retry_delay,
                java_heap_size_initial=args.java_heap_size_initial,
                java_heap_size_max=args.java_heap_size_max,
            )

            logger.debug(f"Result location for manifest {manifest_for_submission}: {output_dir}")
            result_status = check_result(
                result_location=output_dir, context=args.context, assembly_name=assembly_name, mode=args.mode, test=args.test
            )

            if result_status:
                if args.mode == "submit":
                    report_text = (Path(output_dir) / REPORT_FILE).read_text()
                    found = re.findall(ENA_ASSEMBLY_ACCESSION_REGEX, report_text)
                    logger.debug(f"Accession matches for assembly {assembly_name}: {found}")
                    if found and len(found) == 1:
                        accessions_to_write[assembly_name] = found[0]
                    else:
                        logger.warning(f"No accession found in {REPORT_FILE} for {assembly_name}, skipping.")
                logger.info(f"Submission/validation done for {manifest_for_submission}")
            else:
                logger.error(f"Submission/validation failed for {manifest_for_submission}")
                failed_manifests.append(manifest_for_submission)

        if args.mode == "submit":
            write_assigned_accessions(accessions_to_write, args.output_accessions)
            logger.debug(f"Submission mode complete with {len(accessions_to_write)} accession(s) captured")

        if failed_manifests:
            logger.error(f"There are failed submissions for manifests: {','.join(failed_manifests)}")
            return 1
        return 0

    except WebinCredentialsError as e:
        logger.error(f"Credentials error: {e}")
        return 2
    except ManifestValidationError as e:
        logger.error(f"Manifest validation error: {e}")
        return 3
    except WebinCLINotFoundError as e:
        logger.error(f"Webin-CLI not found: {e}")
        return 4
    except WebinCLIExecutionError as e:
        logger.error(f"Webin-CLI execution error: {e}")
        return 5
    except DownloadError as e:
        logger.error(f"Download error: {e}")
        return 6
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 99


if __name__ == "__main__":
    sys.exit(main())
