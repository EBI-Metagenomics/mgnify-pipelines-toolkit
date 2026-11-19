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


class WebinAbsConsentError(Exception):
    """Raised when the Webin account has not completed the Access and Benefit Sharing (ABS) acknowledgement."""

    pass


class WebinReportValidationError(Exception):
    """Raised when the webin-cli report is missing, does not confirm a successful submission or does not contain accession(s)."""

    pass


ENA_WEBIN = "ENA_WEBIN"
ENA_WEBIN_PASSWORD = "ENA_WEBIN_PASSWORD"
REPORT_FILE = "webin-cli.report"
RETRIES = 3
RETRY_DELAY = 5
INSDC_CENTRE_PREFIXES = "EDS"
ENA_ASSEMBLY_ACCESSION_REGEX = f"([{INSDC_CENTRE_PREFIXES}]RZ[0-9]{{6,}})"
ENA_RUN_ACCESSION_REGEX = f"([{INSDC_CENTRE_PREFIXES}]RR[0-9]{{6,}})"
ENA_EXPERIMENT_ACCESSION_REGEX = f"([{INSDC_CENTRE_PREFIXES}]RX[0-9]{{6,}})"

WEBIN_CLI_EXIT_SYSTEM_ERROR = 1  # server/network problems, receipt errors, ABS consent - may be transient
WEBIN_CLI_EXIT_USER_ERROR = 2  # invalid credentials, invalid manifest - retrying will not help
WEBIN_CLI_EXIT_VALIDATION_ERROR = 3  # data validation failed - retrying will not help
WEBIN_CLI_NON_RETRYABLE_EXIT_CODES = {WEBIN_CLI_EXIT_USER_ERROR, WEBIN_CLI_EXIT_VALIDATION_ERROR}

# webin-cli (>= 10.0.0) messages
MSG_VALIDATE_SUCCESS = "Submission(s) validated successfully."
MSG_SUBMIT_SUCCESS = "The submission has been completed successfully."
MSG_ALREADY_EXISTS = "The object being added already exists in the submission account with accession:"
MSG_INVALID_CREDENTIALS = "Invalid submission account user name or password."
MSG_ABS_CONSENT = "Incomplete Access and Benefit Sharing (ABS) acknowledgement"


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
            - resume (bool): Whether to resume previous run by skipping already processed aliases
            - test (bool): Whether to use Webin test server
            - input_dir (str | None): Path to the data files (FASTA, FASTQ or BAM) declared in the manifest file(s),
              required if manifest(s) reference file(s) that cannot be accessed from the execution directory
            - outdir (str | None): Base output directory for webin-cli files
            - output_accessions (str): File to write assigned accessions to (TSV, default: ena_accessions.tsv)
            - download_webin_cli (bool): Whether to download Webin-CLI jar
            - download_webin_cli_directory (str): Path to store downloaded Webin-CLI
            - download_webin_cli_version (str | None): Specific Webin-CLI version
            - webin_cli_jar (str | None): Path to pre-downloaded jar file
            - retries (int): Number of retry attempts
            - retry_delay (int): Initial retry delay in seconds
            - java_heap_size_initial (int | None): Java initial heap size in GB
            - java_heap_size_max (int | None): Java maximum heap size in GB
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
    parser.add_argument(
        "--input-dir",
        required=False,
        help="Path to the data files (FASTA, FASTQ or BAM) declared in the manifest file(s), required if manifest(s) reference file(s) that cannot be accessed from the execution directory",
    )
    parser.add_argument("--outdir", required=False, help="Base output directory for webin-cli files")
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
    parser.add_argument(
        "--resume",
        required=False,
        action="store_true",
        help="Resume previous webin_cli_handler run: skip aliases already present in accessions output file (submit mode only)",
    )
    args = parser.parse_args()

    if args.resume and args.mode != "submit":
        parser.error("--resume can only be used with --mode submit")

    if args.mode == "submit" and not args.resume and Path(args.output_accessions).exists():
        parser.error(
            f"Accessions file {args.output_accessions} already exists; use --resume to continue the previous run, "
            "or pass a different --output-accessions"
        )

    return args


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
        raise DownloadError(f"Failed to fetch latest webin-cli version: {e}")


def download_webin_cli(version: str | None = None, dest_dir: str = ".") -> str:
    """
    Download the Webin-CLI `.jar` file.

    Args:
        version (str, optional): Version to download. If None, fetch latest.
        dest_dir (str): Directory to save the jar file.

    Returns:
        str: Path to the downloaded jar file.

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
        raise DownloadError(f"Failed to download webin-cli.jar: {e}")
    return str(dest_path)


def get_webin_credentials() -> str:
    """
    Get the Webin username from the environment and check that the password is set.

    Returns:
        str: Webin username.

    Raises:
        WebinCredentialsError: If the username or password environment variable is missing or empty.
    """
    for variable in (ENA_WEBIN, ENA_WEBIN_PASSWORD):
        if not os.environ.get(variable):
            raise WebinCredentialsError(f"The variable {variable} is missing from the environment.")
    return os.environ[ENA_WEBIN]


def parse_manifest(manifest_path: Path) -> dict[str, str | list[str]]:
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
        dict[str, str | list[str]]:
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
            # For reads manifest, we can have multiple FASTQ lines, so we store them in a list
            if field == "FASTQ":
                manifest_dict.setdefault("FASTQ", []).append(value)
            else:
                manifest_dict[field] = value
    return manifest_dict


def check_manifest(manifest_file: str, input_dir: str | None, context: str) -> str:
    """
    Validate a manifest file and ensure referenced data file(s) are accessible.
    Args:
        manifest_file (str): Path to the original manifest.
        input_dir (str or None): Optional root directory for relative data file paths.
        context (str): Submission context (e.g. genome, reads).
    Returns:
        str: Alias (ASSEMBLYNAME or NAME) used for output folder naming and accessions tracking.
    Raises:
        ManifestValidationError: If manifest validation fails or data file(s) are not found.
    """
    manifest_path = Path(manifest_file)
    logger.debug(f"Checking manifest file {manifest_path} with input_dir={input_dir}, context={context}")
    if not manifest_path.exists():
        raise ManifestValidationError(f"Manifest file does not exist: {manifest_path}")
    manifest_dict = parse_manifest(manifest_path)
    alias = manifest_dict.get("ASSEMBLYNAME") or manifest_dict.get("NAME")
    if not alias:
        raise ManifestValidationError(f"Neither ASSEMBLYNAME, nor NAME field found in {manifest_path}")

    if context == "reads":
        # Reads can be submitted either as FASTQ file(s) or as a single BAM file (e.g. long reads), not both
        fastq_values = manifest_dict.get("FASTQ", [])
        bam_value = manifest_dict.get("BAM")
        if fastq_values and bam_value:
            raise ManifestValidationError(f"Both FASTQ and BAM fields found in {manifest_path}, only one file type is allowed")
        if bam_value:
            data_files = [Path(bam_value)]
        elif fastq_values:
            data_files = [Path(value) for value in fastq_values]
        else:
            raise ManifestValidationError(f"Neither FASTQ nor BAM field found in {manifest_path}")
    else:
        fasta_value = manifest_dict.get("FASTA")
        if not fasta_value:
            raise ManifestValidationError(f"FASTA field not found in {manifest_path}")
        data_files = [Path(fasta_value)]

    # Validate files only when absolute paths are used or a root directory is provided
    for data_file in data_files:
        if data_file.is_absolute():
            logger.debug(f"Manifest data file path is absolute: {data_file}")
            if not data_file.exists():
                raise ManifestValidationError(f"Data file not found: {data_file}")
        elif input_dir:
            full_data_path = Path(input_dir) / data_file
            logger.debug(f"Resolved data file path using input_dir: {full_data_path}")
            if not full_data_path.exists():
                raise ManifestValidationError(f"Data file not found: {full_data_path}")

    logger.debug(f"Manifest check complete: {manifest_file}")
    return alias


def get_webin_cli_command(
    manifest: str,
    context: str,
    webin: str,
    mode: str,
    test: bool,
    output_dir: str,
    input_dir: str | None = None,
    jar: str | None = None,
    java_heap_size_initial: int | None = None,
    java_heap_size_max: int | None = None,
) -> list[str]:
    """
    Build the webin-cli command list based on execution method.

    Args:
        manifest (str): Path to manifest file.
        context (str): Submission context (genome, assembly, etc.).
        webin (str): ENA Webin username.
        mode (str): Execution mode ('submit' or 'validate').
        test (bool): Whether to use test server.
        output_dir (str): Directory where webin-cli writes output files.
        input_dir (str | None): Root directory for relative data file paths declared in the manifest.
        jar (str | None): Path to webin-cli jar file, if using jar execution.
        java_heap_size_initial (int | None): Java initial heap size in GB (-Xms).
        java_heap_size_max (int | None): Java maximum heap size in GB (-Xmx).

    Returns:
        list[str]: Complete command list ready for execution.

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
        if not webin_cli_path:
            raise WebinCLINotFoundError(
                "ena-webin-cli was not found in PATH. Install it with mamba/conda or use --download-webin-cli / --webin-cli-jar options."
            )
        logger.info(f"Using ena-webin-cli from: {webin_cli_path}")
        cmd = ["ena-webin-cli"]

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

    if input_dir:
        cmd.append(f"-inputDir={input_dir}")

    return cmd


def prepare_output_dir(manifest: str, assembly_name: str, outdir: str | None = None) -> str:
    """
    Creates a timestamped output directory for each submission made by webin-cli.

    Args:
        manifest (str): Path to manifest file.
        assembly_name (str): Assembly name from manifest ASSEMBLYNAME field.
        outdir (str | None): Base output directory. Defaults to manifest parent directory.

    Returns:
        str: Path to the created output directory.
    """
    # If outdir is not provided, use the parent directory of the manifest because it is a default behaviour of webin-cli
    output_root = Path(outdir) if outdir else Path(manifest).parent

    # We need to have separate folder for each submission, otherwise webin-cli may produce no output files
    # when resubmitting an existing assembly and if output folder already exists and contains files
    timestamp = time.strftime("%d-%b-%y_%H%M%S")  # Format: 10-Apr-26_142315
    output_dir = output_root / f"{assembly_name}_{timestamp}"

    if output_dir.exists() and output_dir.is_dir():
        logger.warning(f"Output directory already exists: {output_dir}")

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
    input_dir: str | None = None,
    jar: str | None = None,
    retries: int = RETRIES,
    retry_delay: int = RETRY_DELAY,
    java_heap_size_initial: int | None = None,
    java_heap_size_max: int | None = None,
) -> None:
    """
    Execute webin-cli with retry logic and error handling.

    Args:
        manifest (str): Path to manifest file.
        context (str): Submission context (genome, assembly, etc.).
        webin (str): ENA Webin username.
        mode (str): Execution mode ('submit' or 'validate').
        test (bool): Whether to use test server.
        output_dir (str): Directory where webin-cli writes output files.
        input_dir (str | None): Root directory for relative data file paths declared in the manifest.
        jar (str | None): Path to webin-cli jar file, if using jar execution.
        retries (int): Number of retry attempts.
        retry_delay (int): Initial retry delay in seconds.
        java_heap_size_initial (int | None): Java initial heap size in GB (-Xms).
        java_heap_size_max (int | None): Java maximum heap size in GB (-Xmx).

    Returns normally when webin-cli succeeds, or when the submitted object already exists.
    Only webin-cli SYSTEM_ERROR failures (exit code 1) are retried. USER_ERROR (2) and
    VALIDATION_ERROR (3) failures are deterministic, so they fail immediately.

    Raises:
        WebinCredentialsError: If credentials are invalid.
        WebinAbsConsentError: If the account has not completed the ABS acknowledgement.
        WebinCLIExecutionError: If webin-cli fails with a non-retryable error or all retry attempts fail.
    """

    cmd = get_webin_cli_command(manifest, context, webin, mode, test, output_dir, input_dir, jar, java_heap_size_initial, java_heap_size_max)
    logger.info(f"Command: {' '.join(cmd)}")

    for attempt in range(1, retries + 1):
        logger.info(f"🚀 Running ena-webin-cli (attempt {attempt}/{retries})...")

        result = subprocess.run(cmd, capture_output=True, text=True)
        output_parts = [part.strip() for part in (result.stdout, result.stderr) if part and part.strip()]
        combined_log = "\n".join(output_parts) or "webin-cli produced no stdout/stderr output"
        logger.debug(
            f"webin-cli process finished with return code {result.returncode}, stdout_len={len(result.stdout or '')}, stderr_len={len(result.stderr or '')}"
        )
        logger.debug("ena-webin-cli output:\n" + combined_log)

        if result.returncode == 0:
            logger.info("✅ Submission validation succeeded" if mode == "validate" else "✅ Successfully submitted object for the first time")
            return

        # Resubmission: webin-cli exits with a system error, but the object is already in the account
        if MSG_ALREADY_EXISTS in combined_log:
            logger.info("✅ Submitted object already exists")
            return

        logger.warning(f"❌ ena-webin-cli exited with non-zero code {result.returncode}")

        if MSG_INVALID_CREDENTIALS in combined_log:
            raise WebinCredentialsError("Invalid credentials for Webin account")

        if MSG_ABS_CONSENT in combined_log:
            raise WebinAbsConsentError(
                "Incomplete Access and Benefit Sharing (ABS) acknowledgement for Webin account. "
                "Log in to the Webin Submissions Portal to complete it."
            )

        logger.error("ena-webin-cli output:\n" + combined_log)

        if result.returncode in WEBIN_CLI_NON_RETRYABLE_EXIT_CODES:
            error_type = "user" if result.returncode == WEBIN_CLI_EXIT_USER_ERROR else "validation"
            raise WebinCLIExecutionError(f"Webin-CLI failed with a {error_type} error (exit code {result.returncode}). See logs for details.")

        if attempt == retries:
            raise WebinCLIExecutionError(f"Webin-CLI execution failed after {retries} attempts. See logs for details.")
        sleep_time = retry_delay * (2 ** (attempt - 1))  # exponential backoff
        logger.warning(f"🔁 Retrying in {sleep_time} seconds...")
        logger.debug(f"Retry backoff details: attempt={attempt}, max_retries={retries}, sleep_time={sleep_time}")
        time.sleep(sleep_time)


def parse_accessions_from_report(report_text: str, context: str) -> str | None:
    """
    Get the accession assigned to the submitted object from a webin-cli report.

    Reads submissions need both an experiment and a run accession, returned as "ERX...,ERR...".
    Other contexts need an assembly (analysis) accession.

    Args:
        report_text (str): Contents of a webin-cli report file.
        context (str): Submission context (genome, reads, etc.).

    Returns:
        str | None: The accession string, or None if a required accession is missing or ambiguous.
    """
    if context == "reads":
        patterns = {"experiment": ENA_EXPERIMENT_ACCESSION_REGEX, "run": ENA_RUN_ACCESSION_REGEX}
    else:
        patterns = {"assembly": ENA_ASSEMBLY_ACCESSION_REGEX}

    accessions = []
    for kind, pattern in patterns.items():
        matches = sorted(set(re.findall(pattern, report_text)))
        if len(matches) > 1:
            logger.warning(f"Multiple {kind} accessions found in report: {','.join(matches)}")
        if len(matches) != 1:
            return None
        accessions.append(matches[0])
    return ",".join(accessions)


def check_report(output_dir: str, mode: str, context: str) -> str | None:
    """
    Parse the webin-cli report file to confirm the operation outcome and get the assigned accession.

    Args:
        output_dir (str): The webin-cli output directory, containing webin-cli.report.
        mode (str): Operation mode - either 'submit' or 'validate'.
        context (str): Submission context used to choose accession parsing logic.

    Returns:
        str | None: Accession assigned to the submitted object (see parse_accessions_from_report()),
            for both first-time submissions and resubmissions. None in validate mode.

    Raises:
        WebinReportValidationError: If the report is missing, has no success message,
            or (in submit mode) does not contain exactly one accession of each expected kind.
            The report content is logged before raising.
    """
    report_path = Path(output_dir) / REPORT_FILE
    logger.info(f"Checking webin report {report_path}")

    if not report_path.exists():
        raise WebinReportValidationError(f"Report file not found: {report_path}")

    report_text = report_path.read_text()

    success_messages = [MSG_VALIDATE_SUCCESS] if mode == "validate" else [MSG_SUBMIT_SUCCESS, MSG_ALREADY_EXISTS]
    accession = None
    error = None
    if not any(msg in report_text for msg in success_messages):
        error = f"{mode.capitalize()} success message not found in report"
        logger.error(error)

    elif mode == "submit":
        accession = parse_accessions_from_report(report_text, context)
        if not accession:
            error = "Expected exactly one accession of each kind in report"
            logger.error(error)

    if error:
        logger.error(f"Report content ({report_path}):\n{report_text}")
        raise WebinReportValidationError(f"{error}: {report_path}")
    return accession


def resolve_manifests(manifest_input: str) -> list[str]:
    """
    Resolve a manifest argument to a list of manifest file paths.

    Accepts either a single manifest file or a directory containing manifest files
    (matched by the ``*.manifest`` glob pattern).

    Args:
        manifest_input (str): Path to a manifest file or a directory of manifests.

    Returns:
        list[str]: Sorted list of manifest file paths.

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


def write_assigned_accessions(accessions: dict[str, str], output_accessions: str) -> None:
    """
    Write assigned accessions for all successful submissions to a TSV file.

    This implementation rewrites the full TSV on each call using an atomic replace pattern (write to .tmp, then replace).
    It is safer than writing in "append" mode and keeps a complete table of successful submissions for resume safety.
    The tradeoff is redundant I/O growth (O(n²) complexity) over long runs (about 5 seconds cumulative write time
    and about 600 MB redundant writes for 5000 submissions).
    TODO: consider switching to buffered checkpoints (flush every N successes) or other more efficient persistence method.

    Args:
        accessions (dict[str, str]): Mapping of assembly name to accession.
        output_accessions (str): Path to output TSV file (overwritten each run).
    """
    logger.debug(f"Preparing to write accessions table to {output_accessions}")
    if not accessions:
        logger.warning("No accessions to write.")
        return
    output_file = Path(output_accessions)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output_file.with_suffix(output_file.suffix + ".tmp")
    with open(tmp_output, "w", newline="") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t")
        writer.writerow(["alias", "accession"])
        for assembly_name, accession in accessions.items():
            writer.writerow([assembly_name, accession])
    tmp_output.replace(output_file)
    logger.debug(f"Assigned accessions written to {output_file}")


def load_assigned_accessions(output_accessions: str) -> dict[str, str]:
    """
    Load previously written accessions TSV, if present.

    Args:
        output_accessions (str): Path to the TSV file containing previously assigned accessions.

    Returns:
        dict[str, str]: Mapping of assembly name to accession.
    """
    output_file = Path(output_accessions)
    if not output_file.exists():
        return {}

    loaded: dict[str, str] = {}
    try:
        with open(output_file, "r", newline="") as tsv_in:
            reader = csv.DictReader(tsv_in, delimiter="\t")
            for row in reader:
                alias = row["alias"]
                accession = row["accession"]
                loaded[alias] = accession
    except Exception as exc:
        logger.warning(f"Failed to read existing accessions file {output_file}: {exc}")
        return {}

    logger.info(f"Loaded {len(loaded)} previously written accession(s) from {output_file}")
    return loaded


def main() -> int:
    """
    Main entry point for the webin-cli handler.

    Returns:
        int: Exit code (0 for success, non-zero for failure).
    """
    try:
        args = parse_arguments()
        configure_logging(debug=args.debug)
        logger.debug(f"Parsed arguments: {vars(args)}")

        webin = get_webin_credentials()

        if args.download_webin_cli:
            execute_jar = download_webin_cli(version=args.download_webin_cli_version, dest_dir=args.download_webin_cli_directory)
        else:
            execute_jar = args.webin_cli_jar

        manifests = resolve_manifests(args.manifest)
        logger.info(f"Resolved {len(manifests)} manifest(s) for processing")

        accessions_to_write: dict[str, str] = {}
        if args.resume:
            accessions_to_write = load_assigned_accessions(args.output_accessions)
            aliases_to_skip = set(accessions_to_write.keys())

        for manifest in manifests:
            logger.debug(f"Starting submission loop for manifest: {manifest}")
            assembly_name = check_manifest(manifest, args.input_dir, args.context)
            if args.resume and assembly_name in aliases_to_skip:
                logger.info(f"Skipping {assembly_name}: already present in accessions file with accession {accessions_to_write[assembly_name]}")
                continue
            output_dir = prepare_output_dir(manifest, assembly_name, args.outdir)
            run_webin_cli(
                manifest=manifest,
                context=args.context,
                webin=webin,
                mode=args.mode,
                test=args.test,
                output_dir=output_dir,
                input_dir=args.input_dir,
                jar=execute_jar,
                retries=args.retries,
                retry_delay=args.retry_delay,
                java_heap_size_initial=args.java_heap_size_initial,
                java_heap_size_max=args.java_heap_size_max,
            )

            logger.debug(f"Result location for manifest {manifest}: {output_dir}")
            accession = check_report(output_dir, args.mode, args.context)

            if args.mode == "submit":
                logger.info(f"Assigned accession for {assembly_name}: {accession}")
                accessions_to_write[assembly_name] = accession
                # Persist progress immediately so long-running runs can resume after interruption
                write_assigned_accessions(accessions_to_write, args.output_accessions)
            logger.info(f"Submission/validation done for {manifest}")

        if args.mode == "submit":
            logger.debug(f"Submission mode complete with {len(accessions_to_write)} accession(s) captured")
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
    except WebinAbsConsentError as e:
        logger.error(f"Webin account error: {e}")
        return 7
    except WebinReportValidationError as e:
        logger.error(f"Webin-CLI report error: {e}")
        return 8
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 99


if __name__ == "__main__":
    sys.exit(main())
