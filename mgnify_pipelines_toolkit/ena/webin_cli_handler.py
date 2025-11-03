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
import logging
import json
import os
import time
import subprocess
import shutil
import urllib.request
from pathlib import Path
import re
from typing import Tuple, Optional

ENA_WEBIN = "ENA_WEBIN"
ENA_WEBIN_PASSWORD = "ENA_WEBIN_PASSWORD"
REPORT_FILE = "webin-cli.report"
WEBIN_SUBMISSION_RESULT_FILES = ["analysis.xml", "receipt.xml", "submission.xml", "webin-submission.xml"]
RETRIES = 3
RETRY_DELAY = 5
JAVA_HEAP_SIZE_INITIAL = 10
JAVA_HEAP_SIZE_MAX = 10
INSDC_CENTRE_PREFIXES = "EDS"
ENA_ASSEMBLY_ACCESSION_REGEX = f"([{INSDC_CENTRE_PREFIXES}]RZ[0-9]{{6,}})"


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_arguments():
    """
    Parse command-line arguments for the Webin-CLI submission tool.

    Returns:
        argparse.Namespace: Parsed arguments object containing attributes:
            - manifest (str): Path to manifest file
            - context (str): Submission context: genome, transcriptome, etc.
            - mode (str): 'submit' or 'validate'
            - test (bool): Whether to use Webin test server
            - workdir (Optional[str]): Working directory for temporary output
            - download_webin_cli (bool): Whether to download Webin-CLI jar
            - download_webin_cli_directory (str): Path to store downloaded Webin-CLI
            - download_webin_cli_version (Optional[str]): Specific Webin-CLI version
            - webin_cli_jar (Optional[str]): Path to pre-downloaded jar file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--manifest",
        required=True,
        type=str,
        help="Manifest text file containing file and metadata fields",
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
    return parser.parse_args()


def get_latest_webin_cli_version() -> str:
    """
    Fetch the latest release version tag of ena-webin-cli from GitHub API.

    Returns:
        str: Latest version string (e.g. '9.0.1')
    """
    api_url = "https://api.github.com/repos/enasequence/webin-cli/releases/latest"
    with urllib.request.urlopen(api_url) as response:
        data = json.load(response)
        return data["tag_name"]


def download_webin_cli(version: Optional[str] = None, dest_dir: str = "") -> None:
    """
    Download the Webin-CLI `.jar` file.

    Args:
        version (str, optional): Version to download. If None, fetch latest.
        dest_dir (str): Directory to save the jar file.

    Returns:
        None
    """
    if version is None:
        version = get_latest_webin_cli_version()

    jar_url = f"https://github.com/enasequence/webin-cli/releases/download/{version}/webin-cli-{version}.jar"
    dest_path = os.path.join(dest_dir, "webin-cli.jar")

    try:
        logging.info(f"Downloading webin-cli.jar version {version} to {dest_path} ...")
        urllib.request.urlretrieve(jar_url, dest_path)
    except Exception as e:
        logging.error(f"Failed to download webin-cli.jar: {e}")


def ensure_webin_credentials_exist() -> None:
    if ENA_WEBIN not in os.environ:
        raise Exception(f"The variable {ENA_WEBIN} is missing from the env.")
    if ENA_WEBIN_PASSWORD not in os.environ:
        raise Exception(f"The variable {ENA_WEBIN_PASSWORD} is missing from the env")


def get_webin_credentials() -> Tuple[str, str]:
    ensure_webin_credentials_exist()
    webin = os.environ.get(ENA_WEBIN)
    password = os.environ.get(ENA_WEBIN_PASSWORD)
    return webin, password


def handle_webin_failures(log: str) -> Tuple[str, float]:
    if "Invalid submission account user name or password." in log:
        return "Invalid credentials for Webin account. Exit 1", 1
    if "The object being added already exists in the submission account with accession:" in log:
        return "Submission already exists. Exit without errors.", 0
    return log, 1


def check_manifest(manifest_file: str, workdir: Optional[str]) -> Tuple[str, str]:
    """
    Validate a manifest file and make paths absolute if needed.

    Args:
        manifest_file (str): Path to the original manifest.
        workdir (str or None): Optional working directory to prefix FASTA paths.

    Returns:
        Tuple[str, str]:
            - assembly_name (str): Parsed value of ASSEMBLYNAME.
            - manifest_for_submission (str): Path to the (possibly updated) manifest file.
    """
    if not os.path.exists(manifest_file):
        logging.error(f"{manifest_file} does not exist. Exit")
        exit(1)

    with open(manifest_file, "r") as file_in:
        for line in file_in:
            if "FASTA" in line:
                fasta_location = "/".join(line.replace("FASTA", "").strip().split("/")[:-1])
            if "ASSEMBLYNAME" in line:
                assembly_name = line.replace("ASSEMBLYNAME", "").strip()

    manifest_for_submission = manifest_file
    # add absolute path to FASTA and create a new updated manifest
    if workdir and not fasta_location.startswith("/"):
        manifest_location = os.path.dirname(manifest_file)
        manifest_name = os.path.basename(manifest_file)
        manifest_for_submission = os.path.join(manifest_location, f"updated_{manifest_name}")
        # add absolute path for FASTA
        with open(manifest_for_submission, "w") as file_out, open(manifest_file, "r") as file_in:
            for line in file_in:
                new_line = line
                if "FASTA" in line:
                    path_fasta = line.split("\t")[1]
                    new_line = f"FASTA\t{workdir}/{path_fasta}"
                file_out.write(new_line)
        logging.info(f"New manifest {manifest_for_submission} was created")
    return assembly_name, manifest_for_submission


def create_execution_command_jar(jar: str) -> str:
    logging.info(f"Using java execution {jar}")
    logback_config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "webincli_logback.xml")
    cmd = ["java", f"-Dlogback.configurationFile={logback_config}", f"-Xms{JAVA_HEAP_SIZE_INITIAL}g", f"-Xmx{JAVA_HEAP_SIZE_MAX}g", "-jar", jar]
    return cmd


def is_webin_cli_installed():
    return shutil.which("ena-webin-cli") is not None


def run_webin_cli(manifest: str, context: str, webin: str, password: str, mode: str, test: bool, jar: bool =False):
    if jar:
        cmd = create_execution_command_jar(jar)  # webin-cli java execution
    else:
        # mamba installation in env
        if is_webin_cli_installed():
            cmd = ["ena-webin-cli"]  # mamba installation
        else:
            logging.error("ena-webin-cli was not found. Install it with mamba or use --download-webin-cli in execution")
            exit(1)
    cmd += [f"-context={context}", f"-manifest={manifest}", f"-userName={webin}", f"-password='{password}'", f"-{mode}"]
    if test:
        cmd.append("-test")
    cmd_str = " ".join(cmd)

    for attempt in range(1, RETRIES + 1):
        logging.info(f"🚀 Running ena-webin-cli (attempt {attempt}/{RETRIES})...")

        result = subprocess.run(
            cmd_str,
            shell=True,
            capture_output=True,
            text=True,
            env=os.environ.copy(),
        )

        if result.returncode == 0:
            logging.info("✅ Command completed successfully")
            return result  # success

        logging.error(f"❌ ena-webin-cli failed with code {result.returncode}")
        message, exit_code = handle_webin_failures(result.stdout)
        logging.info(message)

        if exit_code == 0:
            logging.info("✅ Command completed successfully")
            return result
        elif attempt < RETRIES:
            sleep_time = RETRY_DELAY * (2 ** (attempt - 1))  # exponential backoff
            logging.warning(f"🔁 Retrying in {sleep_time} seconds...")
            time.sleep(sleep_time)
        else:
            logging.error("💥 All retries failed.")
            message, exit_code = handle_webin_failures(result.stdout)
            logging.error(message)
            exit(exit_code)


def check_submission_status_test(report_text: str) -> Tuple[bool, bool]:
    """
    Check the test Webin-CLI submission report and determine submission status.

    Args:
        report_text (str): Contents of the Webin-CLI report file.

    Returns:
        Tuple[bool, bool]: (success, submission_exists)
            - success (bool): True if submission or validation succeeded.
            - submission_exists (bool): True if this is a resubmission
              (object already exists on the test server).
    """
    line_count = report_text.count("\n") + 1 if report_text else 0
    if "This was a TEST submission(s)." not in report_text:
        logging.info("Submission failed on TEST server")
        return False, False
    else:
        if line_count == 2:
            # webin-cli.report has only message 'This was a TEST submission(s).' in case of re-submission
            logging.info("Submitted object already exists on TEST server")
            return True, True
        elif "submission has been completed successfully" in report_text:
            logging.info("Successfully submitted object for the first time on TEST server")
            return True, False
        else:
            logging.info("Submission failed on TEST server")
            return False, False


def check_submission_status_live(report_text: str) -> Tuple[bool, bool]:
    # for assembly webin report doesn't have information about re-submission
    # second attempt just return the same result as first
    submission_exists = False
    if "submission has been completed successfully" in report_text:
        logging.info("Successfully submitted object on MAIN server")
        return True, submission_exists
    elif "object being added already exists in the submission account with accession" in report_text:
        logging.info("Submitted object already exists on MAIN server")
        submission_exists = True
        return True, submission_exists
    else:
        logging.info("Submission failed on MAIN server")
        return False, submission_exists


def check_report(fasta_location: str, mode: str, test: bool) -> Tuple[bool, bool]:
    """
    Function checks webin report and report status of submission and existence
    Returns:
        (success: bool, submission_exists: bool)
    """
    report_file = os.path.join(fasta_location, REPORT_FILE)
    logging.info(f"Checking webin report {report_file}")

    if not Path(report_file).exists():
        logging.warning(f"⚠️ Report file not found: {report_file}")
        return False, False

    with open(report_file) as f:
        report_text = f.read()

        assembly_accession_match_list = re.findall(ENA_ASSEMBLY_ACCESSION_REGEX, report_text)
        if assembly_accession_match_list:
            logging.info(f"Assigned accessions: {','.join(assembly_accession_match_list)}")

    if mode == "submit":
        if test:
            return check_submission_status_test(report_text)
        else:
            return check_submission_status_live(report_text)

    elif mode == "validate":
        if "Submission(s) validated successfully." in report_text:
            logging.info("Submission validation succeeded")
            return True, False
        else:
            logging.info("Submission validation failed")
            return False, False


def check_result(result_location: str, context: str, assembly_name: str, mode: str, test: bool) -> bool:
    report_status, submission_exists = check_report(result_location, mode, test)
    if not report_status:
        logging.info("Command failed. Check logs")
    else:
        # check submission folder
        if mode == "submit" and not submission_exists:
            if not os.path.exists(os.path.join(result_location, context, assembly_name)):
                logging.info(f"No result found for {assembly_name}")
                return False
            else:
                # check mode folder
                submission_path = os.path.join(result_location, context, assembly_name, mode)
                if not os.path.exists(submission_path):
                    logging.info(f"No submission result found for {assembly_name}/{mode}")
                    return False
                else:
                    # check files existence
                    for item in WEBIN_SUBMISSION_RESULT_FILES:
                        if not os.path.exists(os.path.join(submission_path, item)):
                            logging.info(f"No {item} found in {os.path.join(submission_path, item)}")
                            return False
                    logging.info(f"Result folder {submission_path} is OK")
    return True


def main():
    args = parse_arguments()

    if args.download_webin_cli:
        download_webin_cli(version=args.download_webin_cli_version, dest_dir=args.download_webin_cli_directory)

    webin, password = get_webin_credentials()

    assembly_name, manifest_for_submission = check_manifest(args.manifest, args.workdir)

    # pick execution method for webin-cli
    execute_jar = False
    if args.download_webin_cli:
        execute_jar = os.path.join(args.download_webin_cli_directory, "webin-cli.jar")
    elif args.webin_cli_jar:
        execute_jar = args.webin_cli_jar
    run_webin_cli(
        manifest=manifest_for_submission, context=args.context, webin=webin, password=password, mode=args.mode, test=args.test, jar=execute_jar
    )

    result_status = check_result(
        result_location=os.path.dirname(manifest_for_submission), context=args.context, assembly_name=assembly_name, mode=args.mode, test=args.test
    )

    if result_status:
        logging.info(f"Submission/validation done for {manifest_for_submission}")
    else:
        logging.error(f"Submission/validation failed for {manifest_for_submission}")


if __name__ == "__main__":
    main()
