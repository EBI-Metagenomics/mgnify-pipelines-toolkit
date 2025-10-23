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
import urllib.request

from mgnify_pipelines_toolkit.constants.webin_cli import (
    ENA_WEBIN,
    ENA_WEBIN_PASSWORD,
    REPORT_FILE,
    WEBIN_SUBMISSION_RESULT_FILES,
    # RETRIES,
    # RETRY_DELAY,
    # XMS,
    # XMX
)

RETRIES = 3
RETRY_DELAY = 5
XMS = 10
XMX = 10

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_arguments():
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


def get_latest_webin_cli_version():
    api_url = "https://api.github.com/repos/enasequence/webin-cli/releases/latest"
    with urllib.request.urlopen(api_url) as response:
        data = json.load(response)
        return data["tag_name"]


def download_webin_cli(version=None, dest_dir=""):
    if version is None:
        version = get_latest_webin_cli_version()

    jar_url = f"https://github.com/enasequence/webin-cli/releases/download/{version}/webin-cli-{version}.jar"
    dest_path = os.path.join(dest_dir, "webin-cli.jar")

    try:
        logging.info(f"Downloading webin-cli.jar version {version} to {dest_path} ...")
        urllib.request.urlretrieve(jar_url, dest_path)
    except Exception as e:
        logging.error(f"Failed to download webin-cli.jar: {e}")


def ensure_webin_credentials_exist():
    if ENA_WEBIN not in os.environ:
        raise Exception(f"The variable {ENA_WEBIN} is missing from the env.")
    if ENA_WEBIN_PASSWORD not in os.environ:
        raise Exception(f"The variable {ENA_WEBIN_PASSWORD} is missing from the env")


def get_webin_credentials():
    ensure_webin_credentials_exist()
    webin = os.environ.get(ENA_WEBIN)
    password = os.environ.get(ENA_WEBIN_PASSWORD)
    return webin, password


def handle_webin_failures(log):
    message = log
    exit_code = 1
    if "Invalid submission account user name or password." in log:
        message = "Invalid credentials for Webin account. Exit 1"
    if "The object being added already exists in the submission account with accession:" in log:
        message = "Submission already exists. Exit without errors."
        exit_code = 0
    return message, exit_code


def check_manifest(manifest_file, workdir):
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
    return fasta_location, assembly_name, manifest_for_submission


def create_execution_command_jar(jar):
    logging.info(f"Using java execution {jar}")
    logback_config = f"{os.path.abspath(__file__)}/webincli_logback.xml"
    cmd = ["java", f"-Dlogback.configurationFile={logback_config}", f"-Xms{XMS}g", f"-Xmx{XMX}g", "-jar", jar]
    return cmd


def run_webin_cli(manifest, context, webin, password, mode, test, jar=False):

    cmd = ["ena-webin-cli"]  # mamba installation
    if jar:
        cmd = create_execution_command_jar(jar)  # webin-cli java execution

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

        logging.info(f"❌ ena-webin-cli failed with code {result.returncode}")
        message, exit_code = handle_webin_failures(result.stdout)
        logging.info(message)

        if exit_code == 0:
            logging.info("✅ Command completed successfully")
            return result
        elif attempt < RETRIES:
            sleep_time = RETRY_DELAY * (2 ** (attempt - 1))  # exponential backoff
            logging.info(f"🔁 Retrying in {sleep_time} seconds...")
            time.sleep(sleep_time)
        else:
            logging.info("💥 All retries failed.")
            message, exit_code = handle_webin_failures(result.stdout)
            logging.error(message)
            exit(exit_code)


def check_submission_status_test(report_text):
    submission_exists = False
    line_count = report_text.count("\n") + 1 if report_text else 0
    if not "This was a TEST submission(s)." in report_text:
        logging.info("Submission for that object failed on TEST server")
        return False, submission_exists
    else:
        if line_count == 2:
            # webin-cli.report has only message 'This was a TEST submission(s).' in case of re-submission
            logging.info(f"Submission for that object has already exist on TEST server")
            submission_exists = True
            return True, submission_exists
        elif "submission has been completed successfully" in report_text:
            logging.info(f"Submission for that object was done first time on TEST server")
            return True, submission_exists
        else:
            logging.info("Submission for that object failed on TEST server")
            return False, submission_exists


def check_submission_status_live(report_text):
    submission_exists = False
    if "submission has been completed successfully" in report_text:
        logging.info(f"Submission for that object was done first time on MAIN server")
        return True, submission_exists
    elif "object being added already exists in the submission account with accession" in report_text:
        logging.info(f"Submission for that object has already exist on MAIN server")
        submission_exists = True
        return True, submission_exists
    else:
        logging.info(f"Submission for that object failed on MAIN server")
        return False, submission_exists


def check_report(fasta_location, mode, test):
    report_file = os.path.join(fasta_location, REPORT_FILE)
    logging.info(f"Checking webin report {report_file}")
    submission_exists = False
    if os.path.exists(report_file):
        with open(report_file) as f:
            report_text = f.read()

        if mode == "submit":
            if test:
                return check_submission_status_test(report_text)
            else:
                return check_submission_status_live(report_text)

        elif mode == "validate":
            if "Submission(s) validated successfully." in report_text:
                logging.info("Submission validation succeeded")
                return True, submission_exists
            else:
                logging.info("Submission validation failed")
                return False, submission_exists
    else:
        logging.info(f"⚠️ Report file not found: {report_file}")
    return False, False


def check_result(fasta_location, context, assembly_name, mode, test):
    report_status, submission_exists = check_report(fasta_location, mode, test)
    if not report_status:
        logging.info(f"Command failed. Check logs")
    else:
        # check submission folder
        if mode == "submit" and not submission_exists:
            if not os.path.exists(os.path.join(fasta_location, context, assembly_name)):
                logging.info(f"No result found for {assembly_name}")
                return False
            else:
                # check mode folder
                submission_path = os.path.join(fasta_location, context, assembly_name, mode)
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

    manifest_file = args.manifest

    webin, password = get_webin_credentials()

    fasta_location, assembly_name, manifest_for_submission = check_manifest(manifest_file, args.workdir)

    # pick execution method for webin-cli
    execute_jar = False
    if args.download_webin_cli:
        execute_jar = os.path.join(args.download_webin_cli_directory, "webin-cli.jar")
    elif args.webin_cli_jar:
        execute_jar = args.webin_cli_jar
    run_webin_cli(
        manifest=manifest_for_submission, context=args.context, webin=webin, password=password, mode=args.mode, test=args.test, jar=execute_jar
    )

    result_status = check_result(fasta_location=fasta_location, context=args.context, assembly_name=assembly_name, mode=args.mode, test=args.test)

    if result_status:
        logging.info(f"Submission/validation done for {manifest_for_submission}")
    else:
        logging.error(f"Submission/validation failed for {manifest_for_submission}")


if __name__ == "__main__":
    main()
