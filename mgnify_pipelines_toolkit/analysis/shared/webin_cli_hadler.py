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
import os
import subprocess
ENA_WEBIN = "ENA_WEBIN"
ENA_WEBIN_PASSWORD = "ENA_WEBIN_PASSWORD"
REPORT_FILE = "webin-cli.report"
WEBIN_SUBMISSION_RESULT_FILES = ['analysis.xml', 'receipt.xml', 'submission.xml', 'webin-submission.xml']
"""
from mgnify_pipelines_toolkit.constants.webin_cli import \
    (
    ENA_WEBIN,
    ENA_WEBIN_PASSWORD,
    REPORT_FILE,
    WEBIN_SUBMISSION_RESULT_FILES
)
"""

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
        choices=['genome', 'transcriptome', 'sequence', 'polysample', 'reads', 'taxrefset']
    )
    parser.add_argument(
        "--mode", required=True, type=str, help="submit or validate", choices=['submit', 'validate']
    )
    parser.add_argument(
        "--test", required=False, action="store_true", help="Specify to use test server instead of live"
    )
    return parser.parse_args()


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
    if 'Invalid submission account user name or password.' in log:
        return "Invalid submission account user name or password."
    return message


def check_manifest(manifest_file):
    if not os.path.exists(manifest_file):
        print(f"{manifest_file} does not exist. Exit")
        exit(1)

    with open(manifest_file, 'r') as file_in:
        for line in file_in:
            if 'FASTA' in line:
                fasta_location = '/'.join(line.replace('FASTA', '').strip().split('/')[:-1])
            if 'ASSEMBLYNAME' in line:
                assembly_name = line.replace('ASSEMBLYNAME', '').strip()
    return fasta_location, assembly_name


def run_webin_cli(manifest, context, webin, password, mode, test):
    cmd = [
        "ena-webin-cli",
        f"-context={context}",
        f"-manifest={manifest}",
        f"-userName={webin}",
        f"-password=\'{password}\'",
        f"-{mode}"
    ]
    if test:
        cmd.append("-test")

    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, env=os.environ.copy())

    #if result.stderr:
    #    print("STDERR:", result.stderr)

    if result.returncode != 0:
        print(f"❌ Command failed with code {result.returncode}")
        message = handle_webin_failures(result.stdout)
        print(message)
        exit(1)
    else:
        print("✅ Command completed successfully")


def check_report(fasta_location, mode):
    report_file = os.path.join(fasta_location, REPORT_FILE)
    print(f'Checking webin report {report_file}')
    if os.path.exists(report_file):
        with open(report_file) as f:
            report_text = f.read()
        if mode == 'submit':
            if "submission has been completed successfully" in report_text:
                return True
            elif "object being added already exists in the submission account with accession" in report_text:
                return True
            else:
                return False
        elif mode == 'validate':
            if 'Submission(s) validated successfully.':
                return True
            else:
                return False
    else:
        print(f"⚠️ Report file not found: {report_file}")


def check_result(fasta_location, context, assembly_name, mode):
    report_status = check_report(fasta_location, mode)
    if not report_status:
        print(f'Command failed. Check logs')
    else:
        # check submission folder
        if mode == "submit":
            if not os.path.exists(os.path.join(fasta_location, context, assembly_name)):
                print(f'No result found for {assembly_name}')
                return False
            else:
                # check mode folder
                submission_path = os.path.join(fasta_location, context, assembly_name, mode)
                if not os.path.exists(submission_path):
                    print(f'No submission result found for {assembly_name}/{mode}')
                    return False
                else:
                    # check files existence
                    for item in WEBIN_SUBMISSION_RESULT_FILES:
                        if not os.path.exists(os.path.join(submission_path, item)):
                            return False
    return True


def main():
    args = parse_arguments()
    manifest_file = args.manifest

    webin, password = get_webin_credentials()

    fasta_location, assembly_name = check_manifest(manifest_file)

    run_webin_cli(
        manifest=manifest_file,
        context=args.context,
        webin=webin,
        password=password,
        mode=args.mode,
        test=args.test
    )

    result_status = check_result(
        fasta_location=fasta_location,
        context=args.context,
        assembly_name=assembly_name,
        mode=args.mode
    )

    if result_status:
        print(f'Submission/validation done for {manifest_file}')
    else:
        print(f'Submission/validation failed for {manifest_file}')


if __name__ == "__main__":
    main()

