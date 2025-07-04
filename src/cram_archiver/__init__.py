# Copyright (C) 2025 Leiden University Medical Center
# This file is part of cram-archiver
#
# cram-archiver is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# cram-archiver is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with cram-archiver. If not, see <https://www.gnu.org/licenses/

import argparse
import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Dict, Iterator, Optional, Sequence, Set

from ._version import __version__
from .references import ReferenceID

# 3.1 Not supported by some tools currently (2025)
DEFAULT_CRAM_VERSION = "3.0"
DEFAULT_THREADS = 1
DEFAULT_WRITE_INDEX = True
DEFAULT_WRITE_CHECKSUM_FILES = True
DEFAULT_LOG_LEVEL = logging.INFO
DEFAULT_MINIMUM_AGE_DAYS = 0


def convert_to_cram(
        input_file: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = DEFAULT_WRITE_INDEX,
):
    additional_threads = max(0, threads - 1)
    command = [
        "samtools", "view",
        "--output-fmt", f"cram,version={cram_version}",
        "--threads", str(additional_threads),
        "-o", output_file,
        "--reference", reference,
        input_file,
    ]
    if write_index:
        command.append("--write-index")
    logging.debug(f"Conversion command: '{' '.join(command)}'")
    subprocess.run(command, check=True)


def checksum(input_file: str, reference: str, threads: int = 1) -> str:
    additional_threads = max(0, threads - 1)
    result = subprocess.run(
        [
            "samtools",
            "checksum",
            "--all",
            "--threads", str(additional_threads),
            "--reference", reference,
            input_file
        ],
        check=True,
        stdout=subprocess.PIPE)
    return result.stdout.decode("ascii")


def strip_comments_from_checksum(checksum: str) -> str:
    return "".join(
        line for line in checksum.splitlines(keepends=True)
        if not (line.startswith("#") or line == "\n")
    )


def convert_to_cram_and_check(
        input_file: str,
        reference_id_to_path: Dict[ReferenceID, str],
        threads: int = DEFAULT_THREADS,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = DEFAULT_WRITE_INDEX,
        write_checksum_files: bool = DEFAULT_WRITE_CHECKSUM_FILES,
) -> str:
    reference_id = ReferenceID.from_file(input_file)
    reference = reference_id_to_path[reference_id]
    output_file = str(Path(input_file).parent / Path(input_file).stem) + ".cram"
    logging.info(f"Convert '{input_file}' to '{output_file}'.")
    convert_to_cram(input_file, output_file, reference, threads, cram_version,
                    write_index)
    logging.info(f"Checksumming {input_file}.")
    input_checksum = checksum(input_file, reference, threads)
    logging.debug(input_checksum)
    logging.info(f"Checksumming {output_file}.")
    output_checksum = checksum(output_file, reference, threads)
    logging.debug(output_checksum)
    if write_checksum_files:
        with open(input_file + ".checksum", "wt") as f:
            f.write(input_checksum)
        with open(output_file + ".checksum", "wt") as f:
            f.write(output_checksum)
    if (
            strip_comments_from_checksum(input_checksum) !=
            strip_comments_from_checksum(output_checksum)
    ):
        os.unlink(output_file)
        raise RuntimeError(
            f"Input checksum does not match output checksum for {input_file} "
            f"and {output_file}.\n'{input_checksum!r}' != {output_checksum!r}.\n"
            f"{output_file} is removed.")
    return output_file


def handle_file_age(file, file_mtime: float, older_than_timestamp: float
                    ) -> Iterator[str]:
    if file_mtime < older_than_timestamp:
        yield file
    else:
        logging.info(f"Skipping too new file: {file}.")


def _find_bam_files_dirscan(
        input_dir: str,
        older_than_timestamp: float,
        ignore_files: Set[str],
        follow_symlinks: bool,
):
    for entry in os.scandir(input_dir):
        if entry.path in ignore_files:
            logging.info(f"Ignoring {entry.path}")
            continue
        logging.debug(f"Searching: {entry.path}")
        if entry.is_file(follow_symlinks=follow_symlinks):
            if entry.name.endswith(".bam"):
                yield from handle_file_age(
                    entry.path, entry.stat().st_mtime, older_than_timestamp)
        elif entry.is_dir(follow_symlinks=follow_symlinks):
            yield from _find_bam_files_dirscan(
                entry.path, older_than_timestamp, ignore_files, follow_symlinks)


def find_bam_files(
        input_path: str,
        older_than_timestamp: float = time.time(),
        ignore_files: Optional[Sequence[str]] = None,
        follow_symlinks=False,
) -> Iterator[str]:
    # Make input path and ignore files absolute. This also deals with trailing
    # slashes for directories and ../ entries.
    input_path = os.path.abspath(input_path)
    if ignore_files is not None:
        ignore_set = {os.path.abspath(p) for p in ignore_files}
    else:
        ignore_set = set()
    if input_path in ignore_set:
        logging.info(f"Ignoring {input_path}")
        return
    if os.path.islink(input_path) and not follow_symlinks:
        return
    if os.path.isfile(input_path):
        yield from handle_file_age(
            input_path, os.path.getmtime(input_path), older_than_timestamp)
    elif os.path.isdir(input_path):
        yield from _find_bam_files_dirscan(
            input_path, older_than_timestamp, ignore_set, follow_symlinks
        )


def cram_archiver(
        input_path: str,
        reference_files: Sequence[str],
        threads: int = DEFAULT_THREADS,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = DEFAULT_WRITE_INDEX,
        write_checksum_files: bool = DEFAULT_WRITE_CHECKSUM_FILES,
        minimum_age_days: int = 0,
        delete: bool = False,
        dry_run: bool = False,
        ignore_files: Optional[Sequence[str]] = None
):
    if delete and not dry_run:
        logging.warning(
            "WARNING: BAM FILES WILL BE DELETED AFTER SUCCESSFUL CONVERSION!!!"
        )
    older_than_timestamp = time.time() - (minimum_age_days * 24 * 60 * 60)

    ref_dicts: Dict[ReferenceID, str] = {}
    for reference in reference_files:
        fai = reference + ".fai"
        if not os.path.exists(fai):
            raise FileNotFoundError(
                f"Fasta index file for {reference} could not be found.")
        ref_id = ReferenceID.from_file(fai)
        ref_dicts[ref_id] = reference

    bam_files = find_bam_files(input_path, older_than_timestamp, ignore_files)
    number_of_bam_files = 0
    errors = []
    total_bam_size = 0
    total_cram_size = 0
    for number_of_bam_files, bam in enumerate(bam_files, start=1):
        bam_size = os.path.getsize(bam)
        bam_name = os.path.basename(bam)
        total_bam_size += bam_size
        if dry_run:
            print(bam)
            continue
        try:
            cram_file = convert_to_cram_and_check(
                input_file=bam,
                reference_id_to_path=ref_dicts,
                threads=threads,
                cram_version=cram_version,
                write_index=write_index,
                write_checksum_files=write_checksum_files,
            )
            cram_name = os.path.basename(cram_file)
            cram_size = os.path.getsize(cram_file)
            total_cram_size += cram_size
            logging.info(f"{bam_name} size: {bam_size / (1024 ** 3):.2f} GiB")
            logging.info(f"{cram_name} size: {cram_size / (1024 ** 3):.2f} GiB")
            if delete:
                logging.info(
                    f"Conversion successful, deleting BAM file: {bam}. Saved "
                    f"space: {(bam_size - cram_size) / (1024 ** 3):.2f} GiB."
                )
                os.unlink(bam)
        except (FileNotFoundError, RuntimeError) as error:
            logging.error(f"Conversion unsuccessful: {bam}. {str(error)}")
            errors.append(error)
    if number_of_bam_files == 0:
        logging.warning("No BAM files found. Exiting.")
    else:
        logging.info(
            f"Found {number_of_bam_files} BAM files of "
            f"total: {total_bam_size / (1024 ** 3):.2f} GiB.")
        if not dry_run:
            logging.info(
                f"Total generated CRAM size: "
                f"{total_cram_size / (1024 ** 3):.2f} GiB.")
            if delete:
                logging.info(
                    f"Total saved size: "
                    f"{(total_bam_size - total_cram_size) / (1024 ** 3):.2f} "
                    f"GiB."
                )
    if errors:
        raise RuntimeError("Errors occurred during conversions.")


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path", metavar="PATH",
        help="Path to BAM file or directory to be recursively searched."
    )
    parser.add_argument(
        "-r", "--reference", action="append", required=True,
        help="Reference to be used for CRAM conversion. Can be used multiple "
             "times. Reference will be checked with the BAM file."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=DEFAULT_THREADS,
        help=f"The number of threads used for conversion and checksumming."
             f"Default: {DEFAULT_THREADS}."
    )
    parser.add_argument(
        "-d", "--minimum-age-days", type=int,
        default=DEFAULT_MINIMUM_AGE_DAYS,
        help=f"The minimum last modification of the BAM file in days prior. "
             f"This assumes the system clock timezone matches that of the "
             f"file while also assuming that every day has 24x60x60 seconds. "
             f"Default {DEFAULT_MINIMUM_AGE_DAYS}.",
    )
    parser.add_argument(
        "--delete", action="store_true",
        help="Delete BAM files after successful conversion."
    )
    parser.add_argument(
        "--cram-version", default=DEFAULT_CRAM_VERSION,
        help="CRAM version to use for CRAM conversion. "
             f"Default: {DEFAULT_CRAM_VERSION}."
    )
    parser.add_argument(
        "--exclude", action="append",
        help="Exclude file or directory from conversion. "
             "Can be supplied multiple times."
    )
    parser.add_argument(
        "--exclude-list", metavar="PATH",
        help="Supply a newline-separated file with files and directories to "
             "exclude."
    )
    parser.add_argument(
        "--dont-write-checksums", action="store_false", dest="write_checksums",
        help="Do not store samtools checksum output on disk."
    )
    parser.add_argument(
        "--dont-write-index", action="store_false", dest="write_index",
        help="Do not write index files for CRAM files."
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the paths of the to be archived BAM files. Perform no actions."
    )
    parser.add_argument(
        "-v", "--verbose", action="count", default=0,
        help="Display more logging information."
    )
    parser.add_argument(
        "-q", "--quiet", action="count", default=0,
        help="Display less logging information."
    )
    parser.add_argument(
        "--version", action="version", version=__version__
    )
    return parser


def cram_archiver_main(*args):
    arg = argument_parser().parse_args(args or None)
    log_level = (arg.quiet - arg.verbose) * 10 + DEFAULT_LOG_LEVEL
    logger = logging.getLogger()
    logger.name = "cram-archiver"
    logger.setLevel(log_level)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    formatter = logging.Formatter(
        "{name}:{asctime}:{levelname}: {message}",
        style="{")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logging.info(f"This is cram-archiver version: {__version__}.")
    if arg.exclude is not None:
        exclude_list = arg.exclude[:]
    else:
        exclude_list = []
    if arg.exclude_list is not None:
        with open(arg.exclude_list, "rt") as f:
            exclude_list.extend(f.read().splitlines(keepends=False))

    cram_archiver(
        input_path=arg.path,
        reference_files=arg.reference,
        threads=arg.threads,
        cram_version=arg.cram_version,
        write_index=arg.write_index,
        write_checksum_files=arg.write_checksums,
        minimum_age_days=arg.minimum_age_days,
        delete=arg.delete,
        dry_run=arg.dry_run,
        ignore_files=exclude_list,
    )
