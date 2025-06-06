import os
import subprocess
from pathlib import Path
from typing import Dict

from .references import ReferenceID

# 3.1 Not supported by some tools currently (2025)
DEFAULT_CRAM_VERSION = "3.0"


def convert_to_cram(
        input_file: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = True,
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
        if not line.startswith("#")
    )


def convert_to_cram_and_check(
        input_file: str,
        reference_files: Dict[ReferenceID, str],
        threads: int = 1,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = True,
        write_checksum_files: bool = True,
):
    reference_id = ReferenceID.from_file(input_file)
    reference = reference_files[reference_id]
    output_file = str(Path(input_file).stem) + ".cram"
    convert_to_cram(input_file, output_file, reference, threads, cram_version,
                    write_index)
    input_checksum = checksum(input_file, reference, threads)
    output_checksum = checksum(output_file, reference, threads)
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
