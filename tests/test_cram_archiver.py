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
import itertools
import logging
import math
import os.path
import shutil
import time
from pathlib import Path

from cram_archiver import (
    checksum,
    convert_to_cram,
    convert_to_cram_and_check,
    cram_archiver,
    cram_archiver_main,
    find_bam_files,
    handle_file_age,
    strip_comments_from_checksum,
)
from cram_archiver.references import ReferenceID

import pytest

TEST_DATA = Path(__file__).parent / "data"


def get_file_cram_version(cram_file: str):
    with open(cram_file, "rb") as f:
        magic = f.read(6)
    if magic[:4] != b"CRAM":
        raise ValueError("Not a valid CRAM file")
    return f"{magic[4]}.{magic[5]}"


def test_checksum():
    bam = str(TEST_DATA / "GM24385_1.bam")
    sum = checksum(bam, str(TEST_DATA / "NC012920.1.fasta"))
    checksum_expected = (
        "all        all           178  06d90a04  526c1585  3e8f585f  "
        "41c54f43  268b6cbb  1ff11ae8  35f542f6  56a9ac23\n"
        "GM24385_fastq-lib1-fastq all           178  06d90a04  526c1585  "
        "3e8f585f  41c54f43  268b6cbb  1ff11ae8  35f542f6  56a9ac23\n"
    )
    assert checksum_expected in sum
    assert checksum_expected == strip_comments_from_checksum(sum)


def test_convert_to_cram(tmp_path):
    tmp_bam = tmp_path / "GM24385_1.bam"
    tmp_cram = tmp_path / "GM24385_1.cram"
    reference = TEST_DATA / "NC012920.1.fasta"
    shutil.copy(TEST_DATA / "GM24385_1.bam",  tmp_path / "GM24385_1.bam")
    convert_to_cram(str(tmp_bam), str(tmp_cram), str(reference))
    assert tmp_cram.exists()


@pytest.mark.parametrize(
    ["cram_version", "write_index", "write_checksum_files"],
    itertools.product(["3.0", "3.1"], [True, False], [True, False])
)
def test_convert_to_cram_and_check(
        tmp_path, caplog, cram_version, write_index, write_checksum_files):
    caplog.set_level(logging.DEBUG)
    tmp_bam = str(tmp_path / "GM24385_1.bam")
    tmp_cram = str(tmp_path / "GM24385_1.cram")
    tmp_crai = tmp_cram + ".crai"
    tmp_bam_checksum = tmp_bam + ".checksum"
    tmp_cram_checksum = tmp_cram + ".checksum"
    shutil.copy(TEST_DATA / "GM24385_1.bam",  tmp_bam)
    reference = TEST_DATA / "NC012920.1.fasta"
    reference_fai = str(reference) + ".fai"
    output_file = convert_to_cram_and_check(
        str(tmp_bam),
        {ReferenceID.from_file(reference_fai): str(reference)},
        cram_version=cram_version,
        write_index=write_index,
        write_checksum_files=write_checksum_files,
    )
    assert os.path.exists(output_file)
    assert output_file == tmp_cram
    assert "samtools view" in caplog.text
    assert f"version={cram_version}" in caplog.text
    assert os.path.exists(tmp_crai) is write_index
    assert os.path.exists(tmp_cram_checksum) is write_checksum_files
    assert os.path.exists(tmp_bam_checksum) is write_checksum_files
    assert get_file_cram_version(tmp_cram) == cram_version


@pytest.mark.parametrize(
    ["file_mtime", "older_than_timestamp", "success"],
    [
        (10.0, 0.0, False),
        (0.0, 10.0, True),
    ]
)
def test_handle_file_age(file_mtime, older_than_timestamp, success, caplog):
    caplog.set_level(logging.INFO)
    result = list(handle_file_age("test", file_mtime, older_than_timestamp))
    if success:
        assert result == ["test"]
        assert caplog.text == ""
    else:
        assert result == []
        assert "test" in caplog.text
        assert "Skipping" in caplog.text


@pytest.mark.parametrize("debug", [True, False])
def test_find_bam_files(tmp_path, caplog, debug):
    if debug:
        caplog.set_level(logging.DEBUG)
    else:
        caplog.set_level(logging.INFO)
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    decoy1 = subdir / "decoy1.txt"
    decoy2 = tmp_path / "decoy2.txt"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    decoy1.touch()
    decoy2.touch()
    os.utime(bam1, (1000, 300))
    os.utime(bam2, (1000, 200))
    os.utime(bam3, (1000, 100))
    result = list(find_bam_files(str(tmp_path), older_than_timestamp=201))
    assert set(result) == {str(bam2), str(bam3)}
    assert str(bam1) in caplog.text
    assert (str(bam2) in caplog.text) is debug
    assert (str(bam3) in caplog.text) is debug
    assert (str(decoy1) in caplog.text) is debug
    assert (str(decoy2) in caplog.text) is debug


@pytest.mark.parametrize("debug", [True, False])
def test_find_bam_files_no_symlinks(tmp_path, caplog, debug):
    if debug:
        caplog.set_level(logging.DEBUG)
    else:
        caplog.set_level(logging.INFO)
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1_link = tmp_path / "hiding_in_plain_sight.bam"
    bam1_link.symlink_to(bam1)
    subdir_link = tmp_path / "decoy_dir"
    subdir_link.symlink_to(subdir)
    decoy1 = subdir / "decoy1.txt"
    decoy2 = tmp_path / "decoy2.txt"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    decoy1.touch()
    decoy2.touch()
    os.utime(bam1, (1000, 300))
    os.utime(bam2, (1000, 200))
    os.utime(bam3, (1000, 100))
    result = list(find_bam_files(str(tmp_path), older_than_timestamp=201))
    assert set(result) == {str(bam2), str(bam3)}
    assert str(bam1) in caplog.text
    assert (str(bam2) in caplog.text) is debug
    assert (str(bam3) in caplog.text) is debug
    assert (str(decoy1) in caplog.text) is debug
    assert (str(decoy2) in caplog.text) is debug


@pytest.mark.parametrize("debug", [True, False])
def test_find_bam_files_exclude(tmp_path: Path, caplog, debug):
    if debug:
        caplog.set_level(logging.DEBUG)
    else:
        caplog.set_level(logging.INFO)
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    decoy1 = subdir / "decoy1.txt"
    decoy2 = tmp_path / "decoy2.txt"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    decoy1.touch()
    decoy2.touch()
    result = list(find_bam_files(
        str(tmp_path),
        # Make the timestamp very big to avoid testing issues.
        older_than_timestamp=math.inf,
        ignore_files=[str(bam1.absolute())]
    ))
    assert set(result) == {str(bam2), str(bam3)}
    assert str(bam1) in caplog.text
    assert "Ignoring" in caplog.text
    assert (str(bam2) in caplog.text) is debug
    assert (str(bam3) in caplog.text) is debug
    assert (str(decoy1) in caplog.text) is debug
    assert (str(decoy2) in caplog.text) is debug


@pytest.mark.parametrize(
    ["cram_version", "write_index", "write_checksum_files", "delete", "use_cli"],
    itertools.product(
        ["3.0", "3.1"],
        [True, False],
        [True, False],
        [True, False],
        [True, False],
    )
)
def test_cram_archiver(
        cram_version,
        write_index,
        write_checksum_files,
        delete,
        use_cli,
        tmp_path,
        caplog,
):
    caplog.set_level(logging.INFO)
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    cram1 = tmp_path / "bam1.cram"
    cram2 = tmp_path / "bam2.cram"
    cram3 = subdir / "bam3.cram"
    bam1_checksum = tmp_path / "bam1.bam.checksum"
    bam2_checksum = tmp_path / "bam2.bam.checksum"
    bam3_checksum = subdir / "bam3.bam.checksum"
    cram1_checksum = tmp_path / "bam1.cram.checksum"
    cram2_checksum = tmp_path / "bam2.cram.checksum"
    cram3_checksum = subdir / "bam3.cram.checksum"
    cram1_index = tmp_path / "bam1.cram.crai"
    cram2_index = tmp_path / "bam2.cram.crai"
    cram3_index = subdir / "bam3.cram.crai"
    decoy1 = subdir / "decoy1.txt"
    decoy2 = tmp_path / "decoy2.txt"
    shutil.copy(TEST_DATA / "GM24385_1.bam", bam1)
    shutil.copy(TEST_DATA / "GM24385_1.bam", bam2)
    shutil.copy(TEST_DATA / "GM24385_1.bam", bam3)
    decoy1.touch()
    decoy2.touch()
    current_time = time.time()
    os.utime(bam1, (current_time, current_time - 10_000))
    os.utime(bam2, (current_time, current_time - 100_000))  # More than 1 day
    os.utime(bam3, (current_time, current_time - 200_000))  # More than 2 days
    if use_cli:
        args = [
            str(tmp_path),
            "--reference", str(TEST_DATA / "NC012920.1.fasta"),
            "--cram-version", cram_version,
            "-vvvv",
            "--minimum-age-days", "1"
        ]
        if not write_index:
            args.append("--dont-write-index")
        if not write_checksum_files:
            args.append("--dont-write-checksums")
        if delete:
            args.append("--delete")
        cram_archiver_main(*args)
    else:
        cram_archiver(
            input_path=str(tmp_path),
            reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
            cram_version=cram_version,
            write_index=write_index,
            write_checksum_files=write_checksum_files,
            minimum_age_days=1,
            delete=delete,
        )
    assert ("WILL BE DELETED" in caplog.text) is delete
    assert (f"deleting BAM file: {bam2}" in caplog.text) is delete
    assert (f"deleting BAM file: {bam3}" in caplog.text) is delete
    assert bam1.exists()
    assert bam2.exists() is not delete
    assert bam3.exists() is not delete
    assert not cram1.exists()  # Not old enough
    assert not cram1_index.exists()
    assert not bam1_checksum.exists()
    assert not cram1_checksum.exists()
    assert cram2.exists()
    assert get_file_cram_version(str(cram2)) == cram_version
    assert cram2_index.exists() is write_index
    assert bam2_checksum.exists() is write_checksum_files
    assert cram2_checksum.exists() is write_checksum_files
    assert cram3.exists()
    assert get_file_cram_version(str(cram3)) == cram_version
    assert cram3_index.exists() is write_index
    assert bam3_checksum.exists() is write_checksum_files
    assert cram3_checksum.exists() is write_checksum_files
    assert ("Total saved size" in caplog.text) is delete
    assert "Found 2 BAM files" in caplog.text
    assert "Total generated CRAM size" in caplog.text


@pytest.mark.parametrize(
    ["cram_version", "write_index", "write_checksum_files", "delete", "use_cli"],
    itertools.product(
        ["3.0", "3.1"],
        [True, False],
        [True, False],
        [True, False],
        [True, False],
    )
)
def test_cram_archiver_single_file(
        cram_version,
        write_index,
        write_checksum_files,
        delete,
        use_cli,
        tmp_path,
        caplog,
):
    caplog.set_level(logging.INFO)
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam = tmp_path / "bam1.bam"
    cram = tmp_path / "bam1.cram"
    bam_checksum = tmp_path / "bam1.bam.checksum"
    cram_checksum = tmp_path / "bam1.cram.checksum"
    cram_index = tmp_path / "bam1.cram.crai"
    shutil.copy(TEST_DATA / "GM24385_1.bam", bam)
    if use_cli:
        args = [
            str(bam),
            "--reference", str(TEST_DATA / "NC012920.1.fasta"),
            "--cram-version", cram_version,
            "-vvvv",
        ]
        if not write_index:
            args.append("--dont-write-index")
        if not write_checksum_files:
            args.append("--dont-write-checksums")
        if delete:
            args.append("--delete")
        cram_archiver_main(*args)
    else:
        cram_archiver(
            input_path=str(bam),
            reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
            cram_version=cram_version,
            write_index=write_index,
            write_checksum_files=write_checksum_files,
            delete=delete,
        )
    assert ("WILL BE DELETED" in caplog.text) is delete
    assert bam.exists() is not delete
    assert cram.exists()
    assert cram_index.exists() is write_index
    assert bam_checksum.exists() is write_checksum_files
    assert cram_checksum.exists() is write_checksum_files
    assert get_file_cram_version(str(cram)) == cram_version


@pytest.mark.parametrize(
    ["delete", "write_checksum_files", "write_index"],
    itertools.product(
        [True, False], [True, False], [True, False]
    )
)
def test_cram_archiver_dry_run(tmp_path, capsys, delete,
                               write_checksum_files,
                               write_index):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    current_time = time.time()
    os.utime(bam1, (current_time, current_time - 10_000))
    os.utime(bam2, (current_time, current_time - 100_000))  # More than 1 day
    os.utime(bam3, (current_time, current_time - 200_000))  # More than 2 days
    cram_archiver(
        input_path=str(tmp_path),
        reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
        minimum_age_days=1,
        dry_run=True,
        delete=delete,
        write_checksum_files=write_checksum_files,
        write_index=write_index,
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == {str(bam2), str(bam3)}
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {str(subdir), str(bam1), str(bam2)}
    assert {str(x) for x in subdir.iterdir()} == {str(bam3)}


def test_cram_archiver_dry_run_exclude_files(tmp_path, capsys):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    cram_archiver(
        input_path=str(tmp_path),
        reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
        dry_run=True,
        ignore_files=[str(bam1), str(bam3)],
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == {str(bam2)}
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {str(subdir), str(bam1), str(bam2)}
    assert {str(x) for x in subdir.iterdir()} == {str(bam3)}


def test_cram_archiver_dry_run_root_dir_excluded(tmp_path, capsys):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    cram_archiver(
        input_path=str(tmp_path),
        reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
        dry_run=True,
        ignore_files=[str(tmp_path)],
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == set()
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {str(subdir), str(bam1),
                                                    str(bam2)}
    assert {str(x) for x in subdir.iterdir()} == {str(bam3)}


def test_cram_archiver_dry_run_root_file_excluded(tmp_path, capsys):
    bam1 = tmp_path / "bam1.bam"
    bam1.touch()
    cram_archiver(
        input_path=str(bam1),
        reference_files=[str(TEST_DATA / "NC012920.1.fasta")],
        dry_run=True,
        ignore_files=[str(bam1)],
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == set()
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {str(bam1)}


def test_cram_archiver_main_dry_run_exclude_files(tmp_path, capsys):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    cram_archiver_main(
        "--reference", str(TEST_DATA / "NC012920.1.fasta"),
        "--dry-run",
        "--exclude", str(bam1),
        "--exclude", str(bam3),
        str(tmp_path)
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == {str(bam2)}
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {str(subdir), str(bam1), str(bam2)}
    assert {str(x) for x in subdir.iterdir()} == {str(bam3)}


def test_cram_archiver_main_dry_run_exclude_files_list(tmp_path, capsys):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    bam1 = tmp_path / "bam1.bam"
    bam2 = tmp_path / "bam2.bam"
    bam3 = subdir / "bam3.bam"
    bam1.touch()
    bam2.touch()
    bam3.touch()
    exclude_list = tmp_path / "exclude.txt"
    exclude_list.write_text(f"{bam1}\n{bam3}")
    cram_archiver_main(
        "--reference", str(TEST_DATA / "NC012920.1.fasta"),
        "--dry-run",
        "--exclude-list", str(exclude_list),
        str(tmp_path)
    )
    stdout, stderr = capsys.readouterr()
    assert set(stdout.splitlines()) == {str(bam2)}
    # Check if no new files are created or that anything is deleted
    assert {str(x) for x in tmp_path.iterdir()} == {
        str(exclude_list), str(subdir), str(bam1), str(bam2)}
    assert {str(x) for x in subdir.iterdir()} == {str(bam3)}


def test_cram_archiver_no_reference_fai(tmp_path):
    bam = tmp_path / "my.bam"
    reference = tmp_path / "reference.fa"
    with pytest.raises(FileNotFoundError) as error:
        cram_archiver(str(bam), [str(reference)])
    error.match("Fasta index")
    error.match(str(reference))
    error.match("not be found")


def test_cram_archiver_no_bam_files(tmp_path, caplog):
    cram_archiver(str(tmp_path), [str(TEST_DATA / "NC012920.1.fasta")])
    assert "No BAM files found." in caplog.text
