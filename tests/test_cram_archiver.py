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
import os.path
import shutil
from pathlib import Path

from cram_archiver import (
    checksum,
    convert_to_cram,
    convert_to_cram_and_check,
    handle_file_age,
    strip_comments_from_checksum,
)
from cram_archiver.references import ReferenceID

import pytest

TEST_DATA = Path(__file__).parent / "data"


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
