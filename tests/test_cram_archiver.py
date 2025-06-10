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
from pathlib import Path

from cram_archiver import checksum

TEST_DATA = Path(__file__).parent / "data"


def test_checksum():
    bam = str(TEST_DATA / "GM24385_1.bam")
    sum = checksum(bam, str(TEST_DATA / "NC012920.1.fasta"))
    assert (
        "all        all           178  06d90a04  526c1585  3e8f585f  "
        "41c54f43  268b6cbb  1ff11ae8  35f542f6  56a9ac23\n"
        "GM24385_fastq-lib1-fastq all           178  06d90a04  526c1585  "
        "3e8f585f  41c54f43  268b6cbb  1ff11ae8  35f542f6  56a9ac23\n"
    ) in sum
