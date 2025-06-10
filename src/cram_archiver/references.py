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

"""
This module contains utilities to handle different references within one
application.
"""
import gzip
import io
import struct
import subprocess
import sys
from typing import TextIO


class ReferenceID:
    """
    A wrapper class for a string which contains contigs in tabular format
    with the first column the contig name and the second the contig length.

    The order of the contigs and their lengths makes for a unique string which
    can be used as a unique ID for different versions of genome builds.
    """
    _id: str

    def __init__(self, reference_id):
        self._id = reference_id

    def __eq__(self, other):
        if not isinstance(other, ReferenceID):
            raise TypeError(
                f"Can only compare with instances of ReferenceID, got "
                f"{other.__class__.__name__}")
        return self._id == other._id

    def __hash__(self):
        return hash(self._id)

    @property
    def id(self) -> str:
        return self._id

    @classmethod
    def from_file(cls, file: str):
        """
        Read the reference ID from the file. Auto-detects format based on the
        magic bytes.
        """
        with open(file, "rb") as filehandle:  # type: io.BufferedReader
            if filehandle.peek(2)[:2] == b"\x1f\x8b":
                # Gzip magic detected.
                filehandle = gzip.open(filehandle, "rb")  # type: ignore
            magic = filehandle.peek(100)
            text_handle = io.TextIOWrapper(filehandle)
            if magic.startswith(b"CRAM"):
                return cls._from_alignment_file(file)
            elif magic.startswith(b"@HD"):
                return cls._from_sam_header(text_handle)
            elif magic.startswith(b"BAM\x01"):
                # Skip magic 4 bytes. The following 4 bytes are the length of
                # the header in plaintext.
                l_text, = struct.unpack("<xxxxI", filehandle.read(8))
                # Only ASCII is allowed in the SAM/BAM header.
                text = filehandle.read(l_text).decode("ascii")
                return cls._from_sam_header(io.StringIO(text))
            # FAI detection
            first_line = magic.splitlines()[0]
            if first_line.count(b"\t") == 4:
                first_line_text = first_line.decode("ascii")
                contig, contig_length, start, nucs_per_line, line_width = \
                    first_line_text.split("\t")
                if (contig_length.isdecimal() and
                        start.isdecimal() and
                        nucs_per_line.isdecimal() and
                        line_width.isdecimal()):
                    # This is a fasta index.
                    return cls._from_fasta_index(text_handle)
            raise NotImplementedError(f"file with magic {magic[:10]!r} not "
                                      f"implemented.")

    @classmethod
    def _from_alignment_file(cls, alignment_file):
        result = subprocess.run(["samtools", "view", "-H", alignment_file],
                                stdout=subprocess.PIPE, check=True)
        return cls._from_sam_header(io.StringIO(result.stdout.decode('ascii')))

    @classmethod
    def _from_sam_header(cls, filehandle: TextIO):
        """Parses the contigs and lengths from SAM header @SQ lines."""
        id_build = io.StringIO()
        for line in filehandle:
            if not line.startswith("@"):
                break
            if line.startswith("@SQ"):
                line_parts = line.strip().split("\t")
                # line_parts from index 1 to remove @SQ part.
                info_parts = (part.split(":", maxsplit=1) for part in line_parts[1:])
                info_dict = dict(info_parts)
                contig = info_dict["SN"]
                contig_length = info_dict["LN"]
                id_build.write(f"{contig}\t{contig_length}\n")
        return cls(id_build.getvalue())

    @classmethod
    def _from_fasta_index(cls, filehandle: TextIO):
        """Parses the contigs and lengths from the fasta index tabular format."""
        id_build = io.StringIO()
        for line in filehandle:
            contig, contig_length, *rest = line.split()
            id_build.write(f"{contig}\t{contig_length}\n")
        return cls(id_build.getvalue())


if __name__ == "__main__":
    file = sys.argv[1]
    print(ReferenceID.from_file(file).id)
