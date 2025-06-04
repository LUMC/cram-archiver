"""
This module contains utilities to handle different references within one
application.
"""
import gzip
import io
import struct
from typing import BinaryIO, TextIO


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
                return cls._from_cram(filehandle)
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

    def _from_cram(cls, filehandle: BinaryIO):
        format = filehandle.read(4)
        if format != b"CRAM":
            raise ValueError(f"Invalid CRAM magic {format}")
        major, minor = struct.unpack("BB", filehandle.read(2))
        if major != 3:
            raise ValueError(
                f"Only CRAM versions in the 3.x range are supported, "
                f"got {major}.{minor}")
        file_id = filehandle.read(20)


    @classmethod
    def _from_fasta_index(cls, filehandle: TextIO):
        """Parses the contigs and lengths from the fasta index tabular format."""
        id_build = io.StringIO()
        for line in filehandle:
            contig, contig_length, *rest = line.split()
            id_build.write(f"{contig}\t{contig_length}\n")
        return cls(id_build.getvalue())


def _itf_8_from_stream(stream: BinaryIO) -> int:
    """Get ITF-8 value from stream. Based on the htslib code."""

    start_to_nbytes = [
        0,  # 0b0000xxxx
        0,  # 0b0001xxxx
        0,  # 0b0010xxxx
        0,  # 0b0011xxxx
        0,  # 0b0100xxxx
        0,  # 0b0101xxxx
        0,  # 0b0110xxxx
        0,  # 0b0111xxxx
        1,  # 0b1000xxxx
        1,  # 0b1001xxxx
        1,  # 0b1010xxxx
        1,  # 0b1011xxxx
        2,  # 0b1100xxxx
        2,  # 0b1101xxxx
        3,  # 0b1110xxxx
        4,  # 0b1111xxxx
    ]

    value_masks = [
        0b0111_1111,
        0b0011_1111,
        0b0001_1111,
        0b0000_1111,
        0b0000_1111,
    ]

    b = stream.read(1)[0]
    nbytes = start_to_nbytes[b >> 4]
    value = b & value_masks[nbytes]
    if nbytes == 0:
        return value
    elif nbytes == 1:
        return value << 8 | stream.read(1)[0]
    elif nbytes == 2:
        c = stream.read(2)
        return value << 16 | c[0] << 8 | c[1]
    elif nbytes == 3:
        c = stream.read(3)
        return value << 24 | c[0] << 16 | c[1] << 8 | c[2]
    c = stream.read(4)
    value = value << 28 | c[0] << 20 | c[1] << 12 | c[2] << 2 | c[3] & 0b0000_1111
    if value & 0x80_00_00_00:
        # Perform two's complement
        return - (((~value) + 1) & 0xffff_ffff)
    return value
