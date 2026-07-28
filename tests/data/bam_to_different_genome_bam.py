#!/usr/bin/env python3

import gzip
import io
import struct
import subprocess

if __name__ == "__main__":
    new_bam = io.BytesIO()
    with gzip.open("GM24385_1.bam", "rb") as bam:
        assert b"BAM\01" == bam.read(4)
        new_bam.write(b"BAM\01")
        l_text, = struct.unpack("<I", bam.read(4))
        text = bam.read(l_text)
        lines = text.split(b"\n")
        only_header_line = lines[0]
        assert only_header_line.startswith(b"@HD\t")
        new_header_text = only_header_line + b"\n"
        new_bam.write(struct.pack("<I", len(new_header_text)))
        new_bam.write(new_header_text)

        n_ref_bytes = bam.read(4)
        new_bam.write(n_ref_bytes)
        n_ref, = struct.unpack("<I", n_ref_bytes)
        for i in range(n_ref):
            name_length_bytes = bam.read(4)
            new_bam.write(name_length_bytes)
            name_length, = struct.unpack("<I", name_length_bytes)
            new_bam.write(bam.read(name_length))
            length_bytes = bam.read(4)
            chrom_length, = struct.unpack("<I", length_bytes)
            new_bam.write(struct.pack("<I", chrom_length + 1))

    # BAM has been read. Now convert the bytes into a proper bgzipped bam
    with open("GM24385_1_unknown_reference.bam", "wb") as result_bam:
        process = subprocess.Popen(
            ["bgzip", "--stdout", "--compress-level", "9"],
            stdin=subprocess.PIPE, stdout=result_bam)
    stdout, stderr = process.communicate(new_bam.getvalue())
    if process.returncode != 0:
        raise RuntimeError(stderr)
