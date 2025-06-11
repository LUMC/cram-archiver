======================
cram-archiver
======================

A samtools wrapper for CRAM conversion automation.

Introduction
============
cram-archiver was written to help with an archival task where a substantial
volume of BAM files needed to be converted to CRAM in order to save
disk space.

Features:

+ Automated recursive discovery of all ``.bam`` files in a directory.
+ Performs CRAM conversion using ``samtools view``.
+ Performs ``samtools checksum --all`` on the BAM and CRAM file and checks
  if the checksum matches.
+ On by default: writes checksum files for manual verification.
+ On by default: writes CRAM indexes.
+ Optional: deletes BAM file after conversion.
+ Optional: Set a minimum age in days for the BAM file's last modified time.
  If the file is "older" than the set number of days, the file will be
  converted.

Caveats
=======
CRAM was never intended and built as a "pure" archival format with bit-for-bit
reproducibility. As a result
it is impossible to get the original BAM file back from a pure CRAM file.
There are several reasons for this:

+ BAM files are by definition always bgzip compressed using the DEFLATE
  algorithm. Differences in the DEFLATE algorithm implementation can cause
  different outputs.
+ When converting a BAM and its derived CRAM to SAM the two SAMs can have
  differences too:

  + MD and NM tags are not stored in CRAM files but always calculated on the
    fly when decoding. If the MD and NM flags were not present in the
    original BAM, this can cause differences.
  + The order of tags might be different.
  + ``M``, ``=`` and ``X`` in CIGAR strings. ``=`` means that the nucleotide
    is the same at this position. ``X`` means a mismatch at this position.
    ``M`` means that the position matches (no indels), but gives no information
    whether it is ``X`` or ``=``. Since ``X`` and ``=`` can be derived from
    the sequence, the extra information is redundant and CRAM stores everything
    as ``M``. This can give rise to differences.
  + Redundant information in BAM files such as unaligned reads with MAPQ values
    or CIGAR strings. This does not get stored.
  + Errors, such as wrong mate pair information. Some of it may be fixed during
    the CRAM conversion.

To assure the CRAM file is "functionally the same" as the BAM file, the
``samtools checksum`` tool with the ``--all`` flag is run. For more information
about comparing BAM and CRAM checkout `the discussion here
<https://github.com/samtools/samtools/issues/2212>`_.

Acknowledgements
================
A huge thank you to James Bonfield (`@jkbonfield <https://github.com/jkbonfield>`_)
for providing a lot of information and background about CRAM and its tooling.
This was invaluable for creating this project. James Bonfield has also spent
a lot of effort into making CRAM the very usable format it is today for which
we are very grateful.
