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

+ Automated recursive discovery of all ``.bam`` files in a directory. Symbolic
  links are ignored.
+ Multiple reference support. Cram-archiver loads in fasta indexed references
  and checks that the appropriate BAM file is matched to the appropriate
  reference using the contig and length information in the BAM header.
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

Quickstart
==========

Converting a single BAM file::

    cram-archiver -r my_reference.fasta my.bam

This will create ``my.cram``, ``my.cram.crai``, ``my.cram.checksum`` and
``my.bam.checksum``. Checksum file creation can be turned of with
``--dont-write-checksums``. The checkums will still be checked, just not
written to disk.

Archiving a directory with BAMs, but only BAMs that have a lost modified time
older than 30 days. Also, there are hg19 and hg38 BAM files in the directory.::

    cram-archiver --reference hg19.fasta --reference hg38.fasta --minimum-age-days 30 my_directory

If the ``--delete`` flag is added, all the converted BAM files will be deleted
and just the CRAM files remain. This only happens when the conversion is
successful and the checksums match.

Usage
=====

    usage: cram-archiver [-h] -r REFERENCE [-t THREADS] [-d MINIMUM_AGE_DAYS]
                         [--delete] [--cram-version CRAM_VERSION]
                         [--exclude EXCLUDE] [--exclude-list PATH]
                         [--dont-write-checksums] [--dont-write-index] [--dry-run]
                         [-v] [-q] [--version]
                         PATH

    positional arguments:
      PATH                  Path to BAM file or directory to be recursively
                            searched.

    options:
      -h, --help            show this help message and exit
      -r REFERENCE, --reference REFERENCE
                            Reference to be used for CRAM conversion. Can be used
                            multiple times. Reference will be checked with the BAM
                            file.
      -t THREADS, --threads THREADS
                            The number of threads used for conversion and
                            checksumming.Default: 1.
      -d MINIMUM_AGE_DAYS, --minimum-age-days MINIMUM_AGE_DAYS
                            The minimum last modification of the BAM file in days
                            prior. This assumes the system clock timezone matches
                            that of the file while also assuming that every day
                            has 24x60x60 seconds. Default 0.
      --delete              Delete BAM files after successful conversion.
      --cram-version CRAM_VERSION
                            CRAM version to use for CRAM conversion. Default: 3.0.
      --exclude EXCLUDE     Exclude file or directory from conversion. Can be
                            supplied multiple times.
      --exclude-list PATH   Supply a newline-separated file with files and
                            directories to exclude.
      --dont-write-checksums
                            Do not store samtools checksum output on disk.
      --dont-write-index    Do not write index files for CRAM files.
      --dry-run             Print the paths of the to be archived BAM files.
                            Perform no actions.
      -v, --verbose         Display more logging information.
      -q, --quiet           Display less logging information.
      --version             show program's version number and exit

On CRAM format settings
=======================
Cram-archiver uses version 3.0 of the CRAM standard by default. The reason
for this is that CRAM version 3.0 is better supported than version 3.1.
CRAM version 3.1 comes with newer codecs and is able to achieve smaller
file sizes because of that. For more information checkout the
`article on advances in CRAM by James Bonfield
<https://doi.org/10.1093/bioinformatics/btac010>`_.

Cram archiver uses the CRAM default presets. CRAM has some presets: fast,
normal, small and archive. However, the size differences between normal and
archive are quite small (less than 6% smaller in our tests). On top of that,
the memory requirements rise steeply especially on very long read alignments of
ONT data.

Acknowledgements
================
A huge thank you to James Bonfield (`@jkbonfield <https://github.com/jkbonfield>`_)
for providing a lot of information and background about CRAM and its tooling.
This was invaluable for creating this project. James Bonfield has also spent
a lot of effort into making CRAM the very usable format it is today for which
we are very grateful.
