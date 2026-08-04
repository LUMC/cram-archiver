==========
Changelog
==========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

develop
------------------
+ Several improvements to the ``--exclude-list`` file:
    + Empty lines are now ignored.
    + Comments (starting with ``#``) are ignored.
    + Relative paths in the file are now resolved relative to the directory
      of the file.
+ Add a ``processes`` flag to run multiple samtools processes in parallel.
+ Add a ``--exclude-extension`` option to exclude extensions. Useful for
  example for DRAGEN ``.repeats.bam`` files, which are very small.
+ Fix bug where cram-archiver would crash without proceeding to process other
  BAM files when an unknown reference was encountered.
+ Fix issue with references not being read from BAM's reference binary blob if
  ``@SQ`` header lines are not present.

1.1.0
------------------
+ Filesizes are logged as well as space savings. Total filesizes and space
  savings are logged at the end of the program run.
+ Do not follow symbolic links anymore.
+ Add ``--exclude`` and ``--exclude-list`` command line options to allow
  for automatically excluding files.

1.0.0
------------------
+ Create a samtools wrapper that can convert BAM files into CRAM and delete
  them when the checksum tests pass.
