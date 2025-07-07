==========
Changelog
==========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

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
