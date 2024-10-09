# `retrieve-logan-contigs`

This script downloads contigs from a specified Logan assembly. It includes options to:

- Set a minimum sequence length.
- Filter for circular sequences.
- Fix concatemeric circular sequences.

## Usage

```
usage: retrieve-logan-contigs [-h] [--min-sequence-length MIN_SEQUENCE_LENGTH]
                              [--circles-only] [--fix-circles]
                              [--output OUTPUT]
                              accession

Retrieve and process Logan contig data from S3.

positional arguments:
  accession             Sample accession

options:
  -h, --help            show this help message and exit
  --min-sequence-length MIN_SEQUENCE_LENGTH
                        Minimum sequence length to output (default: 0)
  --circles-only        Output only circular sequences (default: False)
  --fix-circles         Fix circular sequences by removing duplicates
                        (default: False)
  --output OUTPUT       Output file path (if not specified, print to stdout)
                        (default: None)
```