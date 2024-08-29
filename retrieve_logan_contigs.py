#!/usr/bin/env python

import argparse
import io
import sys
from typing import Iterator

import boto3
import zstandard as zstd
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError


class Sequence:
    def __init__(self, header: str, seq: str):
        self._header = header
        self._seq = seq.encode("ascii")

    @property
    def header(self):
        return self._header

    @property
    def accession(self):
        return self._header.split()[0]

    @property
    def seq(self):
        return self._seq.decode()

    @property
    def seq_ascii(self):
        return self.seq.upper().encode("ascii")

    def count(self, substring: str):
        return self.seq.count(substring)

    def rc(self):
        tab = self.seq.maketrans("ACTGNactgn", "TGACNtgacn")
        return Sequence(self.header, self.seq.translate(tab)[::-1])

    def has_dtr(self, min_length: int = 30):
        substring = self.seq.casefold()[:min_length]
        pos = self.seq.casefold().rfind(substring)
        if pos < len(self) / 2:
            return False, 0
        substring = self.seq.casefold()[pos:]
        return self.seq.casefold()[: len(substring)] == substring, len(substring)

    def has_itr(self, min_len: int = 30):
        rev = self.rc().seq
        if self.seq.casefold()[:min_len] == rev.casefold()[:min_len]:
            i = min_len + 1
            while self.seq.casefold()[:i] == rev.casefold()[:i] and i <= len(self) // 2:
                i += 1
            return True, i - 1
        else:
            return False, 0

    def fix_circle(self, k: int = 31):
        kmers = set()
        for i in range(len(self._seq) - k + 1):
            kmer = self._seq[i : i + k]
            if kmer in kmers:
                return self[: i + k - 1]
            kmers.add(kmer)
        return self

    def __str__(self):
        return f">{self.header}\n{self.seq}\n"

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, k: int):
        return Sequence(self.header, self.seq[k])

    def __eq__(self, other: object):
        if other.__class__ is self.__class__:
            return self.seq.casefold() == other.seq.casefold()
        elif other.__class__ is str:
            return self.seq.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self):
        return hash(self.seq.casefold())


def read_fasta_from_s3_stream(stream_reader: io.IOBase) -> Iterator[Sequence]:
    buffer = ""
    current_name = None
    current_seq = []
    chunk_size = 4096
    while True:
        chunk = stream_reader.read(chunk_size)
        if not chunk:
            break
        buffer += chunk.decode("utf-8", errors="ignore")
        lines = buffer.split("\n")
        for line in lines[:-1]:  # Process all complete lines
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    seq = "".join(current_seq)
                    if seq:
                        yield Sequence(current_name, seq)
                current_name = line[1:]
                current_seq = []
            elif current_name:
                current_seq.append(line)
        buffer = lines[-1]  # Keep the last incomplete line in the buffer
    # Process the last sequence if any
    if current_name:
        seq = "".join(current_seq) + buffer
        if seq:
            yield Sequence(current_name, seq)


def process_s3_fasta(
    accession,
    bucket_name="logan-pub",
    min_sequence_length: int = 0,
    circles_only: bool = True,
    fix_circles: bool = True,
    output=None,
):
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    object_key = f"c/{accession}/{accession}.contigs.fa.zst"

    try:
        response = s3.get_object(Bucket=bucket_name, Key=object_key)
        stream = response["Body"]
        dctx = zstd.ZstdDecompressor()
        stream_reader = dctx.stream_reader(stream)
        if output:
            with open(output, "w") as fout:
                for sequence in read_fasta_from_s3_stream(stream_reader):
                    has_dtr = sequence.has_dtr()[0]
                    if circles_only and not has_dtr:
                        continue
                    if fix_circles and has_dtr:
                        sequence = sequence.fix_circle()
                    if len(sequence) >= min_sequence_length:
                        fout.write(str(sequence))
        else:
            for sequence in read_fasta_from_s3_stream(stream_reader):
                has_dtr = sequence.has_dtr()[0]
                if circles_only and not has_dtr:
                    continue
                if fix_circles and has_dtr:
                    sequence = sequence.fix_circle()
                if len(sequence) >= min_sequence_length:
                    print(str(sequence), end="")
    except ClientError as e:
        if e.response["Error"]["Code"] == "NoSuchKey":
            print(
                f"Error: The file for accession {accession} does not exist in the S3 bucket."
            )
        else:
            print(f"Error accessing S3: {str(e)}")
    except zstd.ZstdError as e:
        print(f"Error decompressing the file: {str(e)}")
    except Exception as e:
        print(f"Unexpected error processing {accession}: {str(e)}")
    finally:
        if "stream" in locals():
            stream.close()


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve and process Logan contig data from S3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument("accession", help="Sample accession")
    parser.add_argument(
        "--min-sequence-length",
        type=int,
        default=0,
        help="Minimum sequence length to output",
    )
    parser.add_argument(
        "--circles-only",
        action="store_true",
        help="Output only circular sequences",
    )
    parser.add_argument(
        "--fix-circles",
        action="store_true",
        help="Fix circular sequences by removing duplicates",
    )
    parser.add_argument(
        "--output", help="Output file path (if not specified, print to stdout)"
    )

    # Print help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    process_s3_fasta(
        args.accession,
        min_sequence_length=args.min_sequence_length,
        circles_only=args.circles_only,
        fix_circles=args.fix_circles,
        output=args.output,
    )


if __name__ == "__main__":
    main()
