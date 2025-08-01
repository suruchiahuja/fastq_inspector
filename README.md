# fastq_inspector
A Python script to parse and analyze FASTQ files for GC content, average quality scores, and read-length statistics, using multiprocessing for performance. 
# fastq-inspector

A Python tool for high-level analysis of FASTQ files. Computes GC content, average Phred quality scores, and read length statistics using multiprocessing for speed.

## Features

- Parses `.fastq` and `.fastq.gz` formats
- Calculates per-read GC content and quality
- Outputs aggregate statistics
- Parallelized for performance
- Clean modular design with CLI

## Example Usage

```bash
python fastq_inspector.py sample.fastq.gz --chunk-size 5000 --workers 8
