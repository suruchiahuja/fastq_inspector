#!/usr/bin/env python3

import sys
import argparse
import gzip
import multiprocessing as mp
from statistics import mean
from collections import defaultdict


def parse_fastq(file_path):
    """Generator that yields (sequence, quality) tuples from a FASTQ file."""
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        while True:
            header = f.readline()
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            if not qual:
                break
            yield seq, qual


def gc_content(seq):
    """Calculate GC content as a percentage."""
    gc = sum(1 for base in seq if base in "GCgc")
    return (gc / len(seq)) * 100 if seq else 0


def quality_scores(qual_str):
    """Convert ASCII-encoded quality string to list of Phred scores."""
    return [ord(char) - 33 for char in qual_str]


def analyze_chunk(chunk):
    """Process a chunk of reads."""
    stats = {
        'gc_contents': [],
        'qualities': [],
        'lengths': []
    }
    for seq, qual in chunk:
        stats['gc_contents'].append(gc_content(seq))
        stats['qualities'].append(mean(quality_scores(qual)))
        stats['lengths'].append(len(seq))
    return stats


def aggregate_stats(results):
    """Aggregate statistics from multiple chunks."""
    all_gc = []
    all_qual = []
    all_lens = []
    for r in results:
        all_gc.extend(r['gc_contents'])
        all_qual.extend(r['qualities'])
        all_lens.extend(r['lengths'])
    return {
        'avg_gc': mean(all_gc),
        'avg_quality': mean(all_qual),
        'avg_length': mean(all_lens),
        'total_reads': len(all_gc)
    }


def chunked_iterator(iterable, chunk_size):
    """Split iterable into chunks."""
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def run_analysis(file_path, chunk_size=10000, workers=4):
    chunks = chunked_iterator(parse_fastq(file_path), chunk_size)
    with mp.Pool(processes=workers) as pool:
        results = pool.map(analyze_chunk, chunks)
    return aggregate_stats(results)


def main():
    parser = argparse.ArgumentParser(description="Analyze FASTQ files for GC content, quality, and read length.")
    parser.add_argument("input", help="Input FASTQ file (can be .gz)")
    parser.add_argument("--chunk-size", type=int, default=10000, help="Number of reads per processing chunk")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel processes")
    args = parser.parse_args()

    stats = run_analysis(args.input, args.chunk_size, args.workers)

    print("\nAnalysis Summary:")
    print(f"  Total Reads:     {stats['total_reads']}")
    print(f"  Avg GC Content:  {stats['avg_gc']:.2f}%")
    print(f"  Avg Read Length: {stats['avg_length']:.1f} bp")
    print(f"  Avg Quality:     {stats['avg_quality']:.2f}\n")


if __name__ == "__main__":
    main()
