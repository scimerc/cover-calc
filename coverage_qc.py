#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

"""
Main script for coverage quality control.
Requires bedtools2.20.1 or higher.

- Calculates the percentage of `coverageRegions` (for single exons and whole trancript aggregates)
  that are equal to or above a minimum coverage threshold, or below such threshold. Takes a single
  BAM file as input, a file with one path to a BAM file for each line or a file with one path to a
  pre-calculated per-base BAM coverage report [NOT IMPLEMENTED].

- Imports the read-filtered coverage per-base output into a BedTool object that is split into two
  one containing base pairs with coverage above the threshold and another containing those with
  coverage below that threshold. Counts the base pairs above threshold in all `coverageRegions` and
  writes one file per aggregation level. Merges all contiguous regions below threshold and writes
  them to another file.
"""

import argparse
import glob
import gzip
import os
import pybedtools
import shutil
import subprocess
import sys

from collections import OrderedDict


class Gene(object):
    '''Genetic object container at various levels of aggregation'''

    aggregationLevels={1: "gene", 2: "transcript", 3: "exon"}

    def __init__(level=1):
        self._level = cls.aggregation(level)

    @classmethod
    def aggregation(cls, level):
        try:
            return cls.aggregationLevels[level]
        except KeyError as e:
            print("unknown aggregation level '{}': falling back to 'gene'".format(level))
            return 'gene'

    @classmethod
    def level(cls, aggregation):
        match_levels = set(l for l,a in cls.aggregationLevels.items() if a == aggregation)
        assert len(match_levels) == 1, 'aggregation level ambiguity. aborting..'
        return match_levels.pop()


def as_str(c, encoding='ascii'):
    '''Convert bytes or string objects to string'''
    if isinstance(c, str):
        return c
    elif isinstance(c, bytes):
        return c.decode(encoding)
    return str(c)


def parse_cmd(argv):
    parser = argparse.ArgumentParser(
        description="""
        Writes for each region the percentage of it that is above the minimum coverage threshold,
        and which contiguous regions are below such threshold. When processing multiple bam files
        using --bampaths, the real sample names are replaced with 'sampleNN' where 'NN' is the index
        of the file in bampaths
        """
    )
    argGroup = parser.add_mutually_exclusive_group(required=True)
    argGroup.add_argument(
        "--bam", dest="bamPath", default=None, help="Path to BAM file"
    )
    argGroup.add_argument(
        "--bampaths", dest="bamPaths", default=None, help="Path to list of BAM files"
    )
    parser.add_argument(
        "--regions",
        dest="coverageRegions",
        required=True,
        help="Path to coverage regions BED file with the regions to calculate coverage on",
    )
    parser.add_argument(
        "--minX",
        default=10,
        dest="minReads",
        required=False,
        type=int,
        help="Minimum acceptable number of reads",
    )
    parser.add_argument(
        "--minMQ",
        default=20,
        dest="minMQ",
        required=False,
        type=int,
        help="Minimum mapping quality of reads to be counted",
    )
    parser.add_argument(
        "--outputdir",
        default=".",
        dest="outputDir",
        required=False,
        help="Path to output directory",
    )
    parser.add_argument(
        "--samplename",
        default="nosamplename",
        dest="sampleName",
        required=False,
        help="Name of sample",
    )
    parser.add_argument(
        "--Nthreads",
        default=2,
        dest="Nthreads",
        required=False,
        type=int,
        help="Number of computing threads for mosdepth runs",
    )
    parser.add_argument(
        "--excludedRegions",
        dest="excludedRegions",
        required=False,
        help="A bed file with regions to exclude from coverage check",
    )
    parser.add_argument(
        "--extraFilters",
        default="-F 0",
        dest="extraReadFilters",
        required=False,
        help="Additional read filters to be applied (comma-separated)",
    )
    parser.add_argument(
        "--perBP",
        action="store_true",
        default=False,
        dest="outputCoveragePerBP",
        required=False,
        help="Keep intermediate coverage per base-pair file (large file)",
    )

    parsed_args = parser.parse_args(argv)

    if parsed_args.minReads < 1:
        print("Warning: trivial minimum reads requested [%d], using 1." % parsed_args.minReads)
        parsed_args.minReads = 1

    return parsed_args


def run_mosdepth_on(bamPath, sampleName, outputDir, Nthreads=1, extraReadFilters="", minMQ=20):
    """
    Runs mosdepth

    extraReadFilters must be specified as a comma-separated list.
    """

    coverageOutput = os.path.join(outputDir, sampleName + '_coverage')
    tempcoverageOutput = os.path.join(outputDir, sampleName + '_tmp_coverage')

    mosdepthexec = os.environ.get('MOSDEPTH_PATH', 'mosdepth')
    cmd = mosdepthexec + " --fast-mode -t {N} -Q {Q} {flags} {out} {bam}".format(
        N=Nthreads, Q=minMQ, flags=extraReadFilters, out=tempcoverageOutput, bam=bamPath
    )
    returnCode = subprocess.call(cmd.split(), stdout=sys.stderr)
    assert returnCode == 0, (
        "Something went wrong in the mosdepth call [error %d]; make sure `mosdepth` is in your\n"
        "executable search path or point the environmental variable `MOSDEPTH_PATH` to the\n"
        "location of the `mosdepth` executable."
    ) % returnCode

    for tempname in glob.glob(tempcoverageOutput + '*'):
        finalname = tempname.replace(tempcoverageOutput, coverageOutput, 1)
        print("renaming temporary file '%s' to '%s'.." % (tempname, finalname))
        os.rename(tempname, finalname)

    return coverageOutput


def split_bed_on_score(inFile, aboveFile, belowFile, minPass):
    """Writes BED records with score >= minPass to aboveFile, else belowFile"""
    for line in inFile:
        if line.strip().isspace():
            continue
        if int(as_str(line).split("\t")[3]) >= minPass:
            aboveFile.write(as_str(line))
        else:
            belowFile.write(as_str(line))


def write_aggregated_report(aboveBED, regionBED, level, thresh, sampleName, outputDir):
    """
    Writes a coverage report at the given aggregation level

    Assumes all regions at the lowest level of aggregation are given in the BED object and outputs
    them in the same order. NOTE: the `count` member of the `Interval` objects yielded by iterating
    over `BedTool` objects isn't clearly documented at the time of writing but seems to hold an
    integer coercion of whatever is in the "last" field of the corresponding `BedTool` record.
    """

    coverageStats = {}
    outputFilePath = os.path.join(
        outputDir, "{}_{}_coverage_{}x.tsv".format(sampleName, Gene.aggregation(level), thresh)
    )
    aboveBEDmerged = aboveBED.merge()  # merge contiguous aboveBED intervals
    # Find all regionBED's overlaps with aboveBED and add their widths up
    # NOTE: the assumption here is, reasonably enough, that the latter are _disjoint_ intervals.
    #       regions with 0 overlaps are also reported (see BEDTools '-wao' command line flag).
    Na = regionBED.field_count()
    Nb = aboveBEDmerged.field_count()
    Ntotal = Na + Nb + 1
    grouping = [n+1 for n in range(Na)]
    regionCoverageAboveBED = \
        regionBED.intersect(aboveBEDmerged, wao=True).groupby(g=grouping, c=Ntotal, o='sum')
    with open(outputFilePath, "w") as fh:
        fh.write("#sampleName\t" + "\t".join((Gene.aggregation(l) for l in range(1, level+1))))
        fh.write("\tfraction covered\tnucleotides covered\tnucleotides in total\n")
        covered_lengths = OrderedDict()
        total_lengths = OrderedDict()

        for record in regionCoverageAboveBED:
            # regions are named (fourth BED file field) 'gene__hgnc__transcript__exon'
            name = record.name.split("__")[:level+1]
            name = "__".join(name)
            covered_lengths[name] = covered_lengths.get(name, 0) + record.count
            total_lengths[name] = total_lengths.get(name, 0) + record.length

        namesOutput = set()
        for name in covered_lengths:
            if name in namesOutput:
                continue  # Aggregate already output
            try:
                f = 100. * (covered_lengths[name] / (total_lengths[name] * 1.0))
            except ZeroDivisionError as e:
                f = 0.0
            namesOutput.add(name)
            nameParts = name.split("__")
            geneSymbol = nameParts.pop(0)
            geneId = nameParts.pop(0)
            coverageStats[name] = {
              'name': geneSymbol, 'HGNC-ID': geneId, 'descriptors': nameParts, 'f': f
            }
            fh.write("\t".join([sampleName, geneSymbol, geneId] + nameParts))
            fh.write("\t%.1f%%\t%d\t%d\n" % (f, covered_lengths[name], total_lengths[name]))

    return coverageStats


def write_low_coverage_regions_report(belowBED, regionBED, thresh, sampleName, outputDir):
    """Writes a BED6+3 with any contiguous regions below threshold"""
    # Initialize an empty tuple, to create an empty output file by default
    belowMergedWithNamesBED = ()
    if belowBED.count() > 0:
        # Merge contiguous regions retaining the minimum coverage (field 4)
        belowMergedBED = belowBED.merge(c=4, o="min")
        # Intersect with `regionBED` keeping the region names
        belowMergedWithNamesBED = regionBED.intersect(belowMergedBED, wa=True, wb=True)
    outputFilePath = os.path.join(outputDir, "{}_coverage_below_{}x.bed".format(sampleName, thresh))
    with open(outputFilePath, "w") as output:
        for record in belowMergedWithNamesBED:
            # Pybedtools knows what 'name', 'score', 'strand' are
            gene, hgnc, transcript, exon = as_str(record.name).split("__")
            output.write(
                "{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}\t{sampleName}\t{gte}\n".format(
                    chrom=record.chrom,
                    start=record.start,
                    stop=record.stop,
                    name=record.name,
                    score=record.count,
                    strand=record.strand,
                    sampleName=sampleName,
                    gte="\t".join((gene, transcript, exon))
                )
            )


def write_summary_coverage_report(coverageStats, outputDir, minReads):
    """
    Writes a summary coverage report at three aggregation levels (gene, transcript, exon)
    Max, Mean, Median, Min coverage statistics are computed across all samples (top level keys)
    """

    for a in coverageStats:  # iterate over aggregation levels
        level = Gene.level(a)
        outputFilePath = os.path.join(outputDir, "coverage_summary_%s_%dx.tsv" % (a, minReads))
        with open(outputFilePath, "w") as fh:
            print("aggregation level '%s'.." % level)
            fh.write("\t".join((Gene.aggregation(l) for l in range(1, level+1))))
            fh.write("\tMax\tMin\tMean\tMedian\n")
            geneList = {g for b in coverageStats[a] for g in coverageStats[a][b]}
            for gene in geneList:
                print("  processing gene '%s'.." % gene)
                fList = [
                    coverageStats[a][b][g].get('f', 0) \
                        for b in coverageStats[a] for g in coverageStats[a][b] if g == gene
                ]
                geneSymbol = {
                    coverageStats[a][b][g].get('name', 0) \
                        for b in coverageStats[a] for g in coverageStats[a][b] if g == gene
                }
                geneDescriptors = {
                    tuple(d for d in coverageStats[a][b][g].get('descriptors', [])) \
                        for b in coverageStats[a] for g in coverageStats[a][b] if g == gene
                }
                assert len(geneSymbol) == 1, 'Gene names mismatch: aborting..'
                assert len(geneDescriptors) == 1, 'Gene descriptors mismatch: aborting..'
                max_coverage = max(fList)
                min_coverage = min(fList)
                mean_coverage = sum(fList) / len(fList)
                median_coverage = sorted(fList)[int((len(fList) + 1) / 2)]
                fh.write("\t".join([geneSymbol.pop()] + [d for d in geneDescriptors.pop()]))
                fh.write("\t%.1f%%\t%.1f%%\t%.1f%%\t%.1f%%\n" % (
                    max_coverage, min_coverage, mean_coverage, median_coverage
                ))


def main(argv=None):
    argv = argv or sys.argv[1:]
    args = parse_cmd(argv)

    pybedtools.set_tempdir(args.outputDir)

    filteredBEDFile = os.path.join(args.outputDir, 'filtered_regions.bed')

    if args.excludedRegions:
        # assume a single bed file:
        exclusions = pybedtools.BedTool(args.excludedRegions)
        target = pybedtools.BedTool(args.coverageRegions)
        _ = target.subtract(exclusions, output=filteredBEDFile)
    else:
        shutil.copy(args.coverageRegions, filteredBEDFile)

    regionBED = pybedtools.BedTool(filteredBEDFile)

    bamPaths = [args.bamPath]
    sampleNames = [args.sampleName]
    
    if args.bamPaths:
        bamPaths = []
        sampleNames = []
        with open(args.bamPaths) as bamPathsFile:
            idx = 0
            for bamPath in bamPathsFile:
                if bamPath.strip().startswith("#") or len(bamPath.strip()) < 1:
                    continue
                idx += 1
                bamPaths.append(glob.glob(bamPath.strip()).pop())
                sampleNames.append("sample%02d" % idx)

    assert len(bamPaths) == len(sampleNames), '"That was definitely the same cat."'

    coverageStats = {a: {} for a in Gene.aggregationLevels.values()}

    for k in range(len(bamPaths)):
        print("processing '%s'.." % sampleNames[k])

        # Run mosdepth calculation
        coverageOutput = run_mosdepth_on(
            bamPaths[k],
            sampleNames[k],
            args.outputDir,
            args.Nthreads,
            args.extraReadFilters,
            args.minMQ
        )
        # Split coverageOutput into above/below threshold
        aboveThresholdPath = os.path.join(args.outputDir, sampleNames[k] + '_above.bed')
        belowThresholdPath = os.path.join(args.outputDir, sampleNames[k] + '_below.bed')
        with gzip.open(coverageOutput + ".per-base.bed.gz") as fh:
            with open(aboveThresholdPath, "w") as gh:
                with open(belowThresholdPath, "w") as hh:
                    split_bed_on_score(fh, gh, hh, args.minReads)
        aboveBED = pybedtools.BedTool(aboveThresholdPath)
        belowBED = pybedtools.BedTool(belowThresholdPath)
        # Write aggregated reports for all levels
        for level, aggregation in Gene.aggregationLevels.items():
            coverageStats[aggregation][sampleNames[k]] = \
                write_aggregated_report(
                  aboveBED, regionBED, level, args.minReads, sampleNames[k], args.outputDir
                )
            write_low_coverage_regions_report(
              belowBED, regionBED, args.minReads, sampleNames[k], args.outputDir
            )
        # Clean up temporary BED files
        os.remove(aboveThresholdPath)
        os.remove(belowThresholdPath)
        # Save per-base coverage, if requested
        temppattern = coverageOutput + '.per-base.bed.gz'
        finalpattern = os.path.join(args.outputDir, "%s_coveragePerBP.bed.gz" % sampleNames[k])
        for tempname in glob.glob(temppattern + '*'):
            if args.outputCoveragePerBP:
                finalname = tempname.replace(temppattern, finalpattern, 1)
                print("renaming file '%s' to '%s'.." % (tempname, finalname))
                os.rename(tempname, finalname)
            else:
                os.remove(tempname)

    os.remove(filteredBEDFile)

    if len(bamPaths) > 1:
        write_summary_coverage_report(coverageStats, args.outputDir, args.minReads)

    return 0


if __name__ == "__main__":
    sys.exit(main())

