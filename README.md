Coverage calculator
===================

This repository contains a script to calculate coverage over a BAM file or a list of such files.


# Installing

It is advisable to use `pip` to install the necessary packages. Change into the repository
directory and run

`pip install -r requirements.txt`


# Running

(see `script.sh` for a usage example)


# Output

the script writes a number of files at the given location. these contain sample-wise coverage
statistics as well as cumulative summary statistics. the `example_output.png` in this repository
shows the list of output files for a run of `script.sh` assuming `bamlist` contains three sample
bam files.

the `sample##_` files contain individual coverage statistics, the `coverage_*` files contain
cumulative summary coverage statistics:

- <pre>`*_coverage.mosdepth.*    `</pre> mosdepth output files;
- <pre>`*_coveragePerBP.bed.gz*   `</pre> mosdepth extra per base-pair output statistics;
- <pre>`*_coverageReport*         `</pre> sample-wise coverage stats at three aggregation levels;
- <pre>`*_lowCoverage*            `</pre> list of sample's regions with reads below threshold;
- <pre>`coverage_*                `</pre> summary of the sample-wise `*_coverageReport*` files.

