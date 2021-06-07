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

- `*_coverage.mosdepth.*   &nbsp; &nbsp; &nbsp; `
  mosdepth output files;
- `*_coveragePerBP.bed.gz* &nbsp; &nbsp; &nbsp; `
  mosdepth extra per base-pair output statistics;
- `*_coverageReport*       &nbsp; &nbsp; &nbsp; `
  sample-wise coverage stats at three aggregation levels;
- `*_lowCoverage*          &nbsp; &nbsp; &nbsp; `
  list of sample's regions with reads below threshold;
- `coverage_*              &nbsp; &nbsp; &nbsp; `
  summary of the sample-wise `*_coverageReport*` files.

