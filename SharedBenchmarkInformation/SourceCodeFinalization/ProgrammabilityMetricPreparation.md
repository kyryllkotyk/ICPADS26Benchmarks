# Programmability Metric Preparation

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** 2026-05-26

## Purpose

This document describes how finalized source files were prepared before collecting programmability metrics.

The goal is to explain what happens to the source code when finalizing, and what is included in final metrics collection.

## Scope

This document covers source code preparation before metrics collection. It does not define the metric tools, metric commands, runtime performance collection process, or interpretation of metric results.

Runtime performance results and programmability metrics are collected separately.

## Included Files

Only final implementation source files were included for programmability metrics.

Generally, these files were included for metrics:

- Final `.cpp` implementation files
- Final `.h` header files
- Benchmark specific source files needed to accurately show programmer work for fairness 

Header files were included because helper functions, declarations, constants, and library specific definitions may be there.

## Excluded Files

The following files were excluded from programmability metrics:

- `main.cpp` files
- Runner scripts
- Temporary files
- Raw logs
- Result CSV files
- Executables
- Generated files
- All non-final implementation files

`main.cpp` files were excluded because the metrics are meant to focus on the final benchmark implementation files.

## Finalization Before Metrics

Formatting normalization, file header cleanup, helper function organization, dead code removal, and final file selection were completed before metrics were collected.

This preparation was done so the metrics reflect the final implementation rather than inconsistent formatting, old code, or meaningless code left over from development.

## Fairness Notes

PM2/MPI, MASS, and HPX implementations may require different file structures because the libraries express benchmark logic differently. Equal preparation does not require identical file counts or identical organization.

The goal is to apply the same cleanup standard where reasonable, not to force all implementations into the same structure.

## Comment Treatment

NLOC (Non-Comment Lines of Code) metric was used instead of LOC (Lines of Code), so comments were not counted as source lines in the metrics collected by Lizard.

Comments were still normalized before metrics collection to improve consistency and readability.
