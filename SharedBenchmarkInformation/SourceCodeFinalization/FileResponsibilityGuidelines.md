# File Responsibility Guidelines

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** 2026-05-26

## Purpose

This document defines the expected responsibility of each source file used in the final implementations.

The goal is to keep entry points, logic, headers, and library specific code separated clearly enough for fair programmability metric collection.

## Main File Responsibility

`main.cpp` files should be responsible only for launching the benchmark.

A main file may handle:

- Argument parsing
- Basic usage messages
- MPI or library initialization when required
- Benchmark object creation
- Calling the benchmark runner function
- Finalization or cleanup required by the runtime

A main file should not contain the core benchmark algorithm, benchmark phase logic, major helper logic, or library specific implementation details unless the library requires it.

## Benchmark Implementation File Responsibility

Benchmark implementation files should contain the main benchmark logic.

This includes:

- Benchmark setup
- Main benchmark phases
- Computation 
- Communication
- Timing 
- Validation and checksum
- Library specific code

The implementation files may look different for each library because they express logic differently. The goal is to have clear separation 
and to not include important logic in a file that won't be used when collecting metrics (`main.cpp`).

## Header File Responsibility

Header files should contain shared declarations and setup needed by the benchmark implementation.

Headers should not contain unrelated benchmark logic or old code.

For more information about header file helper placement, see `HelperFunctionGuidelines.md`, specifically the sections from **Header Placement Rule** through **Inline Helper Rules**.

## Library Specific File Differences

Some libraries may need additional files or different organization.

For example, MASS Places or Agents naturally belong in separate files. HPX components could also require structure that does not match PM2/MPI.

These differences are allowed when they come from the programming model or library design. Implementations aren't forced into identical file layouts when doing so would make the code less natural or less correct for that library.