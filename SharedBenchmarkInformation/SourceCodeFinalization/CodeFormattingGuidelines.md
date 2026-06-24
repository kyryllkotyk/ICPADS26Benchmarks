# Code Formatting Guidelines

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** 2026-05-25  


## Purpose and Scope
These rules define the formatting and cleanup standard used for the final MPI, MASS, and HPX benchmark implementations by Kyryll Kotyk and Ahmed Bera Pay before collecting programmability metrics.

The goal is to make the implementations stylistically consistent and fair to compare without changing their behavior, logic, or semantics. Because the benchmarks were developed by two authors, the original code had differences in indentation, naming, line length, comment style, and general formatting. This document defines the shared style-normalization rules used to reduce those differences before (programmability) metrics collection.

This document is only concerned with formatting and basic source-code cleanup. Its purpose is to make the implementations follow the same style standard before metrics are collected, so NLOC and basic readability comparisons are based on implementation differences rather than inconsistent author formatting.

It is also not a helper function policy. It does not decide when helpers should be created, removed, merged, split, or preserved. 

Any changes beyond formatting and basic style cleanup are not of concern for this document. If you are looking for that information, refer to other documents (as available).

## Rules

### 1. Braces and Newlines

All control flow blocks must use braces, even if the body only has one statement.

Do not write:

```cpp
if (condition) doWork();
```

Write:

```cpp
if (condition) {
    doWork();
}
```

This applies to `if`, `else`, `for`, `while`, and similar.

### 2. Function Parameters

Any function call or function signature with at least one parameter must write each parameter on its own line.

Do not write:

```cpp
runBenchmark(runs, rows, cols, seed);
```

Write:

```cpp
runBenchmark(
    runs,
    rows,
    cols,
    seed
);
```

This applies to both function definitions and function calls.

This rule also applies to calls such as MPI functions:

```cpp
MPI_Comm_rank(
    MPI_COMM_WORLD,
    &rank
);
```

### 3. Stream Formatting

`cout` and `cin` chains should be split across multiple lines instead of being written as one long line.

Do not write:

```cpp
cout << "Ranks: " << ranks << " Runs: " << runs << "\n";
```

Write:

```cpp
cout
    << "Ranks: "
    << ranks
    << " Runs: "
    << runs
    << "\n";
```

### 4. Ternary Expressions

Ternary expressions should be split to multiple lines, similarly to stream formatting, when they contain important logic or would make the line too long.

Do not write:

```cpp
int count = base + (index < extra ? 1 : 0);
```

Write:

```cpp
const int count = base + (
    index < extra
        ? 1
        : 0
);
```

### 5. Line Length

Lines should be at or below 80 characters when reasonably possible.

When splitting a long line, try to split at cleaner points, around 50 to 60 characters, so the split isn't awkward to read.

### 6. Constants and Macros

Use `#define` for named constants where possible.

Example:

```cpp
#define MICRO_KERNEL_ROWS 4
```

### 7. Naming

Variable and function names should be clear and descriptive.

Avoid short or vague names when the meaning is not very obvious.

Do not write:

```cpp
int p;
int proc;
```

Write:

```cpp
int processes;
int processRank;
int processGridRows;
```

Short loop variables such as `i`, `j`, and `k` should still be used if the loop meaning is obvious / adding a descriptive name wouldn't improve clarity much

### 8. Const Correctness

Use `const` for variables, parameters, pointers, and function inputs that are not modified.

Example:

```cpp
const int processGridSize = processGridRows;
```

Function parameters should also be marked `const` when the function does not modify them:

```cpp
void runBenchmark(
    const short runs,
    const int rows,
    const int cols
);
```

### 9. Integer Data Types

Use the smallest feasible integer type for large integer arrays or large storages when the range is known and doing so won't hurt performance or clarity.

For example, if a large array stores values that fall within `short` range, it shouldn't use `int`.

This rule does not require changing every small variable. Non-critical variables such as the following may stay as `int` or whatever they were defined as:

```cpp
int mpiSize = 24;
```

This prevents wasting memory on large integer storage, and each individual variable won't meaningfully change anything, so it's fine to leave as is.

### 10. Comments

Comments should be placed above the line or block they describe, not beside or below it.

Do not write:

```cpp
MPI_Barrier(MPI_COMM_WORLD); // Synchronize ranks
```

Write:

```cpp
// Synchronize ranks before timed section
MPI_Barrier(
    MPI_COMM_WORLD
);
```

Comments should usually be short, around one or two lines. They should explain the purpose of the code, not restate obvious syntax.

### 11. Dead or Disabled Code

Dead code, unused code, and old commented-out implementations should be removed during cleanup unless there is a specific reason to preserve them.

Do not leave old benchmark branches, unused algorithm calls, or large commented code blocks in the final cleaned implementation.

### 12. Include Cleanup

Unused includes should be removed when they are clearly unnecessary.

For example, if a file does not use anything from `<cstring>`, then it should not include it.

Required includes should be kept explicit so each file is understandable on its own.

### 13. Header and Include Organization

Project-level includes, shared `#define` constants, type definitions, and reusable declarations should be placed in the matching header file whenever possible.

Implementation files should avoid adding extra includes, macros, or shared declarations unless they are only needed privately inside that specific `.cpp` file.

For example, benchmark-wide constants such as this should go in the header:

```cpp
#define MICRO_KERNEL_ROWS 4
```

### 13. Behavior Preservation

Formatting cleanup must not change benchmark behavior, benchmark logic, algorithm semantics, input/output meaning, timing regions, checksum or validation behavior.

The cleaned implementation should produce the same benchmark behavior as the original implementation.