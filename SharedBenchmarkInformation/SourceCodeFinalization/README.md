# Source Finalization Overview

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** 2026-05-26

## Purpose

This document contains a brief overview of the source code finalization documents.

It is meant to help understand what each document is for quickly.

## What Source Finalization Means

Source finalization is the cleanup step applied last before programmability metrics are collected.

This includes standardizing formatting, organizing helper functions, adding file headers, clarifying file responsibilities, removing unused code, and selecting the right final files for metrics.

Source finalization doesn't change benchmark behavior, algorithm logic, timing, validation, or runtime performance methodology.

## Documents in This Folder

### `FormattingRulesUsed.md`

Defines the source formatting and cleanup style used for all implementations.

Use this document for rules about braces, line breaks, parameter formatting, comments, naming, constants, and basic cleanup.

### `HelperFunctionGuidelines.md`

Defines when helper functions should be created, named, placed, and marked as `static` or `inline`.

Use this document when deciding whether logic should stay in place, become a helper, or remain library specific.

### `FileHeaderGuidelines.md`

Defines the standard header format for Markdown documents and C/C++ header files.

Use this document for adding metadata (such as author and benchmark naming).

### `FileResponsibilityGuidelines.md`

Defines what's in `main.cpp`, benchmark implementation files, and header files.

Use this document when deciding whether code belongs in the benchmark logic, main, or header.

### `ProgrammabilityMetricPreparation.md`

Defines how final source files were prepared before programmability metrics were collected.

Use this document for what files were included or excluded from metric collection and what happens before metrics are collected.

## Recommended Reading Order

1. `FormattingRulesUsed.md`
2. `FileResponsibilityGuidelines.md`
3. `HelperFunctionGuidelines.md`
4. `FileHeaderGuidelines.md`
5. `ProgrammabilityMetricPreparation.md`

## Notes

These documents are for finalization guidelines, not algorithm specifications.

When a change goes beyond just formatting, helpers, file responsibility, or metric preparation, it's out of the scope of these documents.