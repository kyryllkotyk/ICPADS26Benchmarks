# File Header Guidelines

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** 2026-05-26

## Purpose

This document defines the file header format used for Markdown files and C/C++ header files.

The goal is to make project files easy to identify without adding unnecessary explanation or changing implementation logic.

## Markdown File Headers

Markdown documents should start with the document title, followed by the project metadata.

Use this format:

```markdown
# Document Title

**Prepared by:** Kyryll Kotyk

**Faculty Advisor:** Munehiro Fukuda

**Last Updated:** YYYY-MM-DD
```

There should be a blank line between each line so Markdown renders them as separate lines, instead of combining together.

## Header File Headers

C/C++ header files should start with a boxed file comment before the include guard.

Use this format:

```cpp
/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Benchmark Name                                                  *
 * Library: Library Name                                                      *
 *                                                                            *
 * Author: Author Name                                                        *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Author of Cleanup Name                                  *
 *                                                                            *
 *****************************************************************************/

#ifndef HEADER_GUARD_NAME_
#define HEADER_GUARD_NAME_

...
#endif
```

The file header should identify the benchmark, library, author, faculty advisor, and final formatting and cleanup author.

The include guard should come immediately after the file header.