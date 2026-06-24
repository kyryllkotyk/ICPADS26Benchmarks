# Helper Function Guidelines

**Prepared by:** Kyryll Kotyk  

**Faculty Advisor:** Munehiro Fukuda  

**Last Updated:** 2026-05-25

## Purpose and Scope

These guidelines define how helper functions should be created, named, placed, and applied in the final implementations before programmability metrics are collected.

The goal is to use helpers consistently enough to make the code easier to understand and maintain, while avoiding changes that affect benchmark behavior, algorithm structure, timing boundaries, or performance.

Helper functions should support clarity and reduce repeated source code where reasonable. They shouldn't be used to hide important benchmark logic or restructure one implementation more aggressively than another.

## When Helper Functions Should Be Created

Helper functions should be created when they represent a clear operation that is repeated or separate from the surrounding code.

Helpers are encouraged when they:

- Remove repeated logic
- Reduce total source code length without hurting performance
- Give a meaningful name to a clear operation
- Separate setup, validation, indexing, printing, or utility logic from core algorithm logic
- Make a large function easier to follow without hiding the benchmark algorithm
- Preserve the original behavior, output, and timing meaning

Examples of reasonable helper uses include:

- Input validation
- Argument parsing
- Decomposition
- Index conversion
- Local buffer initialization
- Checksum calculation
- Repeated output formatting
- Repeated communication setup

## When Helper Functions shouldn't Be Created

Helpers should generally not be created for logic that only happens once unless there is a clear readability advantage.

A helper shouldn't be created when it:

- Wraps only one trivial line without adding meaning
- Makes the code harder to follow
- Hides benchmark control flow
- Hides the main algorithm inside a vague helper 
- Adds measurable overhead and is used frequently enough to be noticable in the final results
- Introduces extra allocation, copying, synchronization, or communication
- Changes where benchmark timing starts or stops
- Changes algorithm behavior, output, validation, or checksum
- Restructures one implementation much more than the others without justification

Single use helpers are allowed only when the operation is important enough to name clearly or when keeping it inline would make the code near it more difficult to understand/follow.

## Reasonable Use Across Libraries

Helper usage should be applied with fairness in mind, but not forced identically across MPI, MASS, and HPX.

These libraries may benefit differently from specific helpers due to the difference in programming models and control flow. A helper in MPI might not make sense for MASS or HPX, and the reverse may be true.

The main standard is reasonability:

- Similar repeated logic should receive similar cleanup when practical
- Library specific structure should be respected
- Helpers shouldn't make one implementation appear artificially better in terms of programmability
- Helpers shouldn't force another implementation into an unnatural structure
- Any major helper based refactoring should be justified separately

## Algorithm Structure Preservation

Helpers may organize code, but they shouldn't change the algorithm being benchmarked.

The cleaned implementation should still make it clear where the major benchmark phases occur, like initialization, communication, computation, synchronization, timing, and validation.

Helpers shouldn't hide enough logic that the benchmark structure becomes difficult to understand.

## Timing Boundary Preservation

Helper functions mustn't change benchmark timing boundaries.

Code that was originally inside a timed region should stay inside the timed region. Code that was originally outside a timed region shouldn't be moved inside the timed region by cleanup.

When a helper is called inside a timed region, it should be clear that the helper contributes to that timed measurement.

When a helper is called outside a timed region, it shouldn't perform hidden benchmark work that should have been measured.

## Helper Naming

Helper names should be clear and descriptive.

Names should explain what the helper does without requiring the reader to check the code inside of it.

Use names such as:

```cpp
split1D
validateInputs
computeChecksum
fillPanelForStep
startBroadcastForStep
```

Avoid vague names such as:

```cpp
validate
fill
process
handleData
```

Short names are acceptable only when the meaning is obvious.

## Header Placement Rule

All helper functions should be placed in the matching header file.

This applies to:

- Helper declarations
- Helper definitions
- Inline helper functions
- Static member helper functions
- Private class helper functions
- Small utility helpers 

For this project, helper placement should follow the same rules across all implementations. These rules are meant for benchmark consistency, and may differ
from the separation practices normally preferred in larger C++ projects.

Implementation files should focus on benchmark implementation logic. Headers should contain shared setup, declarations, helper definitions, constants, and reusable helper logic.

## Public and Private Helpers

Benchmark entry point functions should be public.

Helper functions should generally be private because they exist to support the benchmark implementation.

Use `public` for functions that external runner code needs to call.

For others, use `private`.

A helper should only be public if there is a clear reason for code outside the class to call it.

## Static Helper Rules

Use `static` for class helper functions that do not depend on object state.

The `static` keyword should be written in the class declaration in the header.

Example:

```cpp
class PM2_DGEMM
{
private:
    static pair<int, int> split1D(
        const int totalDimensionSize,
        const int processesInDimension,
        const int processIndex
    );
};
```

When defining a static member helper outside the class body, do not repeat the `static` keyword.

Example:

```cpp
inline pair<int, int> PM2_DGEMM::split1D(
    const int totalDimensionSize,
    const int processesInDimension,
    const int processIndex
) {
    ...
}
```

Static helpers should generally follow the header placement rule unless it creates a build, linkage, or design problem.

Avoid `static` free functions in headers unless there is a specific reason to give each translation unit its own private copy.

## Inline Helper Rules

Use `inline` for helper definitions placed in the header when the helper is short.

For these benchmarks, using inline helpers in the header is acceptable and helps readability via explicit definition of all used functions.

Correct pattern:

```cpp
class PM2_DGEMM
{
private:
    static pair<int, int> split1D(
        const int totalDimensionSize,
        const int processesInDimension,
        const int processIndex
    );
};

inline pair<int, int> PM2_DGEMM::split1D(
    const int totalDimensionSize,
    const int processesInDimension,
    const int processIndex
) {
    ...
}
```

## Lambda Rules

Short local lambdas are allowed when they make nearby code clearer and are only used in a small local scope.

Long lambdas should be avoided. If a lambda becomes long, important, repeated, or difficult to understand quickly, it should become a named helper function in the header.

A lambda shouldn't be used as a way to hide helper logic inside a large function.

## Performance Requirements

Helpers mustn't hurt benchmark performance.

Helpers should avoid:

- Extra memory allocation
- Extra data copying
- Extra synchronization
- Extra MPI communication
- Extra virtual dispatch
- Extra unnecessary abstraction in tight loops

Performance sensitive helpers should use appropriate parameter passing, such as references, pointers, or `const` references where useful.

Small helpers in hot paths may be marked `inline` when placed in the header.

Helpers shouldn't be introduced if they make the benchmark slower in a meaningful way.

## Source Code Length Consideration

Helpers should be used where reasonable to reduce repeated code and lower total source code length.

However, helpers not all helpers affect the results positively. A helper that increases boilerplate, adds extra comments, or splits a simple operation into unnecessary pieces may increase NLOC instead of reducing it. So, helpers only should be added when they actually improve the implementation, not just to remove a few lines total.

The goal is to reduce duplicated logic and improve structure, not to split code blindly.

## Commenting Helpers

Helper functions should have short comments when their purpose is not immediately obvious.

Comments should explain what the helper is responsible for, not restate every line of code.

Example:

```cpp
// Split one global dimension into balanced chunks
```

Comments should be consistent with the project formatting rules.

## Behavior Preservation

Creating, removing, or changing helper functions mustn't change benchmark behavior.

Helper cleanup must preserve:

- Algorithm semantics
- Input meaning
- Output meaning
- Timing region meaning
- Checksum behavior
- Validation behavior
- Communication behavior
- Performance expectations