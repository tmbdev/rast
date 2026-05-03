# C++ Coding Guidelines

## Memory Management

- **No owning raw pointers.** Ownership is always expressed through a smart pointer.
- **Smart pointer choice.** Use `std::unique_ptr<T>` for unique (linear) ownership and `std::shared_ptr<T>` for shared (reference-counted) ownership. Use `std::weak_ptr<T>` to break reference cycles.
- **Allocation.** In place of calling a constructor directly with `new T(args...)`, use one of the standard library's factory function templates, which invoke `T`'s constructor internally and return an owning smart pointer:
  - `std::make_unique<T>(args...)` (C++14) — returns `std::unique_ptr<T>`.
  - `std::make_shared<T>(args...)` (C++11) — returns `std::shared_ptr<T>`; also co-allocates the control block with the object.
  - `std::make_unique_for_overwrite<T>` / `std::make_shared_for_overwrite<T>` (C++20) — same, but with default (non-value) initialization.
  
  These templates *are* the constructor call from the caller's point of view: `args...` are perfect-forwarded to `T`'s constructor, and the result is wrapped in the appropriate smart pointer. Direct use of `new` is disallowed.
- **Returning heap-allocated objects.** Either:
  - *Preferred:* return the smart pointer by value, typically as `return std::make_unique<T>(...);` or `return std::make_shared<T>(...);`. NRVO and move semantics make this free, it composes naturally with expressions and `auto`, and it requires no out-parameter ceremony.
  - *Permitted:* pass a smart pointer by reference as an out-parameter and assign to it.
- **Non-owning observers.**
  - References (`T&`, `const T&`) are strongly preferred when the referent is non-null and non-rebinding.
  - Non-owning raw pointers (`T*`) are permitted but discouraged; reserve them for cases that genuinely require nullability or rebinding.
- **Prefer automatic storage and standard containers.** Use stack objects and standard containers (`std::vector`, `std::string`, `std::array`, etc.) in preference to heap allocation whenever lifetime and size permit. They already encapsulate ownership and obviate the smart-pointer question entirely.
- **Smart-pointer parameter passing.** Do not take a smart pointer as a parameter merely to observe the pointee. To observe, take `T&` (or `T*` if nullability is required). Take `std::unique_ptr<T>` by value to *transfer* ownership into the callee. Take `std::shared_ptr<T>` by value only when the callee will *share* ownership; otherwise observe via reference.
- **Ownership transfer.** `std::unique_ptr` is move-only: assignment between two named `unique_ptr` variables requires `p1 = std::move(p2);`. Assignment from an rvalue (e.g. `p = std::make_unique<T>();` or returning a local) needs no explicit `std::move`.
- **Custom deleters for non-memory resources.** `std::unique_ptr<T, Deleter>` (and `std::shared_ptr` with a deleter) generalize cleanly to file handles, OS handles, and other resources, providing RAII without bespoke wrapper classes.

## Iterators and Containers

- **Container choice.** `std::vector<T>` is strongly preferred as the default container. Reach for other standard containers (`std::array`, `std::map`, `std::unordered_map`, `std::deque`, etc.) only when their specific semantics are required.
- **Iteration by integer index.** Where the container supports random access, iteration by integer index is strongly preferred:
  ```cpp
  for (std::size_t i = 0; i < v.size(); ++i) { /* use v[i] */ }
  ```
  Indexing keeps the cursor as a plain integer, sidesteps iterator-invalidation hazards, and stays legible under debugging.
- **Typed sizes.** Use the container's own size type — `std::size_t`, or equivalently `Container::size_type` — for index variables, loop bounds, and any value derived from `.size()`. Do not use `int`, `unsigned`, or other ad-hoc integer types: doing so invites signed/unsigned comparison bugs and narrowing on 64-bit platforms. When a *signed* index is genuinely required (e.g. reverse iteration with `i >= 0`), use `std::ssize(c)` (C++20), which returns a signed type wide enough to hold any container size.
- **Ranges as the fallback.** When integer indexing is not applicable (associative containers, non-random-access containers, lazy or composed sequences), use the C++20 ranges library — range-based `for`, `std::ranges` algorithms, and `std::views` pipelines — in preference to raw iterator pairs.
- **Raw iterators.** Direct manipulation of `begin()` / `end()` iterator pairs is discouraged outside generic library code; prefer indices or ranges at call sites.
- **Bounds-checked access.** Use `.at(i)` rather than `operator[]` for element access. `.at()` performs a bounds check and throws `std::out_of_range` on violation; `operator[]` is undefined behavior on out-of-range indices. `operator[]` is permitted only inside loops whose bounds have already been checked against `.size()`, and only where profiling demonstrates the bounds-check cost matters.

## Views and Borrowed Data

- **Value semantics by default.** Pass and return owning types (`std::string`, `std::vector<T>`, etc.) by value or by `const&`. Copying is the default; views are not. Modern compilers elide and move aggressively, and the resulting code is far easier to reason about with respect to lifetime.
- **For strings, copy.** Pass and return `std::string` (or `const std::string&`). Do not hand out `std::string_view` as a casual substitute — the lifetime hazard is rarely worth the saved copy, and small-string optimization makes most copies cheap in practice.
- **`std::span<T>` and `std::string_view` are discouraged.** Their non-owning nature reintroduces the dangling-reference class of bugs that the rest of these guidelines are designed to eliminate. Permit them only when:
  - an external API or ABI mandates a `(pointer, length)` interface, or
  - profiling has demonstrated that copying is a serious, otherwise-unavoidable bottleneck.
  In either case, document the lifetime contract at the use site.
- **Never store a view as a data member**, never return a view from a function whose result might outlive the source, and never bind a view to a temporary.
- **Other borrowed/shared structures** (e.g. shared mutable state across components, ad-hoc reference-into-container patterns) are discouraged on the same grounds: prefer ownership transfer or copying to shared borrowing.

## Initialization

- **Always initialize.** Every variable — member, local, or static — is initialized at the point of declaration. Default member initializers in class definitions are preferred over assigning inside constructors.
- **Prefer brace initialization (`T x{…};`)** as the default form. It rejects narrowing conversions at compile time and applies uniformly to fundamental types, aggregates, and class types.
- **Use parentheses for container constructions whose semantics depend on the sized form.** Brace-init prefers a container's `std::initializer_list` constructor when one exists, which silently picks the wrong overload:
  ```cpp
  std::vector<int> a{3, 4};   // two elements: 3, 4
  std::vector<int> b(3, 4);   // three elements: 4, 4, 4
  ```
- **`auto` declarations** use the form `auto x = T{…};` or `auto x = expr;`. Do not write `auto x{expr};` — its meaning has shifted between standards and the result is type-ambiguous to readers.
- **Constructor initializer lists** initialize members in the order of *declaration in the class body*, not the order written in the initializer list. List members in declaration order to make this visible.

## Const-Correctness

- **Mark member functions `const`** whenever they do not modify the observable state of the object. This is required for them to be callable on `const` instances and through `const` references.
- **Prefer `const` parameters and locals.** Read-only parameters of non-trivial type are `const T&`; trivially-copyable types are passed by value. Local variables that are not reassigned should be declared `const`.
- **Style: west const.** Write `const T&` and `const T*`, not `T const&` / `T const*`. Both forms are valid C++; the project uses west const consistently.
- **`const` is shallow.** A `const` member function may still mutate state reached through pointers or references it owns. Reserve `mutable` for genuinely orthogonal state (caches, mutexes); do not use it to defeat const-correctness on logical state.
- **Do not cast away `const`.** `const_cast` is reserved for interoperating with legacy APIs that are observably const-correct but not declared so. Mutating an object originally declared `const` is undefined behavior.

## Enums

- **`enum class` only.** Bare unscoped `enum` is disallowed: it leaks names into the enclosing scope and implicitly converts to integer.
- **Specify the underlying type explicitly** when storage size or ABI matters: `enum class Color : std::uint8_t { Red, Green, Blue };`.
- **Switch over enums exhaustively.** Do not include a `default:` branch when switching on an enum whose values are all handled — let the compiler warn (`-Wswitch`) when a new value is added.

## `constexpr` and `consteval`

- **Mark functions `constexpr`** whenever they can plausibly be evaluated at compile time. This is also a useful reader signal that the function has no side effects.
- **Use `consteval` (C++20)** for functions that *must* be evaluated at compile time.
- **Prefer `constexpr` variables** to preprocessor `#define` for compile-time constants.
- **`if constexpr`** is preferred to SFINAE for compile-time branching inside templates.

## Lambdas and Captures

- **Prefer explicit captures** (`[x, &y]`) over default captures (`[=]`, `[&]`). Explicit captures document intent and make dangling-reference hazards visible at the capture list.
- **Never capture by reference in a lambda whose lifetime exceeds the captured object's** — for example, lambdas stored in a member, posted to a thread, or returned from the function.
- **Capturing `this`.** Plain `[this]` is acceptable only when the lambda's lifetime is bounded by the object's. For asynchronous lambdas on `shared_ptr`-managed objects, capture `[self = shared_from_this()]` instead.
- **Mark lambdas `noexcept`** when they are intended for `noexcept` contexts (e.g. destructors or move operations).

## Globals, Statics, and Initialization Order

- **No mutable namespace-scope state by default.** `constexpr` and `inline const` constants are permitted; in general, mutable globals and namespace-scope `static` variables are not.
- **Static initialization order across translation units is undefined.** Do not rely on it. For most global-like state, use the *Construct-On-First-Use* idiom — a function returning a reference to a function-local `static`:
  ```cpp
  Logger& logger() {
      static Logger instance;
      return instance;
  }
  ```
  Function-local statics are initialized on first call, and their initialization is thread-safe (C++11 magic statics).
- **Class-scope `static` data members** with non-trivial initializers carry the same hazard; use the same idiom (a static member function returning a reference).
- **Permitted exceptions: cross-cutting infrastructure.** Mutable namespace-scope variables are permitted for the following categories, where the alternative would be excessive ceremony at every call site:
  - **Configuration options** (parsed once at startup, then read).
  - **Debug and diagnostic settings** (verbosity levels, feature flags, tracing toggles).
  - **Loggers** and similar process-wide sinks.
  - **Type and plugin registries** populated at startup via static initializers.
  
  These are subject to the following conditions:
  - The variable's constructor must perform *constant initialization* or depend only on built-in types and trivially-constructible state — never on other namespace-scope variables in different translation units. This eliminates the static-initialization-order fiasco.
  - Concurrent access requires the variable to be `std::atomic<T>`, `const` after startup, or guarded by a documented mutex.
  - The declaration carries a comment stating the category and the thread-safety contract.
  - When in doubt, fall back to the function-local-static idiom above.

## Error Handling

For portability across ecosystems — including those that disable C++ exceptions (embedded, kernel, game engines, certain library boundaries, `extern "C"` interfaces) — the following conservative discipline applies:

- **Prefer `std::expected<T, E>` (C++23)** as the return type for any operation that can fail in a recoverable way. It expresses success-or-error in the type system, requires no exception machinery, and crosses translation-unit and ABI boundaries cleanly. Where C++23 is unavailable, an equivalent `Result<T, E>` / `tl::expected` polyfill is acceptable.
- **Reserve exceptions for genuinely exceptional, non-recoverable conditions** — out-of-memory, broken invariants, corrupted state — where unwinding to a top-level handler is the appropriate response. Do not use exceptions for ordinary control flow or for expected failure modes (file not found, parse error, validation failure); those return `std::expected`.
- **Programmer errors and contract violations** (precondition failures, unreachable code, internal invariant breaks) use `assert` in debug builds and a documented termination path (`std::abort`, `std::terminate`, or a project-defined fatal handler) in release builds — not exceptions.
- **C-style error codes with out-parameters are discouraged** in new code: they are easy to ignore, do not compose, and require ad-hoc conventions for the success value. Use `std::expected` instead.
- **Mark functions that cannot throw with `noexcept`**, particularly destructors, move operations, and swap. This is required for correct interaction with standard containers and enables optimizations.

## Concurrency

- **Avoid shared mutable state.** Prefer message passing, immutable data, or copying over locking. Most memory-safety bugs in modern C++ codebases are races, not use-after-free.
- **Locking.** Use `std::scoped_lock` (C++17) for any lock acquisition (single or multiple mutexes). Do not call `mutex.lock()` / `mutex.unlock()` directly. `std::lock_guard` is acceptable for single-mutex cases but offers no benefit over `std::scoped_lock`.
- **Mutex types.** `std::mutex` is the default. Use `std::shared_mutex` only when reader-writer semantics are genuine and measured to matter. `std::recursive_mutex` is disallowed in new code — recursive locking usually signals a design defect.
- **Atomics.** Use `std::atomic<T>` for individual fields requiring lock-free access. Restrict to trivially-copyable types and prefer the default `memory_order_seq_cst` unless a relaxed ordering is justified by profiling.
- **Threads.** Prefer `std::jthread` (C++20) over `std::thread`: it joins on destruction and integrates with `std::stop_token` for cooperative cancellation.
- **Document thread-safety on every interface.** Each interface header states which methods are thread-safe, which require external synchronization, and which are single-threaded.
- **Use Clang's thread-safety analysis** (`-Wthread-safety` with `GUARDED_BY(mutex_)` and related annotations) where the toolchain supports it.

## Classes

Three categories of class are recognized. Every class in the codebase belongs to exactly one.

- **Interfaces.** An interface is fully defined in a header. It contains only virtual functions (typically all pure) and, where needed, common instance variables shared by all implementations. Every interface declares a `virtual` destructor. Copy and move operations are deleted at the interface level *to prevent object slicing*: copying through a base reference would copy only the base subobject and silently truncate derived state. Ownership of interface objects is always through `std::unique_ptr<Interface>` or `std::shared_ptr<Interface>` per the memory-management rules.
- **Implementations.** An implementation class derives from one (or more) interfaces and is fully defined in a `.cc` file — declaration, member functions, and any private state. Its constructor is *not* called directly by clients. Instead, the `.cc` file exposes a `make_…` factory function declared in the corresponding interface header, returning a smart pointer typed as the interface (e.g. `std::unique_ptr<IFoo> make_foo(args…);`). The concrete type is complete only within that `.cc` file, where `std::make_unique<Concrete>(...)` can construct it; everywhere else, only the interface is visible. This is what confines all link-time dependencies on the concrete class to the single implementation translation unit.
- **Template classes.** Templates live entirely in headers, with no `.cc` counterpart. Member functions are defined inline in the header (or in a separately-included `.tpp` / `-inl.h` file included from the header). Explicit instantiation in a `.cc` file is permitted only when a template has a small, fixed set of instantiations and hiding the body materially improves compile time.

This three-bucket discipline yields a stable ABI at interface boundaries, confines implementation churn to single translation units, and makes link-time composition predictable.

### Special member functions: rule of zero, rule of five

- **Rule of zero (default).** Design value types so that none of the special member functions — destructor, copy constructor, copy assignment, move constructor, move assignment — need to be written. Compose the class out of types that already manage their resources correctly (standard containers, smart pointers, RAII wrappers) and let the compiler synthesize correct behavior.
- **Rule of five.** If you define, default, or delete *any* of the five special members, you must consider all five and declare each one explicitly (`= default`, `= delete`, or a hand-written body). Partial definitions cause subtle bugs: synthesized members may be implicitly deleted, may violate the class's invariants, or may surprise readers.
- **Interfaces delete all five** (per the Classes section above): copy/move operations are explicitly deleted to prevent slicing, and the destructor is `virtual`.
- **Move operations should be `noexcept`** wherever possible; standard containers select between copy and move based on this annotation.

### Permitted exception: `std::variant`

When a set of alternatives is *closed* (known and fixed at compile time) and the runtime cost of virtual dispatch plus heap allocation is a measured concern, a `std::variant<A, B, C, …>` with `std::visit` is permitted in place of an interface hierarchy. This buys value semantics, no heap allocation, and exhaustiveness checking, at the cost of forcing all alternatives to be visible in the header. Use sparingly and only where the closed-set property is genuinely stable; an open or growing set of alternatives belongs in the interface/implementation pattern.

## Naming Conventions

The codebase uses the LLVM / Qt / Microsoft naming family:

- **Types** — classes, structs, enums, type aliases, template parameters: `PascalCase`. Examples: `FileReader`, `HttpRequest`, `IFoo` (interfaces may carry an `I` prefix at the project's discretion).
- **Functions and methods**: `camelCase` preferred for new code. `snake_case` is also acceptable, particularly in code that integrates closely with the standard library or with existing modules that use it; consistency within a single class or library is required.
- **Variables** — locals, parameters, and globals: `camelCase`. Examples: `byteCount`, `currentLine`.
- **Non-static data members**: `camelCase`. An `m_` prefix or trailing underscore (`m_count` or `count_`) is *optional* and may be applied where it improves readability or disambiguates from parameters; consistency within a single class is required, but the project does not mandate one form across the codebase.
- **Static data members and global constants**: `camelCase`, optionally with a leading `k` for compile-time constants if the project elects that convention (`kMaxRetries`); again, consistency within a class or file is required, project-wide uniformity is not.
- **Macros**: `ALL_CAPS_WITH_UNDERSCORES`. Reserve macros for include guards, conditional compilation, and cases where no language feature suffices; do not use macros to substitute for `constexpr`, `inline`, or templates.
- **Namespaces**: short `lowercase` (`std`-style). Project-wide namespace is the project name in lowercase.
- **File names**: always all-lowercase, for portability across case-sensitive (Linux), case-preserving-but-insensitive (default macOS, Windows), and mixed filesystems. Use `snake_case` (`file_reader.h`, `file_reader.cc`) to keep multi-word names readable. The file name need not match the primary type's case (`FileReader` lives in `file_reader.h`).
- **Reserved identifiers.** Identifiers beginning with an underscore followed by a capital letter, or containing a double underscore, are reserved by the standard at any scope and must not be used. Identifiers beginning with a single underscore are reserved at namespace scope.

## Header Files

The codebase uses classical headers and `.cc` translation units. C++20 modules are not adopted at this time; revisit when toolchain and build-system support is non-experimental across the supported platforms.

- **File extensions.** `.h` for headers, `.cc` for translation units. Inline template helper bodies, when split out for readability, use `.tpp` or `-inl.h` and are included from the corresponding `.h` at end of file.
- **Header organization within a library.** A library directory's non-template declarations may be grouped into a single header for that library; this is the preferred default. If that header becomes too long or unwieldy, the next preference is one `.h` per `.cc`, accompanied by a library-wide umbrella header that `#include`s each of the per-file headers and is the canonical entry point for clients of the library.
- **Template and significant header-only classes get their own headers.** Every template class generally lives in its own dedicated header file, regardless of the library's grouping scheme above. The same applies to non-template classes that are nonetheless header-only — for example, classes whose implementation is fully `inline` or `constexpr` for performance reasons, or classes large enough that bundling them with the library's umbrella header would bloat compile times for unrelated translation units. Trivial inline helpers do not need their own headers; the rule applies to substantive classes whose definition would meaningfully grow the umbrella header.
- **Include guards.** Every header begins with traditional `#ifndef` / `#define` / `#endif` guards. `#pragma once` is *not* used: it is non-standard, and the traditional form is universally portable.
- **Guard naming.** The guard macro mirrors the header's path under the project root, in upper case, with non-alphanumeric characters replaced by underscores, and a trailing underscore. For `project/foo/bar/baz.h`:
  ```cpp
  #ifndef PROJECT_FOO_BAR_BAZ_H_
  #define PROJECT_FOO_BAR_BAZ_H_
  // ...
  #endif  // PROJECT_FOO_BAR_BAZ_H_
  ```
  Guard macros must not begin with an underscore or contain a double underscore — those identifiers are reserved by the standard.
- **Include order.** Within a `.cc` file, includes appear in the following order, each group separated by a blank line:
  1. The corresponding header (`foo.cc` includes `foo.h` first). This ensures `foo.h` is self-contained.
  2. C system headers (`<sys/types.h>`, etc.).
  3. C++ standard library headers (`<vector>`, `<string>`, etc.).
  4. Other-library headers.
  5. This project's headers.
  Within each group, sort alphabetically.
- **Self-contained headers.** Every header must compile on its own when included first. It must include or forward-declare everything it uses; it must not rely on the includer to have included its dependencies.
- **Forward declarations.** Prefer forward declarations to `#include` when only a pointer or reference to a type is needed. This reduces compile-time coupling and shortens rebuild fan-out.
- **No `using namespace` directives in headers.** They leak into every translation unit that includes the header. Namespace aliases (`namespace fs = std::filesystem;`) and targeted `using` *declarations* inside a class or function scope are permitted; `using namespace` at namespace scope is not.
- **No definitions with external linkage in headers** other than `inline` functions, `inline` (or `constexpr`) variables, templates, and class member functions defined in-class. Non-inline free functions and non-inline namespace-scope variables belong in `.cc` files.
- **`#pragma` and compiler-specific extensions** are confined to clearly marked, narrowly scoped blocks, ideally guarded by compiler-detection macros.

## Linking

- **Static libraries (`.a`) by default** for internal modules: no runtime dependency, no symbol-visibility hazards, no ABI surface beyond the project's own headers.
- **Shared libraries (`.so` / `.dylib`)** are reserved for OS-level boundaries, plugin interfaces, and components large enough that duplication across executables is a real cost.
- **Default symbol visibility is hidden.** Build with `-fvisibility=hidden -fvisibility-inlines-hidden`. Exported symbols are annotated explicitly through a project-wide macro (`PROJECT_EXPORT` expanding to `__attribute__((visibility("default")))` on Linux/macOS and `__declspec(dllexport)` / `dllimport` on Windows). This keeps shared-library ABI surfaces small and intentional.
- **One Definition Rule.** Every non-inline function and non-inline namespace-scope variable is defined exactly once in the program. `inline` functions, templates, and `inline` variables may have multiple definitions only if all are token-identical — meaning all participating translation units must include the same header.
- **Cross-DSO discipline.** Do not pass standard-library types whose layout depends on build configuration (`std::shared_ptr`, `std::function`, `std::string` with debug iterators, etc.) across shared-library boundaries unless both sides are built with identical compiler, standard library, and flags. For genuinely stable boundaries, expose `extern "C"` C-callable functions plus opaque handles, and provide C++ factory wrappers on the consuming side.
- **Static initialization across translation units.** Order is undefined; do not rely on it. Use the Construct-On-First-Use idiom from the previous section. Likewise, do not perform non-trivial work in destructors of namespace-scope objects.
- **The interface/implementation/factory pattern from the Classes section is the project's primary linking discipline:** the interface header is the only ABI surface, the concrete type is private to its `.cc`, and link-time fan-out is bounded.

## Build Hygiene

- **Standard.** C++20 (`-std=c++20`). C++23 features (`std::expected`, `std::print`, etc.) are used through polyfills until the compiler baseline reaches C++23 across the supported platforms.
- **Required warnings.** `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wdouble-promotion`.
- **Warnings are errors.** `-Werror`. Local suppressions for known-benign cases use `#pragma` blocks scoped as narrowly as possible.
- **Sanitizers in test builds.** Address sanitizer and undefined-behavior sanitizer are mandatory (`-fsanitize=address,undefined`). Thread sanitizer (`-fsanitize=thread`) runs on builds that exercise concurrency.
- **Standard-library debug modes.** Debug builds enable `_GLIBCXX_DEBUG` (libstdc++) or `_LIBCPP_HARDENING_MODE=debug` (libc++) to catch iterator-invalidation and bounds bugs.
- **No reliance on undefined behavior.** Code must not depend on signed-integer overflow wrapping, strict-aliasing violations, or other UB-based optimizer assumptions.
- **Reproducible builds.** Builds are deterministic given identical inputs; no embedded build timestamps or hostnames in release artifacts.
- **Build targets.** Three configurations are standard:
  - **`debug`** — no optimization (`-O0 -g`), assertions enabled, sanitizers (ASan, UBSan, and TSan where applicable) on, standard-library debug modes enabled. This is the default for development and CI tests.
  - **`production`** — moderate optimization (`-O2 -g`), assertions retained, no sanitizers, standard-library debug modes off. This is the configuration shipped to users and used for most performance work.
  - **`ultra`** — maximum optimization (`-O3 -DNDEBUG`, link-time optimization `-flto` where supported, machine-specific tuning `-march=…` where the deployment target allows it), assertions stripped. Reserved for performance-critical deployment artifacts; never used for development or for default CI runs because it disables runtime checks.
  
  All three configurations are buildable from the same source tree without source changes.
- **Build system.** See the Recommended Tools section. Whatever the choice, the same three targets (`debug`, `production`, `ultra`) are produced from the same source tree.

## Comments and API Documentation

- **Public interface headers** carry Doxygen-style comments on every declared function, type, and parameter. The minimal contract documents preconditions, postconditions, ownership transfer, error conditions, and threading expectations.
- **Implementation comments** explain *why*, not *what*. Code that requires extensive what-comments should be refactored. Inline comments on non-obvious lines (`// FNV-1a constant` next to a magic number) are encouraged.
- **`TODO` / `FIXME` / `HACK`** markers include an owner or ticket reference (`// TODO(jdoe, #1234): …`).

## Testing

- **Framework.** See the Recommended Tools section.
- **Test file layout.** Tests live alongside the code under test or in a `tests/` subdirectory of the library, at the project's discretion. Per-file tests are named `foo_test.cc` and accompany `foo.cc`; broader integration or scenario tests grouped under `tests/` use descriptive names. The build system links them into a per-library or per-executable test binary either way.
- **Mocking.** The interface/factory pattern (see Classes) is the project's mocking strategy: tests substitute a hand-written `IFoo` implementation and pass it where the real `make_foo()` would be called. No mock-generation framework is required.
- **Sanitizer coverage.** The CI test job runs the suite under ASan + UBSan; concurrency-touching libraries are also exercised under TSan. A test that passes only without sanitizers does not pass.
- **Property and fuzz testing** are encouraged for parsers, serializers, and protocol code; framework choices are listed in the Recommended Tools section.

## Formatting

- **Code is formatted by `clang-format`** against the `.clang-format` file at the project root. CI rejects unformatted code.
- The configuration encodes the existing house style: two-space indent, K&R / attached braces, right-aligned references and pointers, short function bodies permitted on a single line, constructor initializer lists kept on a single line where they fit. The full file is below.

```yaml
---
Language: Cpp
BasedOnStyle: LLVM
Standard: c++20
ColumnLimit: 100

# Indentation
IndentWidth: 2
TabWidth: 2
UseTab: Never
ContinuationIndentWidth: 4
AccessModifierOffset: -2
NamespaceIndentation: None

# Braces — K&R / attached
BreakBeforeBraces: Attach

# Pointers and references attach to the name (right alignment)
PointerAlignment: Right
ReferenceAlignment: Pointer
DerivePointerAlignment: false

# Short forms permitted on a single line (matches house style)
AllowShortFunctionsOnASingleLine: All
AllowShortBlocksOnASingleLine: Empty
AllowShortLambdasOnASingleLine: All
AllowShortIfStatementsOnASingleLine: Never
AllowShortLoopsOnASingleLine: false
AllowShortCaseLabelsOnASingleLine: false

# Constructor initializer lists kept on one line when they fit
PackConstructorInitializers: CurrentLine
BreakConstructorInitializers: BeforeColon
ConstructorInitializerIndentWidth: 4

# Spacing
SpaceBeforeParens: ControlStatements
SpacesInParentheses: false
SpacesInSquareBrackets: false
SpacesInAngles: false
SpaceAfterTemplateKeyword: true
SpaceBeforeAssignmentOperators: true
SpaceBeforeRangeBasedForLoopColon: true

# Includes
IncludeBlocks: Regroup
SortIncludes: CaseSensitive

# Templates
AlwaysBreakTemplateDeclarations: Yes

# Misc
AlignAfterOpenBracket: Align
BinPackArguments: true
BinPackParameters: true
FixNamespaceComments: true
KeepEmptyLinesAtTheStartOfBlocks: false
MaxEmptyLinesToKeep: 1
ReflowComments: true
```

This produces the existing reference style — for example:

```cpp
struct Segment {
  vec2 u, v;
  Segment() {}
  Segment(const vec2 &u, const vec2 &v) : u(u), v(v) {}
  vec2 sample(float lambda) { return u + (v - u) * lambda; }
  float ldist(vec2 p) {
    vec2 n = normal((v - u).normalized());
    float offset = n * u;
    return fabs(n * p - offset);
  }
};

bool within(float x, float y, float eps) { return fabs(x - y) <= eps; }
```

## Recommended Tools

Candidates for each tool category. The project will narrow each list to a single chosen tool; multiple are listed here so the trade-offs are visible.

### Compilers

- **GCC (GNU Compiler Collection).** The reference open-source C++ compiler. Mature, conservative defaults, broad platform and architecture support, excellent diagnostics, and strong sanitizer coverage. Default on most Linux distributions.
- **Clang (LLVM).** The other reference open-source C++ compiler. Often-faster compile times, often-clearer diagnostics, the most active development of advanced C++ features, and the home of the thread-safety analysis attributes recommended in the Concurrency section. Default on macOS and FreeBSD.

Both compilers are used in CI to surface portability bugs that one alone would miss.

### Build Systems

- **GNU Make.** Universally available, dependency-free, and produces a build with no hidden behavior. The verbosity is the cost of having no abstractions. Preferred for small projects.
- **build2.** Direct (non-meta) build system purpose-built for C++ by Boris Kolpackov of the ISO C++ committee. First-class C++20 modules support, integrated package manager (`bpkg`). Strongest technical fit for a serious C++ codebase that wants to stay non-meta. Costs: short, hard-to-search names; documentation has improved but is still uneven in places. (Not available via pixi — build2 is its own build and package system, installed out-of-band through its own script and `bpkg`.)
- **xmake.** Lua-based, single binary under 3 MB; combines a build engine, project generator, and package manager. Modern C++ support including modules. Costs: target names are globally unique (no namespacing), build phases are under-specified in places, and the design is more feature-driven than principled.
- **Tup.** Direct build system that uses a filesystem monitor for dependency tracking. Very fast incremental builds, especially for projects with many small inputs. Niche but well-engineered; best fit for asset pipelines or build graphs with high churn at the leaves. (Not available via pixi — install through OS package managers or from upstream releases.)
- **Pixi.** Rust-based package and environment manager (built on conda-forge) by Prefix Dev. Not itself a build system, but provides reproducible toolchains — compiler, build tool, formatters, analyzers — pinned in a `pixi.toml` so all developers and CI runners use identical versions. Useful as a meta-layer above any of the build tools above.

### Python Bindings

- **pybind11.** The de-facto standard. Header-only C++ library, highly polished, large user base, excellent documentation. The default choice for new bindings.
- **SWIG.** Multi-language code generator. Recommended only for legacy bindings already built on it, or when bindings are needed for Python *and* additional languages (Ruby, Java, etc.) from a single source.

### Command-Line Parsing

- **CLI11.** Declarative, header-only, available on conda-forge and therefore through pixi (`pixi add cli11`). The default choice for any tool with non-trivial options, subcommands, configuration files, or environment-variable bindings. Options bind directly to typed variables; `--help` and `--version` are auto-generated; validators and constraints (`->required()`, `->needs()`, `->excludes()`, `->envname()`) compose declaratively.
- **`std::getenv` with lowercase variable names for simple programs.** When a tool only needs a handful of values, reading them from the environment is acceptable in place of a parser. The program reads each value via `std::getenv`:
  ```cpp
  std::string something = std::getenv("lowercase_name");
  ```
  and the user supplies the values inline at invocation in the standard shell form:
  ```
  lowercase_name=value other_name=value program
  ```
  Lowercase variable names distinguish program-specific options from conventional `ALL_CAPS` system and shell variables. `std::getenv` lives in `<cstdlib>`, returns `const char*`, and returns `nullptr` when unset — wrap with a null check and a default. Skip CLI11 entirely when the tool truly is that simple; reach for it as soon as flags, subcommands, help text, or validation become useful.

### Code Formatting

- **clang-format.** The de-facto standard. Configurable via `.clang-format`, integrates with virtually every editor and CI system, and produces deterministic output. The tool this document already encodes.

### Style Checking and Static Analysis

- **clang-tidy.** Companion to clang-format. Extensive built-in checks across modernization, bug-prone patterns, performance, readability, and full C++ Core Guidelines coverage. The default modern choice; runs as part of CI.
- **cppcheck.** Independent static analyzer with no LLVM dependency. Catches a different class of bugs (bounds, dead code, suspicious patterns); often used alongside clang-tidy rather than as a substitute.

### Documentation

- **mkdocs / mkdocs-material.** Markdown-based documentation generator; simple, clean output, fast to set up. The default choice. API reference can be generated from Doxygen-style comments via the `mkdoxy` plugin (or by running Doxygen separately and linking from the mkdocs site) when needed.
- **Doxygen.** The long-standing standard for C++ API documentation, generating HTML and LaTeX/PDF from specially-formatted comments in headers. An acceptable alternative — and the right primary choice for projects whose deliverable is large, heavily-cross-referenced API reference rather than narrative documentation.

### Testing

- **doctest.** Single-header, dependency-free, fast compile times, expressive assertion macros. The default unit-test framework.
- **rapidcheck.** Property-based testing in the QuickCheck tradition; pairs with doctest rather than replacing it. Used for parsers, serializers, and other code with broad input domains.
