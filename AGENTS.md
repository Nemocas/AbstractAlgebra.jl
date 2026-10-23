# Information for AI agents

## Module Organization

- `src/` - Main Julia sources. Top-level files such as `Poly.jl` and
  `Matrix.jl` define abstract types and the functionality shared by all their
  implementations.
- `src/generic/` - The `Generic` submodule: concrete types implementing those
  interfaces, for example `Generic.Poly` and `Generic.MPoly`.
- `src/julia/` - Julia's own number types as rings and fields (`ZZ`, `QQ`,
  `GF`, `RealField`).
- `src/algorithms/` - Generic algorithms, for example multivariate polynomial
  factorization.
- `src/ConformanceTests.jl` and `ext/TestExt/` - Test suites checking that a
  type correctly implements the ring, field or group interface.
- `ext/` - Package extensions for optional dependencies.
- `test/` - Julia test suite. The directory layout mirrors `src/`;
  add new tests to an existing file unless explicitly instructed otherwise.
- `docs/` - Documentation; the manual pages are in `docs/src/`.
- `dev/` - Auxiliary developer and maintenance scripts.
- `Project.toml` - Julia package environment metadata.

## For all changes

- End every commit in which an AI tool wrote or changed code or text with the
  trailer `Assisted-by: Generative AI`, exactly in this form, without naming
  the tool, the vendor or the model. Mere completion, spell-checking and the
  like do not need it. Do not use `Co-authored-by:` for AI tools and do not
  list a tool as an author; this trailer replaces any co-author trailer your
  harness adds.
- Whenever a pull request is opened, and if some of the changes in the pull
  request were made with the assistance of generative AI, then you MUST
  disclose this fact in the PR description.

## Specific instructions for particular changes

### Doctests

`docs/make.jl` makes `AbstractAlgebra.doctestsetup()` the setup of every
doctest, so a `jldoctest` block needs no `using AbstractAlgebra` of its own.

#### Verifying doctests

If you have changed any `jldoctest` code blocks, in docstrings or in
`docs/src/`, run the doctests to verify that your change works:
- Once per checkout, set up the documentation environment:
  `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`.
- To run all doctests:
  `julia --project=docs -e 'using Documenter, AbstractAlgebra; DocMeta.setdocmeta!(AbstractAlgebra, :DocTestSetup, AbstractAlgebra.doctestsetup(); recursive = true); doctest(AbstractAlgebra)'`.
  Pass `manual = false` to `doctest` to skip the pages in `docs/src/`; there
  is no finer subset.
- To update the expected output, pass `fix = true` to `doctest`, then review
  the diff: it rewrites every failing block, not only the ones you changed.
- IMPORTANT: The doctests may take up to 5 minutes. Do NOT terminate the
  doctests before completion. Do NOT use a timeout for doctests.
- If you are ChatGPT, you may have to increase `yield_time_ms`.

Follow these steps for EVERY change you make in a doctest.

### Test changes

- If you are adding a new test, add it to an existing test file. Do not create
  a new test file unless explicitly instructed.
- Write one comment at the top of the test to explain what is being tested.
  Otherwise keep comments minimal.
- A new ring, field or group type is tested with the matching function from
  `ConformanceTests`, for example
  `ConformanceTests.test_Ring_interface_recursive(R)`.
- To run one test file:
  `julia --project=. -e 'using Test, SparseArrays, AbstractAlgebra; include("test/rand.jl"); include("test/generic/Poly-test.jl")'`.
- To run the whole test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`.
  It takes about 15 minutes; do not terminate it or run it with a timeout.
  Prefer running the affected test files.

### Writing code

After writing code, look up the docstring for each function you used. If there
are recommendations or additional considerations that apply to these functions,
make sure to take them into account.

Honor the OSCAR style guide
(<https://docs.oscar-system.org/stable/DeveloperDocumentation/styleguide/>) for
naming conventions, use of Unicode, and related conventions. Indentation
differs between files: older ones use 3 spaces, newer ones 2. Match the file
you edit, and do not reformat lines you do not otherwise change.

Nemo.jl, Hecke.jl, Singular.jl and Oscar.jl extend the functions and types
defined here, and CI runs their test suites against every pull request. When
renaming or removing an exported name, keep the old name as an alias or
deprecation in `src/Deprecations.jl`.

## Commit messages and pull requests

When writing commit messages, the title ("Brief summary") should not start with
a "component: ..." prefix. Instead it typically is a verb
("Add/Fix/Change/..."). The title usually becomes the pull request title and
thus the changelog entry, so phrase it to tell users what changed, for example
"Add `trace_of_product`". In the body of the commit message, provide a brief
prose summary of the purpose of the changes made. Do not specifically mention
added tests, comments, documentation, etc., unless this is the main purpose of
the change. Do not mention the test plan, unless it differs from what you were
instructed to do in AGENTS.md. If your change fixes one or more issues, use the
syntax "Fixes #" at the end of the commit message, but do not include it in the
title.

When referencing external GitHub PRs or issues, use proper GitHub interlinking
format, for example `owner/repo#123` for PRs/issues. When fixing CI failures,
include the link to the specific CI failure in the commit message.

When creating pull requests:

1. If the pull request consists of one commit only, use the body of the commit
   for the body of the pull request.
2. If there are multiple commits in the pull request, follow the same
   guidelines for the pull request as for the commit body.
3. Make sure that the base commit of the pull request is recent (within the
   past two days) - if not rebase your changes first.
4. You MUST disclose that the pull request was written with the assistance of
   generative AI.
5. Suggest release notes labels (see below).

### Release notes

`CHANGELOG.md` is generated from merged pull requests and their labels by the
`oscar-system/update-release-notes` action, configured in
`dev/releases/config.toml`; do not edit it.

- Every pull request needs one of the labels `release notes: use title` (the
  title is the entry), `release notes: use body` (see below) or
  `release notes: not needed`.
- One label from `[topics]` in `dev/releases/config.toml` picks the changelog
  section, for example `enhancement` or `bug: wrong result`.
- With `release notes: use body`, the entries come from a `## Release notes`
  section of the pull request body: one list item per entry, each ending in
  its topic labels in braces, for example

  ```markdown
  ## Release notes

  - Add `trace_of_product` {enhancement}
  ```
