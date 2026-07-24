# multiOmicsSuite

## Description

multiOmicsSuite (`MOSuite`) is an R package for differential multi-omic analysis.
It defines an S7 class called `multiOmicDataSet` to store data and analyses from
multi-omic experiments.

Development of MOSuite follows the R packages 2nd edition (https://r-pkgs.org/),
with a few minor exceptions noted below.
Helper functions from `usethis` and `devtools` are used extensively for development tasks.

## Package conventions

- **Internal functions should have roxygen2 documentation** (with `@keywords internal`). Do not strip roxygen docs from internal functions.
- R code should pass `lintr` and `air format` (run `air format .` from the package root).
- Tests should be written with `testthat`.
- The package must pass `devtools::check()`.
- R code should adhere to the tidyverse style guide. https://style.tidyverse.org/
- Only include one return statement at the end of a function, if a return statement is used at all. Explicit returns are preferred but not required for R functions.

## Commit messages

- Commit messages must follow [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) (as enforced in `CONTRIBUTING.md`).
- Generate messages from staged changes only (`git diff --staged`); do not include unrelated work.
- Commits should be atomic: one logical change per commit.
- If mixed changes are present, split into multiple logical commits; the number of commits does not need to equal the number of files changed.
- Subject format must be: `<type>(optional-scope): short imperative summary` (<=72 chars), e.g., `fix(profile): update release table parser`.
- Add a body only when needed to explain **why** and notable impact; never include secrets, tokens, PHI, or large diffs.
- For AI-assisted commits, add this final italicized footer line in the commit message body: _commit message is ai-generated_

## Pull request (PR) process

- When opening a PR, use the request template (`.github/PULL_REQUEST_TEMPLATE.md`) and fill out all sections of the template in the PR description.
- Do not allow the developer to proceed with opening a PR if it does not fill out all sections of the template.
- Before a PR can be moved from draft to "ready for review", all of the relevant checklist items must be checked, and any
irrelevant checklist items should be crossed out.
- If code is AI-generated, the PR should be labeled `generated-by-AI`. There should be a brief, concise statement in the PR description of how AI was used in creating the PR (model used, high-level prompt intent, manual review confirmation).
- When new features, bug fixes, or other behavioral changes are introduced to the code,
unit tests must be added or updated to cover the new or changed functionality.
- If there are any API or other user-facing changes, the documentation must be updated both inline via roxygen comments and long-form docs in the `vignettes/` directory as R Markdown files.
- The `R-CMD-check` github actions workflow must pass before the PR can be approved.

### Changelog

The changelog for the repository is maintained in `NEWS.md`  at the root of the repository.
Each pull request that introduces user-facing changes must include a concise
entry with the PR number and author username tagged.
Developer-only changes (i.e. updates to CI workflows, development notes, etc.)
should never be included in the changelog.

Example:

```
## development version

- Fix bug in `detect_absolute_paths()` to ignore comments. (#123, @username)
```
