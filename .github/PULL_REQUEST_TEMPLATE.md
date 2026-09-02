## Changes

<!--
Provide a brief, concise summary of what is included in this Pull Request (PR).
If the PR contains complex changes, also include a brief description of the rationale
behind the implementation decisions.
Focus this section on the actual user-facing behavior change, bug fix, or design decision.
Do not include descriptions of obvious, redundant, automated, or required changes that are already covered by checklist items.
Do not add housekeeping-only bullets. In particular, do not list the following in the Changes section unless one of them is itself the primary behavior change being proposed:
- regenerated `man/*.Rd` or `NAMESPACE` files
- `NEWS.md` bookkeeping updates made to satisfy release/changelog requirements
- routine roxygen2 or test file regeneration that only reflects the underlying feature/fix
- any bullet that only says files were regenerated/updated (for example, "Regenerate `man/*.Rd`" or "Update `NEWS.md`")
-->

## Issues

<!--
Reference any issues related to this PR.
If this PR fixes any issues,
[use a keyword](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue#linking-a-pull-request-to-an-issue-using-a-keyword)
when referring to the issue so it will be closed automatically when the PR is merged.
If there are no relevant issues, uncomment the below line:
None
-->

<!-- If generative AI was used to draft or generate any part of this PR:
1. Uncomment the header "## Generative AI Usage Statement" below
2. Fill in a brief description of how generative AI was used
3. Disclose the tool and model version used

If generative AI was NOT used in any way, leave this entire section commented out.

## Generative AI Usage Statement

Include a brief description of how generative AI assistance was used to generate any of the code or content included in this PR,
e.g. writing code, writing unit tests, troubleshooting problems, software design discussion, commit messages, or preparing the PR.
Disclose the tool and model version used, e.g. Claude Haiku 4.5, GPT-5.6, Copilot w/ Claude Sonnet 5, etc.
-->

## PR Checklist

(~Strikethrough~ any points that are not applicable. Check boxes only, do not append descriptions.)

- [ ] This comment contains a description of changes with justifications, with any relevant issues linked.
- [ ] Write unit tests for any new features, bug fixes, or other code changes.
- [ ] Update the docs if there are any API changes (roxygen2 comments, vignettes, readme, etc.).
- [ ] Update `NEWS.md` with a short description of any user-facing changes and reference the PR number. Follow the style described in <https://style.tidyverse.org/news.html>
- [ ] Run `devtools::check()` locally and fix all notes, warnings, and errors.
