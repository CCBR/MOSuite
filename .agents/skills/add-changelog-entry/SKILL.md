---
name: add-changelog-entry
description: >-
  Add or update a project's changelog entry in NEWS.md or CHANGELOG.md. Use for user-facing features and bug fixes, or when asked to edit, draft, or review a changelog. Do not use for developer-only changes.
metadata:
  quiver-version: "1.0.0"
  quiver-authors: "kelly-sovacool"
  quiver-reviewers: "kopardev"
---

# Add Changelog Entry

## When to use this skill

Use this skill when a project introduces a user-facing feature or bug fix that needs a changelog entry. Also use it when the user explicitly asks to add, update, draft, or review a project's changelog.

## Instructions

1. Confirm the change is user-facing. Do not add entries for developer-only work such as CI updates, development notes, or internal documentation changes.
2. Identify the changelog file. Follow repository instructions and documentation first, then use an existing changelog file and its history. If no convention or file exists, use root `NEWS.md` for an R package (for example, one with a `DESCRIPTION` file) and root `CHANGELOG.md` for other project types. If both files exist or the signals conflict, follow explicit repository guidance or ask which file to use; do not create a second changelog by default.
3. Inspect the repository's current state with `git status` and `git diff` (plus `git diff --staged` when changes are staged), then read the selected changelog before editing. Preserve existing user changes, follow its section placement, and check whether the same change is already listed. If no changelog exists, stop and ask the user whether to create one; only create the file after they approve, using a minimal structure and development section consistent with repository conventions.
4. Determine the pull request number and author from the current branch. The PR number is not known up front, and the entry is often written before the PR exists, so look it up rather than asking first:
   - When the GitHub CLI is available and authenticated (`gh auth status`), run `gh pr list --head "$(git branch --show-current)" --state all --json number,author` to find the branch's PR, including merged and closed ones. (`gh pr view --json number,author` is equivalent for the common case of one open PR on the current branch.) Disable the pager, for example with `GH_PAGER=cat`, so output is captured rather than opened interactively.
   - `gh api user --jq .login` returns the authenticated user's username, which is the author when they are opening the PR.
   - If more than one PR matches the branch, ask the user which one the entry belongs to.
   - If the branch has no PR yet, `gh pr list` returns an empty array and `gh pr view` fails with `no pull requests found for branch`. Ask the user to add the entry after the PR is opened. Do not infer the next number from recent issues or PRs, and do not leave a placeholder in a finished entry.
5. Add or update one entry per user-facing change under the development/unreleased section, never under a released section. A pull request that both adds a feature and fixes a bug gets one entry for each. Use the repository's heading convention; if none exists, use `## development version`, unless repository instructions specify otherwise. Write every new entry in the format shown under Examples, even when older entries in the file use a different format; leave those older entries alone rather than reformatting them. Attribute each entry with the actual pull request number, not a related issue number, and the author's GitHub username.
6. Write each entry as a single line describing the user-visible effect, not the implementation. Use present-tense imperative phrasing, wrap code identifiers in backticks, and avoid internal details such as refactored helper names or file paths.
7. Review the resulting diff to confirm the entry is accurate, not duplicated, and limited to the intended changelog change. Report the file edited and the entry added. Do not commit changes unless the user explicitly asks and approves.

## Examples

For an R package, add the entry to `NEWS.md`; for other project types, use `CHANGELOG.md` unless repository conventions specify otherwise. Do not add an entry for a developer-only change, such as updating a CI workflow.

### Entry format

Write every new entry in this format, replacing each placeholder with a real value:

```md
- <user-visible change>. (#<pr-number>, @<author-username>)
```

### Good examples

A concise entry describes the user-visible effect and carries PR and author attribution:

```md
- Fix `filter_counts()` so selected items remain in the results. (#456, @octocat)
```

A pull request with two user-facing changes gets two entries, both citing that pull request:

```md
- Add `plot_volcano()` for visualizing differential expression results. (#457, @octocat)
- Fix `filter_counts()` so selected items remain in the results. (#457, @octocat)
```

### Bad examples

Do not write entries like these:

```md
- Refactored the internal `.apply_mask()` helper in R/filter.R.   # implementation detail, not user-facing
- Fix filtering bug. (#455, @octocat)                             # #455 is the issue, not the PR
- Fix filtering bug. (#TBD, @octocat)                             # placeholder instead of a real number
```

## Safety and limitations

- Add entries only for user-facing changes. When the changelog file or change classification is unclear, ask before editing.
- Do not create a changelog file on your own initiative. If the repository has none, ask the user first and create one only with their approval.
- `gh` may be missing, unauthenticated, or pointed at a different remote. Treat its failure as "unknown" and ask the user rather than guessing a PR number or username.
- Do not replace, reorder, or rewrite unrelated changelog content. Preserve existing local edits and update a matching entry rather than creating a duplicate.
- This skill edits a project's changelog; it does not publish a release.
