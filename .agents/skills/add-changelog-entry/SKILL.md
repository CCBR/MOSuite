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
3. Inspect the repository's current state and read the selected changelog before editing. Preserve existing user changes, follow its format and section placement, and check whether the same change is already listed. If no changelog exists, create the appropriate file with a minimal structure and development/unreleased section consistent with repository conventions.
4. Determine the pull request number and author from the current branch before writing the entry. The PR number is not known up front, so look it up with the GitHub CLI when it is available and authenticated (`gh auth status`) rather than asking first:
   - Get the branch name with `git branch --show-current`.
   - Look up the PR for that branch: `gh pr list --head "$(git branch --show-current)" --state all --json number,author`. This matches the branch explicitly and also finds merged or closed PRs. Equivalently, `gh pr view --json number,author --jq '"#\(.number), @\(.author.login)"'` infers the current branch's PR when run inside the repository.
   - If more than one PR matches the branch, ask the user which one the entry belongs to.
   - `gh api user --jq .login` returns the authenticated user's username, which is the author when they are opening the PR.
   - Disable the pager (for example, `GH_PAGER=cat`) so output is captured rather than opened interactively.
5. If the branch has no PR yet, `gh pr list` returns an empty array and `gh pr view` fails with `no pull requests found for branch`. In that case, ask the user for the number, or add the entry after the PR is opened. Do not infer the next number from recent issues or PRs.
6. Add or update one concise entry under the development/unreleased section, never under a released section. Use the repository's heading convention; if none exists, use `## development version`, unless repository instructions specify otherwise. Include the actual pull request number, not a related issue number, and the author's GitHub username in the repository's attribution format, defaulting to `(#123, @username)`. If either value is unknown, ask the user; do not guess or use placeholders in a finished entry.
7. Review the resulting diff to confirm the entry is accurate, not duplicated, and limited to the intended changelog change. Report the file edited and the entry added. Do not commit changes unless the user explicitly asks and approves.

## Examples

For example, a concise entry with PR and author attribution might read:

```md
- Fix filtering so selected items remain in the results. (#123, @username)
```

For an R package, add the entry to `NEWS.md`; for other project types, use `CHANGELOG.md` unless repository conventions specify otherwise. Do not add an entry for a developer-only change, such as updating a CI workflow.

## Safety and limitations

- Add entries only for user-facing changes. When the changelog file or change classification is unclear, ask before editing.
- Never invent PR numbers or usernames. Resolve them with `gh` when possible; otherwise ask for either missing value before completing the entry.
- `gh` may be missing, unauthenticated, or pointed at a different remote. Treat its failure as "unknown" and ask the user rather than falling back to a guess.
- Do not replace, reorder, or rewrite unrelated changelog content. Preserve existing local edits and update a matching entry rather than creating a duplicate.
- This skill edits a project's changelog; it does not publish a release.
