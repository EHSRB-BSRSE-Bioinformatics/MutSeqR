# MutSeqR: Development, Release & Versioning Guide

How MutSeqR is developed, versioned, and released. MutSeqR is a
**Bioconductor-managed package**: it lives in two repositories and the
release process is partly run by the Bioconductor team, so a few rules
here are Bioconductor policy, not house style. Read §3 and §5 before your
first push to Bioconductor.

## 1. The two repositories

| Repo | Remote | URL | Role |
|---|---|---|---|
| GitHub | `origin` | `github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR` | Development: branches, PRs, CI, pkgdown, GitHub Releases |
| Bioconductor | `upstream` | `git@git.bioconductor.org:packages/MutSeqR.git` (SSH) | Canonical package; source of `BiocManager::install()` |

```
feature branch ─PR─> github devel ─PR─> github main ─push─> Bioc devel
                                                              │
                                                              │ twice a year, done BY
                                                              │ THE BIOCONDUCTOR TEAM
                                                              ▼
                                                    Bioc RELEASE_3_XX (even y)
```

Key facts about the Bioconductor side (see
[Git version control](https://bioconductor.org/developers/package-guidelines/git-version-control.html)):

- Only `devel` and the **current** `RELEASE_X_Y` branch exist and are
  writable. No new branches, no tags, no force-pushes on that server.
- **Releases are not git tags.** A release = the `RELEASE_X_Y` branch plus
  the Bioc team's version-bump commits. (Our GitHub tags are for local
  traceability only; they do nothing for the Bioc process.)
- Pushing is only accepted as a fast-forward. If you're behind, merge
  `upstream/devel` back first (§7).

## 2. Branch model

| Branch | Where | Version | Purpose |
|---|---|---|---|
| `devel` | GitHub | tracks `main` | Feature staging: feature branches PR here, then here → `main`. Keep it current with `main` (`git merge main` into it after merges) |
| `main` | GitHub **and** Bioc `devel` (kept in sync) | odd `y` (e.g. `1.1.z`) | The push line: 1:1 with Bioc `devel`. Everything user-facing-next-release flows through here |
| `RELEASE_3_XX` | Bioc | even `y` (e.g. `1.0.z`) | Current release. **Cherry-picks only** (§6). No feature work, no merges |
| `v1.0.0-Bioc3.23` | GitHub tag + Release | — | Traceability marker per Bioc release; triggers the pkgdown workflow |

Consequences:

- End users install from Bioconductor, so "clean `main`" is for *us*: it is
  exactly the state we push upstream.
- There is no GitHub mirror of `RELEASE_X_Y`; the release state is visible
  via the tag and on the Bioc server.

## 3. Versioning (Bioconductor policy)

From the [version numbering rules](https://bioconductor.org/developers/package-guidelines/versionnum.html):

- Versions are `x.y.z`. **`y` is even in release, odd in `devel`.**
  `0.99.z` was the pre-1.0 signal; MutSeqR's first Bioc release (3.23)
  converted `0.99.11` → **1.0.0** (release) / **1.1.0** (devel).
- **Authors increment `z` by 1** for work pushed to the devel line.
  **Commits without a `z` bump do not propagate** to
  `BiocManager::install()`.
- **`y` (and `x`) are bumped mechanically by the Bioc team at every
  release, for every package** — release → next even `y`, devel → next odd
  `y`, `z` reset to 0. There is no per-package judgment about whether a
  bump is "deserved": `x.y` encodes the release *cycle*, not the scale of
  change. (Bioc versions are not SemVer; `NEWS.md` + our release tags are
  the user-facing "what changed" channel.)
- The only author-driven bump: setting `y = 99` in devel signals "major
  changes" — the next release then converts `x.99.z` → `(x+1).0.0`. Use
  only when you actually ship breaking changes.
- We never bump `x`/`y` ourselves, except the one-time `y = 99` choice
  above. Our invariant to maintain: GitHub `main` == Bioc `devel` (they
  differ from the *release* by design — odd vs even `y`).
- ⚠️ The `x.x.9000` devel convention is a **CRAN** convention. Do **not**
  use it here — it violates the even/odd-`y` policy.

## 4. Everyday development

Environment: we use **pixi** for the R runtime, tooling, and all development
and test dependencies. Bioconductor packages are supplied by Bioconda, so the
native Pixi development environment supports Linux and Intel macOS. Tests
against Bioconductor devel use Bioconductor's official Docker image and work on
any platform where Pixi, Bash, and Docker are available.

```
pixi install        # one-time: installs R, tooling, package deps, and test data
pixi run test       # devtools::test()
pixi run document        # regenerate roxygen output
pixi run check-docs      # verify generated Rd, NAMESPACE, and metadata are current
pixi run check-readme    # verify README.md matches README.Rmd
pixi run check-reference # verify exported topics are indexed and build reference
pixi run check           # full R CMD check (fails on errors, warnings, or notes)
pixi run bioccheck       # BiocCheck (fails the task when checks fail)
pixi run validate        # run all local documentation and package validation gates

# Official Bioconductor devel environment
pixi run bioc-devel-build     # build/update the local image and its dependencies
pixi run bioc-devel-info      # show the image's R and Bioconductor versions
pixi run bioc-devel-test      # run the testthat suite against Bioconductor devel
pixi run bioc-devel-check     # run R CMD check against Bioconductor devel
pixi run bioc-devel-bioccheck # run BiocCheck against Bioconductor devel
pixi run bioc-devel-validate  # run all three devel validation commands
pixi run bioc-devel-shell     # open an interactive shell in the devel image
```

The first devel command builds `mutseqr-bioc-devel:local` automatically. Source
files are mounted into the container, so normal edits do not require a rebuild.
Run `pixi run bioc-devel-build` after changing dependencies in `DESCRIPTION` or
when you want to pull the newest `bioconductor/bioconductor_docker:devel` base
image. The native Pixi and container environments remain independent; switching
between their commands does not create or destroy either one.

`check-docs`, `check-readme`, and `check-reference` are also run on every pull
request by the `Documentation validation` workflow (`docs-and-readme` and
`reference` jobs). CI checks generated files from a clean checkout; locally, each
freshness task compares the relevant output before and after that run, so unrelated
worktree edits do not affect its result. If a check reports changed generated
files, review and commit those outputs, then rerun the task. Repository branch
protection can require these workflow checks alongside the existing `R-CMD-check`
matrix.

Feature flow:

1. `git checkout -b feature/... main` (or off `devel` — either is fine, PRs
   land on `devel`)
2. Develop; PR → `devel` (CI: R-CMD-check on all PRs).
3. PR `devel` → `main` when the batch is ready to ship.
4. After the merge: bump `devel` current with `main`.
5. Push the batch to Bioconductor (§5).

`NEWS.md`: keep the top entry at or ahead of the current devel version.
One entry per release at minimum; the release ritual (§5) enforces it.

## 5. Pushing work to Bioconductor

Do this after merging to `main`, when the batch is ready:

```
pixi run prepare-release        # documentation + package checks -> z-bump -> NEWS check -> commit
git push origin main
git push upstream main:devel    # 'upstream' = git.bioconductor.org; writes to its 'devel'
```

Verify it landed:

- Next day's [devel build report](https://bioconductor.org/checkResults/devel/bioc-LATEST/)
  shows your commit; `BiocManager::install("MutSeqR", version = "devel")`
  gives you the new `1.1.z`.
- If the push is rejected (non-fast-forward): `git fetch upstream`,
  `git merge upstream/devel`, resolve, then push again.

`prepare-release` intentionally stops short of the `upstream` push —
pushing to Bioconductor is a deliberate human act.

## 6. Hotfixes to the released version

Users of the current release (e.g. 1.0.0 on Bioc 3.23) do **not** see the
devel line. A user-facing hotfix means cherry-picking to the release
branch:

1. Fix on the devel line as usual (PR → `main`, `z` bump, push to
   `upstream devel`) so the fix also rides the next release.
2. Cherry-pick the fix commit (not the version-bump commit) onto
   `upstream/RELEASE_3_XX`, bump its `z` separately (e.g. `1.0.0` →
   `1.0.1`), commit, `git push upstream RELEASE_3_XX`.
3. Only the **current** release branch is writable; for older releases,
   email <maintainer@bioconductor.org>.
4. Never merge `devel` → release; cherry-picks only.
5. Note the fix in `NEWS.md` under the release version.

## 7. When Bioconductor cuts a release (≈ twice a year)

The Bioc team branches `RELEASE_3_XX` and bumps versions on both their
branches (devel: even→odd `y`, e.g. `1.1.0` → `1.3.0` for the cycle after
3.24… see the table in the Bioc docs). Our ritual afterwards:

1. `git fetch upstream`
2. Merge `upstream/devel` into local `main` (usually fast-forward) →
   `git push origin main`.
3. Tag + GitHub Release on the *released* commit:
   `v1.2.0-Bioc3.24` → release notes in NEWS → triggers pkgdown.
4. Update pins for the new Bioc cycle:
   - `pixi.toml`: update `r-base` when Bioconda publishes packages for the new
     Bioconductor/R cycle
   - `pixi run bioc-devel-build`: pull and rebuild the moving devel image
5. Bump `devel` staging branch current with `main`.

## 8. Quick reference / troubleshooting

| Symptom | Cause / fix |
|---|---|
| `git push upstream main:devel` rejected | You're behind Bioc: `git fetch upstream && git merge upstream/devel`, then push |
| Change pushed but `BiocManager::install()` doesn't see it | No `z` bump in DESCRIPTION (Bioc only serves bumped commits) |
| `ssh -T git@git.bioconductor.org` → permission denied | Register your SSH key in the [BiocCredentials app](https://bioconductor.org/help/bioc-credentials/); a correct key prints your package list with `R W` next to MutSeqR |
| Pushing to `upstream` a branch other than `devel`/`RELEASE_3_XX` | Not possible by policy — do it via GitHub + the flows above |
| "Should I use 1.2.9000 on devel?" | No — that's CRAN. Here: odd `y`, `z` bumps (§3) |
| R version / Bioc version mismatch in pixi | `pixi.toml` pins `r-base` to the current Bioc cycle's R; update at each release (§7) |

## 9. Open items

- Pull requests run R CMD check, generated documentation and README freshness,
  and pkgdown reference coverage. The pkgdown site continues to publish on
  release; pull requests validate the reference without publishing it.
- Optional CI gate: PR touching `R/` must change `DESCRIPTION` Version and
  `NEWS.md`.
