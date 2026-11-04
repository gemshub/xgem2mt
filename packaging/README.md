# Prerelease packages

The workflow `.github/workflows/prerelease.yml` builds conda packages of xgem2mt.
It also builds Optima and GEMS3K from their development branches.

## Channels

| Channel | Updated by |
|---|---|
| https://gemshub.github.io/xgem2mt/prerelease | a push to a branch listed under `on.push.branches` in the workflow (now `master` and `prerelease_ci`), a GitHub pre-release, or a manual run (Actions > prerelease > Run workflow) |
| https://gemshub.github.io/xgem2mt/stable | a published GitHub release |

Each run builds the packages for Linux (`linux-64`), macOS on Apple Silicon (`osx-arm64`) and Windows
(`win-64`), and replaces the packages in its channel; the other channel stays as it is.
A push that changes only `*.md` files or `docs/` does not start a build.

## Make a release

1. Set the version in `CMakeLists.txt` (`project(xGEM2MT VERSION ...)`) and merge it into `master`.
2. On GitHub, open Releases > Draft a new release. Create a tag `v<version>` (for example `v1.0.4`) on
   `master`, write the notes and publish. A tag that does not match the version in `CMakeLists.txt`
   stops the workflow; for a pre-release, `v<version>-<suffix>` (for example `v1.0.4-rc1`) is accepted.
3. The workflow updates the stable channel, attaches the packages to the release and adds the install
   notes below your text.

Once, after the first run: Settings > Pages > Source = branch `gh-pages`.

## Install

```bash
conda create -n xgem2mt -c https://gemshub.github.io/xgem2mt/stable -c conda-forge xgem2mt      # release
conda create -n xgem2mt -c https://gemshub.github.io/xgem2mt/prerelease -c conda-forge xgem2mt  # test build
```

## Recipes

| Folder | Package |
|---|---|
| `conda/optima` | Optima, static library. It is linked into GEMS3K and released with it, because the GEMS3K CMake targets need it to build against GEMS3K. |
| `conda/gems3k` | GEMS3K with the Optima solver. It is released, because xgem2mt uses the shared GEMS3K library (Linux, macOS). On Windows the Python module contains a static copy of GEMS3K. |
| `conda/xgem2mt` | xgem2mt and its Python module, one package for each Python version. |

The source branches of Optima and GEMS3K are set at the top of the workflow. Each run takes the head
commit of each branch once, and all platforms build that commit; the release notes list both commits.
