# Prerelease packages

The workflow `.github/workflows/prerelease.yml` builds conda packages of xgem2mt.
It also builds Optima and GEMS3K from their development branches.

## Make a prerelease

1. Set the version in `CMakeLists.txt` (`project(xGEM2MT VERSION ...)`).
2. Push a tag that starts with `v`, for example `git tag v1.0.4 && git push origin v1.0.4`,
   or start the workflow by hand (Actions > prerelease > Run workflow).
3. The workflow builds the packages for Linux, macOS and Windows and publishes them as a conda channel
   on the `gh-pages` branch: https://gemshub.github.io/xgem2mt/prerelease. Each run replaces the channel.
4. For a tag, the workflow also makes a GitHub pre-release on the tag with the packages attached.

Once, after the first run: Settings > Pages > Source = branch `gh-pages`.

## Install a prerelease

```bash
conda create -n xgem2mt -c https://gemshub.github.io/xgem2mt/prerelease -c conda-forge xgem2mt
```

## Recipes

| Folder | Package |
|---|---|
| `conda/optima` | Optima, static library. It is linked into GEMS3K and is not released. |
| `conda/gems3k` | GEMS3K with the Optima solver. It is released, because xgem2mt uses the shared GEMS3K library (Linux, macOS). On Windows the Python module contains a static copy of GEMS3K. |
| `conda/xgem2mt` | xgem2mt and its Python module, one package for each Python version. |

The source branches of Optima and GEMS3K are set at the top of the workflow.
