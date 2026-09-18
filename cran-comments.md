# cran-comments

## Test environments

* local Ubuntu 22.04, R 4.1.2
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macOS-latest (R-release), windows-latest (R-release)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michael Whiteside <mickwh@msn.com>'

New submission
```

This is a new submission.

The DOI in the Description field
(<doi:10.21203/rs.3.rs-9283100/v1>) points to the preprint describing the
method. It resolves, but the automatic checker may report it as "possibly
invalid" because the Research Square landing page is served behind a bot
filter.

## Notes for the reviewer

Examples, tests and the vignette use only data shipped in `inst/extdata`
(a 300-gene subset of the public GEO series GSE11923) and simulated series.
Nothing downloads at build, check or run time, and no examples are wrapped
in \donttest{} or \dontrun{}.

The package has no compiled code and imports only `stats`.
