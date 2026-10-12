## Test environments

* local: macOS 26.5.2, aarch64-apple-darwin20, R 4.4.1

## R CMD check results

0 errors | 0 warnings | 3 notes

The three notes are:

* `New submission` — this is the first release of survHC.
* `checking for future file timestamps ... unable to verify current time` — the
  local machine could not reach a time server. Unrelated to the package.
* `checking HTML version of manual ... Skipping checking HTML validation: no
  command 'tidy' found` — HTML Tidy is not installed locally.

## Notes

survHC implements the method of Kipnis, Galili and Yakhini (2026, Biometrika,
doi:10.1093/biomet/asaf075), with a companion application note at
doi:10.48550/arXiv.2610.07268.

The package depends only on base R, 'stats' and 'graphics'. 'survival' is
suggested rather than imported: a 'Surv' object is recognised structurally, so
the formula interface works without it, and the examples and tests that use
'survival' are guarded with requireNamespace().

Any spellings flagged in DESCRIPTION are author surnames or standard
statistical terminology, such as Berk-Jones and hypergeometric.
