# ROADMAP Simulations

Simulation implementation and simulation results for the trial operating characteristics.

To build this site, you need to have a recent version of quarto install plus there is a dependency on the `roadmap.data` R package, see [https://github.com/maj-biostat/roadmap.data].

Also note, the site makes use of the development version of `data.table`.
Use `data.table::update_dev_pkg()` to install.

To render an individual notebook use, for example:

```
quarto render notebooks/sim-design5-results.qmd --to html
```

## Running simulations

Simulations are configured via Makefile. 
Run `make list` to see all targets. 

At the time of writing `make sim05` has the simulations for the current design.

## R, Rstudio, STAN, Quarto

Rstudio: Version 2024.04.2+764 (2024.04.2+764)

cmdstanr::cmdstan_version()
[1] "2.35.0"

To determine what version of quarto is being bundle with rstudio, open the terminal (inside rstudio) and type:

```
quarto check
```

As of the time of the last commit of this document, the version used is (all paths deleted):

```
Quarto 1.5.52
[✓] Checking versions of quarto binary dependencies...
      Pandoc version 3.2.0: OK
      Dart Sass version 1.70.0: OK
      Deno version 1.41.0: OK
      Typst version 0.11.0: OK
[✓] Checking versions of quarto dependencies......OK
[✓] Checking Quarto installation......OK
      Version: 1.5.52
      Path: -

[✓] Checking tools....................OK
      TinyTeX: v2024.03
      Chromium: (not installed)

[✓] Checking LaTeX....................OK
      Using: TinyTex
      Path: -
      Version: 2024

[✓] Checking basic markdown render....OK

[✓] Checking Python 3 installation....OK
      Version: 3.10.1
      Path: -
      Jupyter: (None)

      Jupyter is not available in this Python installation.
      Install with python3 -m pip install jupyter

[✓] Checking R installation...........OK
      Version: 4.4.1
      Path: -
      LibPaths:
        - -
      knitr: 1.47
      rmarkdown: 2.27

[✓] Checking Knitr engine render......OK
```

In theory, this should pick up whatever is on your system.
I think polyfill.io affected versions of quarto < 1.5 so you should probably update if you are not on the latest.

