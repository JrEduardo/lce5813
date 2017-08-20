# Teaching and Learning Bayesian Inference and Analysis #

> [Eduardo Elias Ribeiro Junior][eduardo] (LCE-ESALQ USP / LEG-UFPR)

## Description ##

The `lce5813` package contains functions to support teaching and
learning Bayesian inference and analysis. The package name (`lce5813`)
comes from code for the course _Introduction to Bayesian Inference_ in
[Statistics and Agricultural Experimentation][ppgeea] program at
University of São Paulo (USP), Brazil.

## Download and install ##

You can install automatically from the GitHub repository using
[`devtools`][devtools]. Just run the code below in a R session.

```r

# install.packages("devtools")
devtools::install_git("git@github.com:jreduardo/lce5813.git")
devtools::install_github("jreduardo/lce5813") # For short.

```

After installing the package, you can load and explore its contents with

```r

# Load and explore the package.
library(lce5813)
ls("package:lce5813")
help(package = "lce5813")

```

## License ##

The `lce5813` package is licensed under the
[GNU General Public License, version 3], see the [`LICENSE`](./LICENSE)
file.

<!--------------------------------------------- -->
<!-- Links -->

[eduardo]: https://jreduardo.github.io/
[ppgeea]: http://www.esalq.usp.br/pg/programas/estatistica/
[devtools]: https://github.com/hadley/devtools
