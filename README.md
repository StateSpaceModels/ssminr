# R interface to the SSM library

[SSM](https://github.com/JDureau/ssm) is a C library to perform inference for time series analysis with State Space Models, like playing with duplo blocks. `SSMinR` is a [R](http://cran.r-project.org) package that provides an interface to SSM.

# Installation

Before installing `SSMinR`, make sure you have installed the SSM library by following these [instructions](https://github.com/JDureau/ssm). Note that SSM only works on Unix machines. If you run Windows, you will need to set up a virtual machine.

Assuming that SSM and R are installed on your machine, the easiest way to install `SSMinR` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("ntncmch/SSMinR")
```

To work with `SSMinR`: fire up `R`, and type 

```r
library(help=SSMinR)
```
