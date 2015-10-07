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
library(SSMinR)
```

# Before compiling your model

`SSMinR` compiles your model by calling `SSM` from your `R` session using the function `system()`. As such, when compilating, `SSM` asks `R` where to find the header files as well as the librairy dependencies needed.

Outside of `R`, these information are generally stored in the environment variables `CPATH` (path to header files) and `LIBRARY_PATH` (path to librairy dependencies). You can check this in your terminal by typing:

```sh
echo $CPATH
echo $LIBRARY_PATH
```

If these environment variables don't exist it's probably a good idea to [set them](http://unix.stackexchange.com/questions/117467/how-to-permanently-set-environmental-variables).

However, even if they exist in your system, the environment variables `LIBRARY_PATH` and `CPATH` might not be imported in your `R` session, which can result in a compilation failure. For instance, if `SSM` doesn't have the right path to the header files when called within `R` you will obtain the following error at compilation:

```sh
Building the model...
gcc  -std=gnu99 -O3 -DGSL_RANGE_CHECK_OFF -I /Users/ntncmch/.ssm/include -o Ht.o -c Ht.c
[91mFAIL[0m: In file included from Ht.c:19:
/Users/ntncmch/.ssm/include/ssm.h:36:10: fatal error: 'gsl/gsl_math.h' file not found
[91mFAIL[0m: #include <gsl/gsl_math.h>
         ^
[91mFAIL[0m: 1 error generated.
[91mFAIL[0m: make: *** [Ht.o] Error 1
[91mFAIL[0m: could not build the model (2).
```

If you experience a similar issue, check whether the environment variables are imported in your `R` session by typing 

```r
Sys.getenv("CPATH")
Sys.getenv("LIBRARY_PATH")
```
If one or both variables are missing or doesn't indicate the appropriate path, you will need to set them permanently in your `~/.Rprofile` (read this [post](http://www.r-bloggers.com/fun-with-rprofile-and-customizing-r-startup/) if you don't have one yet) by adding the following lines:

```r
Sys.setenv(CPATH="/usr/local/include") # header path
Sys.setenv(LIBRARY_PATH="/usr/local/lib") # library path
```
_Note that I have set the paths to my default ones and it might not be the same on your system.
_




