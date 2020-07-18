---
title: "CRAN comments"
---

## Platform
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* win-builder (devel and release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results
There were no ERRORs or WARNINGs. 

## Downstream dependencies
There were no ERRORs or WARNINGs. 

## Responses to previous issues (CRAN submission LAWBL 1.1.0 - 07/18 from Martina Schmirl):
1. *\\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \\dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please unwrap the examples if they are executable in < 5 sec, or replace \\dontrun{} with \\donttest{}.*

Thanks for reminding. \\dontrun{} have been replaced with \\donttest{} (the examples are executable in > 5 sec)

2. *You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions)*

All texts to the console have been removed (except for print, summary, interactive fundtions).

3. *Please do not modify the global environment in your functions. This is not allowed by the CRAN policies.*


I didn't modify the global environment in pcfa(). It's recorded at the beginning, as:

```
    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
```

Then recovered at the end, as:

```
    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)
```

4. *Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos. e.g.:*
```
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)
e.g. vignette.
```

Thanks for reminding. The following commands are added in vignette

```
oldmar <- par("mar")
...
par(mar = oldmar) #reset to old mar
```
     


