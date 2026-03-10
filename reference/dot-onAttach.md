# Behavior when attaching the mbg package

Behavior when attaching the mbg package

## Usage

``` r
.onAttach(libname, pkgname)
```

## Arguments

- libname:

  (character(1)) A character string giving the library directory where
  the package defining the namespace was found.

- pkgname:

  (character(1)) A character string giving the name of the package.

## Value

(invisible) A message may be printed to the console.

## Details

Yields a message if the INLA package namespace is not available.
