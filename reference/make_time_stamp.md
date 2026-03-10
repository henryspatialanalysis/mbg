# Make time stamp

Create a string time stamp based on current detailed date/time

## Usage

``` r
make_time_stamp(suffix = NULL, milliseconds = TRUE)
```

## Arguments

- suffix:

  (`character(1)`, default NULL) suffix to append to the time stamp.
  Useful when running batches of related models

- milliseconds:

  (`logical(1)`, default TRUE) Should milliseconds be appended to the
  timestamp? Useful when launching many models in quick succession.

## Value

A string formatted as
`'YYYYMMDD_HH_MM_SS(_optional MS)(_optional suffix)'`
