# End logging timer

End a nested timer

## Usage

``` r
logging_stop_timer(echo = TRUE)
```

## Arguments

- echo:

  (`logical(1)`, default = TRUE) Should the message be written to
  screen?

## Examples

``` r
mbg::logging_start_timer(msg = 'Test logging')
#> Test logging
Sys.sleep(0.1)
mbg::logging_stop_timer()
#> Test logging: 0.102 sec elapsed
log_results <- mbg::logging_get_timer_log()
print(log_results)
#>       tic    toc          msg                    callback_msg elapsed
#>     <num>  <num>       <char>                          <char>   <num>
#> 1: 16.987 17.089 Test logging Test logging: 0.102 sec elapsed   0.102
#> 2: 17.248 17.350 Test logging Test logging: 0.102 sec elapsed   0.102
#> 3: 17.505 17.607 Test logging Test logging: 0.102 sec elapsed   0.102
```
