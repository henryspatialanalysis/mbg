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
#> Test logging: 0.101 sec elapsed
log_results <- mbg::logging_get_timer_log()
print(log_results)
#>       tic    toc          msg                    callback_msg elapsed
#>     <num>  <num>       <char>                          <char>   <num>
#> 1: 17.869 17.971 Test logging Test logging: 0.102 sec elapsed   0.102
#> 2: 18.143 18.245 Test logging Test logging: 0.102 sec elapsed   0.102
#> 3: 18.412 18.513 Test logging Test logging: 0.101 sec elapsed   0.101
```
