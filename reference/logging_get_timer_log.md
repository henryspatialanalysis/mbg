# Get timer log

Return a log of all timed events as a data.table

## Usage

``` r
logging_get_timer_log(clear_log = FALSE, deindent = TRUE)
```

## Arguments

- clear_log:

  (`logical(1)`, default FALSE) Should the log be cleared afterwards?

- deindent:

  (`logical(1)`, default TRUE) Should leading whitespace be removed from
  timer messages?

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
#> 1: 17.219 17.321 Test logging Test logging: 0.102 sec elapsed   0.102
```
