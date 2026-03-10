# Start logging timer

Start a nested timer with an optional message

## Usage

``` r
logging_start_timer(msg, echo = TRUE, indentation_text = "  ")
```

## Arguments

- msg:

  (`character(1)`) Logging message

- echo:

  (`logical(1)`, default TRUE) Should the message be written to screen?

- indentation_text:

  (`character(1)`, default " ") Text that will be repeated at the
  beginning of the message for each layer of indentation

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
#> 1: 17.869 17.971 Test logging Test logging: 0.102 sec elapsed   0.102
#> 2: 18.143 18.245 Test logging Test logging: 0.102 sec elapsed   0.102
```
