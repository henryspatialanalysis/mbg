test_that("logging_start_timer works correctly", {
  # Test basic functionality
  expect_silent(logging_start_timer("Test message"))
  expect_silent(logging_stop_timer())

  # Test with echo = FALSE
  expect_silent(logging_start_timer("Silent message", echo = FALSE))
  expect_silent(logging_stop_timer())

  # Test indentation
  expect_silent(logging_start_timer("First level"))
  expect_silent(logging_start_timer("Second level"))
  expect_silent(logging_stop_timer())
  expect_silent(logging_stop_timer())
})

test_that("logging_stop_timer works correctly", {
  # Start a timer
  logging_start_timer("Test timer")

  # Test basic stop
  expect_silent(logging_stop_timer())

  # Test with echo = FALSE
  logging_start_timer("Silent stop test")
  expect_silent(logging_stop_timer(echo = FALSE))
})

test_that("logging_get_timer_log works correctly", {
  # Clear any existing logs
  logging_get_timer_log(clear_log = TRUE)

  # Create some test timers
  logging_start_timer("Test timer 1")
  Sys.sleep(0.1)
  logging_start_timer("Test timer 2")
  Sys.sleep(0.1)
  logging_stop_timer()
  logging_stop_timer()

  # Get the log
  log <- logging_get_timer_log()

  # Check log structure
  expect_true(inherits(log, "data.table"))
  expect_true("msg" %in% names(log))
  expect_true("elapsed" %in% names(log))
  expect_true(nrow(log) >= 2)

  # Test clear_log parameter
  log2 <- logging_get_timer_log(clear_log = TRUE)
  expect_true(nrow(log2) >= 2)
  log3 <- logging_get_timer_log()
  expect_true(nrow(log3) == 0)

  # Test deindent parameter
  logging_start_timer("  Indented message")
  logging_stop_timer()
  log4 <- logging_get_timer_log(deindent = TRUE)
  expect_false(any(grepl("^  ", log4$msg)))
  logging_get_timer_log(clear_log = TRUE)
})

test_that("logging functions handle errors gracefully", {
  # Test stopping timer without starting
  expect_silent(logging_stop_timer())

  # Test getting log with no timers
  expect_true(inherits(logging_get_timer_log(), "data.table"))
  expect_true(nrow(logging_get_timer_log()) == 0)

  # Test invalid parameters
  expect_error(logging_start_timer(msg = NULL))
  expect_error(logging_start_timer(echo = "not logical"))
  expect_error(logging_start_timer(indentation_text = NULL))

  expect_error(logging_stop_timer(echo = "not logical"))

  expect_error(logging_get_timer_log(clear_log = "not logical"))
  expect_error(logging_get_timer_log(deindent = "not logical"))
})