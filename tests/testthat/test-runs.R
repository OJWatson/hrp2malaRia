test_that("sim runs", {
  expect_output(
  r <- hrp2malaRia::hrp2_Simulation(EIR = 1, N = 1000, years = 0.4)
  )
  expect_true("Buffer" %in% names(r))
  
  expect_output(
    r <- hrp2malaRia::hrp2_Simulation(EIR = 1, N = 1000, years = 0.4, 
                                      just_storage_results = TRUE)
  )
  expect_false("Buffer" %in% names(r))
})
