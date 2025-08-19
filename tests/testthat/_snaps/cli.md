# mosuite cli

    Code
      system(command)

# mosuite --help

    Code
      cli_exec("--help")

---

    Code
      system(paste(system.file("exec", "mosuite", package = "MOSuite"), "--help"))

---

    Code
      cli_exec(c("filter_counts", "--help"))
    Output
      No documentation for 'filter_counts' in specified packages and libraries:
      you could try '??filter_counts'

