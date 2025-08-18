# mosuite cli

    Code
      system(command)

# cli_exec --json --debug

    Code
      cli_exec(c("create_multiOmicDataSet_from_files", paste0("--json=\"",
        system.file("extdata", "example.json", package = "MOSuite"), "\""), "--debug"))

# cli_exec --help

    Code
      cli_exec("--help")

