# mosuite cli

    Code
      system(command)

# cli_exec --json --debug

    Code
      print(cli_exec(c("create_multiOmicDataSet_from_files", paste0("--json=\"",
        system.file("extdata", "example.json", package = "MOSuite"), "\""), "--debug")))
    Output
      MOSuite::create_multiOmicDataSet_from_files(feature_counts_filepath = "inst/extdata/RSEM.genes.expected_count.all_samples.txt.gz", 
          sample_meta_filepath = "inst/extdata/sample_metadata.tsv.gz")

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

