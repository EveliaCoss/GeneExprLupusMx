# process_go skips files without go_id or with few GO terms

    Code
      res <- process_go(f_no_id, "UP", dirs$outdir, dirs$figdir)
    Message
      
      ------------------------------------
      Dataset: noid
      ------------------------------------
      No 'go_id' column found. Skipping file: <indir>/UP_GO_noid.RData

---

    Code
      process_go(f_few, "DOWN", dirs$outdir, dirs$figdir)
    Message
      
      ------------------------------------
      Dataset: few
      ------------------------------------
      Dataset has fewer than 40 unique GO terms. Skipping.

# process_go skips GO IDs that are not in org.Hs.eg.db

    Code
      suppressWarnings(process_go(f, "UP", dirs$outdir, dirs$figdir))
    Message
      
      ------------------------------------
      Dataset: bad
      ------------------------------------
      Create simMatrix
      preparing gene to GO mapping data...
      preparing IC data...
      Invalid simMatrix structure. Skipping dataset: bad

# reduce_go_terms processes every UP_GO_ and DOWN_GO_ file in a folder

    Code
      res <- reduce_go_terms(dirs$indir, dirs$outdir, dirs$figdir)
    Message
      
      ------------------------------------
      Dataset: first
      ------------------------------------
      Dataset has fewer than 40 unique GO terms. Skipping.
      
      ------------------------------------
      Dataset: second
      ------------------------------------
      Dataset has fewer than 40 unique GO terms. Skipping.

# reduce_go_similarity returns NULL with a message if rrvgo fails

    Code
      res <- reduce_go_similarity(bar_data, "test")
    Message
      Create simMatrix
      Failed to calculate simMatrix: boom
      Invalid simMatrix structure. Skipping dataset: test

# reduce_go_similarity returns NULL with a message if the reduction fails

    Code
      res <- reduce_go_similarity(bar_data, "test")
    Message
      Create simMatrix
      Reduce semantic similarity
      reduceSimMatrix() failed: boom
      No reduced terms could be computed.

