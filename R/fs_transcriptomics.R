#' Feature Selection for transcriptomics
#'
#' @param IDs
#' @param data_IDs
#' @param transcriptomics_data
#'
#' @return
#' @export
#'
#' @examples
fs_transcriptomics <-
  function(train_IDs,
           test_IDs,
           data_IDs,
           phenotype_IDs,
           transcriptomics_data,
           gene_anotation) {
    # To add: lowest_level_pathways
    # Working with gene_IDs, file from Phong, with annotation for the transcriptomics
    #
    # Getting the ID sets from the transcriptomics
    train_transcriptomics_IDs <-
      train_IDs$`Feature Selection IDs`$Transcriptomics
    test_transcriptomics_IDs <-
      test_IDs$`Feature Selection IDs`$Transcriptomics

    # Frame of phenotypes, with all used transcriptomic IDs and their inc3 value
    samples <- unlist(train_transcriptomics_IDs) %>% unique()
    tran_IDs <-
      data_IDs[match(samples, data_IDs$Clinical),] %>% select("Transcriptomics") %>% unlist()

    info <-
      phenotype_IDs[as.character(samples), "inc3", drop = FALSE]
    rownames(info) <- tran_IDs

    # return(info)

    # Creating a model matrix
    tt <- factor(c(0, 1))
    modmatrix = info %>% transmute(inc30 = ifelse(inc3 <= 0, 1, 0),
                                   inc31 = ifelse(inc3 <= 1, 1, 0))

    # return(modmatrix)

    # Creating a list for the gsea
    gsea = list(
      edge = list(),
      ilmn = list(),
      auc = list(),
      pathway = list()
    )

    # Setting the seed for the methode
    set.seed(123)

    # Going through each set of IDs

    for (i in 1:length(train_transcriptomics_IDs)) {
      # Selecting the data_IDs with one set of the data_partitioning
      clin_IDs <- unlist(train_transcriptomics_IDs[[i]])

      tmp_transcriptomics_IDs <-
        data_IDs %>% filter(Clinical %in% clin_IDs) %>%
        select(Transcriptomics) %>% as.list() %>% unlist()

      # return(tmp_transcriptomics_IDs)

      # Showcase for missing values
      #
      # tmp_data <-
      #   transcriptomics_data %>% as.data.frame() %>%
      #   select(suppressWarnings(one_of(
      #     as.character(tmp_transcriptomics_IDs)
      #   )))

      tmp_data <-
        transcriptomics_data[,!is.na(match(colnames(transcriptomics_data), tmp_transcriptomics_IDs))]

      tmp_mod <- modmatrix[colnames(tmp_data), , drop = FALSE]


      # if (any(table(tmp_mod[, 1]) < 2) |
      #     any(table(tmp_mod[, 2]) < 2)) {
      #   print("error")
      #   next
      # }

      # return(tmp_data)

      # return(tmp_mod)

      fit <- lmFit(tmp_data, tmp_mod)


      # return(fit)

      # return(colnames(coef(fit)))

      contrast <-
        makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

      tmp <- contrasts.fit(fit, contrast)

      tmp <- eBayes(tmp)

      topde <-
        topTable(tmp, sort.by = "P", n = Inf) %>% mutate(Probe = rownames(.)) %>%
        mutate(Name = gene_anotation$symbol[match(.$Probe, gene_anotation$Probe_Id)],
               EntrezID = gene_anotation$EntrezID[match(.$Probe, gene_anotation$Probe_Id)])

      return(topde)


      # # Sort the column names (will be removed eventually)
      # tmp_transcriptomics <-
      #   tmp_transcriptomics[, order(names(tmp_transcriptomics))]
    }



  }
