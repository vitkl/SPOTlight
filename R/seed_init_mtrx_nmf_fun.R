#' This functions seeds initialization matrices H and W to perform NMF
#'
#' @param counts_sc Object of class Matrix of shape GENESxCELL, se_obj\@assays$RNA@counts.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @return This function returns a list with the initialized matrices H and W.
#' @export
#' @examples
#'


seed_init_mtrx_nmf <- function(cluster_markers,
                               counts_sc,
                               clust_vr,
                               ntop = NULL) {

  #### Get dataset ready ####
  # Comment because the matrix is already sparse
  # se_nmf_ready <- prep_seobj_topic_fun(se_obj = se_obj)

  # Rank of the model equals the number of cell types
  k <- length(unique(clust_vr))

  # Select all marker genes for each cluster AND compute their Z score
  if (is.null(ntop)) ntop <- max(table(cluster_markers$cluster))
  cluster_markers_cut <- suppressMessages(cut_markers2(markers = cluster_markers,
                                                       ntop = ntop))

  # Select unique markers from each cluster, if there are common markers between clusters lda model gets confused and classifies very different clusters as belonging to the same topic just because the seeding induced it!
  cluster_markers_uniq <- lapply(unique(cluster_markers_cut$cluster), function(clust) {
    ls1 <- cluster_markers_cut[cluster_markers_cut$cluster == clust, "gene"]
    ls2 <- cluster_markers_cut[cluster_markers_cut$cluster != clust, "gene"]
    ls1_unique <- ls1[! ls1 %in% ls2]

    return(cluster_markers_cut[cluster_markers_cut$cluster == clust & cluster_markers_cut$gene %in% ls1_unique, ])
  }) %>%
    bind_rows()

  # Set seedwords from top markers. Here we are setting the weights for each topic, the words that are weighted positively are those belonging to the list of top markers for a cluster.
  # In the seedgenes matrix each row represents a topic and each column represents a gene.

  # initialize matrix
  seedgenes <- matrix(nrow = k, ncol = nrow(counts_sc), data = 1e-10)
  colnames(seedgenes) <- rownames(counts_sc)

  # From cluster_markers_uniq remove genes not present in seedgenes
  cluster_markers_upd <- cluster_markers_uniq[cluster_markers_uniq$gene %in% colnames(seedgenes), ]

  # Add seeds to model, if a cluster-topic has 0 unique markers its row will be set to all 0
  for (i in seq_len(k)) {
    # print(i)
    clust_row <- cluster_markers_upd$cluster == as.character(unique(clust_vr))[[i]]

    # Genes to initialize
    gene_init <- as.character(cluster_markers_upd[clust_row, "gene"])

    seedgenes[i, gene_init] = cluster_markers_upd[clust_row, "weight"]
  }
  # Transpose to get in GENExTOPIC
  W <- t(seedgenes)

  ###################
  #### Seeding H ####
  ###################
  H <- matrix(data = 1e-10,
              nrow = k,
              ncol = ncol(counts_sc))

  for (i in seq_len(length(clust_vr))) {
    # Extract topic row
    h_row <- which(unique(clust_vr) == clust_vr[[i]])

    # Add 1 to that topic-cell relation
    H[h_row, i] <- 1
  }

  # get matrix H
  colnames(H) <- colnames(counts_sc)

  return(list("W" = W, "H" = H))
}
