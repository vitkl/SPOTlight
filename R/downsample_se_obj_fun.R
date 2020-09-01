#' This function downsamples the number of cells and genes used to train the model
#'
#' @param se_obj Object of class Seurat; with the data of interest.
#' @param clust_vr Object of calss character; Name of the variable containing the cell clustering.
#' @param cluster_markers Object of class data.frame; obtained from the function Seurat::FindAllMarkers().
#' @param cl_n Object of integer indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly, 100 by default.
#' @param hvg Object of class numeric or "uns"; Number of highly variable genes to use on top of the marker genes, if "uns" then it is completely unsupervised and uses top 3000 HVG.
#' @return A downsampled Seurat object from the original
#' @export
#' @examples
#'

downsample_se_obj <- function(counts_sc,
                              clust_vr,
                              cluster_markers,
                              cl_n = 10,
                              hvg = 0) {

  # Check variables
  # if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.data.frame(cluster_markers)) stop("ERROR: cluster_markers must be a data frame object returned from Seurat::FindAllMarkers()!")
  if (!is.numeric(cl_n)) stop("ERROR: cl_n must be an object of class integer!")
  if (! (is.numeric(hvg) | hvg == "uns")) stop("ERROR: hvg must be an object of class integer or a string 'uns'!")

  # load required packages
  suppressMessages(require(Seurat))
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

  if (is.null(hvg)) {
    #### Keep marker genes only ####
    keep_genes <- unique(cluster_markers$gene)

  } else if (length(hvg) > 0) {
    #### Union of marker genes and highest variable genes and subset genes ####
    keep_genes <- unique(hvg, cluster_markers$gene)
  }

  #### Get cell IDs to subset by cluster ####
  keep_ids <- lapply(unique(clust_vr), function(ct) {
    ct_sub <- clust_vr[clust_vr == ct]

    # Determine n_sample, if the size of the group is < cl_n use all the group, if not just use cl_n
    n_sample <- dplyr::if_else(length(ct_sub) < cl_n,
                               as.numeric(length(ct_sub)),
                               as.numeric(cl_n))

    # Subset a random selection of that group and get the identifiers
    id_pos <- which(clust_vr == ct)
    id_sel <- sample(x = id_pos,
                     size = n_sample,
                     replace = FALSE)
    return(id_sel)
    }) %>% unlist()


  #### Subset seurat object ####
  counts_sc <- counts_sc[keep_genes, keep_ids]

  return(counts_sc)
}
