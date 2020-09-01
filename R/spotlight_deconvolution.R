#' This functions takes in a Seurat object with the scRNAseq data and the count matrix of the spatial transcriptomics data and returns the deconvoluted spots.
#'
#' @param counts_sc Object of class Matrix of shape GENESxCELL, se_obj\@assays$RNA@counts.
#' @param counts_spatial Object of class Matrix of shape GENESxSPOT, se_obj\@assays$Spatial@counts.
#' @param clust_vr Object of class list/vector containing the cell identity/state of each cell (column) in counts_sc.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param cl_n Object of integer indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly, 100 by default.
#' @param hvg Object of class list/vector containing the HVG gene names to use in the analysis.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), raw (no transformation applied). By default CPM.
#' @param method Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.09, by default, means that we will accept those cell types whose weight coefficient is at least 0.09 of the total.
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples

spotlight_deconvolution <- function(counts_sc,
                                    counts_spatial,
                                    clust_vr,
                                    cluster_markers,
                                    cl_n = 100,
                                    hvg = NULL,
                                    ntop = NULL,
                                    transf = "uv",
                                    method = "nsNMF",
                                    min_cont = 0.09) {

  # Downsample scRNAseq to select gene set and number of cells to train the model
  sc_down_ls <- downsample_se_obj(counts_sc = counts_sc,
                                  clust_vr = clust_vr,
                                  cluster_markers = cluster_markers,
                                  cl_n = cl_n,
                                  hvg = hvg)

  # Extract downsampled counts_sc
  counts_sc_down <- sc_down_ls[[1]]

  # Update cluster_vr to downsampled dataset
  clust_vr_mod <- sc_down_ls[[2]]

  # Train the NMF model
  nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers,
                          counts_sc = counts_sc_down,
                          counts_sp = counts_spatial,
                          ntop = ntop,
                          transf = transf,
                          clust_vr = clust_vr_mod,
                          method = method)

  # Get cell type specific topif profiles
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])

  # Perform deconvolution of the capture location mixtures
  decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmf_mod_ls[[1]],
                                          mixture_transcriptome = counts_spatial,
                                          transf = transf,
                                          reference_profiles = ct_topic_profiles,
                                          min_cont = min_cont)

  return(list(nmf_mod_ls, decon_mtrx))
}
