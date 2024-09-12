mat_chunk_max_length <- function(mat, margin, chunk_max_bytes, cell_bytes = 8) {
  fixed_dim_length <- if (margin == 1) ncol(mat) else nrow(mat)
  fixed_vec_bytes <- fixed_dim_length * cell_bytes
  floor(chunk_max_bytes / fixed_vec_bytes)
}


#' @title Matrix chunking indices.
#' @description Chunking a matrix means preserving at least one of the two dimensions.
#' In the pathological extreme case the smallest chunk is either a single col or a single row.
#' No subset of _both_ rows and cols are supported.
#' @value A list of indices to use for chunking corresponding to the specified `margin`.
#' @export
mat_chunk_iixs <- function(mat, margin, chunk_count = NULL, chunk_max_length = NULL, chunk_max_bytes = NULL, cell_bytes = 8) {
  tryCatch({
    if (!is.null(chunk_count) || !is.null(chunk_max_length)) {
      return(chunk_int(dim(mat)[margin], chunk_count = chunk_count, chunk_max_length = chunk_max_length, method = "seq"))
    }
    if (!is.null(chunk_max_bytes)) {
      return(chunk_int(dim(mat)[margin], chunk_count = NULL, chunk_max_length = mat_chunk_max_length(mat, margin, chunk_max_bytes, cell_bytes), method = "seq"))
    }
    stop("one of (chunk_count, chunk_max_length, chunk_max_bytes) must be non-NULL")
  }, error = function(e) {
    stop(sprintf("failed to chunk matrix. maybe chunk_max_bytes is too small for any matrix slice? errmsg = %s", conditionMessage(e)))
  })
}


#' @title Margin norms computed chunk-wise
#' @description Computes margin norms for a matrix but in chunks.
#' Notably this allows for parallel computation, which is supported out-of-the-box using `parallel::mclapply` and the `mc.cores` parameter.
mat_margin_lpnorms_chunked <- function(mat, margin = 2, p = 1,
                                       chunk_count = NULL,
                                       chunk_max_length = NULL,
                                       chunk_max_bytes = NULL,
                                       cell_bytes = 8,
                                       mc.cores = 1) {
  chunk_iixs <-
    mat_chunk_iixs(mat, margin, chunk_count = chunk_count, chunk_max_length = chunk_max_length, chunk_max_bytes = chunk_max_bytes, cell_bytes = cell_bytes)
  chunked_lpnorms <-
    parallel::mclapply(chunk_iixs, \(chunk_iix) {
      if (margin == 1) {
        return(mat_margin_lpnorms(mat[chunk_iix,], margin, p))
      }
      if (margin == 2) {
        return(mat_margin_lpnorms(mat[,chunk_iix], margin, p))
      }
      stop("margin must be 1 or 2")
    }, mc.cores = mc.cores)
  chunked_lpnorms
}


