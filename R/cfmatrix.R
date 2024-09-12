pow <- `^`


##----------------------------------------
## Always yields a dense base::matrix.
##----------------------------------------

#' @title Basic matrix ops
#' @name basic_matrix_ops
#' @description These will always yield a (dense) `base::matrix`.
#' For `as_row()` and `as_col()`, that'll be a one-row or one-col matrix, but still a matrix.
NULL

#' @rdname basic_matrix_ops
#' @export
as_mat <- function(x) {
  base::as.matrix(x)
}

#' @rdname basic_matrix_ops
#' @export
as_row <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(nrow(x) == 1)) return(x)
  base::matrix(x, nrow = 1)
}

#' @rdname basic_matrix_ops
#' @export
as_col <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(ncol(x) == 1)) return(x)
  base::matrix(x, ncol = 1)
}


##----------------------------------------
## Lp-norms
##----------------------------------------

#' @title Lp-norms
#' @name lp_norms
NULL

#' @rdname lp_norms
#' @export
vec_lpnorm <- function(vec, p = 1) {
  vec |>
    abs() |>
    pow(p) |>
    sum() |>
    pow(1/p)
}

#' @rdname lp_norms
#' @export
mat_margin_lpnorms <- function(mat, p = 1, margin = 2) {
  f_sums <-
    if (margin == 1) {
      rowSums
    } else if (margin == 2) {
      colSums
    } else {
      stop("not a valid margin")
    }
  mat |> abs() |> pow(p) |> f_sums() |> pow(1/p)
}

#' @rdname lp_norms
#' @export
mat_col_lpnorms <- function(mat, p = 1) {
  mat_margin_lpnorms(mat, p, 2)
}

## Earlier version of mat_row_lpnorms before deciding to just use the col method on t(mat).
## mat_row_lpnorms <- function(mat, p = 1) {
##   mat_margin_lpnorms(mat, p, 1)
## }

#' @rdname lp_norms
#' @export
mat_row_lpnorms <- function(mat, p = 1) {
  mat_col_lpnorms(t(mat), p = 1)
}


##----------------------------------------
## Lp-normalizers.
##----------------------------------------

#' @title Lp-normalizers
#' @name lp_normalizers
NULL

#' @rdname lp_normalizers
#' @export
vec_lpnormalize <- function(vec, p = 1) {
  vec / vec_lpnorm(vec, p)
}

#' @rdname lp_normalizers
#' @export
mat_row_lpnormalize <- function(mat, p = 1) {
  Matrix::Diagonal(x = 1 / mat_row_lpnorms(mat, p)) %*% mat
}

#' @rdname lp_normalizers
#' @export
mat_col_lpnormalize <- function(mat, p = 1) {
  mat %*% Matrix::Diagonal(x = 1 / mat_col_lpnorms(mat, p))
}

#' @rdname lp_normalizers
#' @export
mat_margin_lpnormalize <- function(mat, p = 1, margin = 2) {
  if (margin == 1) {
    return(mat_row_lpnormalize(mat, p))
  }
  if (margin == 2) {
    return(mat_col_lpnormalize(mat, p))
  }
  stop("margin must be 1 or 2")
}


##----------------------------------------
## Vector binary product operations.
##----------------------------------------

#' @title Vector binary product ops
#' @name vector_binary_product_ops
NULL

## x & y because the inner product is symmetric.
#' @rdname vector_binary_product_ops
#' @export
vec_inner_prod <- function(vec_x, vec_y = vec_x) {
  crossprod(vec_x, vec_y)
}

## l and r because the outer product is not symmetric.
#' @rdname vector_binary_product_ops
#' @export
vec_outer_prod <- function(vec_l, vec_r = vec_l) {
  tcrossprod(vec_l, vec_r)
}


##----------------------------------------
## Simple cosine similarities.
##----------------------------------------

#' @title Cosine similarities
#' @name cosine_sim

#' @rdname cosine_sim
#' @export
vec_cos_sim <- function(vec_x, vec_y) {
  numer <- vec_inner_prod(vec_x, vec_y)
  denom <- vec_lpnorm(vec_x, 2) * vec_lpnorm(vec_y, 2)
  numer / denom
}

#' @rdname cosine_sim
#' @export
mat_row_cos_sim <- function(mat_l, mat_r = mat_l) {
  numer <- tcrossprod(mat_l, mat_r)
  mat_l_row_lpnorms <- mat_row_lpnorms(mat_l, 2)
  mat_r_row_lpnorms <- mat_row_lpnorms(mat_r, 2)
  denom <- vec_outer_prod(mat_l_row_lpnorms, mat_r_row_lpnorms)
  numer / denom
}

#' @rdname cosine_sim
#' @export
mat_col_cos_sim <- function(mat_l, mat_r = mat_l) {
  numer <- crossprod(mat_l, mat_r)
  mat_l_col_lpnorms <- mat_col_lpnorms(mat_l, 2)
  mat_r_col_lpnorms <- mat_col_lpnorms(mat_r, 2)
  denom <- vec_outer_prod(mat_l_col_lpnorms, mat_r_col_lpnorms)
  numer / denom
}

#' @rdname cosine_sim
#' @export
mat_margin_cos_sim <- function(mat_l, mat_r = mat_l, margin = 2) {
  if (margin == 1) return(mat_row_cos_sim(mat_l, mat_r))
  if (margin == 2) return(mat_col_cos_sim(mat_l, mat_r))
  stop("margin must be either 1 or 2")
}
