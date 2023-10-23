pow <- `^`


##----------------------------------------
## Always yields a dense base::matrix.
##----------------------------------------

#' @export
as_mat <- function(x) {
  base::as.matrix(x)
}

#' @export
as_row <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(nrow(x) == 1)) return(x)
  base::matrix(x, nrow = 1)
}

#' @export
as_col <- function(x) {
  if (isTRUE(is.matrix(x)) && isTRUE(ncol(x) == 1)) return(x)
  base::matrix(x, ncol = 1)
}


##----------------------------------------
## Lp-norms
##----------------------------------------

#' @export
vec_lpnorm <- function(vec, p = 1) {
  vec |>
    abs() |>
    pow(p) |>
    sum() |>
    pow(1/p)
}

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

#' @export
mat_col_lpnorms <- function(mat, p = 1) {
  mat_margin_lpnorms(mat, p, 2)
}

## Earlier version of mat_row_lpnorms before deciding to just use the col method on t(mat).
## mat_row_lpnorms <- function(mat, p = 1) {
##   mat_margin_lpnorms(mat, p, 1)
## }

#' @export
mat_row_lpnorms <- function(mat, p = 1) {
  mat_col_lpnorms(t(mat), p = 1)
}


##----------------------------------------
## Lp-normalizers.
##----------------------------------------

#' @export
vec_lpnormalize <- function(vec, p = 1) {
  vec / vec_lpnorm(vec, p)
}

#' @export
mat_row_lpnormalize <- function(mat, p = 1) {
  Matrix::Diagonal(x = 1 / mat_row_lpnorms(mat, p)) %*% mat
}

#' @export
mat_col_lpnormalize <- function(mat, p = 1) {
  mat %*% Matrix::Diagonal(x = 1 / mat_col_lpnorms(mat, p))
}

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

## x & y because the inner product is symmetric.
#' @export
vec_inner_prod <- function(vec_x, vec_y = vec_x) {
  crossprod(vec_x, vec_y)
}

## l and r because the outer product is not symmetric.
#' @export
vec_outer_prod <- function(vec_l, vec_r = vec_l) {
  tcrossprod(vec_l, vec_r)
}


##----------------------------------------
## Simple cosine simialrities.
##----------------------------------------

#' @export
vec_cos_sim <- function(vec_x, vec_y) {
  numer <- vec_inner_prod(vec_x, vec_y)
  denom <- vec_lpnorm(vec_x, 2) * vec_lpnorm(vec_y, 2)
  numer / denom
}

#' @export
mat_row_cos_sim <- function(mat_l, mat_r = mat_l) {
  numer <- tcrossprod(mat_l, mat_r)
  mat_l_row_lpnorms <- mat_row_lpnorms(mat_l, 2)
  mat_r_row_lpnorms <- mat_row_lpnorms(mat_r, 2)
  denom <- vec_outer_prod(mat_l_row_lpnorms, mat_r_row_lpnorms)
  numer / denom
}

#' @export
mat_col_cos_sim <- function(mat_l, mat_r = mat_l) {
  numer <- crossprod(mat_l, mat_r)
  mat_l_col_lpnorms <- mat_col_lpnorms(mat_l, 2)
  mat_r_col_lpnorms <- mat_col_lpnorms(mat_r, 2)
  denom <- vec_outer_prod(mat_l_col_lpnorms, mat_r_col_lpnorms)
  numer / denom
}

#' @export
mat_margin_cos_sim <- function(mat_l, mat_r = mat_l, margin = 2) {
  if (margin == 1) return(mat_row_cos_sim(mat_l, mat_r))
  if (margin == 2) return(mat_col_cos_sim(mat_l, mat_r))
  stop("margin must be either 1 or 2")
}
