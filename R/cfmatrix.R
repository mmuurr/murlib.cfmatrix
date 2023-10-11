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
mat_dim_lpnorms <- function(mat, dim = 2, p = 1) {
  f_sums <-
    if (dim == 1) {
      rowSums
    } else if (dim == 2) {
      colSums
    } else {
      stop("not a valid dim")
    }
  mat |> abs() |> pow(p) |> f_sums() |> pow(1/p)
}

#' @export
mat_col_lpnorms <- function(mat, p = 1) {
  mat_dim_lpnorms(mat, dim = 2, p = p)
}

## Earlier version of mat_row_lpnorms before deciding to just use the col method on t(mat).
## mat_row_lpnorms <- function(mat, p = 1) {
##   mat_dim_lpnorms(mat, dim = 1, p = p)
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

