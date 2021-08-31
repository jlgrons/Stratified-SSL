# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param S Vector with stratum id.
#' @param num_knots Number of knots.
#' @param basis_type Type of basis to use.
#' @export
#' @return Matrix containing basis.
#'

AlternativeBasis <- function(X, S, num_knots, basis_type = 'interact'){

  if (basis.type == 'interact'){

    basis.X <- InteractionBasis(X)
    basis.S <- OneHotEncoding(S)

    basis <- cbind(basis.X, basis.S)

  }

  # Will change name later.
  if (basis.type == 'IC1'){

    basis.X <- InteractionBasis(X)
    basis.S <- TwoWayInteractionBasis(basis.X, S)

    basis <- cbind(basis.X, basis.S)
  }

  # Will change name later.
  if (basis.type == 'II1'){

    natural_spline_basis <- NaturalSplineBasis(X)
    interaction_basis <- TwoWayInteractionBasis(X, S)

    basis <- cbind(natural_spline_basis, S)
  }


  return(basis)
}
