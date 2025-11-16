
#' Example Diffusion Model Data
#'
#' A dataset containing reaction times and responses from a study on wishful predictions,
#' used as an example for fitting the Diffusion Decision Model (DDM) with \code{\link{driftd}}.
#'
#' @format A data frame with \eqn{N} rows and \eqn{K} columns:
#' \describe{
#'   \item{ID}{Unique participant identifier.}
#'   \item{val}{Experimental condition indicating trial type (critical-outcome type: "pos" vs "neg").}
#'   \item{prob}{Experimental condition indicating probability of critical outcome.}
#'   \item{R}{The binary response (1 = critical outcome; 0 = neutral outcome).}
#'   \item{RT}{Reaction time in seconds (s).}
#' }
#'
#' @details
#' This dataset is derived from the study cited below. Note that at the time of
#' this package release, the citation refers to a preprint that has not yet
#' undergone peer review. Participants in this study were presented with bi-colored
#' stimuli (grids composed of squares) that varied in color discriminability. On each trial,
#' one square in the grid was randomly selected by the computer and participants had to predict
#' the color of the randomly chosen square. Each trial had a critical outcome color (associated with a gain or a loss)
#' and a neutral outcome color (associated with no points).
#'
#' @source
#' Strueder, J. D., Park, I., Mordkoff, J. T., Pleskac, T. J., & Windschitl, P. D. (2025).
#' The Cognitive Mechanisms Behind Wishful Predictions: A Diffusion Model Decomposition. PsyArXiv.
#' \href{https://doi.org/10.31234/osf.io/bt7m5_v1}{DOI:10.31234/osf.io/bt7m5_v1}
#'
#' @usage data(example_data)
"example_data"
