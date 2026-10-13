
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

#' Speed-Accuracy Lexical Decision Data (Wagenmakers et al., 2008)
#'
#' Responses and response times of 17 participants in a lexical decision task (deciding
#' whether a letter string is a word) in which instructions emphasized either speed or
#' accuracy. This is the main example dataset of driftd. It is a prepared version of
#' \code{\link[rtdists]{speed_acc}} from the rtdists package.
#'
#' @format A data frame with 31,351 rows (trials) and 5 columns:
#' \describe{
#'   \item{id}{Participant identifier (factor, 17 levels).}
#'   \item{condition}{Instruction in the block: \code{"accuracy"} or \code{"speed"}.}
#'   \item{stimulus}{Stimulus type: a \code{"high"}, \code{"low"} or \code{"very_low"}
#'     frequency word, or a \code{"nonword"}.}
#'   \item{response}{Response: 1 = "word", 0 = "nonword".}
#'   \item{rt}{Response time in seconds.}
#' }
#'
#' @details
#' Prepared from \code{rtdists::speed_acc} by removing the trials flagged in its
#' \code{censor} column (uninterpretable responses, RTs below 180 ms or above 3 s) and
#' combining the three nonword types into one category. The script is in the
#' \code{data-raw} folder of the package's source repository. Participant 2 contributed
#' fewer trials than the others.
#'
#' The data were collected by Wagenmakers et al. (2008, Experiment 1) and, with two
#' additional participants, analyzed by Heathcote and Love (2012). They are distributed with
#' the rtdists package (Singmann et al.) under the GPL (>= 3). Please cite the original
#' study when using these data.
#'
#' @source
#' Wagenmakers, E.-J., Ratcliff, R., Gomez, P., & McKoon, G. (2008). A diffusion model
#' account of criterion shifts in the lexical decision task. \emph{Journal of Memory and
#' Language, 58}(1), 140-159. \doi{10.1016/j.jml.2007.04.006}
#'
#' Singmann, H., Brown, S., Gretton, M., & Heathcote, A. rtdists: Response Time
#' Distributions. R package. \url{https://CRAN.R-project.org/package=rtdists}
#'
#' @references
#' Heathcote, A., & Love, J. (2012). Linear deterministic accumulator models of simple
#' choice. \emph{Frontiers in Psychology, 3}, 292. \doi{10.3389/fpsyg.2012.00292}
#'
#' @examples
#' data(speed_acc)
#' head(speed_acc)
#' # proportion of "word" responses per stimulus type and instruction
#' with(speed_acc, round(tapply(response, list(stimulus, condition), mean), 2))
#'
#' @usage data(speed_acc)
"speed_acc"
