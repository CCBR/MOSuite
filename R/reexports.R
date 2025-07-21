#' dplyr pipe
#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' magrittr assignment pipe
#' @importFrom magrittr %<>%
#' @export
magrittr::`%<>%`

#' walrus operator
#' @importFrom rlang :=
#' @export
rlang::`:=`

#' bang-bang
#' @importFrom rlang !!
#' @export
rlang::`!!`

#' rlang data pronoun
#' @importFrom rlang .data
#' @export
rlang::.data

# Suppress R CMD check note 'All declared Imports should be used'.
# These are used in S7 methods.
#' @importFrom dendextend rotate
#' @importFrom matrixStats rowVars
NULL
