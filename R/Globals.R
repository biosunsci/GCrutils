#' @name HEIGHT
#' @title use to control the default output fig height
#' @description This variable defines the default height of the output figs, commonly used by
#'   function ysave, ydumpto. like cairo_pdf(..., height = HEIGHT, width = WIDTH)
#' @export
HEIGHT <<- 8

#' @name WIDTH
#' @title use to control the default output fig width
#' @description This variable defines the default width of the output figs, commonly used by
#'   function ysave, ydumpto. like cairo_pdf(..., height = HEIGHT, width = WIDTH)
#' @export
WIDTH <<- 10


#' @name OUTPUTROOT
#' @title Output root directory
#' @description This variable defines the root directory for output files.
#' @export
OUTPUTROOT <<- 'report'


#' @name INPUTROOT
#' @title Input root directory
#' @description This variable defines the root directory for input files.
#' @export
INPUTROOT <<- 'export'


warning('Global Vars set: INPUTROOT OUTPUTROOT WIDTH HEIGHT')
