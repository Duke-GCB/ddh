## A funciton to generate color palettes

generate_colors <- function(hex) {
  #' @param hex Hex code for main color as string
  #' @return A set of three colors
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .5), .3),
    hex,
    colorspace::darken(hex, .5, space = "HLS")
  )
  
  return(color_set)
    
  #assign("color_pal", color_pal, envir = .GlobalEnv)
}


#' @param color_set Number of unqiue colors
#' @param n Number of unqiue colors
#' @return A function to generate color palette
color_pal <- function(n) {
  pal <- grDevices::colorRampPalette(color_set)
  return(pal)
}













## A function to generate color palettes

generate_colors <- function(hex) {
  #' @param hex Hex code for main color as string
  #' @return A color palette function with argument `n` representing number of colors to be returned
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .5), .3),
    hex,
    colorspace::darken(hex, .5, space = "HLS")
  )
  color_pal <- grDevices::colorRampPalette(color_set)
  
  assign("color_pal", color_pal, envir = .GlobalEnv)
  return(color_set)
}






