# When deployed on the server we intermitently get a bug.
# After specifying a valid disease model and pedigree the user presses the calc risk button
# However risk calculation fails because
# can't write/access Rplots.pdf
# Given the code we shouldn't be plotting to this file.
# see
# https://groups.google.com/forum/#!topic/shiny-discuss/gWuY9OY6RIs
# http://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
# I think this is a problem with library(gplots)
# I have put the below function which uses library(gplots) in a separate source file in order to isolate it.
# Different programs can then choose whether or not to source this code.
# I will ensure server.R doesn't
#
# D Campbell 2.12.15



#### adapted from
# A short tutorial for decent heat maps in R, Sebastian Raschka
# http://sebastianraschka.com/Articles/heatmaps_in_r.html
# don't call library() in a package
suppressPackageStartupMessages(library(gplots))

#' Plots a heatmap representation of a covariance matrix
#'
#' Plots a heatmap representation of a covariance matrix
#' @param mat_data matrix to be plotted
#' @param nofSignifDigits number of significant digits for cell contents
#' @return NULL
#' @export
fnHeatMapCovMatrix <- function(mat_data, nofSignifDigits=2)
{
  rnames <-colnames(mat_data)
  # creates a own color palette from red to green
  nofBreaks <- 3*15
  my_palette <- grDevices::cm.colors(nofBreaks-3)
  my_palette <- grDevices::colorRampPalette(c("red", "orange", "yellow", "white", "light green", "light blue","violet"))(n = nofBreaks-3)

  # (optional) defines the color breaks manually for a "skewed" color transition
  myTrans <- 0.4
  col_breaks = c(seq(-1,-myTrans,length=nofBreaks/3),  # for red
                 seq(-myTrans,myTrans,length=nofBreaks/3)[-1],      # for white
                 seq(myTrans,1,length=nofBreaks/3)[-1]              # for dark green
  )

  mat_data_signif <- signif(mat_data,nofSignifDigits)
  suppressWarnings(
    gplots::heatmap.2(mat_data
              ,Rowv=NULL
              ,symm=T
              ,cellnote = mat_data_signif  # same data set for cell labels
              ,main = "Correlation" # heat map title
              ,notecol="black"      # change font color of cell labels to black
              ,density.info="none"  # turns off density plot inside color legend
              ,trace="none"         # turns off trace lines inside the heat map
              ,margins =c(12,9)     # widens margins around plot
              ,col=my_palette       # use on color palette defined earlier
              ,breaks=col_breaks    # enable color transition at specified limits
              ,dendrogram="none"     # only draw a row dendrogram
    )
  )
}
