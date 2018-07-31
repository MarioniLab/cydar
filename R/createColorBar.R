#' @export
#' @importFrom utils head tail
#' @importFrom graphics rect text
createColorBar <- function(colors, top.text=NULL, bottom.text=NULL, lower=-0.5, upper=0.5, x.pos=0, width=1, cex=1.5) 
# Creates a color bar.
{
    start.loc <- seq(lower, upper, length.out=length(colors)+1)
    interval <- diff(start.loc)[1]
    start.loc <- head(start.loc, -1L)
    half.width <- width/2
    rect(x.pos - half.width, start.loc, x.pos + half.width, start.loc + interval, col=colors, border=colors)
    
    if (is.null(top.text)) {
        top.text <- tail(names(colors), 1L)
    } 
    if (!is.null(top.text)) {
        text(x.pos, upper, pos=3, top.text, cex=cex)
    }

    if (is.null(bottom.text)) {
        bottom.text <- head(names(colors), 1L)
    } 
    if (!is.null(bottom.text)) {
        text(x.pos, lower, pos=1, bottom.text, cex=cex)
    }
    return(invisible(NULL))
}
