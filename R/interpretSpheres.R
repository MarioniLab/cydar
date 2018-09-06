#' @export
#' @importFrom shiny plotOutput reactiveValues fluidRow column textInput actionButton hr HTML htmlOutput selectInput
#' pageWithSidebar headerPanel sidebarPanel h4 tableOutput sliderInput mainPanel renderPlot reactive renderTable
#' observeEvent updateTextInput nearPoints updateSelectInput renderUI shinyApp runApp br observe stopApp
#' @importFrom flowCore markernames
#' @importFrom utils tail
#' @importFrom grDevices col2rgb
interpretSpheres <- function(x, markers=NULL, labels=NULL, select=NULL, 
    metrics=NULL, num.per.row=6, plot.height=100, xlim=NULL, p=0.01, 
    red.coords=NULL, red.highlight=NULL, red.plot.height=500, 
    add.plot=NULL, add.plot.height=500, run=TRUE, ...)   
# Creates a Shiny app to assist interpretation of the hyperspheres.
#
# written by Aaron Lun
# created 1 November 2016    
{
    nrows <- ceiling(length(markernames(x, mode="all"))/num.per.row)
    plot.height <- plot.height*nrows
    collim <- intensityRanges(x, p=p)
    all.dens <- .prepareDensity(x, ...)

    # Checking marker order.
    all.markers <- markernames(x, mode="all")
    if (is.null(markers)) { 
        markers <- all.markers 
    } else if (!all(markers %in% all.markers)) {
        stop("requested markers not present in available markers in 'x'")
    }

    # Checking hypersphere order.
    if (!is.null(select)) {
       select <- .subsetToIndex(select, nrow(x), "select") 
    } else {
        select <- seq_len(nrow(x))
    }
    if (length(select)==0L){ 
        stop("empty 'select' provided")
    }

    # Should we add a nav plot?
    if (!is.null(red.coords)) { 
        add.nav <- TRUE
        stopifnot(identical(nrow(red.coords), nrow(x)))
        stopifnot(identical(ncol(red.coords), 2L))
        if (!is.null(red.highlight)) {
            red.highlight <- .subsetToIndex(red.highlight, nrow(x), "red.highlight")
        }
        red.coords <- data.frame(x=red.coords[,1], y=red.coords[,2])
    } else {
        add.nav <- FALSE
    }

    # Using existing labels if provided.
    if (is.null(labels)) {
        labels <- character(nrow(x))
    } else {
        stopifnot(identical(length(labels), nrow(x)))
    } 

    # Checking metrics if provided.
    if (!is.null(metrics)) {
        stopifnot(identical(nrow(metrics), nrow(x)))
    }

    # Setting up the internal plotting and storage.
    collected <- reactiveValues(current=select[1], labels=labels, history=rep(NA_integer_, 5))
   
    # Main panel arguments.
    main.args <- list(plotOutput("histograms", height = plot.height), hr())
    main.args <- append(main.args, list(
        fluidRow(
            column(
                textInput("label", "Label sphere:"),
                width=4
            ),
            column(
                textInput("gotonum", "Go to sphere:", select[1]),
                actionButton("previous", "Previous"),
                actionButton("continue", "Next"),
                width=4
            ),
            column( 
                actionButton("finish", "Save to R", style="width:120%; line-height: 50px; background-color: #FA8072"),
                width=2
            )
        )
    ))

    if (add.nav) { 
        main.args <- append(main.args, list(hr(), 
            plotOutput("navplot", height = red.plot.height, click = "nav_click"),
            hr(),
            fluidRow(
                column(
                    actionButton("updatelabels", "Update labels"), br(), br(),
                    selectInput("labeltouse", "Labels to plot:", choices="", multiple=TRUE, selectize=FALSE, size=10),
                    HTML("<b>Legend:</b>"),
                    htmlOutput("lablegend"),
                    width=3
                ), 
                column(
                    plotOutput("labplot", height = red.plot.height),
                    width=8
                )
            )
        ))
    }    
    if (!is.null(add.plot)) {
        main.args <- append(main.args, list(hr(), plotOutput("addplots", height = add.plot.height)))
    }

    # Generating the page layout.
    ui <- pageWithSidebar(
        headerPanel("Interpreting hypersphere coordinates"),
        sidebarPanel(
            h4("Metrics"),
            tableOutput("metrics"),
            hr(size=30),
            h4("History:"),
            tableOutput("history"), 
            hr(size=30),
            h4("Closest labelled:"),
            tableOutput("closest"), 
            hr(size=30),
            sliderInput("intbar", "Intensity interval for current hypersphere (%):",
                min = 0, max = 100, value = 95),
            textInput("extraplot", "Add more hyperspheres:", value=""),
            actionButton("clearplot", "Clear")
        ),
        do.call(mainPanel, main.args)
    )

    # Setting up the server actions.
    server <- function(input, output, session) {
        reactiveHistPlot <- makeHistograms(input, mfrow=c(nrows, num.per.row), density.data=all.dens, 
            xlim=xlim, collim=collim, collected=collected, markers=markers, x=x)
        output$histograms <- renderPlot({ reactiveHistPlot() })

        if (add.nav) {
            reactiveNavPlot <- makeNavPlot(input, red.highlight=red.highlight, red.coords=red.coords, collected=collected)
            output$navplot <- renderPlot({ reactiveNavPlot() })
        }

        if (!is.null(add.plot)) { 
            reactiveAddPlot <- reactive({ add.plot(as.integer(input$gotonum), x) })
            output$addplots <- renderPlot({ reactiveAddPlot() })
        }

        reactiveMetricTable <- makeMetricTable(input, metrics=metrics, collected=collected)
        output$metrics <- renderTable({
            reactiveMetricTable()
        }, colnames=FALSE, sanitize.text.function=identity)

        reactiveHistoryTable <- makeHistoryTable(input, collected=collected)
        output$history <- renderTable({ reactiveHistoryTable() })

        reactiveClosestTable <- makeClosestTable(input, x=x, collected=collected)
        output$closest <- renderTable({ reactiveClosestTable() })

        # Setting up events to observe.
        observeEvent(input$gotonum, {
            current <- as.integer(input$gotonum)
            if (is.na(current) || current < 1L || current > nrow(x)) {
                warning("specified index is not a number within range")
            } else {
                collected$current <- current
                updateTextInput(session, "label", value=collected$labels[current])
            }
        }) 

        observeEvent(input$continue, {
            current <- collected$current
            attempt <- select[which(select > current)[1]]
            if (is.na(attempt)) {
                warning("no index in 'select' is larger than the current index")
            } else {
                collected$current <- attempt
                updateTextInput(session, "gotonum", value=attempt)
                updateTextInput(session, "label", value=collected$labels[current])
            }
        })                  

        observeEvent(input$previous, {                     
            current <- collected$current
            attempt <- select[tail(which(select < current), 1)]
            if (length(attempt)==0) {
                warning("no index in 'select' is smaller than the current index")
            } else { 
                collected$current <- attempt
                updateTextInput(session, "gotonum", value=attempt)
                updateTextInput(session, "label", value=collected$labels[current])
            }
        })

        observeEvent(input$label, {
            collected$labels[collected$current] <- input$label 
        }, priority=1)


        if (add.nav) {
            observe({
                all.near <- nearPoints(red.coords, input$nav_click, xvar = "x", yvar = "y", threshold=Inf, maxpoints=1)
                if (nrow(all.near)) {
                    current <- as.integer(rownames(all.near)[1])
                    collected$current <- current 
                    updateTextInput(session, "gotonum", value=current)
                    updateTextInput(session, "label", value=collected$labels[current])
                }
            })

            observeEvent(input$updatelabels, {
                available <- setdiff(unique(collected$labels), "")
                new.select <- intersect(input$labeltouse, available)
                updateSelectInput(session, "labeltouse", choices=available, selected=new.select)
                collected$more.labels <- labelSpheres(x, collected$labels)
            })

            reactiveLabPlot <- makeLabPlot(input, red.coords, collected)
            output$labplot <- renderPlot({ reactiveLabPlot() }) 

            output$lablegend <- renderUI({
                cols <- .obtainColours(input$labeltouse)
                as.rgb <- col2rgb(cols)
                leg.box <- sprintf('<div style="display:inline-block; background:rgb(%i, %i, %i); margin-right:5px; width:10px; height:10px"></div>', 
                                   as.rgb[1,], as.rgb[2,], as.rgb[3,]) 
                leg.names <- sprintf('<span>%s</span><br />', names(cols))
                stuff <- paste(paste0(leg.box, leg.names), collapse="")
                HTML(stuff)
            })
        }

        observeEvent(input$clearplot, {
            updateTextInput(session, "extraplot", value="")
        })

        observeEvent(input$finish, {
            stopApp(collected$labels)
        })
    }
    
    app <- shinyApp(ui, server)
    if (run) {
        return(runApp(app))
    } else {
        return(app)
    }
}

#######################################################
# Setting up functions to return histograms.

#' @importFrom graphics par
makeHistograms <- function(input, mfrow, ..., x, collected) {
    reactive({
        extras <- as.integer(unlist(strsplit(input$extraplot, " ") ))    
        extras <- extras[!is.na(extras)]
        invalid <- !is.finite(extras) | extras < 1L | extras > nrow(x)
        if (any(invalid)) {    
            warning("indices out of range of number of hyperspheres")
            extras <- extras[!invalid]
        }
        par(mfrow=mfrow, mar=c(2.1, 1.1, 2.1, 1.1))
        .generateDensity(..., x=x, current=collected$current, extras=extras, interval=input$intbar)
    })
}

#' @importFrom shiny reactive
#' @importFrom graphics plot par legend points
makeNavPlot <- function(input, red.highlight, red.coords, collected) { 
    reactive({
        par(mar=c(5.1, 4.1, 1.1, 10.1))
        col <- rep("grey", nrow(red.coords))
        if (!is.null(red.highlight)) {
            col[red.highlight] <- "orange"
        }
        plot(red.coords$x, red.coords$y, xlab="Dimension 1", ylab="Dimension 2", col=col, pch=16, cex.lab=1.4)

        has.label <- collected$labels!=""
        points(red.coords$x[has.label], red.coords$y[has.label], col="black", pch=16, cex=1.5)
        current <- collected$current
        points(red.coords$x[current], red.coords$y[current], col="red", pch=16, cex=2)

        uout <- par()$usr
        par(xpd=TRUE)
        legend(uout[2] + (uout[2]-uout[1])*0.01, uout[4], col=c("red", "black", "orange", "grey"), 
            pch=16, legend=c("Current", "Labelled", "Highlighted", "Other"))
    })
}

#' @importFrom shiny reactive
makeMetricTable <- function(input, metrics, collected) {
    reactive({ 
        current <- collected$current
        labels <- collected$labels
        out.metrics <- c("Number", "Label")
        out.vals <- c(current, labels[current])

        if (!is.null(metrics)) { 
            extra.metrics <- colnames(metrics)
            extra.vals <- vector("list", length(extra.metrics)) 
            for (met in seq_along(extra.metrics)) {
                my.val <- metrics[collected$current,met]
                if (is.double(my.val)) {
                    my.val <- format(my.val, digits=5)
                } else {
                    my.val <- as.character(my.val)
                }
                extra.vals[[met]] <- my.val
            }
            out.metrics <- c(out.metrics, extra.metrics)
            out.vals <- c(out.vals, unlist(extra.vals))
        }
        data.frame(paste0("<b>", out.metrics, "</b>"), out.vals)
    })
}

history_len <- 5

#' @importFrom shiny reactive
makeHistoryTable <- function(input, collected) { 
    reactive({ 
        current <- collected$current
        if (is.na(collected$history[1]) || current!=collected$history[1]) { 
            collected$history <- c(current, collected$history)[seq_len(history_len)]
        }
        data.frame(Number=as.character(collected$history),
            Label=collected$labels[collected$history])
    })
}

#' @importFrom shiny reactive
makeClosestTable <- function(input, x, collected) {
    reactive({
        labels <- collected$labels
        current <- collected$current
        is.anno <- setdiff(which(labels!=""), current)

        if (length(is.anno)) { 
            intvals <- .raw_intensities(x)
            all.dist <- sqrt(colSums((t(intvals[is.anno,,drop=FALSE]) - intvals[current,])^2))
        } else {
            all.dist <- numeric(0)
        }

        o <- order(all.dist)[seq_len(history_len)]
        closest <- is.anno[o]
        data.frame(Distance=all.dist[o], Number=as.character(closest), Label=labels[closest])
    })
}

#' @importFrom shiny reactive
#' @importFrom graphics plot plot.new points
makeLabPlot <- function(input, red.coords, collected) {
    reactive({
        if (length(input$labeltouse)==0 || identical(input$labeltouse, "")) { 
            plot.new()
        } else {
            plot(red.coords$x, red.coords$y, xlab="Dimension 1", ylab="Dimension 2", col="grey80", pch=16, cex.lab=1.4)
            all.colours <- .obtainColours(input$labeltouse)
            for (lab in input$labeltouse) { 
                autolab <- which(collected$more.labels==lab)
                points(red.coords$x[autolab], red.coords$y[autolab], col=all.colours[lab], pch=16, cex=1.1)
            }
        }
    })
}

#######################################################
# Utility functions.

.subsetToIndex <- function(subset, N, argname) {
    if (is.logical(subset)) { 
        stopifnot(identical(length(subset), N))
        subset <- which(subset)
    } else if (is.numeric(subset)) { 
        subset <- sort(as.integer(subset))
        stopifnot(all(subset >= 1L & subset <= N))
    } else {
        stop(sprintf("unknown type for '%s'", argname))
    }
    return(subset)
}

#' @importFrom stats density 
#' @importFrom flowCore markernames
.prepareDensity <- function(x, ...) 
# Computing the density once to speed up plotting.
# Done separately for both the used and unused markers.
{ 
    used.markers <- markernames(x)
    used.int <- .raw_cellIntensities(x)
    used.collected <- vector("list", length(used.markers))
    for (m in seq_along(used.markers)) {
        cur.intensities <- used.int[m,]
        used.collected[[m]] <- density(cur.intensities, ...)
    }
    names(used.collected) <- used.markers

    unused.markers <- markernames(x, mode="unused")
    unused.int <- .raw_unusedIntensities(x)
    unused.collected <- vector("list", length(unused.markers))
    for (m in seq_along(unused.markers)) {
        cur.intensities <- unused.int[m,]
        unused.collected[[m]] <- density(cur.intensities, ...)
    }
    names(unused.collected) <- unused.markers

    c(used.collected, unused.collected)
}

#' @importFrom viridis viridis
#' @importFrom stats approx quantile median
#' @importFrom graphics plot polygon lines points par text segments
.generateDensity <- function(density.data, markers, collim, xlim, current, extras, interval, x, ...) 
# Plotting the densities and adding the point corresponding to the current coordinates.
{
    all.cols <- viridis(100)
    coords <- intensities(x)
    used.markers <- markernames(x)

    unused.markers <- markernames(x, mode="unused")
    used.int <- .raw_cellIntensities(x)
    unused.int <- .raw_unusedIntensities(x)
    sid <- .raw_sample_id(x)
    cur.assign <- .raw_cellAssignments(x)[[current]]

    for (m in markers) {
        was.used <- m %in% used.markers

        # Defining the colour and position for (un)used markers.
        # For unused markers, we take the mean of sample-wise medians, reflecting the linear modelling of the medians.
        if (was.used) { 
            curpos <- coords[current, m]
            curlim <- collim[,m]
            if (diff(curlim) > 0) {
                curdex <- round(approx(collim[,m], c(1, 100), xout=curpos, rule=2)$y)
            } else {
                curdex <- 1L
            }
            curcol <- all.cols[curdex]
            cur.cell.int <- used.int[match(m, used.markers), cur.assign, drop=FALSE]
        } else {
            cur.cell.int <- unused.int[match(m, unused.markers), cur.assign, drop=FALSE]
            by.sample <- split(cur.cell.int, sid[cur.assign])
            curpos <- mean(vapply(by.sample, FUN=median, FUN.VALUE=0))
            curcol <- "grey80"
        }

        if (is.null(xlim)) { 
            xlim2 <- collim[,m]
        } else {
            xlim2 <- xlim
        }
    
        curdens <- density.data[[m]]
        plot(0, 0, type="n", xlab="", ylab="", yaxt="n", bty="n", main=m, xlim=xlim2, ylim=c(0, max(curdens$y)), ...)
        my.x <- c(curdens$x[1]-10, curdens$x, curdens$x[length(curdens$x)]+10)
        my.y <- c(0, curdens$y, 0)
        polygon(my.x, my.y, col=curcol, border=NA)
        lines(my.x, my.y)
   
        # Adding the point of the median intensity.
        curpos <- pmin(xlim2[2], pmax(xlim2[1], curpos))
        cury <- approx(my.x, my.y, curpos, rule=2)$y
        par(xpd=TRUE)
        points(curpos, cury, pch=16, col="red", cex=1.5)
        
        # Adding intervals specifying the spread of cells within this hypersphere.
        if (interval > 0) {
            half.left <- (1-interval/100)/2
            prob.interval <- c(half.left, interval/100 + half.left)
            q.int <- quantile(cur.cell.int[m,], prob.interval)
            segments(q.int[1], cury, q.int[2], col="red")
        }

        # Adding previous points for comparison.
        if (length(extras)) {
            text(curpos, cury, pos=3, current, col="red")
            for (ex in extras) {
                expos <- coords[ex,m]
                expos <- pmin(xlim2[2], pmax(xlim2[1], expos))
                exy <- approx(my.x, my.y, expos, rule=2)$y
                points(expos, exy, pch=16, col="red", cex=1.5)
                text(expos, exy, pos=3, ex, col="red")
            }
        }
        par(xpd=FALSE)
    }

    return(invisible())
}

#' @importFrom grDevices rainbow
.obtainColours <- function(labels) 
# Generate colours.
{
    all.colours <- rainbow(length(labels))
    names(all.colours) <- labels
    all.colours
}

