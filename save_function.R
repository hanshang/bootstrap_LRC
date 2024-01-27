savefig = function (filename, height=10, width = (1 + sqrt(5))/2*height, type=c("eps","pdf","jpg","png"), pointsize = 10, family = "Helvetica", sublines = 0, toplines = 0, leftlines = 0, res=300) 
{
    type <- match.arg(type)
    filename <- paste(filename, ".", type, sep = "")
    if(type=="eps")
    {
        postscript(file = filename, horizontal = FALSE, 
                width = width/2.54, height = height/2.54, pointsize = pointsize, 
                family = family, onefile = FALSE, print.it = FALSE)
    }
    else if(type=="pdf")
    {
        pdf(file = filename, width=width/2.54, height=height/2.54, pointsize=pointsize,
            family=family, onefile=TRUE)
        par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
    }
    else if(type=="jpg")
    {
        jpeg(filename=filename, width=width, height=height, res=res,quality=100, units="cm")#, pointsize=pointsize*50)
    }
    else if(type=="png")
    {
        png(filename=filename, width=width, height=height, res=res, units="cm")#, pointsize=pointsize*50)
    }
    else
        stop("Unknown file type")
    par(mgp = c(2.2, 0.45, 0), tcl = -0.4, mar = c(3.2 + sublines + 0.25 * (sublines > 0), 
         3.5 + leftlines, 1 + toplines, 1) + 0.1)
    par(pch = 1)
    invisible()
}

savepdf = function(...)
{
    savefig(...,type="pdf")
}
