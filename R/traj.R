### Generic function to extract time and population size(s) from
### population and interaction models.  The method functions are in
### the same files with the main functions.

`traj` <-
    function(x, ...)
{
    UseMethod("traj")
}

`plot.traj` <-
    function(x, xlab = "Time (t)", ylab = "Population size", lwd =2, lty = 1,
             col =c(2,4), type = "l", ...)
{
    matplot(x[,1], x[,-1], xlab = xlab, ylab = ylab, lwd = lwd, lty = lty,
            col = col, type = type, ...) 
}

`lines.traj` <-
    function(x, ...)
{
    matlines(x[,1], x[,-1], ...)
}
