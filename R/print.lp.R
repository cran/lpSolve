print.lp <- function(x, ...)
{
	if(x$status == 0)
		cat("Success: the objective function is", x$objval, "\n")
	else if(x$status == 2)
		cat("Error: no feasible solution found")
	else cat("Error: status", x$status, "\n")
}
