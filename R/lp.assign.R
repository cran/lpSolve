lp.assign <- function(cost.mat)
{
	#
	# lp.assign: use lpsolve.dll to solve an assignment problem. This is a linear program
	# with an ixj matrix of decision variables, and i+j constraints: that the rows and
	# columns all add up to one.
	#
	# Arguments: matrix or data.frame of costs
	#
	# Return value: list from lpsolve, including objective and assignments.
	#
	# Check for the lpslink function, dyn.open if needed. (It should have been loaded
	# by the library() function, though.)
	#
	if(!is.loaded(symbol.C("lpslink"))) {
		base <- "d:/sam/students/lpsolve/lp_solve_4.0/lpsolve.dll"
		if(any(names(version) == "language")) {
			options(show.error.messages = FALSE)
			load.ret <- try(dyn.load(base))
			options(show.error.messages = TRUE)
			if(inherits(load.ret, "try-error"))
				stop("Sorry, error loading the lpsolve.dll")
		}
		else load.ret <- try(dyn.open(base))
		if(inherits(load.ret, "Error"))
			stop("Sorry, error loading the lpsolve.dll")
		if(!is.loaded(symbol.C("lpslink")))
			stop("Sorry, lpsolve.dll not loaded")
	}
	#
	# Check that the cost matrix is in fact a matrix; convert from data.frame if needed.
	#
	if(!is.matrix(cost.mat)) stop("Matrix of costs required.")
	if(is.data.frame(cost.mat))
		cost.mat <- as.matrix(cost.mat)
	#
	# Set up the stuff. The direction is 0, for minimization.
	#
	nr <- nrow(cost.mat)
	nc <- ncol(cost.mat)
	if(nr != nc)
		stop("Cost matrix must be square.")
	direction <- as.integer(0)
	varcount <- as.integer(nr * nc)
	objective <- as.double(c(0, c(t(cost.mat))))
	#
	# Set up the row and column constraints. Each
	#
	constcount <- as.integer(nr + nc)
	row.constraints <- array(0, c(nr, nc, nr))
	for(i in 1:nr)
		row.constraints[i,  , i] <- rep(1, nc)
	row.constraints <- matrix(c(row.constraints), nrow = nr)
	row.constraints <- cbind(row.constraints, rep(1, nr), rep(1, nr))
	#
	col.constraints <- array(0, c(nr, nc, nc))
	for(i in 1:nc)
		col.constraints[, i, i] <- rep(1, nr)
	col.constraints <- matrix(c(apply(col.constraints, c(1, 2), t)), nrow = nc, byrow
		 = TRUE)
	col.constraints <- cbind(col.constraints, rep(1, nc), rep(1, nc))
	all.constraints <- rbind(row.constraints, col.constraints)
	all.constraints <- t(all.constraints)
	constvec <- as.double(c(0, c(all.constraints)))
	intcount <- as.integer(varcount)
	intvec <- as.integer(1:varcount)
	#
	# Prepare objective value, solution, and status
	#
	objval <- as.double(0.)
	solution <- as.double(numeric(nc * nr))
	status <- as.double(0.)
	lps.out <- .C("lpslink",
		direction = direction,
		varcount = varcount,
		objective = objective,
		constcount = constcount,
		constvec = constvec,
		intcount = intcount,
		intvec = intvec,
		objval = objval,
		solution = solution,
		status = status)
	lps.out$solution = matrix(lps.out$solution, nr, nc)
	if(any(names(version) == "language"))
		class(lps.out) <- "lp"
	else oldClass(lps.out) <- "lp"
	lps.out
}
