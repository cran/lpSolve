lp.assign <- function(cost.mat, presolve=0, compute.sens=0)
{
	#
	# lp.assign: use lpsolve.dll to solve an assignment problem. This is a linear program
	# with an ixj matrix of decision variables, and i+j constraints: that the rows and
	# columns all add up to one.
	#
	# Arguments: matrix or data.frame of costs
	#  presolve: numeric. Presolve? Default 0. Ignored.
	#  compute.sens: numeric. Compute sensitivities? Default 0 (no).
	#                Any non-zero number means "yes"
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
	# Set up the row and column constraints. Each is of the
	# "=1" type, represented by 3 (for "equals") 1.
	#
	constcount <- as.integer(nr + nc)
	row.constraints <- array(0, c(nr, nc, nr))
	for(i in 1:nr)
		row.constraints[i,  , i] <- rep(1, nc)
	row.constraints <- matrix(c(row.constraints), nrow = nr)
	row.constraints <- cbind(row.constraints, rep(3, nr), rep(1, nr))
	#
	col.constraints <- array(0, c(nr, nc, nc))
	for(i in 1:nc)
		col.constraints[, i, i] <- rep(1, nr)
	col.constraints <- matrix(c(apply(col.constraints, c(1, 2), t)), nrow = nc, byrow
		 = TRUE)
	col.constraints <- cbind(col.constraints, rep(3, nc), rep(1, nc))
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
	#
	# Set up sensitivity stuff
	#
	sens.coef.from <- sens.coef.to <- 0
	duals <- duals.from <- duals.to <- 0
	if(compute.sens) {
		sens.coef.from <- sens.coef.to <- numeric(x.count)
		duals <- duals.from <- duals.to <- numeric(x.count + const.count)
	}
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
		presolve = as.integer(presolve),
		compute.sens = as.integer(compute.sens),
		sens.coef.from = as.double(sens.coef.from),
		sens.coef.to = as.double(sens.coef.to),
		duals = as.double(duals),
		duals.from = as.double(duals.from),
		duals.to = as.double(duals.to),
		status = status, PACKAGE="lpSolve")
	#
	# Make stuff that should be matices into matrices
	#
	lps.out$solution = matrix(lps.out$solution, nr, nc, byrow=TRUE)
	if(length(duals) > 0) {
		lps.out$sens.coef.from <- matrix(lps.out$sens.coef.from, nr, nc, byrow =
			TRUE)
		lps.out$sens.coef.to <- matrix(lps.out$sens.coef.to, nr, nc, byrow = TRUE
			)
		lps.out$duals <- matrix(lps.out$duals, nr, nc, byrow = TRUE)
		lps.out$duals.from <- matrix(lps.out$duals.from, nr, nc, byrow = TRUE)
		lps.out$duals.to <- matrix(lps.out$duals.to, nr, nc, byrow = TRUE)
	}
	if(any(names(version) == "language"))
		class(lps.out) <- "lp"
	else oldClass(lps.out) <- "lp"
	lps.out
}
