lp.transport <- function(cost.mat, row.signs, row.rhs, col.signs, col.rhs, presolve = 0, compute.sens = 0)
{
	#
	# lp.transport: use lpsolve.dll to solve a transportation problem. This is a linear program
	# with an ixj matrix of decision variables, and constraints on the row and column sums.
	#
	# Arguments: cost.mat: matrix or data.frame of costs
	#           row.signs: signs for row constraints
	#             row.rhs: values for row constraints
	#           col.signs: signs for column constraints
	#             col.rhs: values for column constraints
	#            presolve: Numeric: should we presolve? Default 0 (no); non-0
	#                      values means "yes." Currently mostly ignored.
	#        compute.sens: Numeric: compute sensitivities? Default 0 (no);
	#                      non-zero value means "yes."
	#
	# Return value: list from lpsolve, including objective and optimal values.
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
	# Set up the stuff.
	#
	nr <- nrow(cost.mat)
	nc <- ncol(cost.mat)
	#
	# Ensure that row and column stuff is of the correct size.
	#
	if(is.matrix(row.signs)) row.signs <- as.vector(row.signs)
	if(length(row.signs) != nr)
		stop(paste("Error: We have", length(row.signs), "signs, but", nr, "rows"))
	if(is.matrix(row.rhs))
		row.rhs <- as.vector(row.rhs)
	if(length(row.rhs) != nr)
		stop(paste("Error: We have", length(row.rhs), "rhs's, but", nr, "rows"))
	if(is.matrix(col.signs))
		col.signs <- as.vector(col.signs)
	if(length(col.signs) != nc)
		stop(paste("Error: We have", length(col.signs), "signs, but", nc,
			"columns"))
	if(is.matrix(col.rhs))
		col.rhs <- as.vector(col.rhs)
	if(length(col.rhs) != nc)
		stop(paste("Error: We have", length(col.rhs), "rhs's, but", nc, "rows"))
	#
	# The direction is 0, for minimization.
	#
	direction <- as.integer(0)
	varcount <- as.integer(nr * nc)
	objective <- as.double(c(0, c(t(cost.mat))))
	#
	# Set up the row and column constraints.
	#
	constcount <- as.integer(nr + nc)
	row.constraints <- array(0, c(nr, nc, nr))
	for(i in 1:nr)
		row.constraints[i,  , i] <- rep(1, nc)
	row.constraints <- matrix(c(row.constraints), nrow = nr)
	num.signs <- rep(-1, nr)
	num.signs[row.signs == "<" | row.signs == "<="] <- 1
	num.signs[row.signs == "=" | row.signs == "=="] <- 3
	num.signs[row.signs == ">" | row.signs == ">="] <- 2
	if(any(num.signs == -1))
		stop(paste("Unknown row sign in position ", which(num.signs == -1)[1]))
	row.constraints <- cbind(row.constraints, num.signs, row.rhs)
	#
	col.constraints <- array(0, c(nr, nc, nc))
	for(i in 1:nc)
		col.constraints[, i, i] <- rep(1, nr)
	col.constraints <- matrix(c(apply(col.constraints, c(1, 2), t)), nrow = nc, byrow
		 = TRUE)
	num.signs <- rep(-1, nc)
	num.signs[col.signs == "<" | col.signs == "<="] <- 1
	num.signs[col.signs == "=" | col.signs == "=="] <- 3
	num.signs[col.signs == ">" | col.signs == ">="] <- 2
	if(any(num.signs == -1))
		stop(paste("Unknown column sign in position ", which(num.signs == -1)[
			1]))
	col.constraints <- cbind(col.constraints, num.signs, col.rhs)
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
	lps.out$solution = matrix(lps.out$solution, nr, nc, byrow = TRUE)
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
