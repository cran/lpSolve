lp <- function(direction = "min", objective.in, const.mat, const.dir, const.rhs,
	transpose.constraints = TRUE, int.vec, presolve = 0, compute.sens = 0)
{
	#
	# lp: solve a general linear program
	#
	# Arguments:
	#     direction: Character: direction of optimization: "min" (default) or "max."
	#  objective.in: Numeric vector (or one-column data frame) of coefficients
#                of objective function
	#     const.mat: Matrix of numeric constraint coefficients, one row  per
#                constraint, one column per variable (unless
#                transpose.constraints =  FALSE; see below).
	#     const.dir: Vector of character strings giving the direction of the
#                constraints: each value should be one of "<," "<=," "=," "==,"
#                ">," or ">=."
	#     const.rhs: Vector of numeric values for the right-hand sides of  the
#                constraints.
	# transpose.constraints: By default each constraint occupies a row  of
#                const.mat, and that matrix needs to be transposed before
#                being passed  to the optimizing code.  For very large
#                constraint matrices it may be wiser  to construct the
#                constraints in a matrix column-by-column. In that case set
#                transpose.constraints to FALSE.
	#       int.vec: Numeric vector giving the indices of variables that are
#                required to be integer. The length of this vector will
#                therefore be the  number of integer variables.
	#	 presolve: Numeric: Should presolve be done (in lp_solve)? Default: 0 (no).
	#                A non-zero value means "yes." Currently mostly ignored.
	#  compute.sens: Numeric: compute sensitivities? Default 0 (no). Any non-zero
	#                value means "yes."
	#
	# Set up the direction.
	#
	if(direction == "min")
		direction <- 0
	else if (direction == "max")
                direction <- 1
             else stop ("Direction must be 'max' or 'min'")
	#
	# Convert one-column data frame objective to vector. Add leading 0 to obejctive.
	#
	if(is.data.frame(objective.in)) {
		if(ncol(objective.in) > 1)
			stop("Objective vector has more than one column")
		objective.in <- unlist(objective.in)
		names(objective.in) <- NULL
	}
	#
	# Set up solution, status, x.count (= number of variables)
	#
	objective <- c(0, objective.in)
	solution <- numeric(length(objective.in))
	status <- objval <- 0
	x.count <- length(objective.in)
	#
	# Convert "constraints" to a matrix if necessary; set NAs to 0.
	#
	if(is.data.frame(const.mat)) {
		cm <- as.numeric(unlist(const.mat))
		names(cm) <- NULL
		const.mat <- matrix(cm, nrow = nrow(const.mat))
	}
	const.mat[is.na(const.mat)] <- 0
	#
	# Transpose if necessary.
	#
	if(transpose.constraints)
		const.mat <- t(const.mat)
	#
	# Set up constraint signs...
	#
	const.dir.num <- rep(-1, length(const.dir))
	const.dir.num[const.dir == "<" | const.dir == "<="] <- 1
	const.dir.num[const.dir == "=" | const.dir == "=="] <- 3
	const.dir.num[const.dir == ">" | const.dir == ">="] <- 2
	if(any(const.dir.num == -1))
		stop("Unknown constraint direction found\n")
	#
	# ...constraint count, and right-hand sides.
	#
	const.count <- ncol(const.mat)
	if(is.data.frame(const.rhs))
		const.rhs <- as.matrix(const.rhs)
	const.rhs <- c(const.rhs)
	names(const.rhs) <- NULL
	#
	# Set up big matrix of constraint info; add a 0 on the front.
	#
	big.const.mat <- rbind(const.mat, const.dir.num, const.rhs)
	constraints <- c(0, c(big.const.mat))
	#
	# Set up int.vec.
	#
	if(missing(int.vec)) {
		int.count <- 0
		int.vec <- 0
	}
	else {
		int.count <- length(int.vec)
	}
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
	# Set up sensitivity stuff.
	#
	sens.coef.from <- sens.coef.to <- 0
	duals <- duals.from <- duals.to <- 0
	if(compute.sens != 0) {
		sens.coef.from <- sens.coef.to <- numeric(x.count)
		duals <- duals.from <- duals.to <- numeric(x.count + const.count)
	}
	#
	lp.out <- .C("lpslink",
		direction = as.integer(direction),
		x.count = as.integer(x.count),
		objective = as.double(objective),
		const.count = as.integer(const.count),
		constraints = as.double(constraints),
		int.count = as.integer(int.count),
		int.vec = as.integer(int.vec),
		objval = as.double(objval),
		solution = as.double(solution),
		presolve = as.integer(presolve),
		compute.sens = as.integer(compute.sens),
		sens.coef.from = as.double(sens.coef.from),
		sens.coef.to = as.double(sens.coef.to),
		duals = as.double(duals),
		duals.from = as.double(duals.from),
		duals.to = as.double(duals.to),
		status = as.integer(status), PACKAGE="lpSolve")
        lp.out$objective <- objective.in
        lp.out$constraints <- big.const.mat
	if(any(names(version) == "language"))
		class(lp.out) <- "lp"
	else oldClass(lp.out) <- "lp"
	return(lp.out)
}
