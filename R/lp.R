lp <- function(direction = "min", objective.in, const.mat, const.dir, const.rhs,
	transpose.constraints = TRUE, int.vec)
{
	if(direction == "min")
		direction <- 0
	else direction <- 1
	if(is.data.frame(objective.in)) {
		if(ncol(objective.in) > 1)
			stop("Objective vector has more than one column")
		objective.in <- unlist(objective.in)
		names(objective.in) <- NULL
	}
	objective <- c(0, objective.in)
	solution <- numeric(length(objective.in))
	status <- objval <- 0
	x.count <- length(objective.in)
	#
	# "objective" and "constraints" get extra 0's on the front.
	#
	if(is.data.frame(const.mat)) {
		cm <- as.numeric(unlist(const.mat))
		names(cm) <- NULL
		const.mat <- matrix(cm, nrow = nrow(const.mat))
	}
	const.mat[is.na(const.mat)] <- 0
	if(transpose.constraints)
		const.mat <- t(const.mat)
	const.dir.num <- rep(-1, length(const.dir))
	const.dir.num[const.dir == "<" | const.dir == "<="] <- 0
	const.dir.num[const.dir == "=" | const.dir == "=="] <- 1
	const.dir.num[const.dir == ">" | const.dir == ">="] <- 2
	if(any(const.dir.num == -1))
		stop("Unknown constraint direction found\n")
	const.count <- ncol(const.mat)
	if(is.data.frame(const.rhs))
		const.rhs <- as.matrix(const.rhs)
	const.rhs <- c(const.rhs)
	names(const.rhs) <- NULL
	big.const.mat <- rbind(const.mat, const.dir.num, const.rhs)
	constraints <- c(0, c(big.const.mat))
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
		status = as.integer(status))
	if(any(names(version) == "language"))
		class(lp.out) <- "lp"
	else oldClass(lp.out) <- "lp"
	return(lp.out)
}
