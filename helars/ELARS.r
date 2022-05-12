library(Matrix)
library(nleqslv)

Diabetes <- read.table("diabetes.data", header = T)

# Scale
Diabetes[-11] <- scale(Diabetes[-11])[ , ]
Diabetes[11] <- scale(Diabetes[11], scale = FALSE)[ , ]
attach(Diabetes)


# Useful global variables
n <- nrow(Diabetes) # Sample size
p <- ncol(Diabetes) - 1 # Parameters
r <- 1 # Depends on the model


# Design matrix
X <- as.matrix(Diabetes[-11])
X.tilde <- cbind(1,X)
X.B <- as.matrix(bdiag(X.tilde, diag(rep(1,r))))

# Response
y <- Diabetes$Y


# Convert theta to xi
theta2xi <- function(theta) {
	if (length(theta) != p+r+1) {
		print(paste("theta2xi error, argument needs size ", p+r+1))
		return(0)
	}

	return(X.B %*% theta)
}

# Convert xi to mu
# MODIFY FOR EACH MODEL
xi2mu <- function(xi) {
	if (length(xi) != n+r) {
		print(paste("xi2mu error, argument needs size ", n+r))
		return(0)
	}
	mu <- -xi[1:442]/(2*xi[443])
	mu <- c(mu, sum(xi[1:442]^2)/(4*xi[443]^2) - 442/(2*xi[443]))
	return(mu)
}

# Convert mu to xi
# MODIFY FOR EACH MODEL
mu2xi <- function(mu) {
	if (length(mu) != n+r) {
		print(paste("mu2xi error, argument needs size ", n+r))
		return(0)
	}

	xi <- 442*mu[1:442]/(mu[443] - sum(mu[1:442]^2))
	xi <- c(xi, -442/(2*(mu[443] - sum(mu[1:442]^2))))
	return(xi)
}


# Convert from mu to eta
mu2eta <- function(mu) {
	if (length(mu) != n+r) {
		print(paste("mu2eta error, argument needs size ", n+r))
		return(0)
	}

	return(t(X.B)%*%mu)
}


theta2eta <- function(theta) {
	if (length(theta) != p+r+1) {
		print(paste("theta2eta error, argument needs size ", p+r+1))
		return(0)
	}

	return (mu2eta(xi2mu(theta2xi(theta))))
}

eta2theta <- function(eta) {
	if (length(eta) != p+r+1) {
		print(paste("eta2theta error, argument needs size ", p+r+1))
		return(0)
	}

	Y <- ginv(X.B) # Y is the generalized inverse of X.B

	return (Y %*% mu2xi(t(Y) %*% eta))
}

# MODIFY FOR EACH MODEL
psi.star <- function(xi) {
	if (length(xi) != n+r) {
		print(paste("psi.star error, argument needs size ", p+r+1))
		return(0)
	}

	xi <- -sum(xi[1:442]^2)/(4*xi[443]) - 442/2*log(-xi[443]) + 442/2*log(pi)
	return(xi)
}

phi.star <- function(mu) {
	if (length(mu) != n+r) {
		print(paste("phi.star error, argument needs size ", n+r))
		return(0)
	}

	return (sum(mu2xi(mu)*mu) - psi.star(mu2xi(mu)))
}

psi <- function(theta) {
	if (length(theta) != p+r+1) {
		print(paste("psi error, argument needs size ", p+r+1))
		return(0)
	}

	return (psi.star(theta2xi(theta)))
}

psi.I <- function(theta,I) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("psi.I error, argument needs size ", p+r+1))
		return(0)
	}

	return (psi(theta*I))
}

phi <- function(theta) {
	if (length(theta) != p+r+1) {
		print(paste("phi error, argument needs size ", p+r+1))
		return(0)
	}

	return (sum(theta2eta(theta) * theta) - psi(theta))
}

phi.I <- function(theta, I) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("phi.I error, argument needs size ", p+r+1))
		return(0)
	}

	return (sum(theta2eta(theta) * theta * I) - psi.I(theta, I))
}

# Divergence in M(I)
div.I <- function(theta1, theta2, I=!logical(p+r+1)) {
	return(phi.I(theta1, I) + psi.I(theta2, I) - sum(theta2eta(theta1) * theta2 * I))
}

# Fisher information matrix, i.e. second derivatives of psi.star
# MODIFY FOR EACH MODEL
fisher <- function(xi) {
	if (length(xi) != n+r) {
		print(paste("fisher error, argument needs size ", n+r))
		return(0)
	}

	inf <- matrix(nrow = 443, ncol = 443)
	for (i in 1:442) {
		for (j in 1:442) {
			inf[i,j] = 0
		}
		inf[i,i] = -1/(2*xi[443])
		inf[i,443] = xi[i]/(2*xi[443]^2)
		inf[443,i] = inf[i,443]
	}
	inf[443,443] = -sum(xi[1:442]^2)/(2*xi[443]^3) + 442/(2*xi[443]^2)
	return(inf)
}

# Compute the Jacobian of the transformation eta to theta
jacob <- function(theta) {
	return(t(X.B) %*% fisher(theta2xi(theta)) %*% X.B)
}

# Convert mixed coordinates to theta. idx.eta is a logical vector, where TRUE position
# corresponds to a eta-coordinate in mix, and false corresponds to 
# an theta-coordinate
mixed2theta <- function(mix, idx.eta, theta.guess = theta1) {
	# Function whose root we want to find. The argument vars refers to the
	# unknown thetas in mix (i.e. theta[idx.eta])
	froot <- function(vars) {
		new.theta <- mix
		new.theta[idx.eta] <- vars
		return(mix[idx.eta] - theta2eta(new.theta)[idx.eta])
	}

	# Jacobian of froot. Again, vars contain the unknown thetas
	fjac <- function(vars) {
		new.theta <- mix
		new.theta[idx.eta] <- vars
		G <- jacob(new.theta)
		# Return a submatrix of G. The element (i,j) belongs to the submatrix
		# if both idx.eta[i]==TRUE and idx.eta[j]==TRUE
		return(-G[which(idx.eta), which(idx.eta)])
	}
	# First guess at the MLE of zero model TODO change maybe
	guess <- theta.guess[idx.eta]

	# Numerical solver may need tweaking
	res <- nleqslv(guess,froot, fjac, xscalm="auto", global = "none", control=list(xtol=1e-16))
	print(res$message)
	if(res$termcd != 1) browser()
	mix[idx.eta] <- res$x
	return(mix)
}


# Compute theta coordinates of the m-projection of theta to M(i,a,I)
m.projection <- function(theta, i, I, a=0) {
	eta <- theta2eta(theta)
	mix <- eta
	mix[!I] = 0
	mix[i] = a
	I[i] = FALSE
	return(mixed2theta(mix,I, theta))
}

# Get the value of component i such that the divergence from cur.theta is
# equal to t.star
ger.component <- function(t.star, i, I, cur.theta = theta1) { 
	r <- c(0, cur.theta[i])
	for(k in 1:10) {
		tmp = mean(r)
		theta.hat <- m.projection(cur.theta, i, I, tmp)
		if (div.I(cur.theta, theta.hat, I) - t.star < 0) {
			r[2] <- tmp
		} else {
			r[1] <- tmp
		}
	}
	return (mean(r))
}


# MLE of the full model
fit <- glm(y ~ X, family = gaussian(link = "identity"))
s2 <- sum(fit$residuals^2)/fit$df.residual # Estimate of the variance
theta1 <- c(fit$coefficients/s2, -1/(2*s2))
names(theta1)[12] <- "S2"


# MLE of the zero model
fit0 <- glm(y ~ 1, family = gaussian(link = "identity"))
s20 <- sum(fit0$residuals^2)/fit0$df.residual # Estimate of the variance
theta0 <- c(coefficients(fit0)/s20, rep(0,10), -1/(2*s20))
names(theta0) <- names(theta1)
eta0 <- theta2eta(theta0)




main <- function() {
	theta <- theta1
	df <- data.frame(t(theta))
	maxdiv <- div.I(theta,theta0)
	df <- cbind(df, 1)
	names(df)[length(theta)+1] <- "Div/max(Div)"

	# Main loop
	I <- !logical(p+r+1) # Variables present
	for (k in 1:p) {
		I.small <- I[(1:p)+1] # Active theta 1 to d
		nn <- sum(I.small) # Number of active covariates
		t <- rep(0,nn)
		ind <- 1
		for (i in which(I.small)+1) { # For every active covariate, compute m-proj and divergence
			theta.bar <- m.projection(theta, i, I)
			t[ind] <- div.I(theta,theta.bar,I)
			ind <- ind+1

			print(paste(k,i))
		}
		t.star <- min(t)
		i.star <- which(I.small)[which.min(t)]

		I.small[i.star] <- FALSE

		print("Get rest of components")
		theta.next <- theta
		for (i in which(I.small)+1) {
			theta.next[i] <- ger.component(t.star, i, I, theta)
		}
		theta.next[i.star + 1] <- 0


		print("Wrapup step")
		# Find theta[0] and theta[12] given eta[0], eta[12]
		mix<-c(eta0[1], theta.next[(1:p)+1])
		if (r>0) mix <- c(mix, eta0[(1:r)+p+1])
		II <- logical(p+r+1)
		II[1] <- TRUE
		if(r > 0) II[(1:r)+p+1] <- TRUE
		theta <- mixed2theta(mix,II,theta)

		#theta <- theta.next
		I[i.star+1] <- FALSE
		df <- rbind(c(theta, div.I(theta,theta0)/maxdiv),df)
	}
	
	return(df)
}
	### PLOTTING
	n.thetas <- p+r+1
	plot(df[,n.thetas+1], df[,2], ylim=c(min(df[,2:n.thetas]), max(df[,2:n.thetas])),
		xlab = "div/max(div)", ylab="value of parameter", main="ELARS, Diabetes data")
	lines(df[,n.thetas+1], df[,2])
	for (k in (2:p)+1) {
		lines(df[,n.thetas+1], df[,k], type="o")
	}
	abline(v=df[ ,n.thetas+1])
	text(1.02, df[p+1, 1:p + 1], 1:p)
	
	############

get.order <- function(df) {
	dims <- dim(df)
	I <- !logical(dims[2])
	I[c(1, dims[2])] <- FALSE
	order <- NULL
	for(i in (dims[1]-1):1) {
		zero <- which(df[i,] == 0 & I)
		order <- c(order, zero)
		I[zero] <- FALSE
	}
	try(print(colnames(X)[order]))
	return(order - 1)
}




