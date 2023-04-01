require(fastmatrix)

# Returns a n by n matrix with nondiagonalelements from -100 to 100 where the matrix has proportion of 
# zero elements sparsity
vals <- (-10^8:10^8)/10^6
generateSddMatrix <- function(sparsity, n)
{
  matrix <- sample(vals, n^2, replace = TRUE)
  
  sparseVals <- sample(n^2 - n, min(round(sparsity*n^2), n^2 - n), replace = FALSE)
  for (val in sparseVals)
  {
    val2 <- val + 1 + floor((val - 1)/n)
    matrix[val2] <- 0
  }
  for (i in (1: n))
  {
    multiplier <- rexp(1, rate = 1) + 1
    sign <- sample(c(-1, 1), 1)
    rowSum <- 0
    for (j in (1:n))
    {
      rowSum <- rowSum + abs(matrix[(i - 1)*n + j])
    }
    matrix[1 + (i - 1)*(n + 1)] <- rowSum*multiplier*sign  
  }
  matrix(matrix, ncol = n)
}

#Creates a vector u with distance d from vector v
generateVector <- function(v, d, n)
{
  u <- runif(n, -1, 1)
  u <- u / sqrt(sum(u^2))
  
  v + u*d
}

# Does SOR and returns err at each step
SOR <- function(w, A, b, x0, tol, max_iter) {
  n <- nrow(A)
  
  k <- 0
  error <- tol + 1
  x <- x0
  errors <- c(0)
  empty <- TRUE
  while (error > tol && k < max_iter) {
    k <- k + 1
    
    for (i in 1:n) {
      x[i] <- (1 - w) * x[i] + (w / A[i,i]) * (b[i] - A[i,] %*% x + A[i,i] * x[i])
    }
    
    
    # Calculate the error
    error <- max(abs(A %*% x - b))
    if (empty)
    {
      errors <- c(error)
      empty <- FALSE
    }
    else
    {
      errors <- c(errors, error)
    }
  }
  
  errors
}




gaussian_elimination <- function(A, b, n) {

  for (k in 1:(n-1)) {
    for (i in (k+1):n) {
      factor <- A[i, k] / A[k, k]
      A[i, (k+1):n] <- A[i, (k+1):n] - factor * A[k, (k+1):n]
      b[i] <- b[i] - factor * b[k]
    }
  }
  
  # Backward sub
  x <- (1:n)
  x[n] <- b[n] / A[n, n]
  for (i in (n-1):1) {
    x[i] <- (b[i] - sum(A[i, (i+1):n] * x[(i+1):n])) / A[i, i]
  }
  x
}

