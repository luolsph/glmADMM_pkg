Choleski_factors <-
function(E,rho)
{
  m=nrow(E);
  n=ncol(E);
  if (m>=n)
    U=chol(t(E)%*%E+rho*diag(n)) else 
      U=chol(diag(m)+1/rho*(E%*%t(E)))
  
  L=t(U);
  return(list(L,U))
}
