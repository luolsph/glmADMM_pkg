shrinkage <-
function(x,kappa){
	z=pmax(0,x-kappa)-pmax(0,-x-kappa);
}
