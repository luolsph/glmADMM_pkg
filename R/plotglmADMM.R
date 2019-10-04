plotglmADMM<-function(fit){
	if (fit$intercept==0){
		matplot(log(fit$lambda),t(fit$solution),type="l",lty=1,
	xlab="log lambda",ylab="solution of coefficients");
		}else{
			#do not plot intercept since it is not shrinked.
        matplot(log(fit$lambda),t(fit$solution)[,-1],type="l",lty=1,
        	xlab="log lambda",ylab="solution of coefficients");
		}
}
