NR_update_x <-
function(A,F,b,x,u1,u2,z,rho,C,d,family){
    tolerance=1e-5;
    max_iter=50;
    m=nrow(A);
    n=ncol(A);
    Ctc=t(C)%*%C;
    #anonymous function
    #f=(function(w) sum(log(1+exp(F%*%w)))+(rho/2)*vector.2.norm(w-z+u2)^2+(rho/2)*(as.numeric(C%*%w)+u1)^2);
    for (j in 1:max_iter) {
        if (family=="binomial"){
            Fx=exp(F%*%x);
            g=t(F)%*%(Fx/(1+Fx))+rho*(x-z+u2)+rho*(as.numeric(C%*%x)+u1-d);
            H=t(F)%*%diag((Fx/(1+Fx)^2)[,1])%*%F+rho*diag(n)+rho*Ctc;
        } else {
            Ax=exp(A%*%x);
            g=t(A)%*%(Ax-b)+rho*(x-z+u2)+rho*(as.numeric(C%*%x)+u1-d);
            H=t(A)%*%diag(Ax[,1])%*%A+rho*diag(n)+rho*Ctc;
        }
        dx=-solve(H,g); #Newton step (n)*1
        dfx=t(g)%*%dx; 
        if (abs(dfx)<tolerance){
          break
        }
        x=x+dx;
    }
    return(x);
}
