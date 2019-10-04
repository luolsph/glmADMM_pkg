glmADMM <-
function(A,b,lambda_seq,family,rho=1,intercept=TRUE,equality=FALSE,C,d,inexact=FALSE){


#start the clock
ptm<-proc.time()

if (inexact==TRUE) {
  ABSTOL=1e-8;
  RELTOL=1e-6;
  MAX_ITER=2000;
} else {
  ABSTOL=1e-6;
  RELTOL=1e-4;
  MAX_ITER=1000;
}

m=nrow(A);
n=ncol(A);

if (equality==TRUE) {
    if (missing(C)) { #C is not specified, set it to a row of 1.
    C=matrix(1,nrow=1,ncol=n);
    }
} else {
    # no equality constraints
    C=matrix(0,nrow=1,ncol=n); 
}

if (missing(d)) d=0;

if (intercept==TRUE) {
    A=cbind(1,A);
    C=cbind(0,C);
}

p=ncol(A);
#initial value of x and z
x=matrix(0,nrow=p,ncol=1); 
z=matrix(0,nrow=p,ncol=1);
u1=0;
u2=matrix(0,nrow=p,ncol=1); 

#B=rbind(C,diag(p));

if (family=="binomial"){
    F=-diag(b[,1])%*%A; #pre-calculate the -b_i%*%a_i in logistic loss
    if (missing(lambda_seq)) {
        b_hat=matrix(0,m);
        ratio=sum(b==1)/m;
        for (i in 1:m){
            if (b[i]==1) b_hat[i]=1-ratio else 
            b_hat[i]=-ratio
        }
        lambda_max=norm(t(A)%*%b_hat,'i');  
    }
    if (inexact==TRUE) {
        J=t(F)%*%diag((exp(F%*%x)/(1+exp(F%*%x))^2)[,1])%*%F;
    }
} else if (family=="poisson") {
    if (missing(lambda_seq)) {
        lambda_max=norm(t(A)%*%(b-exp(A%*%x)),'i');
    }
    if (inexact==TRUE) {
        J=t(A)%*%diag(exp(A%*%x)[,1])%*%A;
    }
    } else {
    Atb=t(A)%*%b; 
    if (missing(lambda_seq)){
        lambda_max=norm(Atb,'i');
    }
    #cache the factorization
    if (inexact==TRUE) {
        J=t(A)%*%A;
        } else {
            E=rbind(A,sqrt(rho)*C); 
            LU=Choleski_factors(E,rho);
            L=LU[[1]]; #m*m lower triangular
            U=LU[[2]];
        }
}

if (missing(lambda_seq)){
lambda_min=0.01*lambda_max;
step=(log(lambda_max)-log(lambda_min))/99;
log_lambda_seq=seq(log(lambda_min),log(lambda_max),step)
log_lambda_seq=sort(log_lambda_seq,decreasing=TRUE)
lambda_seq=exp(log_lambda_seq)
}

if (inexact==TRUE) {
    eigenvalues=eigen(J)$values;
    h=max(eigenvalues);
    P=1/(h+rho)*(diag(p)-rho/(h+rho*(p+1))*t(C)%*%C);
}

iter_num=c()
x_result=c()

for (i in 1:length(lambda_seq)){
    lambda=lambda_seq[i]; 
  
    for (k in 1:MAX_ITER){
        if (family=="gaussian") {
            if (inexact==TRUE){
                    if (p>50*m) {
                    x=P%*%(rho*(z-u1-u2+d)+h*x-t(A)%*%(A%*%x)+Atb);
                    } else {
                    x=P%*%(rho*(z-u1-u2+d)+h*x-J%*%x+Atb);
                    }
            } else {
                q=Atb+rho*(z-u1-u2+d); #n*1 vector
                if (m>=p){
                    x=backsolve(U,forwardsolve(L,q));
                }else{
                    x=q/rho-(t(E)%*%(backsolve(U,forwardsolve(L,(E%*%q)))))/rho^2;
                }  
            }
          
        } else {
            if (inexact==TRUE) {
                if (family=="binomial") {
                 x=P%*%(h*x+rho*(z-u1-u2+d)-t(F)%*%(exp(F%*%x)/(1+exp(F%*%x))));
                } else {
                x=P%*%(h*x+rho*(z-u1-u2+d)+t(A)%*%b-t(A)%*%exp(A%*%x));
                }
            } else {
            x=NR_update_x(A,F,b,x,u1,u2,z,rho,C,d,family) 
            }
        }

    # z-update 
    zold=z;
    z=x+u2; 
    if (intercept==TRUE) {
    z2=shrinkage(z[-1],lambda/rho);
    z=c(z[1],z2);
    } else {
    z=shrinkage(z,lambda/rho); 
    }
    
    #u-update
    Cx=as.numeric(C%*%x);
    u1=u1+Cx-d;
    u2=u2+(x-z);

    #diagnostics, reporting, termination checks    
    r_norm=vector.2.norm(rbind(Cx-d,x-z));    #r_norm
    s_norm=vector.2.norm(-rho*(z-zold));  #s_norm
    eps_pri=sqrt(p+1)*ABSTOL+RELTOL*max(vector.2.norm(rbind(Cx,x)),vector.2.norm(-z),abs(d)); #eps_pri
    eps_dual=sqrt(p)*ABSTOL+RELTOL*vector.2.norm(rho*(u1+u2)); #eps_dual

    #termination criterion
    if (r_norm<eps_pri && s_norm<eps_dual){
      break
    }
  }
  iter_num=c(iter_num,k)
  x_result=cbind(x_result,t(t(x)))
 }
 time=proc.time()-ptm
 return(list(solution=x_result,time=time,iter=iter_num,lambda=lambda_seq,intercept=p-n))
}
