library(shiny)
library(shinydashboard)

library(fda)
library(plotly)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)

v = function(ell){
   coef = rep(0,D)
   coef[ell] = 1
   fd(coef, BASIS)
}

change = function(k, permutation){
   m = permutation[1:k]
   sum(v(m))
}

rarh1 = function(n , Sigma=NULL, theta, k, SNR){
   #burnin
   l_burnin = n/2
   perm = sample(c(1:D), D)
   data = matrix(0, 21, (n+l_burnin))
   for (i in 1:(n+l_burnin)){
      data[,i] = rnorm(21, 0, Sigma[order(perm)])
   }
   # Psi operator
   Psi = matrix(0, D, D)
   for (i in 1:D){
      for (j in 1:D){
         Psi[i, j] = rnorm(1, 0, t(Sigma[i])%*%Sigma[j])
      }
   }
   #adjust the norm
   Psi = Psi/(2*norm(Psi))
   coef = matrix(0, D, (n+l_burnin))
   coef[ ,1] = data[,1]
   # recursion
   for (i in 2:(n+l_burnin)){
      coef[,i] = Psi %*% coef[ ,i-1] + data[,i]
   }
   dat = coef[, (l_burnin+1):(n+l_burnin)]
   tr = sum(diag(cov(t(dat))))
   v_hat = change(k, perm)
   # pick c based on SNR
   c = SNR*tr/(theta*(1-theta)*sqrt(D))
   Change = v_hat$coefs*sqrt(c/k)
   newdata = dat
   x = theta*n
   for (i in (x+1):n){
      newdata[ ,i] = dat[ ,i] + Change
   }
   fd(newdata, BASIS)
}

iid = function(n , Sigma=NULL, theta, k, SNR){
   #burnin
   l_burnin = n/2
   perm = sample(c(1:D), D)
   data = matrix(0, 21, (n+l_burnin))
   for (i in 1:(n+l_burnin)){
      data[,i] = rnorm(21, 0, Sigma[order(perm)])
   }
   dat = data[, (l_burnin+1):(n+l_burnin)]
   tr = sum(diag(cov(t(dat))))
   v_hat = change(k, perm)
   # pick c based on SNR
   c = SNR*tr/(theta*(1-theta)*sqrt(D))
   Change = v_hat$coefs*sqrt(c/k)
   newdata = dat
   x = theta*n
   for (i in (x+1):n){
      newdata[ ,i] = dat[ ,i] + Change
   }
   fd(newdata, BASIS)
}


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# PCA
k.star1 = function(data, d){
   D = dim(data$coef)[1]
   n = dim(data$coef)[2]
   Fpca=pca.fd(data, nharm=D, centerfns=TRUE)
   eta = Fpca$scores[,1:d]
   Sigma = matrix(0,d,d)
   for (i in 1:d){
      Sigma[i,i] = Fpca$values[i]
   }
   if (d==1){
      Kappa = function(k){ 
         if (k==1){
            eta[1]-mean(eta)
         } else {
            sum(eta[1:k])-k*mean(eta)
         }
      }
   } else {
      Kappa = function(k){ 
         if (k==1){
            eta[1,]-colMeans(eta)
         } else {
            colSums(eta[1:k,])-k*colMeans(eta)
         }
      }
   }
   Q = function(k){
      1/n*t(Kappa(k))%*%solve(Sigma)%*%Kappa(k)
   }
   val = sapply(1:n, function(k) Q(k))
   return(which.max(val))
}
#fully functional (faster)
k.star2 = function(fdata){
   samp = fdata$coefs
   N = ncol(samp)
   Sn=(1:N)
   Sn[1]=0
   for(j in (2:N)){
      Sn[j]= sum((rowSums(samp[,1:j]) - (j/N)*rowSums(samp[,1:N]))^2)
   }
   min(which(Sn==max(Sn)))
}
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
f.fun = function(t, t1, t2){
   2*t1*(t1+2*t2)*exp(2*t2*(t1+t2)*t)*pnorm(-(t1+2*t2)*sqrt(t))-2*t1^2*pnorm(-t1*sqrt(t))
}
max.density = function(t, c1, c2, s1, s2){
   if(t<0){
      t.p = -t
      pdf = f.fun(t.p, c1/s1, c2*s1/s2^2)
   } else{
      pdf = f.fun(t, c2/s2, c1*s2/s1^2)
   }
   return(pdf)
}
int.x = function(x, t1, t2){
   t2/(t1+t2) + 2*t1*(2*pi)^(-1/2)*sqrt(x)*exp(-t1^2*x/2)+
      ((t1*(t1+2*t2))/(t2*(t1+t2)))*exp(2*t2*(t1+t2)*x)*pnorm(-(t1+2*t2)*sqrt(x))-
      (2*t1^2*x+(t1^2+2*t2^2+2*t1*t2)/(t2*(t1+t2)))*pnorm(-t1*sqrt(x))
}
CdF = function(x, c1, c2, s1, s2){
   t1 = c1/s1
   t2 = c2*s1/s2^2
   t1p = c2/s2
   t2p = c1*s2/s1^2  
   if(x>0 & !(is.na(x))){
      cdf = t2/(t1+t2) + int.x(x, t1p, t2p)
      if (is.na(cdf)){out = 1}else(out=cdf)
   } else{
      cdf = -int.x(-x, t1, t2) + t2/(t1+t2)
      if (is.na(cdf)){out = 0}else(out=cdf)
   }
   return(out)
}
newton <- function(f, tol=1e-7, x0=0, N=1000){
   h <- 1e-7
   i <- 1
   x1 <- x0
   p <- numeric(N)
   while(i <= N){
      df.dx <- (f(x0+h)-f(x0))/h
      x1 <- (x0 - (f(x0)/df.dx))
      p[i] <- x1
      i = i+1
      if(abs(x1-x0)<tol & !(is.na(abs(x1-x0)))) break
      x0=x1
   }
   return(p[(i-1)])
}
quant <- function(c.1, c.2, s.1, s.2, p){
   ftn = function(x){
      CdF(x, c1=c.1, c2=c.2, s1=s.1, s2=s.2) - p
   }
   newton(ftn)
}

LongRun <- function(fdobj, h, basis, kerneltype = "bartlett"){
   N = ncol(fdobj$coefs)
   D = nrow(fdobj$coefs) 
   Kernel <- function(i, h) {
      x = i/h
      if (kerneltype == "flat") {
         return(1)
      }
      if (kerneltype == "simple") {
         return(0)
      }
      if (kerneltype == "bartlett") {
         return(1 - x)
      }
      if (kerneltype == "flat_top") {
         if (x < 0.1) {
            return(1)
         } else {
            if (x >= 0.1 & x < 1.1) {
               return(1.1 - x)
            } else {
               return(0)
            }
         }
      }
      if (kerneltype == "parzen") {
         if (x < 1/2) {
            return(1 - 6 * x^2 + 6 * abs(x)^3)
         } else {
            return(2 * (1 - abs(x))^3)
         }
      }
   }
   D_mat = matrix(0, D, D)
   fdobj_centered = center.fd(fdobj)
   # Long Run Cov Est
   for (k in 1:D) {
      for (r in k:D) {
         s = fdobj_centered$coefs[k, 1:N] %*% fdobj_centered$coefs[r, 1:N]
         if (h > 0) {
            for (i in 1:h) {
               a = fdobj_centered$coefs[k, 1:(N - i)] %*% fdobj_centered$coefs[r, 
                                                                               (i + 1):N]
               a = a + fdobj_centered$coefs[r, 1:(N - i)] %*% fdobj_centered$coefs[k, 
                                                                                   (i + 1):N]
               s = s + Kernel(i, h) * a
            }
         }
         D_mat[k, r] = s
         D_mat[r, k] = D_mat[k, r]
      }
   }
   D_mat = D_mat/N
   # compute eigenvalues of the coefficient matrix
   eigen_struct = eigen(D_mat, symmetric = TRUE)
   # compute eigenfunctions using the eigenvalues
   eigenfunc = fd(eigen_struct$vectors, basisobj = basis)
   # using abs(...) is an ad-hoc regularization solution
   list(ef = eigenfunc, ev = abs(eigen_struct$values), covm = D_mat)
}

Change.point.q = function(fundata, h=0, L, alpha, basis){
   N = ncol(fundata$coefs)
   D =  nrow(fundata$coefs)
   # change point estimation (fully functional)
   k.star = k.star2(fundata)
   theta.hat = k.star/N
   # size of the change
   delta.h = mean(fundata[(k.star+1):N]) - mean(fundata[1:k.star])
   norm.d = inprod(delta.h, delta.h)
   LongRunC = LongRun(fdobj=fundata, h=h, basis = basis)
   lambda.hat = LongRunC$ev
   phi.hat = LongRunC$ef
   # trunctaion parameter L, to estimate sigma
   sigma.h = (norm.d)^-1*sum(sapply(1:L, function(j) lambda.hat[j]*(inprod(delta.h,phi.hat[j]))^2))
   sigma.hat = sqrt(sigma.h)
   q1 = quant((1-theta.hat), theta.hat, sigma.hat, sigma.hat, alpha/2)
   q2 = quant((1-theta.hat), theta.hat, sigma.hat, sigma.hat, 1-alpha/2)
   # upper and lower 1-alpha CI
   upper = k.star + q2/norm.d
   lower = k.star + q1/norm.d
   out = c(lower, k.star, upper)
   return(out)
}
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# -- picking d based on TVE
pick.d = function(fdata, c.t){
   D = nrow(fdata$coefs)
   Fpca <- pca.fd(fdata, nharm=D, centerfns=TRUE)
   Lambda = Fpca$values
   tves = sapply(1:D, function(j) sum(Lambda[1:j])/sum(Lambda))
   a = length(which(tves<c.t))
   if(a==0){out=1}else{out=a}
   return(out)
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

sim = function(N, DataK, Setting, theta, k, SNR, TVE){
   if (Setting=="fast"){
      Sigma = 3^-(1:21)
   }
   if (Setting=="slow"){
      Sigma = (1:21)^-1
   }
   if (Setting=="custom"){
      Sigma = c(1,1,1,rep(0,18))
   }
   if (DataK=="IID"){
      fdata = iid(n=N, Sigma=Sigma, theta=theta, k=k, SNR=SNR)
   }
   if (DataK=="FAR"){
      fdata = rarh1(n=N, Sigma=Sigma, theta=theta, k=k, SNR=SNR)
   }
   ff = Change.point.q(fdata, h=0, L=21, alpha=0.1, basis=BASIS)
   dd = pick.d(fdata, TVE)
   fpca = k.star1(fdata,dd)
   
   l = ff[1]
   u = ff[3]
   fff = ff[2] 
   
   datt = data.frame(method=c("FF", "fPCA"),Lower = c(fff-l,NA), Estimate=c(fff,fpca), Upper=c(u-fff, NA))
   p = plot_ly(data = datt, y = method, x = Estimate, mode = "markers", color = method,
               error_x = list(type = "data", symmetric = FALSE,
                              array = Upper, arrayminus = Lower))
   
   dat = data.frame(method=c("FF", "fPCA"),Lower = c(l,NA), Estimate=c(fff,fpca), Upper=c(u, NA))
   list(fdata=fdata, plot=p, table=dat)
}



server <- function(input, output) {
   
   out <- reactive({
      out = sim(as.numeric(input$N), input$DataK, input$Sigma, input$theta, input$k, input$SNR, input$TVE)
   })
   
   output$ff <- renderTable({
      out()$table
   })
   
   output$plot <- renderPlotly({
      out()$plot
   })
   
   output$plotd <- renderPlot({
      plot(out()$fdata, main="Plot of Simulated Data")
   })
}