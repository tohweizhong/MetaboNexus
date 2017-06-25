#House Blend recipe for generating all plots & functions

#R script for generating 3D plots of scores and circle of correlation

#Functions available
#=======================================================================================
# plot_3d(x,type=c("scores","variables")) : Plots 3D models for PCA
# plot_pls3D(x,type=c("scores","variables")) : Plots 3D models for PLS
# scree(x,line=TRUE,cum=TRUE) : Plots scree plot for PCA
# plot_R2Q2(x,type=c("Q2","Q2cum","PRESS","R2","R2cum")) : Plots PLS quality parameters
# impute(df) : imputes missing value with half of the minimum in the dataframe
# nearZeroVar(x, freqCut = 95/5, uniqueCut = 10, saveMetrics = FALSE) 
#            borrowed from caret package to detect variables with no/little variance
#=======================================================================================

#For drawing circle of correlation

plot3d_var<-function(x,title=NULL, col.score=NULL){
  require(rgl)
  open3d()
    plot3d(x$cor.xt[,1],x$cor.xt[,2],x$cor.xt[,3],type="s",col="lightblue",
           size=0.5,axes=TRUE,xlab="Component 1",ylab="Component 2",zlab="Component 3")
    abclines3d(0,0,0, a=diag(3), col="gray")
    rad.in=1
    spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", 
          emission = gray(0.7), alpha = 0.3,smooth=TRUE,lwd=0.5)
#spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", 
#emission = gray(0.9),smooth=TRUE)
#For creating triangles (borrowed from mixOmics)
x = c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 
      1.2, 1.09, 1.09, 0, 0, 0, 0, 0.035, -0.035, 0, 0.035 * 
        sin(pi/4), -0.035 * sin(pi/4), 0, 0.035 * sin(pi/4), 
      -0.035 * sin(pi/4), 0, 0, 0, 0, 0, 0, 0, 0.035, -0.035, 
      0, 0.035 * sin(pi/4), -0.035 * sin(pi/4))
y = c(0, 0, 0, 0, 0, 0, 0, 0.035, -0.035, 0, 0.035 * 
        sin(pi/4), -0.035 * sin(pi/4), 1.2, 1.09, 1.09, 1.2, 
      1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 0, 
      0.035, -0.035, 0, 0, 0, 0, 0.035 * sin(pi/4), -0.035 * 
        sin(pi/4), 0, -0.035 * sin(pi/4), 0.035 * sin(pi/4))
z = c(0, 0.035, -0.035, 0, 0.035, -0.035, 0, 0, 0, 0, 
      0.035 * sin(pi/4), -0.035 * sin(pi/4), 0, 0.035, 
      -0.035, 0, 0, 0, 0, 0.035 * sin(pi/4), -0.035 * sin(pi/4), 
      0, -0.035 * sin(pi/4), 0.035 * sin(pi/4), 1.2, 1.09, 
      1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 
      1.09)
triangles3d(x = x, y = y, z = z, col = "black")
    play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
 }


plot3d_obs<-function(x,title=NULL, col.score=NULL){
  require(rgl)
  open3d()
  if(length(col.score)>0){
    plot3d(x$scores[,1],x$scores[,2],x$scores[,3],xlab="Component 1",ylab="Component 2",zlab="Component 3")
    text3d(x$scores[,1],x$scores[,2],x$scores[,3],text=rownames(x$scores),col=as.numeric(col.score)+1)
    play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
  }
  else{
    plot3d(x$scores[,1],x$scores[,2],x$scores[,3],xlab="Component 1",ylab="Component 2",zlab="Component 3")
    text3d(x$scores[,1],x$scores[,2],x$scores[,3],text=rownames(x$scores))
    play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
  }
}  



#3D plot for PLS
pls3d_var<-function(x, title=NULL, col.score=NULL){
  require(rgl)
  open3d()
  
    plot3d(x$cor.xyt[,1],x$cor.xyt[,2],x$cor.xyt[,3],type="s",col="lightblue",
           size=0.5,axes=TRUE,xlab="Component 1",ylab="Component 2",zlab="Component 3")
    abclines3d(0,0,0, a=diag(3), col="gray")
    rad.in=1
    spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", 
              emission = gray(0.7), alpha = 0.3,smooth=TRUE,lwd=0.5)
    #spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", 
    #emission = gray(0.9),smooth=TRUE)
    #For creating triangles (borrowed from mixOmics)
    x = c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 
          1.2, 1.09, 1.09, 0, 0, 0, 0, 0.035, -0.035, 0, 0.035 * 
            sin(pi/4), -0.035 * sin(pi/4), 0, 0.035 * sin(pi/4), 
          -0.035 * sin(pi/4), 0, 0, 0, 0, 0, 0, 0, 0.035, -0.035, 
          0, 0.035 * sin(pi/4), -0.035 * sin(pi/4))
    y = c(0, 0, 0, 0, 0, 0, 0, 0.035, -0.035, 0, 0.035 * 
            sin(pi/4), -0.035 * sin(pi/4), 1.2, 1.09, 1.09, 1.2, 
          1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 0, 
          0.035, -0.035, 0, 0, 0, 0, 0.035 * sin(pi/4), -0.035 * 
            sin(pi/4), 0, -0.035 * sin(pi/4), 0.035 * sin(pi/4))
    z = c(0, 0.035, -0.035, 0, 0.035, -0.035, 0, 0, 0, 0, 
          0.035 * sin(pi/4), -0.035 * sin(pi/4), 0, 0.035, 
          -0.035, 0, 0, 0, 0, 0.035 * sin(pi/4), -0.035 * sin(pi/4), 
          0, -0.035 * sin(pi/4), 0.035 * sin(pi/4), 1.2, 1.09, 
          1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 
          1.09)
    triangles3d(x = x, y = y, z = z, col = "black")
    play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
}

pls3d_obs<-function(x, title=NULL, col.score=NULL){
      require(rgl)
      open3d()
      if(length(col.score)>0){
        plot3d(x$x.scores[,1],x$x.scores[,2],x$x.scores[,3],xlab="Component 1",ylab="Component 2",zlab="Component 3")
        text3d(x$x.scores[,1],x$x.scores[,2],x$x.scores[,3],text=rownames(x$x.scores),col=as.numeric(col.score)+1)
        play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
      }
      else{
        plot3d(x$x.scores[,1],x$x.scores[,2],x$x.scores[,3],xlab="Component 1",ylab="Component 2",zlab="Component 3")
        text3d(x$x.scores[,1],x$x.scores[,2],x$x.scores[,3],text=rownames(x$x.scores))
        play3d(spin3d(axis=c(0,0,1), rpm=8), duration=7.5)
      }
    }


#Screeplot for plsdepot NIPALS
scree<-function(x,line=TRUE,cum=TRUE){
  if(cum==TRUE){
  comps<-nrow(x$values)
  labs<-paste(1:comps)
  barplot(x$values[,3],ylim=c(0,100),xlab="Principal Components",ylab="Cumulative % of Explained Variance",
        names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Principal Components Screeplot",border=NA)
  axis(side=2,col="gray75")
  if(line==TRUE){
  abline(h=80,col="red",lty=3)
  }
  }
  if(cum==FALSE){
    comps<-nrow(x$values)
    labs<-paste(1:comps)
    barplot(x$values[,2],ylim=c(0,100),xlab="Principal Components",ylab="% of Explained Variance",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Principal Components Screeplot",border=NA)
    axis(side=2,col="gray75")
    
    }
  }

plot_R2Q2<-function(x,type=c("Q2","Q2cum","PRESS","R2","R2cum")){
  if(type=="Q2"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,3],ylim=c(0,1),xlab="PLS Components",ylab="Q2 value",
           names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Q2 value",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="Q2cum"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,5],ylim=c(0,1),xlab="PLS Components",ylab="Cumulative Q2 value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Cumulative Q2 value",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="PRESS"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    max_h<-round(max(x[,1]),0)
    barplot(x[,1],ylim=c(0,max_h), xlab="PLS Components",ylab="PRESS",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Predicted residual sums of squares",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="RSS"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    max_h<-round(max(x[,2]),0)
    while(max_h%%5!=0){
      max_h<-max_h+1
    }
    barplot(x[,2],ylim=c(0,max_h),xlab="PLS Components",ylab="RSS",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Residual sums of squares",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="R2"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x,ylim=c(0,1),xlab="PLS Components",ylab="R2 value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="R2",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="R2Xcum"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,2],ylim=c(0,1),xlab="PLS Components",ylab="Cumulative R2X value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Cumulative R2X",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="R2X"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,1],ylim=c(0,1),xlab="PLS Components",ylab="R2X value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="R2X",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="R2Y"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,3],ylim=c(0,1),xlab="PLS Components",ylab="R2Y value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="R2Y",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="R2Ycum"){
    x1<-data.frame(x)
    comps<-nrow(x1)
    labs<-paste(1:comps)
    barplot(x[,4],ylim=c(0,1),xlab="PLS Components",ylab="R2Ycum value",
            names.arg=labs,width=1.1,col="lightblue",cex.names=0.8,main="Cumulative R2Y",border=NA)
    axis(side=2,col="gray75")
  }
  if(type=="y.pred"){
    
    plot(x,ylab="Predicted Y",
            col="skyblue",pch=16,cex=2,main="Predicted Y")
    box(col="gray75")
  }
  if(type=="resid"){
    
    plot(x,ylab="Residuals",
         col="skyblue",pch=16,cex=2,main="Residuals Plot")
    box(col="gray75")
  }
}

###Data control & processing
impute<-function(df){
  
  imp<-min(df[df!=min(df,na.rm=TRUE)],na.rm=TRUE) #inception
  imp<-imp/2
  #less observations than variables
  for(i in 1:nrow(df)){
    
    zeroes<-which(df[i,]==0)
    nas<-which(is.na(df[i,]))
    
    df[i,zeroes]<-imp
    df[i,nas]<-imp
  }
  return (df)
}

###Modification of plsreg1 from plsdepot (Gaston Sanchez) to create model diagnostics

mod_plsreg1<-function (predictors, response, comps = 2, crosval = TRUE) 
{
  X = as.matrix(predictors)
  n = nrow(X)
  p = ncol(X)
  if (p < 2) 
    stop("\npredictors must contain more than one column")
  if (is.null(colnames(X))) 
    colnames(X) = paste(rep("X", p), 1:p, sep = "")
  if (is.null(rownames(X))) 
    rownames(X) = 1:n
  Y = as.matrix(response)
  if (ncol(Y) != 1) 
    stop("\nresponse must be a single variable")
  if (any(is.na(response))) 
    stop("\nresponse must not contain missing values")
  if (nrow(X) != nrow(Y)) 
    stop("\npredictors and response have different number of rows")
  if (is.null(colnames(Y))) 
    colnames(Y) = "Y"
  if (is.null(rownames(Y))) 
    rownames(Y) = 1:n
  if (any(is.na(X))) 
    na.miss = TRUE
  else na.miss = FALSE
  if (!is.null(comps)) {
    nc = comps
    if (mode(nc) != "numeric" || length(nc) != 1 || nc <= 
          1 || (nc%%1) != 0 || nc > min(n, p)) 
      nc = min(n, p)
    if (nc == n) 
      nc = n - 1
  }
  else {
    if (na.miss) {
      crosval = FALSE
      nc = 2
    }
    else {
      if (n >= 10) 
        crosval = TRUE
      else crosval = FALSE
      nc = min(n, p)
    }
  }
  if (!is.logical(crosval)) 
    crosval = FALSE
  Xx = scale(X)
  Yy = scale(Y)
  X.old = Xx
  Y.old = Yy
  Th = matrix(NA, n, nc)
  Ph = matrix(NA, p, nc)
  Wh = matrix(NA, p, nc)
  Uh = matrix(NA, n, nc)
  ch = rep(NA, nc)
  Hot = matrix(NA, n, nc)
  hlim = rep(NA, nc)
  if (crosval) {
    RSS = c(n - 1, rep(NA, nc))
    PRESS = rep(NA, nc)
    Q2 = rep(NA, nc)
    sets_size = c(rep(n%/%10, 9), n - 9 * (n%/%10))
    obs = sample(1:n, size = n)
    segments = vector("list", length = 10)
    ini = cumsum(sets_size) - sets_size + 1
    fin = cumsum(sets_size)
    for (k in 1:10) segments[[k]] = obs[ini[k]:fin[k]]
  }
  w.old = rep(1, p)
  t.new = rep(1, n)
  p.new = rep(NA, p)
  h = 1
  repeat {
    if (na.miss) {
      for (j in 1:p) {
        i.exist = which(complete.cases(X[, j]))
        w.old[j] = sum(X.old[i.exist, j] * Y.old[i.exist])
      }
      w.new = w.old/sqrt(sum(w.old^2))
      for (i in 1:n) {
        j.exist = which(complete.cases(X[i, ]))
        t.new[i] = sum(X.old[i, j.exist] * w.new[j.exist])
      }
      for (j in 1:p) {
        i.exist = intersect(which(complete.cases(X[, 
                                                   j])), which(complete.cases(t.new)))
        p.new[j] = sum(X.old[i.exist, j] * t.new[i.exist])/sum(t.new[i.exist]^2)
      }
      c.new = t(Y.old) %*% t.new/sum(t.new^2)
      u.new = Y.old/as.vector(c.new)
    }
    if (!na.miss) {
      w.old = t(X.old) %*% Y.old/sum(Y.old^2)
      w.new = w.old/sqrt(sum(w.old^2))
      t.new = X.old %*% w.new
      p.new = t(X.old) %*% t.new/sum(t.new^2)
      c.new = t(Y.old) %*% t.new/sum(t.new^2)
      u.new = Y.old/as.vector(c.new)
      if (crosval) {
        RSS[h + 1] = sum((Y.old - t.new %*% c.new)^2)
        press = rep(0, 10)
        for (i in 1:10) {
          aux = segments[[i]]
          Xy.aux = t(X.old[-aux, ]) %*% Y.old[-aux]
          wh.si = Xy.aux %*% sqrt(solve(t(Xy.aux) %*% 
                                          Xy.aux))
          th.si = X.old[-aux, ] %*% wh.si
          ch.si = t(Y.old[-aux]) %*% th.si %*% solve(t(th.si) %*% 
                                                       th.si)
          ch.si = as.vector(ch.si)
          Yhat.si = ch.si * X.old[aux, ] %*% wh.si
          press[i] = sum((Y.old[aux] - Yhat.si)^2)
        }
        PRESS[h] = sum(press)
        Q2[h] = 1 - PRESS[h]/RSS[h]
      }
    }
    Y.old = Y.old - (t.new %*% c.new)
    X.old = X.old - (t.new %*% t(p.new))
    Th[, h] = t.new
    Ph[, h] = p.new
    Wh[, h] = w.new
    Uh[, h] = u.new
    ch[h] = c.new
    Hot[, h] = (n/(n - 1)) * t.new^2/(sum(t.new^2)/(n - 1))
    hlim[h] = qf(0.95, h, n - h) * (h * (n^2 - 1))/(n * (n - 
                                                           h))
    if (is.null(comps) && crosval) {
      if (Q2[h] < 0.0975 || h == nc) 
        break
    }
    else {
      if (h == nc) 
        break
    }
    h = h + 1
  }
  if (crosval) {
    q2cum = rep(NA, h)
    for (k in 1:h) q2cum[k] = prod(PRESS[1:k])/prod(RSS[1:k])
    Q2cum = 1 - q2cum
    Q2cv = cbind(PRESS[1:h], RSS[1:h], Q2[1:h], rep(0.0975, 
                                                    h), Q2cum)
    dimnames(Q2cv) = list(1:h, c("PRESS", "RSS", "Q2", "LimQ2", 
                                 "Q2cum"))
    if (is.null(comps)) 
      h = h - 1
  }
  if (!crosval) 
    Q2cv = NULL
  Th = Th[, 1:h]
  Ph = Ph[, 1:h]
  Wh = Wh[, 1:h]
  Uh = Uh[, 1:h]
  ch = ch[1:h]
  Ws = Wh %*% solve(t(Ph) %*% Wh)
  Bs = as.vector(Ws %*% ch)
  if (!na.miss) {
    Br = Bs * (rep(apply(Y, 2, sd), p)/apply(X, 2, sd))
    cte = as.vector(colMeans(Y) - Br %*% apply(X, 2, mean))
    y.hat = as.vector(X %*% Br + cte)
    cor.xyt = cor(cbind(Xx, y = Yy), Th)
  }
  else {
    mu.x <- attributes(Xx)$"scaled:center"
    sd.x <- attributes(Xx)$"scaled:scale"
    X.hat = Th %*% t(Ph) %*% diag(sd.x, p, p) + matrix(rep(mu.x, 
                                                           each = n), n, p)
    Br = Bs * (rep(apply(Y, 2, sd), p)/sd.x)
    cte = as.vector(colMeans(response) - Br %*% mu.x)
    y.hat = as.vector(X.hat %*% Br + cte)
    cor.xyt = matrix(NA, p + 1, h)
    for (j in 1:p) {
      i.exist <- which(complete.cases(X[, j]))
      cor.xyt[j, ] = cor(Xx[i.exist, j], Th[i.exist, ])
    }
    cor.xyt[p + 1, ] = cor(Yy, Th)
  }
  resid = as.vector(Y - y.hat)
  R2 = as.vector(cor(Th, Yy))^2
  R2x = cor(X, Th)^2   #added
  R2y = cor(Y, Th)^2   #added
  Rdx = colMeans(R2x)  #added
  Rdy = colMeans(R2y)  #added
  EV = cbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy))  #added
  Rd.mat = matrix(0, h, h) #added
  for (j in 1:h) Rd.mat[1:j, j] = Rdy[1:j]  #added
  VIP = sqrt((Wh^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), h, h)) #added
  #   VIP = round(VIP,2) #added
  R2Xy = t(apply(cor.xyt^2, 1, cumsum))
  T2hot = rbind(hlim[1:h], t(apply(Hot[, 1:h], 1, cumsum)))
  dimnames(Wh) = list(colnames(X), paste(rep("w", h), 1:h, 
                                         sep = ""))
  dimnames(Ws) = list(colnames(X), paste(rep("w*", h), 1:h, 
                                         sep = ""))
  dimnames(Th) = list(rownames(X), paste(rep("t", h), 1:h, 
                                         sep = ""))
  dimnames(Ph) = list(colnames(X), paste(rep("p", h), 1:h, 
                                         sep = ""))
  dimnames(Uh) = list(rownames(Y), paste(rep("u", h), 1:h, 
                                         sep = ""))
  names(ch) = paste(rep("c", h), 1:h, sep = "")
  dimnames(T2hot) = list(c("T2", rownames(X)), paste(rep("H", 
                                                         h), 1:h, sep = ""))
  names(Bs) = colnames(X)
  names(Br) = colnames(X)
  names(resid) = rownames(Y)
  names(y.hat) = rownames(Y)
  names(R2) = paste(rep("t", h), 1:h, sep = "")
  dimnames(EV) = list(paste(rep("t", h), 1:h, sep = ""), c("R2X", 
                                                           "R2Xcum", "R2Y", "R2Ycum"))
  dimnames(VIP) = list(colnames(X), paste(rep("t", h), 1:h, 
                                          sep = ""))
  colnames(R2Xy) = paste(rep("t", h), 1:h, sep = "")
  dimnames(cor.xyt) = list(c(colnames(X), colnames(Y)), colnames(Th))
  res = list(x.scores = Th, x.loads = Ph, y.scores = Uh, y.loads = ch, 
             cor.xyt = cor.xyt, raw.wgs = Wh, mod.wgs = Ws, std.coefs = Bs, 
             reg.coefs = c(Intercept = cte, Br),expvar = EV, VIP=VIP, R2 = R2, R2Xy = R2Xy, 
             y.pred = y.hat, resid = resid, T2 = T2hot, Q2 = Q2cv, 
             y = response)
  class(res) = "plsreg1"
  return(res)
}

#system.time is 2.69sec for 1000X10000 matrix
# impute<-function(df){
#   if(min(df)%in%0||min(df)%in%NA){
#     imp<-min(df[df!=min(df,na.rm=TRUE)],na.rm=TRUE) #inception
#     imp<-imp/2
#   }
#   else{
#     imp<-min(df)
#     imp<-imp/2
#   }
#   
#   #less observations than variables
#   for(i in 1:nrow(df)){
#     
#     zeroes<-which(df[i,]==0)
#     nas<-which(is.na(df[i,]))
#     
#     df[i,zeroes]<-imp
#     df[i,nas]<-imp
#   }
#   return (df)
# }

#system.time is 1.33sec for 1000X10000 matrix (faster, need more trials)
impute<-function(df){
  df<-as.matrix(df)
  if(min(df)%in%0||min(df)%in%NA){
    imp<-min(df[df!=min(df,na.rm=TRUE)],na.rm=TRUE) #inception
    imp<-imp/2
  }
  else{
    imp<-min(df)
    imp<-imp/2
  }
  
  #less observations than variables
#   for(i in 1:nrow(df)){
    
    zeroes<-which(df==0)
    nas<-which(is.na(df))
    DF<-as.matrix(df)
    DF[zeroes]<-imp
    DF[nas]<-imp
#   }
  return (data.frame(DF))
}

#Borrowed from caret package

# nearZeroVar<-function (x, freqCut = 95/5, uniqueCut = 10, saveMetrics = FALSE) 
# {
#   if (is.vector(x)) 
#     x <- matrix(x, ncol = 1)
#   freqRatio <- apply(x, 2, function(data) {
#     t <- table(data[!is.na(data)])
#     if (length(t) <= 1) {
#       return(0)
#     }
#     w <- which.max(t)
#     return(max(t, na.rm = TRUE)/max(t[-w], na.rm = TRUE))
#   })
#   lunique <- apply(x, 2, function(data) length(unique(data[!is.na(data)])))
#   percentUnique <- 100 * lunique/apply(x, 2, length)
#   zeroVar <- (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
#   if (saveMetrics) {
#     out <- data.frame(freqRatio = freqRatio, percentUnique = percentUnique, 
#                       zeroVar = zeroVar, nzv = (freqRatio > freqCut & percentUnique <= 
#                                                   uniqueCut) | zeroVar)
#   }
#   else {
#     out <- which((freqRatio > freqCut & percentUnique <= 
#                     uniqueCut) | zeroVar)
#     names(out) <- NULL
#   }
#   out
# }

nearZeroVar<-function (x, freqCut = 95/5, uniqueCut = 10, saveMetrics = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = 1)
  require(snow)
  cl<-makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
  freqRatio <- parApply(cl,x, 2, function(data) {
    t <- table(data[!is.na(data)])
    if (length(t) <= 1) {
      return(0)
    }
    w <- which.max(t)
    return(max(t, na.rm = TRUE)/max(t[-w], na.rm = TRUE))
  })
  stopCluster(cl)
  lunique <- apply(x, 2, function(data) length(unique(data[!is.na(data)])))
  percentUnique <- 100 * lunique/apply(x, 2, length)
  zeroVar <- (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
  if (saveMetrics) {
    out <- data.frame(freqRatio = freqRatio, percentUnique = percentUnique, 
                      zeroVar = zeroVar, nzv = (freqRatio > freqCut & percentUnique <= 
                                                  uniqueCut) | zeroVar)
  }
  else {
    out <- which((freqRatio > freqCut & percentUnique <= 
                    uniqueCut) | zeroVar)
    names(out) <- NULL
  }
  out
}

parLOG10<-function(X){
  require(snow)
  cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
  parLOG<-snow::parSapply(cl,X,log10)
  snow::stopCluster(cl)
  rownames(parLOG)<-rownames(X)
  parLOG<-data.frame(parLOG)
  
}

parLOG2<-function(X){
  require(snow)
  cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
  parLOG<-snow::parSapply(cl,X,log2)
  snow::stopCluster(cl)
  rownames(parLOG)<-rownames(X)
  parLOG<-data.frame(parLOG)
  
}

LN<-function(X){
  out<-log(X,base=exp(1))
}

parLN<-function(X){
  require(snow)
  cl<-snow::makeCluster(as.numeric(Sys.getenv()["NUMBER_OF_PROCESSORS"]),type="SOCK")
  parLOG<-snow::parSapply(cl,X,LN)
  snow::stopCluster(cl)
  rownames(parLOG)<-rownames(X)
  parLOG<-data.frame(parLOG)
  
}