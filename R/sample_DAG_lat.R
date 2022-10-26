sample_DAG_lat <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  Y = DAG$Y
  
  err=matrix(0,nsamps,r)
  for (i in 1:r){
    if (DAG$errors[i] == 1){
      err[,i]=matrix(rt(nsamps,df=5),nsamps)
    } else if (DAG$errors[i]==2){
      err[,i]=matrix(runif(nsamps,-1,1),nsamps)
    } else if (DAG$errors[i]==3){
      err[,i]=matrix(rchisq(nsamps,df=3)-3,nsamps)
    }
  }
  
  err[,Y]=0 #error for diagnosis is zero
  # data=normalizeData(err) # dont need to normalize data here
  data = err
  
  done=which(colSums(G)==0) # variables without parents
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]==1) 
        
        h=intersect(pa,done)
        
        if (setequal(h,pa)){ # if all parents already done
          
          A = (data[,h,drop=FALSE]%*%DAG$weights[h,c,drop=FALSE])
          data[,c]=A+err[,c]
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == r){
      stop=1;
    }
  }
  
  Y0 = data[,Y]
  pY = logistic(Y0)
  for (i in 1:nsamps){
    data[i,Y] = rbinom(n=1,size=1,prob=pY[i]) 
  }
  
  return(list(data=data, E=err, Y0=Y0))
}


logistic <- function(X){
  
  return(1/(1+exp(-X)))
}