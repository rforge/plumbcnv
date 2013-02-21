options(stringsAsFactors=FALSE)


########### DEFINE FUNCTIONS ############

cmp.there <- require(compiler)


narm <- function(X) { return(X[!is.na(X)]) }


matmul <- function(A, x, transpose=FALSE)
{
   # bigalgebra friendly version of matrix multiplier function
   if(transpose) {
    return(t( t(x) %*% A)) 
   } else {
    return (A %*% x) 
   }
}

# mean replacement code not used at this stage
row.rep <- function(X) { X[is.na(X)] <- mean(X,na.rm=T); X }


spc <- function(X,char=" ") { 
  ## print 'X' spaces
  lX <- length(X); out <- rep("",lX)
  for(j in 1:lX) {
    if(X[j]>0) out[j] <- paste(rep(char,times=X[j]),collapse="")
  } 
  return(out) 
}


rmv.spc <- function(str,before=T,after=T,ch=" ") {
  # remove spaces at start and end of string
  kk <- (length(str))
  if(length(str)<1) { return(str) }
  for (cc in 1:kk) {
    if(before){
      while(substr(str[cc],1,1)==ch) {
       if(nchar(str[cc])>1) {
        str[cc] <- substr(str[cc],2,nchar(str[cc])) 
       } else {
        str[cc] <- gsub(ch,"",str[cc])
       }
      }
    }
    if(after) {
      while(substr(str[cc],nchar(str[cc]),nchar(str[cc]))==ch) {
       if(nchar(str[cc])>1) {
         str[cc] <- substr(str[cc],1,nchar(str[cc])-1)
       } else {
         str[cc] <- gsub(ch,"",str[cc])
       }
      }
    }
  }
  return(str)
}
  

printBigSummary <- function(bigMat,name=NULL,dat=T,descr=NULL,bck=NULL,mem=F,rw=3,cl=2,...) {
  # print summary of big matrix contents
  if(!is.null(name)) { 
    cat("Big matrix; '",name,"', with: ",sep="") 
  } else { cat("Big matrix with: ") }
  nC <- ncol(bigMat); nR <- nrow(bigMat)
  cat(nR," rows (SNPs), ",nC," columns (Samples)\n",sep="") 
  if(is.sub.big.matrix(bigMat)) { cat("[a sub.big.matrix object]\n")}
  cat(" - data type:",is(bigMat[1,1])[1],"\n")
  if(!is.null(descr)) { cat(" - descriptor file;",descr,"\n") }
  if(is.filebacked(bigMat)) { 
    if(!is.null(bck)) {
      cat(" - backing file;",bck,"\n") }
  } else {
    cat(" - not filebacked! [only recommended when RAM is high versus datasize]")
  }
  if(dat) {
    print.large(bigMat,rw=rw,cl=cl,rlab="SNP-id",clab="Sample IDs",...)
  } else {
    if(!is.null(colnames(bigMat))) {
      cat(" - columns:",paste(colnames(bigMat)[1:max(rw,cl)],collapse=", "),
          "...",colnames(bigMat)[nC],"\n")
    } else { cat(" - no column names\n") }
    if(!is.null(rownames(bigMat))) {
      cat(" -    rows:",paste(rownames(bigMat)[1:max(rw,cl)],collapse=", "),
          "...",rownames(bigMat)[nR],"\n")  
    } else { cat(" - no row names\n") }
  }
  if(mem) {
    total.datapoints <- nR*nC
    disk.est <- round(estimate.memory.for.dat(bigMat))
    cat("Total of",total.datapoints,"data-points, using",disk.est,"GB estimated disk space\n")
  }
  cat("\n")
}

print.large <- function(largeMat,rw=3,cl=2,dg=4,rL="Row#",rlab="rownames",clab="colnames",rownums=T,ret=F) 
{
  # nicely print a large matrix without overloading the output space
  # can return result as lines of text instead of printing to screen (for printing to file)
  # allows customization of row and column labels
  nC <- ncol(largeMat); nR <- nrow(largeMat)
  if(is.null(colnames(largeMat)) | is.null(rownames(largeMat))) {
    print(largeMat[c(1:rw,nR-1,nR),c(1:cl,nC)],digits=dg)
    if(ret) { return("<matrix not returned from print.large(), no colnames/rownames>") }
  } else {
    rD <- spc(min(2,max(nchar(paste(nR)))),".")
    rnD <- spc(min(4,max(nchar(rownames(largeMat)[c(1:rw,nR)]))),".")
    linez <- vector("list",rw+3) #rw,cl =number of rows,cols to print
    rwn <- max(nchar(paste(nR)),nchar(rL))*as.numeric(rownums)
    hdr <- (nchar(colnames(largeMat)[c(1:cl,nC)]))
    hdr[hdr<7] <- 7
    idln <- max(nchar(rlab),nchar(rownames(largeMat)[c(1:rw,nR)]))
    pad <- function(X,L) { paste(spc(L-nchar(X)),X,sep="") }
    if(!ret) { cat("\n"); cat(spc(rwn),spc(idln),clab,"\n") }
    dotz <- "  ...  "; dotzh <- " ..... "; dotzn <- "..."
    # make adjustments if matrix is small enough to display all rows/cols
    if(nC<=cl) { dotz <- dotzh <- "" ; cl <- cl-1 }
    if(nR<=rw) { lstln <- 1 } else {  lstln <- 3 }
    ## make adjustments if not displaying rownumbers
    if(!rownums) {
      lstR <- "" ; rD <- ""; jstr <- rep("",times=rw); rL=""
    } else {
      lstR <- nR; jstr <- paste(1:rw)
    }
    linez[[1]] <- c(pad(rL,rwn),pad(rlab,idln),pad(colnames(largeMat)[c(1:cl)],hdr[1:cl]),
                    dotzh,pad(colnames(largeMat)[nC],tail(hdr,1)))
    for (j in 1:rw) { 
      linez[[j+1]] <- c(pad(jstr[j],rwn),pad(rownames(largeMat)[j],idln),
                        pad(round(largeMat[j,1:cl],dg),hdr[1:cl]),dotz,
                        pad(round(largeMat[j,nC],dg),tail(hdr,1)))
    }
    linez[[rw+2]] <- c(pad(rD,rwn),pad(rnD,idln),pad(rep(dotzn,times=cl),
                       hdr[1:cl]),dotz,pad(dotzn,tail(hdr,1)))
    linez[[rw+3]] <- c(pad(lstR,rwn),pad(rownames(largeMat)[nR],idln),
                       pad(round(largeMat[nR,1:cl],dg),hdr[1:cl]),
                       dotz,pad(round(largeMat[nR,nC],dg),tail(hdr,1)))
    if(!ret) {
      for (j in 1:(rw+lstln)) {
        cat(paste(linez[[j]],collapse=" "),"\n")
      }
    } else {
      # remove last two lines if all rows are displayed
      if(lstln==1) { for(ii in 1:2) { linez[[length(linez)]] <- NULL }  }
      return(linez)
    }
  }
}


estimate.memory.for.dat <- function(dat)
{
  # based on a numeric object, estimate the minimum memory requirement
  if(!is.null(dim(dat))) { dimz <- dim(dat) } else { dimz <- dat }
  if(length(dimz)==1) { dimz[2] <- 1 }
  if(length(dimz)==2) {
    rws <- dimz[1]; cls <- dimz[2]
    cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
    memory.estimate <- as.double((as.double(rws)*as.double(cls))/cells.per.gb)
    return(memory.estimate)
  } else {
    cat("tried to estimate memory for object which is neither dimensions or a dataframe/matrix\n") 
  }
}


LRR.PCA <- function(subDescr,pcs.to.keep=9,SVD=F,LAP=F,saver=T,pcs.fn="PCsEVsFromPCA.RData") 
{
  # run principle components analysis on the SNP subset of the LRR snp x sample matrix
  # various methods to choose from with pro/cons of speed/memory, etc.
  # must use SNP-subset to avoid LD, destroying main effects, +avoid huge memory requirements
  pcaMat <- getBigMat(subDescr,dir)
  cat("\nRunning Principle Components Analysis (PCA), using LRR-data subset:\n")
  printBigSummary(pcaMat,name="pcaMat")
  est.mem <- estimate.memory.for.dat(pcaMat)
  cat(" estimated memory required for",nrow(pcaMat),"x",ncol(pcaMat),"matrix:",round(est.mem,2),
      "GB. If this exceeds available,\n expect PCA to take a long time or fail!\n")
  subMat <- as.matrix(pcaMat)
  rm(pcaMat)
  # center using row means
  cat(" centering data by row means...\n")
  subMat <- subMat - rowMeans(subMat)  #matrix(rep(rowMeans(subMat),times=ncol(subMat)),ncol=ncol(subMat))
  cat(" means for first 10 snps:\n")
  print(round(head(rowMeans(subMat)),10)) # show that centering has worked
  subMat[is.na(subMat)] <- 0 # replace missing with the mean
  cat(" replaced missing data with mean (PCA cannot handle missing data)\n")
    #subMat <- t(subMat) # transpose
  dimz <- dim(subMat)
  if(!SVD & (dimz[2]>dimz[1])) {
    cat(" PCA using 'princomp' (faster for datasets with more samples than markers)\n")
    print(system.time(result <- princomp(t(subMat))))
    PCs <- result$scores[,1:pcs.to.keep]
    Evalues <- result$sdev
  } else {
    if(!SVD) {
      cat(" PCA by crossproduct and solving eigenvectors\n")
      cat(" obtaining crossproduct of the matrix and transpose XtX...")
      uu <-(system.time(xtx <- crossprod(subMat)))
      cat("took",round(uu[3]/60,1),"minutes\n")
      cat(" obtaining eigen vectors of the crossproduct XtX...")
      uu <-(system.time(result <- eigen((xtx/nrow(subMat)),symmetric=T)))
      cat("took",round(uu[3]/60,1),"minutes\n")
      PCs <- result$vectors[,1:pcs.to.keep]
      Evalues <- result$values
    } else {
      cat(" PCA by singular value decomposition...") # La.svd gives result with reversed dims. (faster?)
      if(!LAP) {
        if((require(irlba) & require(bigalgebra))) {
          uu <-(system.time(result <- irlba(subMat,nv=pcs.to.keep,nu=0,matmul=matmul))) 
        } else {
          uu <-(system.time(result <- svd(subMat,nv=pcs.to.keep,nu=0)))
        }
        cat("took",round(uu[3]/60,1),"minutes\n")
        PCs <- result$v[,1:pcs.to.keep]
        Evalues <- result$d
      } else {
        cat("\n [using LAPACK alternative with La.svd]")
        uu <- (system.time(result<- La.svd(subMat,nv=pcs.to.keep,nu=0)))
        cat("took",round(uu[3]/60,1),"minutes\n")
        PCs <- t(result$vt)[,1:pcs.to.keep]  ##?
        Evalues <- result$d
      }
    }
  }
  rownames(PCs) <- colnames(subMat)
  colnames(PCs) <- paste("PC",1:pcs.to.keep,sep="")
  if(saver) {
    ofn <- paste(dir$pc,pcs.fn,sep="")
    cat(paste("saved PC data to file:",ofn,"\n"))
    save(PCs,Evalues,file=ofn) }
  out.dat <- list(PCs,Evalues)
  names(out.dat) <- c("PCs","Evalues")
  return(out.dat)
}


LRR.PCA.correct <- function(pca.result,descr.fn,num.pcs=9,pref="corrected",write=F)
{
  ## using results of a PCA analysis, run correction for 'num.pcs' PCs on a dataset
  # uncorrected matrix
  origMat <- getBigMat(descr.fn,dir)
  cat("\nRunning Principle Components correction (PC-correction), using LRR-dataset:\n")
  printBigSummary(origMat,name="origMat")
  # get filenames now to add to result later
  rN <- rownames(origMat); cN <- colnames(origMat)
  # run pca.correction using eigenvectors (PCs) and eigenvalues from LRR.PCA
  if(!is.list(pca.result)) {
    if(is.character(pca.result)) {
      ofn <- paste(dir$pc,pca.result,sep="") 
      if(file.exists(ofn))
      {
        pca.file <- get(load(ofn))
        cat(" l PCA eigenvalues and eigenvectors\n")
        PCs <- pca.file$PCs
      } else {
        stop("Error: file",ofn,"does not exist\n")
      }
    } else {
      if(ncol(origMat) %in% dim(pca.result))
      {
        #given dimensions = number of samples, assume PCs entered as matrix
        PCs <- pca.result
      } else {
        stop("Error: expecting file name or PC matrix: pca.result\n")
      }
    } 
  } else {
    PCs <- pca.result$PCs
  }

  # create new matrix same size, ready for corrected values
  nR <- nrow(origMat); nC <- ncol(origMat)
  cat(" creating new file backed big.matrix to store corrected data...")
  pcCorMat <- filebacked.big.matrix(nR,nC, backingfile=paste(pref,"Bck",sep=""),
                         backingpath=dir$big, descriptorfile=paste(pref,"Descr",sep=""))
  cat("done\n")
  if(!is.filebacked(pcCorMat) | !is.filebacked(origMat)) {
    cat("Warning: at least one of the big.matrices is not filebacked, memory problems may be encountered\n")
  }
  # in result$vectors, PCs are the columns, only need first 10 or so
  # rows are subjects / samples
  col.sel <- 1:ncol(origMat)
  nPCs <- PCs[,1:num.pcs]
  cat(" correcting by principle components, taking the LRR lm-residual for each SNP\n")
  jj <- proc.time()
  num.snps <- nrow(origMat); sampz <- 1:ncol(origMat)
  stepz <- round(seq(from=1,to=num.snps+1,by=200))
  if((tail(stepz,1)) != num.snps+1) { stepz <- c(stepz,num.snps+1) }
  split.to <- length(stepz)-1
  big.extras <- T # flush memory every 'n' iterations.
  flush.freq <- 20
  
  # this simple way works (instead of big for-loop) but hogs memory and is no faster
  # [NB: requires transpose of target corrected big matrix dimensions]
  ### pcCorMat <- apply(origMat,1,PC.fn,nPCs=nPCs,col.sel=sampz)
  
  for (dd in 1:split.to)
  {
    x1 <- stepz[dd]; x2 <- stepz[dd+1]-1 #subset row selection

    # use of this 'sub.big.matrix' structure, stops the memory leak behaviour which spirals
    # the memory relating to 'origMat' out of control. 
    next.rows <- sub.big.matrix(origMat, firstRow=x1, lastRow=x2, backingpath=dir$big )
    # next.rows is now a pointer to a matrix subset, must use 'as.matrix' to coerce to a regular R object 
    pcCorMat[x1:x2,] <- PC.fn.mat(as.matrix(next.rows),nPCs)
    loop.tracker(dd,split.to)
    
    ## Every 'flush.freq' iterations, clean up the memory, remove the 
    ##  big.matrix object 'pcCorMat' and re-attach it 
    if(dd %% flush.freq == 0) {    
      fl.suc <- flush(pcCorMat) & flush(next.rows)
      if(!fl.suc) { cat("flush failed\n") } 
      gc()  # garbage collection
      if(big.extras) {
        RR <- describe(pcCorMat)
        rm(pcCorMat)
        pcCorMat <- attach.big.matrix(RR,path=dir$big)
      }
    }
    rm(next.rows) # remove the sub-matrix pointer each iteration or this memory builds up 
  }

  options(bigmemory.allow.dimnames=TRUE)
  rownames(pcCorMat) <- rN;  colnames(pcCorMat) <- cN 
  ll <- proc.time()
  cat(paste(" LRR PC-Correction took",round((ll-jj)[3]/3600,3),"hours\n"))
  
  cat("\nPC-corrected dataset produced:\n")
  printBigSummary(pcCorMat,name="pcCorMat")
  
  if(write) {
    mat.ref <- describe(pcCorMat)
    big.fn <- paste("describePCcorrect",num.pcs,".RData",sep="")
    ofn <- paste(dir$big,big.fn,sep="")
    save(mat.ref,file=ofn)
    cat(paste("wrote PC-corrected data description file to file:",ofn,"\n"))
    return(big.fn)
  } else {
    return(mat.ref)
  }
}



PC.fn.mat <- function(next.rows,nPCs)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time)
  col.sel <- 1:ncol(next.rows)
  for (dd in 1:nrow(next.rows)) {
    # compiled PC.fn should speed up these ops a little
    next.rows[dd,] <- PC.fn(next.rows[dd,],nPCs,col.sel) 
  }  
  return(next.rows)
}

PC.fn.mat.apply <- function(nextrows,nPCs)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time), vectorized version
  # testing shows the for-loop (non-vectorized) to be slightly faster, maybe because of t()
  col.sel <- 1:ncol(nextrows)
  nextrows <- t(apply(nextrows,1,PC.fn,nPCs=nPCs,col.sel=col.sel))
  return(nextrows)
}


lrr.sample.dists <- function(bigMat,snp.info.fn,mode,gc=F,dlrs=T,dir,pref="",med.chunk.fn="",
                              plotfn="DistributionsBoxPlot.pdf",tabfn="StatsPerSample.tab",...) 
{
  #calculates sample distributions for mean, optionally DLRS/StDev, GC-wave.
  #passes arguments to calc.chr.info via get.chr.filt.info and calculate.gc.for.sample
  CHR.INFO <- get.chr.filt.info(rownames(bigMat),dir,ret=T, recalcO=T,mode="map3",snp.info.fn=snp.info.fn,...)
  snp.info <- convert.chr.obj.IRanges(CHR.INFO)
  
  stat.mat <- calc.main.stats.on.big(bigMat,dir,deleteDLRS=T,skip.dlrs=!dlrs,apply.method=T)
  if(gc) {
    cat("\nAdding GC wave to LRR stats, be aware this can be slow\n")
    gc.wave <- calculate.gc.for.samples(bigMat,snp.info,dir,med.chunk.fn=med.chunk.fn)
    ## combine GC and main, then report
    c.nms <- rownames(stat.mat)[row.names(stat.mat) %in% row.names(gc.wave)]
    if(!dlrs) { stat.set <- c("Mean","StDev") } else { stat.set <- c("Mean","DLRS") }
    stat.table <- as.data.frame(cbind(stat.mat[c.nms,stat.set],gc.wave[c.nms,2]))
    colnames(stat.table)[ncol(stat.table)] <- "GCWave"
  } else {
    cat(" calculating LRR stats without GC wave (faster)\n")
    stat.table <- stat.mat
  }
  if(!is.na(plotfn)) {
    ofn <- dir.fn.cmb(dir$qc.lrr,plotfn,pref=pref)
    LRR.boxplot(stat.table,ofn)
  }
  if(!is.na(tabfn)) {
    ofn <- dir.fn.cmb(dir$qc.lrr,tabfn,pref=pref)
    write.table(round(stat.table,5),sep="\t",col.names=T,row.names=T,quote=F,file=ofn)
    cat("produced table\n",ofn,"\n")
  } else {
    cat("Warning: sample-wise statistics table not written to file, no file name given\n")
  }
  outlist <- list(stat.table,snp.info,CHR.INFO)
  names(outlist) <- c("stat.table","snp.info","CHR.INFO")
  return(outlist)
}

do.scree.plots <- function(eigenv,dir,fname="ScreePlotPCA.pdf",elbow=9,n.comp=30,printvar=T,writefig=T,...) 
{
  # do SCREE PLOTS AND calculate EIGENVALUE VARIANCE after a PCA
  dir <- validate.dir.for(dir,"pc")
  cat(" generating scree plots for principal components analyses\n")
  if(writefig) {  ofn <- paste(dir$pc,fname,sep="");   pdf(ofn) }
  plot(eigenv[1:n.comp],bty="l",xlab="number of principle components",ylab="eigenvalues",bg="green",pch=21)
  abline(v=(elbow+.5),lty="dashed",...)
  legend("topright",legend=c("Principle components","scree plot 'elbow' cutoff"),
         pt.bg=c("green",NA),pch=c(21,NA),lty=c(NA,"dashed"),bty="n")
  if(writefig) { dev.off() ;  cat(paste("Wrote file:",ofn,"\n")) }
  psuedo.var.pc.post <- eigenv[1:n.comp]/sum(eigenv)
  if(printvar) {
    cat(" sum of eigen-variance:",round(sum(eigenv),2),"\n")
    cat(" variance % estimates: \n ",round(psuedo.var.pc.post,2),"\n")
  }
  return(psuedo.var.pc.post)
}


get.excl.filter.for <- function(object,dir)
{
  # create sample-wise filter for any object, excluding samples in exclusion files directory
  bad.samps <- paste(get.all.samp.fails(dir))
  if (is.list(object))
  {
    for(cc in 1:length(object)) {  object[[cc]] <- get.excl.filter.for(object[[cc]],dir) }
  } else {
    if(is.character(object))
    {
      object <- !(object %in% bad.samps)
    } else {
      if(!is.null(colnames(object))) {
        object <- !(colnames(object) %in% bad.samps)
      } else {
        if(!is.null(rownames(object))) {
          object <- !(rownames(object) %in% bad.samps)
        } else {
          if(is.factor(object)){
            object <- paste(object)
            object <- !(object %in% bad.samps)
          } else {
            cat(" object",is(object)[1],"seemed unsuitable for exclusion of ids, was left unchanged.\n")
          }
        }
      }
    }
  }
  return(object)
}

  
  
exclude.combine.big.mats <- function(bigList,dir,pref=paste("SampQC_Combined",sep=""),debug=F) {
  ## apply exclusions from sample QC and combine all matrices from different 
  ## datasets into a single big.matrix object
  dir <- validate.dir.for(dir,c("big"),warn=F)
  pulp <- function(vec,ch="-") { paste(vec,collapse=ch) }
  if(length(bigList)==1)
  {
    bigUnified <- big.exclude.sort(bigList[[1]], dir=dir, pref=pref)
  } else {
    # obtain all dimensions, names, etc
    multi.file.list <- load.big.pointer.list(bigList,dir)
    bigMatLst <- multi.file.list$bigMatLst  
    all.samps <- multi.file.list$combined.samples
    samp.list <- multi.file.list$samples
    smp.szs <- multi.file.list$grp.sizes
    n.data.files <- length(bigMatLst)
    if(all((rownames(bigMatLst[[1]])==rownames(bigMatLst[[2]])) & 
      (rownames(bigMatLst[[3]])==rownames(bigMatLst[[2]])))) {
      snp.labels <- rownames(bigMatLst[[1]])          } else {
      stop("Error: Rownames (snp labels) of each big.matrix file were different") 
    }
    ## need mask filters for each input sub-dataset
    ## T/F masks remove bad samples from each subset
    f.SL <- get.excl.filter.for(samp.list,dir) # separate filters per subgroup
    f.AS <- get.excl.filter.for(all.samps,dir) # all samples filter
    smp.szsF <- sapply(f.SL,function(X) { length(which(X)) })
    cat(" original sample sizes in separate matrices:",paste(smp.szs,collapse=","),"\n")
    cat(" new sample sizes in combined matrix:",paste(smp.szsF,collapse=","),"\n")
    ## need to adjust final matrix size accordingly and all refs
    
    if(all((diff(as.numeric(sapply(bigMatLst,nrow))))==0))
    {
      des <- paste(pref,"descrFile",sep="_")
      bck <- paste(pref,"bckFile",sep="_")
      nR <- length(snp.labels)
      nC <- length(all.samps);  nCf <- length(all.samps[f.AS])
      cat(" combing matrices, expect this to take some time\n")
      bigUnified <- big.matrix(nR,nCf, backingfile=bck,
                               backingpath=dir$big, descriptorfile=des)
      d2 <- d1 <- 0
      for (dd in 1:n.data.files)
      {
        split.to <- max(1,(smp.szs[dd]) %/% 300) # save RAM without creating groups too small to process
        stepz <- round(seq(from=1,to=smp.szs[dd]+1,length.out=round((split.to+1))))
        if (debug) { print(stepz) ; print(dim(bigUnified)) }
        #^ check this in FN!
        cat(" copying file",dd,"of",n.data.files,"to combined big.matrix object:\n")
        for (cc in 1:split.to)
        {
          # within submatrix cols
          c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
          # within combined matrix cols
          #d1 <- c1+c(0,(cumsum(smp.szs)))[dd]; d2 <- d1 + (c2-c1)
          # do the copying
          lilColRange <- c(c1:c2)[ (f.SL[[dd]])[c(c1:c2)] ] # filter samps in the old sub-matrix
          nvalid <- length(which((f.SL[[dd]])[c(c1:c2)]))
          d1 <- d2+1; d2 <- (d1 + nvalid - 1)
          bigColRange <- c(d1:d2) # filter samps in the new big matrix
          if(debug) {
            cat(" comb'd range:",pulp(range(bigColRange)),"length:",diff(range(bigColRange))+1,"\n")
            cat(" subset range:",pulp(range(lilColRange)),"length:",diff(range(lilColRange))+1,"valid:",nvalid,"\n\n")
          } else {
            loop.tracker(cc,split.to,freq=1)
          }
          if(is.finite(sum(lilColRange))) {
            bigUnified[1:nR,bigColRange] <- bigMatLst[[dd]][1:nR,lilColRange]
          } else {
            cat(" Warning: empty interval ignored\n")
          }
        }
      }
      cat(" combining complete, converting result to big matrix\n")
    } else {
      stop("Error: number of rows (snps) for each group matrix don't match\n")
    }
    options(bigmemory.allow.dimnames=TRUE)
    colnames(bigUnified) <- all.samps[f.AS]
    cat(" added colnames\n")
    rownames(bigUnified) <- snp.labels  
    cat(" added rownames\n")
    descr <- describe(bigUnified)
    flush(bigUnified) # hopefully this will ensure the row/colnames are added to the file backing
  }
  return(descr)
}


load.big.pointer.list <- function(descr.list,dir)
{
  ## taking a list of big matrix descriptors (of some sort) return
  # list of pointers to each big.matrix, plus information about the comprised samples
  n.data.files <- length(descr.list)
  bigMatLst <- samps <- vector("list", n.data.files)
  samp.list <- NULL
  smp.szs <- numeric()
  for (cc in 1:n.data.files)
  { 
    bigMatLst[[cc]] <- getBigMat(descr.list[[cc]],dir) 
    samps[[cc]] <- colnames(bigMatLst[[cc]])
    samp.list <- c(samp.list,samps[[cc]])
    smp.szs[cc] <- length(samps[[cc]])
  }
  outlist <- list(bigMatLst,samp.list,smp.szs,samps)
  names(outlist) <- c("bigMatLst","combined.samples","grp.sizes","samples")
  return(outlist)
}

do.quick.LRR.snp.assocs <- function(descr.fn,sample.info,dir)
{
  ## use sample.info
  # create list of bigMatrix locations
  # go 'N' snps at time, concatenate into 1 file, run regression assocs
  # for 2 phenotypes gives t values - ordinally equivalent to logistic regression with 2 groups
  cat(" running SNP-wise LRR-intensity tests against phenotype to filter most associated")
  bigMat <- getBigMat(descr.fn,dir)
  samp.list <- colnames(bigMat)
  tot.samps <- length(samp.list)
  # check that sample.info is valid, if not attempt to fix it
  sample.info <- validate.samp.info(sample.info,dir,QC.update=F)
  n.phenos <- length(table(sample.info$phenotype,useNA=NULL))
  t.type <- "single"
  if(n.phenos==2) { t.type <- "t.test"}
  if(n.phenos>2) { t.type <- "anova"}
  cat(" found",n.phenos,"phenotypes, ",t.type," will be used to summarise most associated SNPs for LRR vs Phenotype\n")
  three.test <- function(col,pheno) { return(summary(aov(col~pheno))[[1]][["F value"]][1]) }
  two.test <- function(col,pheno) { return((cor.test(col,pheno)$statistic)^2)  }
  ph.test <- switch(t.type,anova=three.test,t.test=two.test,single=NULL)
  if(is.null(ph.test)) { stop("Error: use option for association test by phenotype but there is only 1 type")}
  
  snp.labels <- rownames(bigMat)
  full.size <- length(snp.labels)
  good.mem.lim <- 10^6
  opt.size <- round(good.mem.lim/tot.samps)
  n.segs <- ceiling(full.size/opt.size)
  seg.starts <- 1+c(0:(n.segs+2))*opt.size
  seg.starts[seg.starts>full.size] <- full.size+1 # +1 ensures final snp included
  seg.starts <- seg.starts[!duplicated(seg.starts)]
  results <- vector("list", n.segs)
  pheno <- sample.info$phenotype[rownames(sample.info) %in% samp.list]
  test.seg <- matrix(numeric(),nrow=opt.size,ncol=tot.samps)
  
  # break analysis into chunks that fit in memory
  # NB: avoided parallel computing - it does not seem to add any speed here
  for (dd in 1:n.segs)
  {
    row.subset <- (seg.starts[dd]):(seg.starts[dd+1]-1)
    nr <- length(row.subset)
    # test subset of snps for each segment in turn
    test.seg[1:nr,] <- bigMat[row.subset,] 
    results[[dd]] <- apply(test.seg[1:nr,],1,ph.test,pheno=pheno)
    loop.tracker(dd,n.segs)
  }
  
  Fvalues <- round((do.call(c,results)),2) # remember t^2 = F
  if(length(snp.labels)==length(Fvalues))
  { names(Fvalues) <- snp.labels } else {
    stop("Error: analysis failure as list of SNPs did not match length of F values")
  }
  return(Fvalues)
}


snp.select.least.assoc <- function(sample.info,descr.fn,pc.to.keep=.05,dir)
{
  ##association version of selecting subset of snps for PCA
  # takes the 'pc.to.keep'% least associated with #sample.info$phenotype'
  myFs <- do.quick.LRR.snp.assocs(descr.fn,sample.info,dir)
  myFs <- sort(myFs,decreasing=T)
  n.to.keep <- round(max(min(1,pc.to.keep),0)*length(myFs))
  kept.snps <- names(myFs)[1:n.to.keep]
  return(kept.snps)
}

random.spacing.snp.select <- function(snp.info,pc.to.keep=.05,dir) {  
  # uses the chr, pos, label of each snp, stored in a genoset RangedData object
  # to yield a % (e.g, 5%) subset of snps evenly spaced throughout the genome
  hcl <- get.chr.lens(dir)
  chr.starts <- c(0,cumsum(hcl)[1:21])
  cat(" choosing spaced",round(pc.to.keep*100,1),"% subset of SNPs\n")
  # to use in PCA [only for Chr 6, exclude the MHC region]
  ## MHC chr6:30,018,000-33,606,563; build 37 only slighty earlier
  must.use.package("genoset",T)
  n.chr <- 22; skip.chr.mhc <- 6
  mhc <- c(29500000,34000000) # ok for build 36/37

  ratio = min(1,(1.15*pc.to.keep))  # upward adjust to counter empty bits and exclusions
  ngPCA <- (nrow(snp.info)*ratio) # number of intervals in genome
  av.int <- max(sum(hcl[1:(n.chr-1)])+max(start(snp.info[n.chr])))/ngPCA # mean interval length
  chr.uniq <- list()
  sumo <- 0 # counter for non-unique mappings
  for (cc in 1:n.chr)
  {
    rr <- range(start(snp.info[cc])) # start and end position of chromosome
    l.out <- (rr[2]-rr[1])/av.int
    # create a sequence of evenly spaced locations in the chromosome
    lociz <- seq(from=((av.int/2)+rr[1]),to=(rr[2]-(av.int/2)),length.out=l.out)
    if (cc==skip.chr.mhc) {
      nskip <- length(which(lociz>mhc[1] & lociz<mhc[2]))
      lociz <- lociz[lociz<mhc[1] | lociz>mhc[2]]   
      cat(" skipped",nskip,"locations in MHC region\n")
    }
    # find closest SNPs to each location (lociz)
    lc <- (matchpt(lociz, start(snp.info[cc])))
    # calculate number of locations not matching uniquely
    not.uniq <- length(lc[,1])-length(unique(lc[,1]))
    ##cat("",paste("chr",cc,": for",not.uniq,"locations the most proximal snp was already used\n"))
    sumo <- sumo + not.uniq
    chr.uniq[[cc]] <- unique(lc[,1])
  }
  kpt <- sum(sapply(chr.uniq,length))
  
  cat(paste(" total number of SNPs kept:",kpt,"[",round(kpt/nrow(snp.info)*100,2),"%]\n"))
  cat(paste(" total of",sumo,"removed because they did not uniquely map to 1 location\n"))
  
  # combine snps from each chromosome get the full list of SNP IDs for PCA
  snp.to.pca <- character()
  for (cc in 1:n.chr) { snp.to.pca <- c( snp.to.pca,rownames(snp.info[cc])[chr.uniq[[cc]]] ) }
  return(snp.to.pca)
}


extract.snp.subset <- function(snp.info,sample.info,pc.to.keep=.05,assoc=F,
                               writeResultsToFile=T,big.fn="combinedBigMat.RData",out.fn="pca.snp.subset.txt",dir)
{
  # extract small percentage of SNPs as a subset to commit to PCA to correct batch effects
  dir <- validate.dir.for(dir,c("big","ano"),warn=F)
  if(assoc) {
    ##association version of selecting subset of snps for PCA
    cat("\nSelecting SNP subset for PCA based on the ",round(pc.to.keep*100,1),"% least associated SNPS\n",sep="")
    keep.snps <- snp.select.least.assoc(sample.info,big.fn,pc.to.keep,dir)
  } else {
    cat("\nSelecting SNP subset for PCA based on a ",round(pc.to.keep*100,1),"% random selection of SNPS\n",sep="")
    #subset of SNPs to use (based on even as possible spacing thru genome)
    keep.snps <- random.spacing.snp.select(snp.info,pc.to.keep,dir)  
  }
  # write to file
  if(writeResultsToFile) {
    ofn <- paste(dir$ano,out.fn,sep="")
    writeLines(keep.snps,con=ofn)
    cat(paste("wrote list of SNPs to use in the PCA to file:\n",ofn,"\n"))
  } else {
    cat(" results not written to file\n")
    print(head(keep.snps)) ; cat("\n\t . . . .\n")
  }
  return(keep.snps)
}


batch.effects.density.plot <- function(plot.stats, grp.dat, grp.lab="Plate", npg=8, 
         ylz=c(25,25), xllz=c(-.2,.0), xlhz=c(.1,.5), dir, extrapref="")
{
  # detailed density plot for each batch (e.g, plate) for mean and DLRS distributions
  dir <- validate.dir.for(dir,c("qc.pl"),warn=F)
  #make the histos for all regions
  grps <- levels(as.factor(grp.dat))
  ng <- ceiling(length(grps)/npg)
  colz <- c("red","green","blue","orange","purple","lightblue","pink","brown")
  refcol <- "black"
  statz <- c("Mean","DLRS")
  spotz <- c("topleft","topright")
  
  ofn <- paste(dir$qc.pl,grp.lab,"DistributionsMnDlrs",extrapref,".pdf",sep="")
  pdf(ofn,height=3*ng,width=10)
  
  par(mfcol= c(ng,length(statz)))
  for (tt in 1:2)
  { 
    stat <- paste(statz[tt])
    bg.ref <- density(plot.stats[[stat]])
    yl <- ylz[tt]; xll <- xllz[tt]; xlh <- xlhz[tt]
    spott <- spotz[tt]
    for (pp in 0:(ng-1))
    {
      plot(bg.ref,xlab=paste(stat,"Log-R Ratio"),
           ylab="Frequency amongst samples", col=refcol, lwd=1, ylim=c(0,yl),xlim=c(xll,xlh),
           main=paste("Distribution of",stat,"Log-R Ratio (Call Rate>=.95)"),bty="l")
      XY <- approx(x=bg.ref$x,y=bg.ref$y,n=length(bg.ref$x)*5)
      x <- XY$x; y <- XY$y
      segments(x,0,x,y,col="grey")
      lines(bg.ref, col="white", lwd=1)
      lines(bg.ref, col=refcol, lwd=1,lty="dashed")
      for (cc in 1:npg)
      {
        subrowz <- which(grp.dat==grps[(pp*npg)+cc])
        if(length(subrowz)>1) {
          lines(density(plot.stats[[stat]][subrowz]), col=colz[cc], lwd=1)
        } else {
          # warn if not on the last row of plots
          if (pp!=(ng-1)) {
          cat(" Warning: skipped plotting",stat,"series",cc,"from row",pp,"as the",grp.lab,"had less than 2 samples\n")
          }
        }
      }
      leg.its <- c(paste(grp.lab,narm(grps[(pp*npg)+(1:npg)])),"All (Ref.)")
      legend(spott,legend=narm(leg.its),
             text.col=c(colz[1:(length(leg.its)-1)],refcol),bty="n",ncol=2)
    }
  }
  dev.off()
  cat(paste("produced file:",ofn,"\n"))
}

lrr.boundary.scatter <- function(plot.stats,pass.tab,dir,fn.pre="",fn.suf="",...)
{
  # makes scatterplot for each combination of stats versus each other
  # requires 'plot.stats'; also uses 'pass.stats' which is the reference (may be the same as plot.stats)
  colorz <- "darkblue";   extxt <- ""
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  ## PLOT SD/MEAN/DLRS SCATTERS ##
  stat1 <- c("DLRS","GCWave","GCWave")
  stat2 <- c("Mean","Mean","DLRS")
  stat2p <- c("LRR-Mean","LRR-Mean","DLRS")
  stat1p <- c("DLRS","GCWave","GCWave")
  titx <- c("Sample LRR Mean Versus DLRS", "Sample GC Wave Factor Versus LRR Mean",
           "Sample GC Wave Factor Versus DLRS" )
  cat("produced file(s):\n")
  for (ss in 1:length(stat1))
  {
    if(all(c(stat1[ss],stat2[ss]) %in% colnames(pass.tab)) & all(c(stat1[ss],stat2[ss]) %in% colnames(plot.stats)))
    { 
      # all req'd data is there for this iteration!
      #SD vs MEAN #SD vs DLRS  # MEAN vs DLRS
      ofn <- dir.fn.cmb(dir$qc.lrr,paste(stat1[ss],"vs",stat2[ss],sep=""),
                                         suf=fn.suf,pref=fn.pre,ext=".pdf")
      pdf(ofn)
      valz <- c(pass.tab[c("UB","LB"),stat2[ss]],pass.tab[c("UB","LB"),stat1[ss]])
      plot(plot.stats[[stat2[ss]]], plot.stats[[stat1[ss]]], pch=".",
           xlab=stat2p[ss], ylab=stat1p[ss], ...)
      abline(v=valz[1],lty="dotted")
      abline(v=valz[2],lty="dotted")
      abline(h=valz[3],lty="dotted")
      abline(h=valz[4],lty="dotted")
      bnd <- c("Upper","Lower")
      stt <- rep(c(stat2[ss],stat1[ss]),each=2)
      leg.txt <- paste(stt,bnd,"Bound:",round(valz,3))
      legend("topright",legend=c(titx[ss],
                                 leg.txt),bg="white",col=c((colorz)[1],rep("black",4)),
             lwd=c(10,1,1,1,1),lty="dotted" ,cex=.75)
      dev.off()
      cat(" ",rmv.dir(ofn),"\n")
    } else {
      cat("Warning: skipped plot for ",stat1[ss],"vs",stat2[ss],"as some data was missing\n")      
    }
  }
  cat("\n")
}


get.plate.lrr.stats <- function(plt,stat.table,samps=plt$id)
{
  ## table of plate-wise LRR stats
  plt <- plt[plt$id %in% samps,]
  if(nrow(plt)<2) { stop("Error: no samples were found in plate index table") }
  plate.lrr.stats <- as.data.frame(matrix(ncol=ncol(stat.table)*2,nrow=length(unique(plt$plate))))
  colnames(plate.lrr.stats) <- paste(rep(colnames(stat.table),each=2),rep(c("(Av)","(SD)"),2))
  rownames(plate.lrr.stats) <- levels(as.factor(plt$plate))
  
  for(j in 1:ncol(stat.table))
  {
    plate.lrr.stats[,1+((j-1)*2)] <- tapply(stat.table[plt$id,j],as.factor(plt$plate),mean,na.rm=T)
    plate.lrr.stats[,(j*2)] <- tapply(stat.table[plt$id,j],as.factor(plt$plate),sd,na.rm=T)
  }
  return(plate.lrr.stats)
}


dual.cat <- function(file="",...,to.scrn=T,to.file=T) {
  # simultaneously apply the cat function to the screen and an open connection
  if(to.scrn){ cat(...,file="") }
  if(to.file){ if(isOpen(file)) { cat(...,file=file) } }
}


print.large.to.file <- function(tab,con="",dg=3,rlab="",rownums=F,...,to.scrn=T) 
{
  #print a matrix to file nicely
  linez <- print.large(tab,rw=nrow(tab),cl=ncol(tab),
                       rlab=rlab,dg=dg,rownums=rownums,...,ret=T)
  for (j in 1:(length(linez))) {
    dual.cat(paste(paste(linez[[j]],collapse=" "),"\n",sep=""),file=con,to.scrn=to.scrn)
  }
}


make.QC.summary.table <- function(sample.list,dir,pass.fn="PassQCSamples.txt",sum.fn="QCsummaryTable.txt",to.scrn=T)
{
  # create exclusion comparison table
  # show overlap between samples excluded by various QC metrics
  # use just annotated exclusion lists
  dir <- validate.dir.for(dir,c("ano","qc.lrr"),warn=F)
  sampwise <- sample.bad.count.table(dir,sample.list,type=2)
  keep.samples <- rownames(sampwise)[sampwise$TOTAL==0]
  ofn <- dir.fn.cmb(dir$ano,pass.fn)
   writeLines(keep.samples,con=ofn)
  cat(paste("wrote file of sample IDs passing all QC to:\n",ofn,"\n"))
  
  # prepare index/names for loop of n x n comparisons
  tab.list <- colnames(sampwise[,-which(colnames(sampwise)=="TOTAL")])
  tl <- length(tab.list)
  # unique row at bottom will show samples excluded only on each criteria
  res.tab <- matrix(nrow=(tl+1),ncol=tl)
  # fill in QC table to show # samples excluded by combinations of metrics
  for (ii in 1:(tl+1))
  {
    for (jj in 1:tl)
    {
      if (ii!=(tl+1)) {
        res.tab[ii,jj] <- with ( sampwise, length(which(get(tab.list[ii]) & get(tab.list[jj]))) )
      } else {
        #last row, uniques counts..
        uniqz <- (rowSums(sampwise)==2)
        res.tab[ii,jj] <- with ( sampwise, length(which(uniqz & get(tab.list[jj]))) )
      }
    }
  }
  # tidy up table for display
  colnames(res.tab) <- tab.list
  row.names(res.tab) <- c(tab.list,"No Others")
  tot.uniq <- sum((sampwise$TOTAL))
  
  # make failure counts table
  failedN <- table(rowSums(sampwise))
  failedN.tab <- cbind(failedN,100*round(failedN/sum(failedN),3))
  rownames(failedN.tab) <- c("None","1","2","3","All") #,"4"
  colnames(failedN.tab) <- c("count","%")
  
  if(nchar(paste(sum.fn))>1) {
    ofn <- dir.fn.cmb(dir$qc.lrr,sum.fn,ext=".txt")
    my.fn <- file(ofn,open="w") # write to file
  } else {
    my.fn <- "" # write to display
  }
  
  ## Print tables to screen / and or file
  dual.cat("\ncounts\n",file=my.fn,to.scrn=to.scrn)
  print.large.to.file(res.tab,con=my.fn,to.scrn=to.scrn)  #print((res.tab))
  dual.cat("\npercent of failers\n",file=my.fn,to.scrn=to.scrn)
  fpc.tab <- ((round((res.tab/tot.uniq),3)))
  print.large.to.file(fpc.tab,con=my.fn,to.scrn=to.scrn) #print(fpc.tab)
  dual.cat("\npercent of pass 95% callrate sample\n",file=my.fn,to.scrn=to.scrn)
  pc.tab <- (round((res.tab/length(sample.list)),3))
  print.large.to.file(pc.tab,con=my.fn,to.scrn=to.scrn) #print(pc.tab)
  # print counts table
  dual.cat("\nFrequency count of number of separate QC indices failed per sample\n",
           file=my.fn,to.scrn=to.scrn)
  print.large.to.file(failedN.tab,con=my.fn,to.scrn=to.scrn)
  
  if(my.fn!="") { close(my.fn); cat("\nwrote QC summary to:\n",ofn,"\n") }
  out.list <- list(keep.samples,res.tab,failedN.tab)
  names(out.list) <- c("pass.samples","qcsummary","failcounts")
  return(out.list)
}

wait <- function(dur,unit="s",silent=T) {
  ## do nothing for a period of time
  jj <- proc.time()[3]
  if(unit=="s") { mm <- 1 }
  if(unit=="m") { mm <- 60 }
  if(unit=="h") { mm <- 3600 }
  if(!silent) { timz <- c("hour","minute","second");
                cat("waiting ",dur," ",timz[grep(unit,timz)],"s...",sep="") }
  while((proc.time()[3]-jj)< (mm*dur)) { NULL  }
  if(!silent) { cat("done\n") }
}


LRR.boxplot <- function(stats.table,f.name=NULL,label="Sample")
{
  # make boxplot of LRR for each statistic in 'stats.table'; if file name present to file, 
  # else to screen; add appropriate labels
  if(!is.null(f.name)) {  pdf(f.name) }
    boxplot(stats.table,pch=".",main=paste(label,"Distributions of LRR Statistics"),
          xlab=c("LRR Distribution Metric"), ylab="Value (LRR units)",yaxp=c(-.5,1.25,7),bty="l")
  if(!is.null(f.name)) { dev.off() ;  cat(paste("produced boxplot:\n",f.name,"\n")) }
}

calc.main.stats.on.big <- function(des.fn,dir,dlrs.pref="DLRS",deleteDLRS=T,skip.dlrs=F,apply.method=T) 
{
  ## do the calculation of columnwise Mean, stdev and DLRS(optionally) for a big.matrix
  dir <- validate.dir.for(dir,c("big"),warn=F)
  # attach datafile bigmemory object
  bigMat2 <- getBigMat(des.fn,dir)
  
  if(!skip.dlrs & !apply.method) {
    # make DLRS bigmat
    # imply DLRS description and backing file names
    bck.fn.dlrs <- paste(dlrs.pref,"bckfile",sep="")
    des.fn.dlrs <- paste(dlrs.pref,"descrFile",sep="")
    cat("\nCreating DLRS big.matrix file for calculations\n")
    start.tm <- proc.time()
    ## note that it is sometimes quicker to do this seemingly convoluted creation of matrices
    ## than to calculate using 'apply' with an explicit DLRS function ['colsd' is really fast]
    lastr <- nrow(bigMat2)
    post.set <- sub.big.matrix(bigMat2, firstRow=2,lastRow=lastr, backingpath=dir$big)
    pre.set <- sub.big.matrix(bigMat2, firstRow=1,lastRow=(lastr-1), backingpath=dir$big)
    DLRSMAT <- big.matrix.operation(post.set,pre.set,"-",bck.fn.dlrs,des.fn.dlrs,dir=dir,low.ram=T)
    cat(" [DLRS matrix creation took",round(proc.time()[3]-start.tm[3]),"seconds]\n")
  }
  # calculate mean, sd, median of the big matrix
  cat("\nCalculating sample-wise statistics for LRR\n")
  cat(" calculating column means...")
  uu <- system.time(mN <- colmean(bigMat2,na.rm=T)) # 100 sec
  cat("took",round(uu[3],1),"seconds\n")
  cat(" calculating column SDs...")
  uu <- system.time(sD <- colsd(bigMat2,na.rm=T)) # 30 sec
  cat("took",round(uu[3],1),"seconds\n")
  
  if(!skip.dlrs)
  {
    if(apply.method)
    {
      cat(" calculating column DLRS (apply method) ...")
      dlrs.fast <- function(X) { sd(diff(X),na.rm=T) }
      cmp.there <- require(compiler); if (cmp.there) { dlrs.fast <- cmpfun(dlrs.fast) }
      uu <- system.time(dlrs <- apply(bigMat2,2,dlrs.fast)) # ? sec
      dlrs <- dlrs/sqrt(2)
      cat("took",round(uu[3],1),"seconds\n")
    } else {
      cat(" calculating column dlrs (big.matrix method ...\n")
      uu <- system.time(dlrs <- colsd(DLRSMAT,na.rm=T)) # 80 sec
      cat("took",round(uu[3],1),"seconds\n")
      dlrs <- dlrs/sqrt(2)
      if(deleteDLRS)
      {
        cat(" deleting DLRS bigmemory objects to save disk space...\n")
        bfn <- paste(dir$big,des.fn.dlrs,sep="")
        unlink(bfn); cat("deleted:",bfn,"\n")
        bfn <- paste(dir$big,bck.fn.dlrs,sep="")  
        unlink(bfn); cat("deleted:",bfn,"\n")
      } else {
        cat(" manually delete DLRS bigmemory objects to save disk space if desired\n")    
      }
    }
    out.mat <- cbind(mN,sD,dlrs)
    colnames(out.mat) <- c("Mean","StDev","DLRS")
  } else {
    out.mat <- cbind(mN,sD)
    colnames(out.mat) <- c("Mean","StDev")
  }
  rownames(out.mat) <- colnames(bigMat2)
  return(out.mat)
}

do.venn.QC <- function(tab,dir,pref="LRR_VennDiagram") 
{
  # for up to 3 columns of sample statistic pass/failures in a table, do venn diagram
  must.use.package("limma",bioC=T)
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  if(ncol(tab)>=2)
  {
    # set to all ids to get in bottom right number NOT excluded:
    vc <- vennCounts(tab[,1:min(ncol(tab),3)])

      ofn <- paste(dir$qc.lrr,pref,paste(colnames(vc)[-ncol(vc)],collapse="+"),".pdf",sep="")
      pdf(ofn)
      vennDiagram(vc,main=c(rep(" ",5),"Degree of overlap between samples","outside bounds for LRR QC Metrics"))
      text(2.2,-2.7,"# not rejected")
      dev.off()
    cat(paste("produced file:",ofn,"\n"))
  } else {
    cat("Error: could not generate Venn diagram, input not a table of >=2 columns")
  }
}


trim.plate.duplicates <- function(plate.lookup,by=c("counts","first")[1])
{
  ## remove duplicate entries in a plate.lookup index dataframe that
  # was generated by 'get.plate.info'
  if(!is.list(plate.lookup) | length(plate.lookup)!=3) {
    cat("Warning: expecting a plate.lookup object, a list of length 3")
    cat(" with 1) index; 2) duplicates list (may be empty); 3) counts per plate\n")
    cat("returning unchanged lookup object, made no attempt to detect or remove duplicates\n")
    return(plate.lookup)
  } 
  if(!is.null(plate.lookup[[2]])) {
    pre.len <- nrow(plate.lookup[[1]])
    uniqz <- unique(plate.lookup[[2]][,1])
    cntr <- 0
    for (cc in 1:length(uniqz))
    {
      sub.sel <- which(plate.lookup[[2]][,1]==uniqz[cc])
      sub.tab <- plate.lookup[[2]][sub.sel,]
      if(by=="counts") {
        # keep plate assignment with highest number of samples on it
        # (idea is that this increases the probability this was the true source plate)
        if(sub.tab$count[2]>sub.tab$count[1]) { lose <- 2 } else { lose  <- 1 }
      } else {
        lose <- 2 # i.e, do it by order, keep first occuring plate assignment for each id
      }
      one.to.kill <- paste(plate.lookup[[2]][sub.sel[lose],1:2])
      to.kill <- which(plate.lookup[[1]][,1]==one.to.kill[1] & plate.lookup[[1]][,2]==one.to.kill[2])[1]
      if(length(to.kill)==1) {
        plate.lookup[[1]] <- plate.lookup[[1]][-to.kill,]
        cntr <- cntr+1
      } else {
        stop("error: duplicated plate not found in master list")
      }
    }
    len.chng <- (pre.len-nrow(plate.lookup[[1]]))
    if(cntr==len.chng) {
      cat("To ensure unique id/plate mapping,",cntr,"duplicates were removed from plate table\n")
      cat("[For best accuracy, you should instead manually remove duplicates from the\n")
      cat("source 'plate.lookup.txt' file based on verified annotation]\n")
    } else {
      cat("warning: number of duplicate plates removed",cntr,
          "doesn't match change in table size",len.chng,"\n")
    }
    plate.lookup[[2]] <- plate.lookup[[3]] <- NULL # clear prev info - order important
    plate.lookup[[3]] <- table(plate.lookup[[1]]$plate)
    return(plate.lookup)
  } else {
    #no duplicates
    return(plate.lookup)
  }
}

get.plate.info <- function(dir, id.col=1, plate.col=2, well.col=NA,
                           fn="plate.lookup.txt",dup.action=c("","trim","print","return"))
{
  ## get a table of plate information from a plate lookup file.
  # fairly flexible but lookup must contain plate ids alongside sample ids as
  # a minimum requirement. may include additional information such as 'well'
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  if(fn %in% list.files(dir$ano))
  {
    plate.lookup <- read.delim(paste(dir$ano,fn,sep=""),header=T,stringsAsFactors=F)
    cat("\nRetrieving Plate/Well information from lookup table:\n")
    print(head(plate.lookup,4))
    select.cols <- c(id.col,plate.col)
    if (!is.na(well.col)) { select.cols <- c(select.cols,well.col) }
    plate.lookup <- plate.lookup[,select.cols]
    colnames(plate.lookup) <- c("id","plate","well")[!is.na(c(T,T,well.col))]
    counte <- table(plate.lookup$plate)
    plate.lookup[["count"]] <- (counte[match(plate.lookup$plate,names(counte))])
    if(anyDuplicated(plate.lookup[,1]))
    {
      cat("Warning: duplicate records found in plate info:\n")
      dupz <- plate.lookup[(plate.lookup[,1] %in% plate.lookup[duplicated(plate.lookup[,1]),1]),]
      if (dup.action=="print") { print(dupz) }
      if (dup.action=="return") { return(list(plate.lookup,dupz,counte)) }
      if (dup.action=="trim") {
        plate.obj <- list(plate.lookup,dupz,counte)
        plate.obj <- trim.plate.duplicates(plate.obj,by="counts")
        return(plate.obj)
      }
    } 
    plate.obj <- list(plate.lookup,NULL,counte)
    return(plate.obj)
  } else {
    cat("\nUsing plate QC in this pipeline requires an annotation ")
    cat("file containing each sample id and the plate id. This can be ")
    cat("generated using something like the script in getPlateSupport.sh. ")
    cat("Then ensure the correct column numbers are entered for the function ")
    cat("'get.plate.info' to process this annotation file.\n")
    stop(paste("Error: could not find plate file",fn,"in",dir$ano))
  }
}


excl.bad.plates <- function(plate.qc.tab,dir,badPlateThresh=0.33,writeExclList=T,append=T)
{
  ## get list of plates with more than 'badPlateThresh' percent of samples failing on
  # at least one QC measure. return list of plates, and samples in these plates, by
  # default write the sample information to an exclusion file 'BadPlates.txt'
  taken.out.list <- plate.excl <- NULL
  dir <- validate.dir.for(dir,c("ano"))
  bad.plts <- which((plate.qc.tab$TOTAL/plate.qc.tab$SIZE)>badPlateThresh)
  n.bad <- length(bad.plts)
  cat("\nTesting for plates where more than ",round(100*badPlateThresh,1),"% of samples fail QC\n",sep="")
  if(n.bad>0)
  {
    plate.excl <- plate.qc.tab$ID[bad.plts]
    cat(" found",n.bad,"bad plates\n")
    for (ii in 1:n.bad)
    {  
      nxt.plt <- plate.excl[ii] 
      in.plt <- (plate.lookup[[1]]$id)[plate.lookup[[1]]$plate %in% nxt.plt]
      taken.out.list <- c(taken.out.list,in.plt)
      cat(paste(" removed",length(in.plt),"samples from plate",nxt.plt,"\n"))
    }
  }
  # re-write to original file name
  if(writeExclList & dir$ano!="") {
    ofn <- paste(dir$ano,"SAMPLE_EXCLUDE/BadPlates.txt",sep="")
    # write sample exclude list to annotation directory
    if(!append) {
      writeLines(taken.out.list,con=ofn) 
      cat("\nwrote bad plate excluded subject list to:\n",ofn,"\n")
    } else {
      if(file.exists(ofn)) {
        exist.lns <- readLines(ofn) ####
      } else {
        exist.lns <- NULL
      }
      to.wrt <- paste(unique(c(exist.lns,taken.out.list)))
      writeLines(to.wrt,con=ofn) 
      cat(paste("appended bad plate excluded subject list:\n",ofn,"\n"))     
    }
  } else {
    cat("plate excl. file not written\n") 
    print(head(sub.left)) ;cat("\n\t . . . .\n")
  }
  outlist <- list(taken.out.list,plate.excl)
  names(outlist) <- c("samples.from.bad.plates","bad.plates")
  return(outlist)
} 


chr.ab.report <- function(chr.stat,chrWarns,excluded.samps,dir,writeExclList=F,makeGraphs=F,pref="",append=F)
{
  # make chromosomal abberation graphs, write sample exclusion file 'ChrAb.txt' to dir$ano
  nC <- 22
  dir <- validate.dir.for(dir,c("ano","qc.cab"))
  # put stats in summary table for display
  chr.stats.table <- data.frame(LRR.Mean.Av=round(colMeans(chr.stat$chr.mean,na.rm=T),3),
                                LRR.Mean.Sd=round(apply(chr.stat$chr.mean,2,sd,na.rm=T),3),
                                LRR.StDev.Av=round(colMeans(chr.stat$chr.sd,na.rm=T),3),
                                LRR.StDev.Sd=round(apply(chr.stat$chr.sd,2,sd,na.rm=T),3))
  rownames(chr.stats.table) <- paste("Chr",1:nC,sep="")
  cat("\nTable of LRR-Mean and LRR-SD for each chromosome\n")
  print(chr.stats.table,digits=3)
  
  # print number of outliers found for each chromosome
  cat("\nNumber of outliers per chromosome:\n")
  print((sapply(chrWarns,length)))
  cat("Number of unique samples with warnings:",length(unique(unlist(chrWarns))),"\n")
  cat("Total number of warnings:",length((unlist(chrWarns))),"\n")
  
  cat("\nHistogram of number of chromosome outliers per sample:\n")
  textogram(excluded.samps)
  
  # create box plot of distributions per chromosome
  if(makeGraphs) {
    ofn <- paste(dir$fig,"ChrMeansBoxPlot",pref,".pdf",sep="")
    pdf(ofn)
    boxplot(chr.stat$chr.mean,pch=".",xlab="Chromosome number",ylab="LRR-Mean",
            main="Boxplot of LRR-Mean distributions for each chromosome",bty="l")
    dev.off()
    cat(paste("produced plot:",ofn,"\n"))
  }
  if(writeExclList & dir$ano!="") {
    ofn <- paste(dir$ano,"SAMPLE_EXCLUDE/ChrAb.txt",sep="")
    # write sample exclude list to annotation directory
    if(!append) {
      writeLines(names(excluded.samps),con=ofn) 
      cat("\nwrote chr.ab excluded subject list to:",ofn,"\n")
    } else {
      if(file.exists(ofn)) {
        exist.lns <- readLines(ofn) ####
      } else {
        exist.lns <- NULL
      }
      to.wrt <- unique(c(exist.lns,names(excluded.samps)))
      writeLines(to.wrt,con=ofn) 
      cat(paste(" appended chr.ab excluded subject list:\n",ofn,"\n"))     
    }
  }
}




get.chr.stats <- function(bigMat,CHR.INFO,dir="")
{
 # get the LRR mean and SD for each chromosome (autosome) separately
 # used for detection of chromosomal abberations
 must.use.package(c("bigmemory","biganalytics"))
 dir <- validate.dir.for(dir,c("big"),warn=F)
 nC <- 22; range.chr <- 1:nC
 if((is(CHR.INFO)[1]!="list") | (is(bigMat)[1]!="big.matrix") )
 { stop("Error: parameter 'CHR.INFO' should be a list, and 'bigMat' should be a big.matrix") }
 # Start making histograms of LRR-Mean for each chromosome
 dns=list(colnames(bigMat),paste("chr",range.chr,sep="")) # get rownames,colnames
 # create matrices for LRR-stats
 chr.mean <- chr.sd <- matrix(nrow=ncol(bigMat),ncol=nC,dimnames=dns) # used to be list()
 cat(" processing chr: ")
 for (dd in range.chr) 
 { 
  cat(dd,"..",sep="")
  LRR.dat <- sub.big.matrix(bigMat, firstRow=CHR.INFO$chrStrtsF[dd], 
   lastRow=CHR.INFO$chrEndsF[dd], backingpath=dir$big)
  chr.mean[,dd] <- colmean(LRR.dat,na.rm=T) # 50 sec
  chr.sd[,dd] <- colsd(LRR.dat,na.rm=T)
 }
 cat("done\n")
 out.list <- list(chr.mean,chr.sd)
 names(out.list) <- c("chr.mean","chr.sd")
 return(out.list)
}
 
 
make.chr.fail.table <- function(chr.mean,pctile.bound=.01)
{
 # using chromosome (autosome) means from 'get.chr.stats', create list of
 # failures per sample for their chr.mean to within confidence limits
 range.chr <- 1:22
 headz <- c("ID","LRR.Mean","LoHi","ChrNum")
 listOfFails <- matrix(ncol=4)
 colnames(listOfFails)<- headz
 for (dd in range.chr) 
 { 
  badcuts <- get.pctile(chr.mean[,dd],pctile.bound)
  lo.peeps <- rownames(chr.mean)[chr.mean[,dd] <= badcuts[1]]
  hi.peeps <- rownames(chr.mean)[chr.mean[,dd] >= badcuts[2]]
  new.chunk <- rbind(cbind(lo.peeps,chr.mean[lo.peeps,dd],
          rep("low",length(lo.peeps)),rep(dd,length(lo.peeps))),
          cbind(hi.peeps,chr.mean[hi.peeps,dd],
          rep("high",length(hi.peeps)),rep(dd,length(hi.peeps))) )
  colnames(new.chunk) <- headz # see def'n above
  listOfFails <- rbind(listOfFails, new.chunk)
  if(is.na(listOfFails[1,1])) { listOfFails <- listOfFails[-1,]}
  # record list of all samples exceeding 'pctile.bound' upper or lower, per chr.
 }
 return(listOfFails)
}

 
do.chr.ab.contrasts <- function(chr.mean,lob=2,hib=2.5,pctile.bound=.01,nC=22)
{
 ## Set up contrast matrices to test for chromosomal aberrations (autosomes) ##
 ## two matrices, one a simple average comparison, target vs rest,
 ## the other a local comparison, target versus adjacent weighted
 ## 2^n
 # generate local chromosomal contrast matrix
 # comparison of target chromosome to weighted mean of adjacent chromosomes
 if(ncol(chr.mean)!=nC) { stop("Object should contain 1 column for each chromosome") }
 chromocontr <- matrix(integer(),nrow=nC,ncol=nC)
 # weightings for the 3 chromosomes to the left/right of target chromosome
 i <- c(1,2,4)
 versuz <- c(sum(i),sum(i)+i[3],(2*sum(i))-i[1],rep(12,16),(2*sum(i))-i[1],sum(i)+i[3],sum(i))
 stepo <- list(NULL,-i[3],c(-i[2],-i[3]),c(-i[1],-i[2],-i[3]),c(-i[3],-i[2],-i[1]),c(-i[3],-i[2]),-i[3],NULL)
 for (jj in 1:nC) {
 chromocontr[jj,] <- c( rep(0,max(0,jj-4)),stepo[[max(1,min(jj,4))]],
                    versuz[jj], stepo[[max(5,min(jj-14,8))]],
                     rep(0,max(0,19-jj)) ) /versuz[jj]
 }
 # generate global chromosomal contrast matrix
 # simple comparison of target chromosome to mean of all others
 chromocontr2 <- matrix(-1/(nC-1),nrow=nC,ncol=nC)
 diag(chromocontr2) <- 1
 
 # calculate contrast matrices using matrix multiplication #
 # LOCAL
 chr.conts <- t(chromocontr %*% t(chr.mean))
 colnames(chr.conts) <- paste("chr",1:nC,sep="")
 # GLOBAL
 chr.conts2 <- t(chromocontr2 %*% t(chr.mean))
 colnames(chr.conts2) <- paste("chr",1:nC,sep="")
 
 listOfFails <- make.chr.fail.table(chr.mean)	

 # standardize contrasts to z scores to help determine outliers
 chr.contsZ <- apply(chr.conts,2,StandardizeX)
 chr.contsZ2 <- apply(chr.conts2,2,StandardizeX)

 cat("\nTable of Exclusion processing per chromosome...\n")
 chrWarns <- list()
 for (cc in 1:nC)
 {
  # get T/F vector for all samples if contrast values are HIGH (above thresh)
  hiz <- (chr.contsZ[,cc] > hib & chr.contsZ2[,cc] > lob) | (chr.contsZ2[,cc] > hib & chr.contsZ[,cc] > lob)  
  # get T/F vector for all samples if contrast values are LOW (below thresh)
  loz <- (chr.contsZ[,cc] < -hib & chr.contsZ2[,cc] < -lob) | (chr.contsZ2[,cc] < -hib & chr.contsZ[,cc] < -lob)
 # lookup whether samples also failed on initial raw outlier (percentile) status
  failz <- rownames(chr.contsZ) %in% listOfFails[as.numeric(listOfFails[,"ChrNum"])==cc,"ID"]
  # select list of bad subjects for current chromosome as those failing
  # either the high (hiz) or low (loz) threshold, plus are an outlier in absolute terms (failz)
  cat(paste("chr:",cc,"too high:",length(which(hiz)),"too low:",length(which(loz)),
              "excluded top/bottom",pctile.bound*100,"%:",length(which(failz))),"\n")
  chrWarns[[cc]] <- rownames(chr.contsZ)[which( (hiz | loz) & failz)]
 } 
 # set names of warning list components to allow determination of which
 # chromosome the warning came from once it's been converted to a vector
 names(chrWarns) <- paste("chr",1:nC,"_",sep="")
 return(chrWarns)
}

calculate.gc.for.samples <- function(bigMat,snp.info,dir,med.chunk.fn="",hg.version=18)
{
 # calculate GC-wave stats for each sample
 must.use.package(c("BiocGenerics","IRanges","Biobase"),bioC=T)
 ## parameters in order are:
 ## bigmatrix, directory with bigmatrix, snp.info object, gc.dat object, file name, regenerate?
 dir <- validate.dir.for(dir,c("gc","big","qc.gc"))
 # get genome wide human GC data at megabase resolution 
 # load average data for 10^6 windows prepared using 'extractGCwindows.R'
 gc.sav <- paste(dir$qc.gc,"GCAv6.RData",sep="")
 if(!file.exists(gc.sav))
 {
   gc.dat <- get.gc.human(10^6,hgV=hg.version)
   save(gc.dat,file=gc.sav)
 } else {
   gc.dat <- get(paste(load(gc.sav)))
 }
 #####print(is(bigMat)); print(dim(bigMat))
 # get list of snps for each valid megabase windows that has at least 10 snps
 med.list <- get.ranges.for.median(snp.info,gc.dat,bigMat)
 rng <- med.list$gc.dat
 if(!(med.chunk.fn=="")) { 
   med.chunk.fn <- dir.fn.cmb(dir$qc.gc,med.chunk.fn) 
   if(!file.exists(med.chunk.fn)) { 
     gotMedian <- F 
   } else {
     ## if median matrix has been calculated previously load it and skip this step
     med.store <- get(paste(load(paste(med.chunk.fn))))
     if(!((ncol(med.store)==ncol(bigMat)) & (nrow(med.store)==nrow(rng))))
     { cat("Warning: can't used loaded object 'med.store'.",(dim(med.store)),"\n") ;
       cat("Object has the wrong dimensions versus other parameters\n") 
       gotMedian <- F
     } else {
       gotMedian <- T
     }
   }
 }
 if(!gotMedian | !exists("med.store")) {
   ## calculate the medians in each range for each sample, about 1hr per 100K-snps*5K-samples
   med.store <- do.median.for.ranges(med.list[[1]],bigMat,dir,cont=!med.list[[3]],
                              use.big.list=med.list[[2]], med.sav.nm=med.chunk.fn)
 } 
 r_GC <- S_WF <- S_GCWF <- numeric(ncol(bigMat))
 cat(" calculating Diskin constants...\n")
 for (j in 1:ncol(bigMat)) {
  r_GC[j] <- cor(rng$gc,med.store[,j],use="pairwise.complete.obs")
  S_WF[j] <- (1-2*(as.numeric(r_GC[j]<0))) *
    median(abs(med.store[,j]-median(med.store[,j])))
 }
 S_GCWF <- S_WF * abs(r_GC)
 dat.out <- cbind(r_GC,S_WF,S_GCWF)
 rownames(dat.out) <- colnames(bigMat)
 colnames(dat.out) <- c("r_GC","S_WF","S_GCWF")
 return(dat.out)
}


do.median.for.ranges <- function(ranges.list,bigMat,dir,cont=T,
  use.big.list=NULL, med.sav.nm="")
{
  # the diskin method to calculate GC requires sample-wise medians for snps
  # in each 1MB window of the genome
  dir <- validate.dir.for(dir,c("big"),warn=F)
  if((is(ranges.list)[1]!="list") | (is(bigMat)[1]!="big.matrix") )
  { 
    if((is(bigMat)[1]=="matrix")) {
      cat(" for GC calc, bigMat was a matrix, attempting to convert to big.matrix...")
      bigMat <- as.big.matrix(bigMat)
      cat("ok\n")
    } else {
      stop("Error: parameter 'ranges.list' should be a list, and 'bigMat' should be a big.matrix") 
    }
  }
  num.ranges <- length(ranges.list)
  ncb <- ncol(bigMat)
  cat(" doing median calculation\n")
  med.store <- matrix(nrow=num.ranges,ncol=ncb)
  jj <- proc.time()
  for (cc in 1:num.ranges) {
    if(cont) {
      # faster method (when each range is continuous)
      r1 <- ranges.list[[cc]][1]; r2 <- ranges.list[[cc]][2]
      if(cc %in% use.big.list) {
    	# this provides large speed increase when chunks have >200 snps
        LRR.chunk <- sub.big.matrix(bigMat, firstRow = r1, lastRow = r2, backingpath=dir$big)
      } else {
        LRR.chunk <- bigMat[c(r1:r2),c(1:ncb)]
      }
    } else {
      # method if intervals can be discontinuous [e.g, if dataset not sorted by chr,pos] #
      LRR.chunk <- bigMat[ranges.list[[cc]],c(1:ncb)]
    }
    med.store[cc,] <- apply(LRR.chunk,2,median,na.rm=T)
    # add a tracker, saving backups as we go
    loop.tracker(cc, num.ranges, jj, mode="time", sav.obj=med.store, sav.fn=med.sav.nm)
  }
  kk <- proc.time()
  cat(" took",round((kk-jj)[3]/3600,2),"hours\n")
  ## save median matrix to pick up later
  if(!(med.sav.nm=="")) {
    save(med.store,file=med.sav.nm) }
  return(med.store)
}


get.ranges.for.median <- function(snp.info,gc.dat,bigMat,min.size=10)
{
 # get list of position ranges for snps in each 1MB window of the genome,
 # ignoring any window with less than 'min.size' snps
 must.use.package(c("genoset","BiocGenerics"),T); must.use.package("bigmemory")
 if((is(gc.dat)[1]!=is(snp.info)[1]) | (is(gc.dat)[1]!="RangedData") )
 { stop("Error: parameters snp.info and gc.dat should be of type 'RangedData'") }

 olp <- findOverlaps(snp.info,gc.dat)

 # get list of snp.info row positions within each gc.dat range
 mbGroups <- tapply(queryHits(olp),subjectHits(olp),c)
 # get number of snps in each of these
 mbGroupsC <- sapply(mbGroups,length)
 too.small <- which(mbGroupsC<min.size)
 if(any(too.small)) { mbGroups <- mbGroups[-too.small]; mbGroupsC <- mbGroupsC[-too.small] }
 if(length(mbGroups)<2) { stop("Error: insufficient snp density for GC calculations. Check data or turn off GC") }
 
 gc.dat.filt <- subsetByOverlaps(gc.dat,snp.info)
 gc.dat.filt <- gc.dat.filt[-too.small,]
 
 # Prepare list of ranges for median calcs (depends on whether ranges will be continuous)
 if( any(rownames(snp.info)!=rownames(bigMat)) )
 { 
  cat("Warning: snp.info object rownames don't match big.matrix rownames.")
  cat("suggest using vectorToBigMatrix.R script to reduce the datafile to only the ")
  cat("needed snp-set. Otherwise this routine will be greatly slowed by using name indexing.")
  # identify any ranges that will cause a large enough block of memory allocation
  # that should use big memory subset instead of standard R object for subset chunk
  cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
  max.gb <- 1.5  #maximum number of gigabytes of chunk before might cause probs
  bad.row.max <- round(max.gb*((cells.per.gb/as.double(ncol(bigMat)))))
  name.ind <- T
  # convert position info to snp labels
  mbGroupsN <- lapply(mbGroups,function(X,rn) { rn[X] },rn=rownames(snp.info))
  ranges.list <- lapply(mbGroups,function(X,rn) { narm(match(X,rn)) },rn=rownames(bigMat))
  use.big.chunk <- which(sapply(ranges.list,length)>bad.row.max)
  if(any(use.big.chunk)) 
  { 
   stop("Error: some ranges are too large for name indexing, regenerate datafile as suggested above")
  }
 } else {
  #snp.info rownames match bigmatrix snps exactly
  cat(" 100% match between snp.info annotation and big.matrix file\n")
  name.ind <- F
  ranges.list <- lapply(mbGroups,range)
  row.max <- 200  
  use.big.chunk <- which(sapply(ranges.list,diff)>row.max)
 }
 outlist <- list(ranges.list,use.big.chunk,name.ind,gc.dat.filt)
 names(outlist) <- c("ranges.list","big.list","name.ind","gc.dat") 
 return(outlist)
}

loop.tracker <- function(cc, max, st.time=NULL, mode=c("time","perc")[1], 
                     sav.obj=NULL, sav.fn=NA,freq=round(max/min(max,50)), 
                     sav.freq=round(max/(min(max,4))), unit=c("m","s","h")[1])
{
 # insert this function into any loop with a counter and it will provide
 # a progress tracker either as a percentage or estimate of time remaining
 ## cc is loop counter, max is the last value the loop should take
 ## st.time is the result of a call to proc.time().
 if(cc==1)
 {
 	if(mode=="perc" | is.null(st.time)) {
   scale <- "0         25%         50%         75%         100%\n"
 	 cat(scale)
 	} else {
 	 cat("Processing: time left (",unit,"):\n",sep="")
 	}
 }
 if (cc %% freq == 0) {
   if(mode=="perc" | is.null(st.time))
   {
  	 intv <- diff(round(seq(from=1,to=51,length.out=(max/freq))))[cc %/% freq]
   	 if(!is.na(intv)) { if(intv>0) { cat(rep(".",intv),sep="") } }
   } else {
     time.now <- proc.time()[3]-st.time[3]; time.per <- time.now/cc
     tm.u <- switch(unit,m=60,s=1,h=3600)
     to.go <- round(((max-cc)*time.per/tm.u))
     cat(to.go,"..",sep="") 
   }
   if((cc+freq)>max) { cat("\n") }
   ## save as we go - in case of crash
   if (cc %% sav.freq == 0)
   {
   	if(!is.null(sav.obj) & !is.na(sav.fn) & ((max-cc)>5) ) {
       save(sav.obj,file=sav.fn)
     }
   }
 }
}


remove.qc.failers <- function(tab,sample.info,fail.col="QCfail") 
{
  # use column 'QCfail' in sample.info to remove failing ids from any frame
  # sample.info could be any dataframe where rownames are sample ids.
  dat <- get.frame.multi.type(tab) # allows flexible input, forces dataframe result
  if(!fail.col %in% colnames(sample.info)) 
    { stop(paste("Error: sample.info object must have column '",fail.col,"'",sep="")) } 
  pass.ids <- rownames(sample.info)[sample.info[[fail.col]]==0]
  index <- find.id.col(dat,ids=pass.ids)$index
  return(dat[index,])
}



add.to.sample.info <- function(sample.info,more.info,to.add=NULL,overwrite=F)
{
  # adds columns 'to.add' or else all valid columns from table 'more.info' to
  # the object sample.info, for matching rownames. autodetects id column in 
  # more.info based on the rownames of sample.info. can optionally overwrite existing columns
  cat("\nAdding batch information to 'sample.info' table\n")
  if(is.null(colnames(more.info))) 
  { 
    fake.names <- T # ie, indicating that no colnames present so using fake ones
    #colnames(more.info) <- paste("col",1:ncol(more.info)) 
  } else {
    fake.names <- F
  }
  more.info <- get.frame.multi.type(more.info)
  mat.fil.list <- find.id.col(more.info,ids=rownames(sample.info))
  cat("",round(mat.fil.list$maxpc*100,1),"% of ids in batch lookup table matched sample.info\n")
  indx <- mat.fil.list$index

  if(length(narm(indx))>0) {
    whichoz <- c(1:ncol(more.info))
    # remove ID column and if overwrite=F, any columns already in sample.info
    if(!overwrite) {
     nn <- which(colnames(sample.info) %in% colnames(more.info))
    } else { nn <- NULL }
    mm <- (mat.fil.list$col); mm <- c(mm[mm!=0],nn[nn!=mm])
    if(length(mm)!=0) { whichoz <- whichoz[-mm] }
    # left with list of all valid columns in more.info
    to.try <- paste(colnames(more.info)[whichoz])
    to.use <- to.add[to.add %in% to.try] # select only requested ones
    if(length(to.use)==0) {
      to.use <- to.try # if no request colnames found, use all possible
    }
    ## decide what column names will be in sample.info
    if(fake.names & length(to.add)==1 & length(to.use)==1) {
      #only 1 var to be added, only 1 in dataframe but unnamed, so use to.add name
      to.name <- to.add
    } else {
      to.name <- to.use
    }
    for(jj in 1:length(to.use)) {
      if(to.use[jj] %in% colnames(more.info))
      {  
        sample.info[[to.name[jj]]][indx] <- more.info[[to.use[jj]]] 
      }
    }
    cat(" added columns",to.name,"[from columns",to.use,"in 'more.info'] to sample.info\n",sep=", ")
  } else {
    stop("Warning: no subject ids matched ids in plate lookup file\n")
  }
  return(sample.info)
}


find.id.col <- function(frame,ids,ret=c("col","maxpc","index","result"))
{
  # looks across each column of a dataframe 'frame' and selects the
  # column that has the most matches to the vector 'ids'
  # can return any/all of:
  #   the number of best match column, the match %, index numbers, matches
  opts <- c("col","maxpc","index","result")
  targ.cols <- ncol(frame)+1 # +1 if for looking at 'rownames'
  pcs <- numeric(targ.cols)
  num.ids <- length(ids)
  coln <- 1; best <- NA
  for (cc in 1:targ.cols)
  {
    if(cc==targ.cols)
    {
      if(!is.null(rownames(frame))) { posit <- match(ids,rownames(frame)) } else { break }
    } else {
      posit <- match(ids,frame[,cc])
    }
    pcs[cc] <- length(which(!is.na(posit)))/num.ids
    if(pcs[cc]>max(pcs[-cc])) { best <- posit ; coln <- cc}
    if(pcs[cc]==1) { break } # exit if found a perfect match
  }
  maxpc <- max(pcs)
  if(cc==targ.cols) {
    result <- rownames(frame)[best]
    coln <- 0
  } else {
    result <- frame[best,coln]
  }
  out <- list(coln,maxpc,best,result)
  names(out) <- opts
  for (cc in length(opts):1)
  {
    if (!(opts[cc] %in% ret)) { out[[cc]] <- NULL }
  }
  return(out)
}

get.frame.multi.type <- function(unknown.data,too.big=10^7)
{
  # returns a dataframe if 'unknown.data' can in anyway relate to such:
  # it can be:
  # - dataframe, matrix, big.matrix, sub.big.matrix, big.matrix.descriptor,
  # a bigmatrix description file, an RData file containing one of these objects,
  # the name of a text or RData file, a named vector (names become rownames),
  # or a list containing a matrix or dataframe. Using this in functions allows
  # flexibility in specification of a datasource
  uis <- is(unknown.data)
  if(uis[1]=="character" & length(unknown.data)==1)
  {
    # filename?
    if(file.exists(unknown.data)) { 
      out.data <- reader(probably.is)
    } else {
      stop("Error: argument seemed to be a filename (char, length 1) but file did not exist")
    }
  } else {
    if(uis[1] %in% c("matrix","data.frame")) {
      out.data <- unknown.data
    } else {
      if(uis[1]=="list") {
        types <- sapply(unknown.data,is)
        wty <- which(types %in% c("matrix","data.frame"))[1]
        if(length(wty)==1) {
          out.data <- unknown.data[[wty]]
        } else {
          stop("Error: object was list and no elements were a matrix or data.frame")
        }
      } else {
        if(is.null(dim(unknown.data)) & !is.null(names(unknown.data)))
        {
          cat(" converting named vector into dataframe\n")
          out.data <- as.matrix(unknown.data)
        } else {
          if(length(grep("big.matrix",uis))>0) {
            if(nrow(unknown.data)*ncol(unknown.data) <= too.big) {
              cat(" converting big.matrix object into a dataframe\n")
              out.data <- as.matrix(unknown.data) 
            } else {
              stop(paste("Error: big matrix object was too big to convert to a dataframe [>",
                         too.big,"cells]\n"))
            }
          } else {
            stop(paste("Error: object type",uis[1],"could not be converted into a dataframe"))
          }
        }
      }
    }
  }
  return(as.data.frame(out.data))
}


column.salvage <- function(frame,desired,testfor) 
{
  ## attempt to find predictable misnaming of column (contents of 'testfor')
  # and change name of the first detected misname in the 'testfor' list to 'desired'
  # e.g, want desired column 'GRP' and look for misnamings: 'group', 'Grp', 'grp', etc
  if(!desired %in% colnames(frame)) {
    if(any(testfor %in% colnames(frame))) {
      colnames(frame)[(narm(match(testfor,colnames(frame)))[1])] <- desired
    } else {
      stop("Error: couldn't find columns",paste(testfor,collapse=", "),"to change to '",desired,"' in frame")
    }
  }
  return(frame)
}


validate.samp.info <- function(sample.info,dir,QC.update=T)
{
  ## make sure current sample.info object conforms to the expected structure ##
  # update QC failure flags based on sample exclusion directory
  # update 'grp' membership based on file.spec.txt column 'GRP'
  dir <- validate.dir.for(dir,c("lrr.dat","ids","ano"))
  if(!is.data.frame(sample.info)) {
    cat(" attempting to coerce sample info to data.frame\n")
    sample.info <- as.data.frame(sample.info)
  }
  if(!( "grp" %in% colnames(sample.info)))
  {
    stop("Error: sample.info must contain a column named 'grp'")
  } else {
    if(is.null(rownames(sample.info)))
    { stop("Error: sample.info must have subject Ids as rownames (unique)") }
  }
  sub.id.list <- get.subIDs(dir,"all")
  ID.lists <- sub.id.list$original
  file.info <- get.file.specs(dir)
  file.info <- column.salvage(file.info,"GRP",c("grp","Grp","Group","group","GROUP"))
  
  if(length(unique(narm(file.info$GRP)))!=length(unique(narm(sample.info$grp))))
  {
    cat("Warning: number of 'grp's in sample info clashes with file.spec.txt in:\n ",dir$ids,"\n")
    cat(" updating sample.info$grp to match 'GRP' coding in file.spec.txt:\n")
    for (jj in 1:length(ID.lists))
    { 
      fl.lst <- paste(file.info[,1])
      nfn <- (sub.id.list$files)[jj]
      next.grp.num <- match(gsub(".ids","",nfn),fl.lst)
      if(is.na(next.grp.num)) {
        next.grp.num <- match(nfn,fl.lst)
        if(is.na(next.grp.num)) { stop(paste("Error: couldn't find",nfn,"in file.spec.txt")) }
      }
      sample.info$grp[rownames(sample.info) %in% ID.lists[[jj]]] <- file.info$GRP[next.grp.num]
    }
  }
  if(!( "QCfail" %in% colnames(sample.info)))
  {
    cat(" adding qc failure flag column to sample.info\n")
    sample.info[["QCfail"]] <- rep(0,times=nrow(sample.info))
    QC.update <- T
  }
  if(QC.update) {
    failed.samps <- get.all.samp.fails(dir)
    cat(" updating QC failure flags to reflect lists in '/SAMPLE_EXCLUDE/'\n")
    sample.info$QCfail[!rownames(sample.info) %in% failed.samps] <- 0
    sample.info$QCfail[rownames(sample.info) %in% failed.samps] <- 1
  }
  return(sample.info)
}


get.subIDs <- function(dir,ret=c("all","lists","combined","files","groups"),verbose=F,combine.grp=F)
{
  # if separate id lists for each cohort are in 'dir' then create list of these
  dir <- validate.dir.for(dir,c("ids","col","ano"))
  if(verbose) { cat("\nSearching for subject IDs\n") }
  file.info <- get.file.specs(dir)
  id.fnz <- paste(file.info[,1],".ids",sep="")
  grpz <- unique(file.info$GRP)
  num.filz <- length(id.fnz); num.grpz <- length(grpz)
  subIDs.list <- vector("list", num.filz)
  names(subIDs.list) <- id.fnz
  path.fn <- dir.fn.cmb(dir$ids,id.fnz,must.exist=T)  # sample id file name
 ## print(path.fn)
  if(verbose) { cat(paste(" reading IDs from: ",id.fnz,"\n",sep="")) }
  ID.list <- lapply(path.fn,readLines)  #readLines(sample.fn)
  cmb.ID.list <- paste(do.call("c",ID.list))
  if(tolower(ret[1]) %in% c("all","groups","group","group.files","group.file","groupfiles","groupfile")) {
    next.grp.files <- next.grp.ids <- vector("list",length(grpz))
    names(next.grp.files) <- names(next.grp.ids) <- grpz
    dd <- 1
    for (cc in grpz) {
      next.fn <- dir.fn.cmb(dir$ids,id.fnz[which(file.info$GRP %in% cc)])
      next.idz <- lapply(next.fn,readLines)
      next.grp.files[[dd]] <- next.fn
      if(combine.grp) {
        next.grp.ids[[dd]] <- paste(do.call("c",next.idz)) 
      } else {
        next.grp.ids[[dd]] <- next.idz
      }
      dd <- dd + 1
    }
  }
  out.obj <- NULL
  if (tolower(ret[1]) %in% c("group","groups")) { out.obj <- next.grp.ids }
  if (tolower(ret[1]) %in% c("group.file","group.files","groupfile","groupfiles")) 
    { out.obj <- next.grp.files }
  if (tolower(ret[1]) %in% c("list","lists")) { out.obj <- ID.list }
  if (tolower(ret[1]) %in% c("combined","combine")) { out.obj <- cmb.ID.list }
  if (tolower(ret[1]) %in% c("file","files")) { out.obj <- id.fnz }
  if (tolower(ret[1]) %in% c("","all","every")) {
    out.obj <- list(cmb.ID.list,next.grp.ids,ID.list,id.fnz,next.grp.files)
    names(out.obj) <- c("combined","group","original","files","group.files")
  }
  return(out.obj)
}


get.subIDs.old <- function(dir,id.fnz=NULL,dir.spec=NULL)
{
  # if separate id lists for each cohort are in 'dir' then create list of these
  dir <- validate.dir.for(dir,c("ids","col","ano"))
  cat("\nSearching for subject IDs\n")
  if(is.null(id.fnz))
  {
    # no names given so try to get the id file names using file.spec.txt
    if(is.null(dir.spec)) { dir.spec <- dir }
    id.fnz <- paste(get.file.specs(dir.spec)[,1],".ids",sep="")
  }
  num.filz <- length(id.fnz)
  subIDs.list <- vector("list", num.filz)
  names(subIDs.list) <- id.fnz
  for (cc in 1:num.filz) {
    if(id.fnz[cc] %in% c(list.files(dir$ids),list.files(dir$ano)))
    {
      cat(paste(" reading IDs from: ",id.fnz[cc],"\n",sep=""))
      if(id.fnz[cc] %in% c(list.files(dir$ids))) {
        ids.fn <- paste(dir$ids,id.fnz[cc],sep="") 
      } else {
        ids.fn <- paste(dir$ano,id.fnz[cc],sep="") 
      }
      subIDs.list[[cc]] <- readLines(ids.fn)
    } else {
      subIDs.list[[cc]] <- NULL
      stop(paste("Error: subject ids file",ids.fn,"required for group",cc))
    }
  }
  return(subIDs.list)
}

get.file.specs <- function(dir,fn="file.spec.txt")
{ 
  # check for and read specifications file 'file.spec.txt' from directory dir$lrr.dat
  # input as 'dir' list or otherwise will try to convert to this format from character()
  dir <- validate.dir.for(dir,c("lrr.dat","ano","raw"),warn=T)
  if (fn %in% c(list.files(dir$lrr.dat),list.files(dir$ano)))
  {
    if (fn %in% c(list.files(dir$lrr.dat))) {
      file.info <- read.table(paste(dir$lrr.dat,fn,sep=""),header=T) 
    } else {
      file.info <- read.table(paste(dir$ano,fn,sep=""),header=T)
    }
    alph.ord <- order(file.info[,1])
    if (any(alph.ord!=1:nrow(file.info)))
    { cat("rearranging ",fn," into alphabetical order by filename")
      file.info <- file.info[alph.ord,]  } #alphabetical order  
  } else {
    linez <- character()
    linez[1] <- paste("The location",dir$lrr.dat)
    linez[2] <- paste("should contain a file '",fn,"' which is tab delimited file with columns: FILE GRP TYPE SAMP SNP A1 A2 LRR BAF")
    linez[3] <- "FILE=file name (no path), GRP = 1:ngrps, TYPE=txt or gzip, then remainder are the column number "
    linez[4] <- "in each raw file that contain data for SAMP(les) SNP(s) (alleles,A1 A2) and LRR, BAF, respectively."
    linez[5] <- "\nThis file should have been used by 'bash_getAColumn.sh to extract the data initially."
    if ("raw" %in% names(dir)) {
      linez[6] <- paste("File locations should be relative to the path:",dir$raw) }
    cat(c(linez))
    stop("Error: critical file missing!")
  }
  return(file.info)
}

get.file.lens <- function(dir,fn="file.lengths.txt")
{
  ## read in file lengths from text file 'file.lengths.txt' (generated during bash import)
  dir <- validate.dir.for(dir,c("lrr.dat","ano"),warn=T)
  if (fn %in% c(list.files(dir$lrr.dat),list.files(dir$ano)))
  {
    if (fn %in% c(list.files(dir$lrr.dat))) {
      lens.fn <- paste(dir$lrr.dat,fn,sep="")
    } else {
      lens.fn <- paste(dir$ano,fn,sep="")
    }
    cat(" reading file lengths of raw datafiles from",fn,"\n")
    len.tab <- read.table(lens.fn)
    alph.ord <- order(len.tab[,2])
    if (any(alph.ord!=1:nrow(len.tab)))
    { cat(" rearranging",fn,"into alphabetical order by filename")
      len.tab <- len.tab[alph.ord,]  }
  } else {
    linez <- character()
    linez[1] <- paste("The location",dir$lrr.dat,"")
    linez[2] <- paste("should also contain a file '",fn,"' specifying the length of each data file")
    linez[3] <- " (alphabetical order); generate using a command like this in the terminal:\n\n"
    linez[4] <- paste("    wc -l *myfiles*.dat | sed '$d' | sort -b -k2 > ",fn,"\n")
    linez[5] <- paste("\nThe file '",fn,"' should have lengths for the following files:\n")
    fnz <- get.file.specs(dir)[,1]
    cat(c(linez,paste("\n",fnz)))
    stop("Error: critical file missing!")
  }
  len.tab <- as.data.frame(len.tab)
  colnames(len.tab) <- c("length","file.name")
  return(len.tab)
}


validate.dir.for <- function(dir,lst,warn=F) {
  # in case the 'dir' input list object is not the standardised list form, convert
  # allows flexible use of list or regular directory specifications in plumbCNV functions
  if(!is.list(dir)) {
    if(warn) { cat(lst[cc],"'dir' object wasn't a list\n")}
    dir <- as.list(dir); names(dir)[1:length(dir)] <- lst[1:length(dir)] 
  }
  for (cc in 1:length(lst)) {
    if(!lst[cc] %in% names(dir)) { 
      dir[[paste(lst[cc])]] <- "" ;
      if(warn) { stop(paste("dir$",lst[cc]," was empty.. set to current\n",sep="")) } 
    }
  }
  return(dir)
}


get.lrr.files <- function(dir,grps=NULL,...) {
  # wrapper for get.lrr.dat.file.names
  all.fn <- get.lrr.dat.file.names(dir,...)
  file.info <- get.file.specs(dir)
  if(!is.null(grps) & length(all.fn)>1) {
    select <- which(paste(file.info$GRP) %in% paste(grps))
    if(length(select)==0) {
      stop(paste("No LRR files found for group(s):",paste(grps,collapse=",")))
    } 
    return(all.fn[select])
  } else {
    return(all.fn)
  }
}


get.lrr.dat.file.names <- function(dir,suffix="LRR.dat")
{
  ## tries really hard to find the right datafiles to read in for LRR raw data
  dir <- validate.dir.for(dir,"col")
  all.in.dir <- list.files(dir$col)  
  file.info <- get.file.specs(dir)
  exp.fns <- paste(file.info[,1],suffix,sep=".")
  cat(" searching for raw datafiles in processed long format (e.g, converted from GenomeStudio)\n")
  if(nrow(file.info)==1) {
    # simplest case - only 1 file
    if(file.exists(dir.fn.cmb(dir$col,suffix))) {
      cat(" found combined data file:",suffix,"\n")
      names.out <- suffix
    } else {
      if(file.exists(dir.fn.cmb(dir$col,exp.fns))) {
        names.out <- exp.fns
        cat(" found sole data file:",exp.fns,"\n")
      } else {
        any.dats <- grep(".DAT",toupper(all.in.dir),fixed=T)
        if(length(any.dats)>0) {
          cat("Warning: was expecting file named ",exp.fns,"or",suffix)
          cat(", was not in",dir$col,", but there were",length(any.dats),"*.dat files\n")
          cat("Will try to proceed with the first file:",all.in.dir[any.dats][1],"\n")
          cat("If this is not the LRR raw datafile then please CTRL-C and review file.spec.txt\n")
          names.out <- dir.fn.cmb(dir$col,all.in.dir[any.dats][1])
        } else {
          if(length(all.in.dir)==1) {
            cat("Warning: was expecting file named ",exp.fns,"or",suffix)
            cat(", these were not in ",dir$col,", but there was only 1 file present\n")
            cat("Will try to proceed with file:",all.in.dir[1],"\n")
            cat("If this is not the LRR raw datafile then please CTRL-C and review file.spec.txt\n")
            names.out <- dir.fn.cmb(dir$col,all.in.dir[1])
          } else {
            cat("Error: found multiple files in",col$dir,", and could not tell which if any were the raw LRR datafile\n")
            stop("Please review file.spec.txt")
          }
        }
      }
    }
    return(names.out)
  } 
  
  if(length(unique(file.info$GRP))<=1) {
    # only 1 grp, so hopefully a combined datafile already, if not need to manage that
    if(file.exists(dir.fn.cmb(dir$col,suffix))) {
      cat(" found combined data file:",suffix,"\n")
      names.out <- suffix
      return(names.out)
    } else {
      any.dats <- grep(".DAT",toupper(all.in.dir),fixed=T)
      if(length(any.dats)==1) {
        cat("Warning: was expecting file named",suffix,". This was not in")
        cat(dir$col,", but there was 1 dat file\n",
            all.in.dir[any.dats][1],"\n which will be used. If this is not")
        cat(" the combined LRR raw datafile then please CTRL-C and review file.spec.txt\n")
        names.out <- dir.fn.cmb(dir$col,all.in.dir[any.dats][1])
        return(names.out)
      } else {
        # move onto next part, hoping now there are multiple files 
        # which are meant to be combined into 1 group
      }
    }
  }
  
  if(nrow(file.info)>1) {
    # if we've made it to here there must be more than one cohort in the files or
    #  there are multiple files that need to be combined into one group   
    ##print(file.exists(dir.fn.cmb(dir$col,file.info[,1])))
    ##print(dir.fn.cmb(dir$col,file.info[,1]))    
    exp.fnp <- dir.fn.cmb(dir$col,exp.fns)
    if (any(file.exists(exp.fnp)))
    {
      if (all(file.exists(exp.fnp)))
      {
        cat(" all expected data files found in",dir$col,"\n")
        names.out <- exp.fnp
      } else {
        which.left <- length(which(file.exists(exp.fnp)))
        cat("Error: only",length(which.left),"/",nrow(file.info),
            "of the expected datafiles were found in",dir$col,"\n")
        stop("please review file.spec.txt")
      }
      return(names.out)
    } else {
      cat("Error: none of the expected data files found in",dir$col,"\n")
      cat("was looking for:\n",paste(exp.fnp,collapse="\n"),"\n")
      stop("Please review file.spec.txt")
    }
  } else {
    cat("No LRR datafiles found with suffix",suffix,"in ",dir$col,"\n")
    stop("Please review file.spec.txt")
  }
}


time.est.snps <- function(tab,lines.per.min=12*10^6)
{
  times <- round(tab[,1]/lines.per.min,1)
  times <- c(times,sum(times))
  files <- c(tab[,2],"All files")
  result <- cbind(times,files)   
  colnames(result) <- c("~ mins to process","Raw file name")
  return(result)
}


import.snp.matrix.list <- function(snp.list,dir,data=NULL,samples=NULL,field.list=NULL,HD=F,multi=F)                       
{
 # Reads in the SNP data for the whole dataset, with a separate snpMatrix entry per cohort
 # setting HD=T uses a lot less RAM, can be used for very large datasets (but slow)
 # option multi allows use of parallel processing, but the main limited is hard disk
 # speed, so for now this is likely to decrease performance if anything.
 num.genes <- length(snp.list)
 # generate locations of the data files:
 dir <- validate.dir.for(dir,c("raw","col","cr","lrr.dat","ano"),warn=T)
 # dir$raw should be the location of the raw data files specified by
 ### locations of files specifying raw data (separate file for zipped/unzipped data)
 if(is.null(data)) {
   file.info <- get.file.specs(dir)
   fns <- file.info[,1]
   cat("",length(fns),"data files found\n")
   data <- dir.fn.cmb(dir$raw,fns,must.exist=T)
 }
 if(is.null(samples)) 
 { 
   subs.list <- get.subIDs(dir,"list") 
 } else {
   if(length(samples)>10) {
     subs.list <- samples
   } else {
     subs.list <- dir.fn.cmb(dir$ano,samples,must.exist=T)
   }
 } 
 num.filz <- length(data)
 # define names and column locations of fields in the data file
 if(is.null(field.list)) {
   # no field.list object passed in so try to read from file.spec.txt
   file.info <- get.file.specs(dir)
   field.list <- list() #defined using columns of 'file.spec.txt'
   for (cc in 1:num.filz) {
     field.list[[cc]] <- c(sample = file.info$SAMP[cc], snp = file.info$SNP[cc], 
                          allele1 = file.info$A1[cc], allele2 = file.info$A2[cc]) }
 } else {
   cat("Error: import requires setup file 'file.spec.txt or else 'field.list': a list for each file,\n")
   cat("with elements: sample,snp,allele1,allele2; giving columns numbers for each in the datafile\n")
   stop()
 }
 ########################
 len.tab <- get.file.lens(dir)
 cat("\nTime estimates for importing SNP data:\n")
 print(time.est.snps(len.tab),quote=F)
 file.lens.nh <- (len.tab[,1])  #lengths only 
 # get header lengths (may not be uniform)
 header.lens <- get.hdr.lens(data)
 file.lens.sub <- (file.lens.nh)/num.genes
 cat("Files contain:",paste(file.lens.sub,collapse=","),"samples, respectively\n")
 # go to /CALLRATE directory
 setwd(dir$cr)
 ######
 #############id.file.names <- paste(data,".ids",sep="")
 if(!is.list(subs.list)) { subs.list <- list(subs.list) } # force list
 if(length(subs.list)!=length(data)) { stop("Must be same number of datafiles as lists of IDs") }
 # Read each set of SNPs into a SnpMatrix object stored in a list element
 snpMatLst <- list()
 #if HD=T only store list of locations of saved objects (in case of memory limitation)i
 cat("\n")
 #else list contains the objects themselves
 if(!multi) {
   options(warn = -1)
   for (tt in 1:num.filz) {
     kk <- proc.time()
     snpMat <- read.snps.long(files = data[tt], sample.id = subs.list[[tt]],    
                              snp.id = snp.list, diploid = NULL, 
                              fields = field.list[[tt]], 
                              codes = "nucleotide", sep = "\t", comment = "#", 
                              skip = header.lens[tt], simplify = c(FALSE,FALSE),
                              verbose = T, in.order = TRUE, every = num.genes)
     if(HD) {
       ofn <- paste("snpMat",tt,".RData",sep="")
       save(snpMat,file=ofn)
       snpMatLst[[tt]] <- ofn
     } else {
       snpMatLst[[tt]] <- snpMat ; snpMat <- NULL
     }
     jj <- proc.time()
     cat(paste("file",tt,"took",round((jj[3]-kk[3])/60),"minutes\n"))
   }
   options(warn = 0)
 } else {
   for (tt in 1:num.filz) {
     snpMatLst[[tt]] <- parallel(read.snps.long(files = data[tt], sample.id = subs.list[[tt]],    
                                                snp.id = snp.list, diploid = NULL, 
                                                fields = field.list[[tt]], 
                                                codes = "nucleotide", sep = "\t", comment = "#", 
                                                skip = header.lens[tt], simplify = c(FALSE,FALSE),
                                                verbose = T, in.order = TRUE, every = num.genes))
   }	
   snpMatLst <- collect(snpMatLst)
 }
 return(snpMatLst)
}


draw.density.plots <- function(fn,sample.info,snp.info,samp.result=NULL,snp.result=NULL,
                            callrate.samp.thr=.95, callrate.snp.thr=.95, anot=T, par=NULL) 
{
 # draw density plots for SNP/sample-wise callrate evaluation
 cat(" generating density plots\n")
 # set text parameters for samples (1) and snps (2); i,j = both
 par.list <- list(a1=1,b1=140,x1=28,y1=14,a2=27,b2=144,x2=32,y2=14,i=0.17,j=0.3)
 # any parameters entered as a named list in par will override these defaults^
 if(!is.null(par)) { par.list <- update.list.with.list(par.list,par) }
 if(!is.null(sample.info$call.rate) & !is.null(snp.info$call.rate))
 {
   pdf(fn)
   par(mfrow=c(2,1))
   #samples
   plot(density(sample.info$call.rate),xlim=c(callrate.samp.thr,1),main="A. Sample-wise call rate distribution")
   if(anot & length(dim(samp.result))>0) {
     if(all(dim(samp.result)==c(8,3))) {
       attach(par.list)
       tt <- paste(as.matrix(samp.result,byrow=T))
       tt <- c(colnames(samp.result),tt)
       tt <- tt[c(1,4:11,2,12:19,3,20:27)]
       spn <- 1-callrate.samp.thr
       mdpt <- callrate.samp.thr+ spn/2
       if(median(sample.info$call.rate,na.rm=T)<mdpt) {
         ofs <- callrate.samp.thr + spn*.6 
       } else {
         ofs <- callrate.samp.thr
       }
       text(rep(c(0,i,j)*spn,each=9)+ofs,a1*rep(seq(b1,x1,by=-y1),times=3),labels=tt,cex=0.6,pos=4)
       detach(par.list)
     }
   }
   #snps
   plot(density(snp.info$call.rate),xlim=c(callrate.snp.thr,1),main="B. SNP-wise call rate distribution")
   if(anot & length(dim(snp.result))>0) {
     if(all(dim(snp.result)==c(8,3))) {  
       attach(par.list)
       tt<-paste(as.matrix(snp.result,byrow=T))
       tt <- c(colnames(snp.result),tt)
       tt <- tt[c(1,4:11,2,12:19,3,20:27)]
       spn <- 1-callrate.snp.thr
       mdpt <- callrate.snp.thr+ spn/2
       if(median(snp.info$call.rate,na.rm=T)<mdpt) {
         ofs <- callrate.snp.thr + spn*.6 
       } else {
         ofs <- callrate.snp.thr
       }
       text(rep(c(0,i,j)*spn,each=9)+ofs,a2*rep(seq(b2,x2,by=-y2),times=3),labels=tt,cex=0.6,pos=4)
       detach(par.list)
     }
   }
   dev.off()
   cat(paste("wrote density figure to file:",fn,"\n"))
 } else {
   cat("Warning: invalid snp.info or sample.info objects, could not produce graph\n")
 }
}


call.rate.summary <- function(s.info,snp=(is(s.info)[1]=="RangedData"),print=T) 
{
 # generate summary of callrate performance for snps/samples
 cat("\nGenerating report")
 if(snp) { cat(" (snp)\n") } else { cat(" (sample)\n") }
 # Report
 thresh.lt <- c(.5,.75,.9,.95)
 thresh.gt <- c(.95,.97,.99,.999)
 counts.lt <- integer(length(thresh.lt))
 counts.gt <- integer(length(thresh.gt))
 # samples or snps (autodetected)
 for (cc in 1:length(thresh.lt))
 { counts.lt[cc] <- length(s.info$call.rate[s.info$call.rate<thresh.lt[cc]]) }
 for (cc in 1:length(thresh.gt))
 { counts.gt[cc] <- length(s.info$call.rate[s.info$call.rate>=thresh.gt[cc]]) }
 denom <- nrow(s.info)  
 counts <- c(counts.lt, counts.gt)
 counts.pc <- round(counts/denom*100,2)
 rn.resul <- c(paste("callrate <",thresh.lt),paste("callrate >=",thresh.gt))
 if(snp) {
   result <- data.frame(CallRate=rn.resul,SNPs=counts,Percent=paste(counts.pc,"%",sep=""))
 } else {
   result <- data.frame(CallRate=rn.resul,Samples=counts,Percent=paste(counts.pc,"%",sep=""))
 }
 if(print) { print(result); cat("\n") }
 return(result)
}

update.list.with.list <- function(orig.list, add.list)
{
 # idea is start with a full list, and if any elements with
 # the same name are in add.list, and they are the same type,
 # then the original list entry will be updated with the add.list
 # entry. the idea is to facilitate easy control of a large number
 # of parameters (which have defaults) while only having to specify
 # those to change, not all of them
 if(is(add.list)[1]=="list" & is(orig.list)[1]=="list")
 {
   for (cc in 1:length(add.list)) {
     nxt.nm <- names(add.list)[cc]
     if(nxt.nm %in% names(orig.list))
     { 
       if(is(add.list[[cc]])[1]==is(orig.list[[nxt.nm]])[1]) { 
         orig.list[[nxt.nm]] <- add.list[[cc]] 
       } 
     }
   }
 }
 return(orig.list)
}

doSampQC <- function(dir, plink=T,callrate.samp.thr=.95, snpMatLst=NULL, subIDs.actual) 
{
 # Do sample callrate quality control on a snpMatrix object (or just use plink files
 # if the QC has already been performed in plink)
 ## get sample call rates ##
 dir <- validate.dir.for(dir,c("ano","cr.plk"),warn=F)
 if(is(snpMatLst[[1]])[1]=="SnpMatrix") {
   group.nums <- rep(c(1:length(snpMatLst)),sapply(snpMatLst,dim)[1,])
 } else {
   if(plink) {
     group.nums <- rep(1,times=length(subIDs.actual))
   } else {
     lsp <- get.snpMat.spec(snpMatLst)
     group.nums <- rep(c(1:length(snpMatLst)),times=lsp[1,])
   }
 }
 sample.info <- data.frame(grp=group.nums,row.names=subIDs.actual,QCfail=rep(0,times=length(group.nums)))
 if (plink | is.null(snpMatLst)){
   if(!plink) { cat(" snpMatLst is NULL but plink=FALSE, trying to see whether plink QC data is available\n") }
   if (!"snpdataout.irem" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink samples-removed-by-QC file: snpdataout.irem, in",dir$cr.plk)) }
   if (!"snpdataout.imiss" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink samples-QC output file: snpdataout.imiss, in",dir$cr.plk)) }
   call.rate.excl.grp <- paste(read.table(paste(dir$cr.plk,"snpdataout.irem",sep=""),header=F)[,2])
   sample.info[call.rate.excl.grp,"QCfail"] <- 1
   cat(paste(" plink sample QC removed",length(call.rate.excl.grp),"samples\n"))
   imisser <- read.table(paste(dir$cr.plk,"snpdataout.imiss",sep=""),header=T)
   sample.info[["call.rate"]] <- 1-(imisser[match(rownames(sample.info),imisser$IID),"F_MISS"])
   sample.info$call.rate[is.na(sample.info$call.rate)] <- 0
   more.fail <- which(sample.info$call.rate[!is.na(sample.info$call.rate)]<callrate.samp.thr)
   if(length(more.fail)>0)
   {
     cat("Warning: threshold entered (",callrate.samp.thr,") may differ from plink threshold (MIND=...) ")
     cat("because",length(more.fail),"additional samples failed the call rate threshold. Removing additional samples.\n")
     cat("Note that this will affect the SNP-qc as some snps with missing calls for these ")
     cat("bad samples might have failed unnecessarily.\n")
     sample.info$QCfail[more.fail] <- 1
     call.rate.excl.grp <- c(call.rate.excl.grp,rownames(sample.info)[more.fail])
   }
 } else {
   sample.qc <- list.rowsummary(snpMatLst)
   sample.info[["call.rate"]] <- sample.qc[rownames(sample.info),"Call.rate"]
   sample.info$call.rate[is.na(sample.info$call.rate)] <- 0
   call.rate.excl.grp <- paste(rownames(sample.info)[sample.info$call.rate<callrate.samp.thr])
   call.rate.excl.grp <- narm(call.rate.excl.grp)
   if(length(call.rate.excl.grp)>0) {
     sample.info[call.rate.excl.grp,"QCfail"] <- 1
   }
   #to.cut <- match(call.rate.excl.grp,rownames(sample.info))
   #if(length(to.cut)>0) { 
   #  sample.info <- sample.info[-to.cut,] 
   #}
 }
 out.list <- list(call.rate.excl.grp,sample.info)
 names(out.list) <- c("CR.EXCL","SAMPLE.INFO")
 return(out.list)
}



doSnpQC <- function(dir, plink=T,
                   callrate.snp.thr=.95, hwe.thr=(10^-7),
                   snpMatLst=NULL, subIDs.actual, snp.info) 
{
 # Do SNP callrate/HWE quality control on a snpMatrix object (or just use plink files
 # if the QC has already been performed in plink)
 dir <- validate.dir.for(dir,c("ano","cr.plk"),warn=F)
 ## get snp call rates ##
 if(is(snpMatLst[[1]])[1]=="SnpMatrix") { HD <- F } else { HD <- T }
 if (plink | is.null(snpMatLst)){
   if (!"snpdataout.lmiss" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink snp-QC output file: snpdataout.lmiss, in",dir$cr.plk)) }
   if (!"snpdataout.hwe" %in% list.files(dir$cr.plk)) {
     stop(paste("Error: expecting plink snp-HWE output file: snpdataout.hwe, in",dir$cr.plk)) }
   lmisser <- read.table(paste(dir$cr.plk,"snpdataout.lmiss",sep=""),header=T)
   snp.info[["call.rate"]] <- 1-(lmisser[match(rownames(snp.info),lmisser$SNP),"F_MISS"])
   hweer <- read.table(paste(dir$cr.plk,"snpdataout.hwe",sep=""),header=T)
   snp.info[["P.hwe"]] <- hweer[match(rownames(snp.info),hweer$SNP),"P"]
 } else {
   if(length(snpMatLst)!=22) {
     snpMatLst <- convert.smp.to.chr22(snpMatLst,snpInfo=snp.info,HD=F,loc=dir$cr.plk) }
   snp.qc <- list.colsummary(snpMatLst)
   snp.info[["call.rate"]] <- snp.qc[rownames(snp.info),"Call.rate"]
   zz <- snp.qc[rownames(snp.info),"z.HWE"]
   pp <- pmin(pnorm(as.numeric(zz)),1-pnorm(as.numeric(zz)))*2  #two-tailed p value
   snp.info[["P.hwe"]] <- pp
 }
 ## now remove samples failing HWE + callrate 95 regardless of source!
 snp.info[["QCfail"]] <- rep(0,times=nrow(snp.info))
 snp.info$call.rate[is.na(snp.info$call.rate)] <- 0
 snp.info$P.hwe[is.na(snp.info$P.hwe)] <- 0
 cr.cond <- snp.info$call.rate < callrate.snp.thr
 call.rate.excl.snps <- row.names(snp.info)[cr.cond]
 HWE.cond <- snp.info$P.hwe < hwe.thr
 HWE.exclude.snps <- row.names(snp.info)[HWE.cond]
 to.remove <- unique(c(call.rate.excl.snps,HWE.exclude.snps))
 if(length(to.remove)>0)
 {
   snp.info[["QCfail"]][match(to.remove, rownames(snp.info))] <- 1
 }
 out.list <- list(call.rate.excl.snps,HWE.exclude.snps,snp.info,cr.cond,HWE.cond)
 names(out.list) <- c("CR.EXCL","HWE.EXCL","SNP.INFO","CR.cond","HWE.cond")
 return(out.list)
}



get.hdr.lens <- function(fnz,max.feas.hdr=500,firstrowhdr=T,sep.chr="\t") 
{
 # get the lengths of headers in genome studio files
 # assuming the header must be less than 'feas.hdr'/500 lines long.
 # starts from line 500 counts the number of tabs/separators and 
 # looks for when the number first changes. if it doesn't assumes 1st row only is header
 # if firstrowhdr = T, else assume no headers at all
 header.lens <- 1
 num.tabs <- function(str,ch="\t") {  return(-1+sapply(strsplit(str,ch,fixed=T),length)) }

 for (cc in 1:length(fnz))
 {
   tabs.per.line <- num.tabs(readLines(fnz[cc],n=max.feas.hdr),ch=sep.chr)
   full.n.tabs <- tail(tabs.per.line,1)
   fl <- 1; notfnd <- T
   while(notfnd) {
     if(tabs.per.line[fl]!=full.n.tabs) { fl <- fl + 1 } else { notfnd <- F }
   }
   if(firstrowhdr) { header.lens[cc] <- fl } else { header.lens[cc] <- fl-1 }
 }
 return(header.lens)
}

colMedianQS <- function(mmm)
{
 # apply 'quicksort' C implemented median function to each column of a matrix
 len <- as.integer(ncol(mmm)) 
 apply(mmm,2,medianQS, len)
}

medianQS <- function(X,len=NULL,na.rm=T)
{
 # calculate median using 'quicksort' algorithm written in C
 if(na.rm) { X <- X[!is.na(X)]; len <- NULL }  
 if(is.null(len)) { len <- length(X) }
 C.out = .C("quickerselect",len,X,0.01,NAOK=T,DUP=F)
 return((C.out[[3]]))
}

textogram <- function(X,...)
{
  # print a text based histogram of data, extra args may be passed to hist()
  hdat <- hist(X,plot=F,...)
	dens <- round(100*hdat$density/sum(hdat$density))
	if(max(dens)>90) {
		cat(" using halved %'s as maximum freq is too big for terminal\n")
		dens <- dens/2
	}
	for (cc in 1:length(dens))
	{
		label <- pad.zeros(hdat$breaks[cc],char=" ")
		cat(label," ",paste(rep(".",times=dens[cc]),collapse=""),"\n")
	}
}
	

preload.bioC <- function(bio.packs=c("IRanges","BiocGenerics","Biobase","GenomicRanges","genoset")) 
{
  ## load bioconductor packages without annoying package loading text ##
  cat(" silently loading bioconductor packages:",paste(bio.packs,collapse=", "),"... ")
  suppressWarnings(suppressMessages(must.use.package(bio.packs,T,quietly=T)))
  cat("done\n")
}


must.use.package <- function(pcknms,bioC=F,...)	
{
 ## like 'base::library()' but can easily work for bioconductor packages, and will
 # automatically attempt to install any required packages not found
 for (cc in 1:length(pcknms))
 {
   nxt.pck <- pcknms[cc]
   got1 <- require(nxt.pck,character.only=T,warn.conflicts=F,...)
   if(!got1) {
     if(bioC) {
       source("http://bioconductor.org/biocLite.R")
       biocLite(nxt.pck)
       library(nxt.pck,character.only=T,warn.conflicts=F,...)
     } else {
        install.packages(nxt.pck,character.only=T)
        library(nxt.pck,character.only=T,warn.conflicts=F,...) 
     }
   } else {
     #
   }
 }
}



list.rowsummary <- function(snpMatLst,mode="row")
{
 # performs 'row.summary' or 'col.summary' snpStats functions
 # on a list of snpMatrix objects
 fail <- F
 typz <- sapply(snpMatLst,is)[1,]
 if(all(typz=="SnpMatrix"))
 { HD <- F } else {
   if (all(typz=="character")) { HD <- T } else {
     cat("Warning: snpMatLst doesn't appear to be a list of SnpMatrix objects",
         "or a list of file locations, row.summary impossible\n")
     fail <- T
   } 
 }
 # this line allows this function to be either row.summary or col.summary (mode)
 if(mode!="col" & mode!="row") { mode <- "row" }
 my.summary <- switch(mode,row=row.summary,col=col.summary)
 if (!fail) {
   must.use.package("snpStats",T)
   rowsum.list <- list()
   for (dd in 1:length(snpMatLst))
   {
     cat(" processing element ",dd,"...",sep="")
     if(HD) {
       rowsum.list[[dd]] <- my.summary(get(paste(load(snpMatLst[[dd]]))))
     } else {
       rowsum.list[[dd]] <- my.summary(snpMatLst[[dd]])
     }
     cat("done\n")
   }
   return(do.call("rbind",rowsum.list))
 } else {
   return(NULL)
 }
}


list.colsummary <- function(snpChrLst,mode="col")
{
 # wrapper to make 'list.rowsummary' work for 'col.summary' too
 if(mode!="col" & mode!="row") { mode <- "col" }
 return(list.rowsummary(snpChrLst,mode=mode))
}

snp.mat.list.type <- function(snpMatLst,fail=T)
{
 # plumbCNV snp-qc can be run on a list of snpMatrix objects in memory ('memory') or
 # on a list of RData files on 'disk' to save RAM, this fn detects which type is in use
 typz <- sapply(snpMatLst,is)[1,]
 if(all(typz=="SnpMatrix"))
 { HDt <- "memory" } else {
   if (all(typz=="character")) { HDt <- "disk" } else {
     HDt <- "error"
     if(fail) {  stop("Error: not a valid snpMatLst!") 
     } else { cat("Warning: not a valid snpMatLst!\n")  }
   }
 }
 return(HDt)
}

get.snpMat.spec <- function(snpMatLst,fail=T)
{
 # get dimensions of snpMatLst (SnpMatrix list) regardless
 # of whether it's a set of SnpMatrix objects or list of file locations
 HD <- switch(snp.mat.list.type(snpMatLst,fail),memory=F,disk=T,error=NULL)
 if(HD) {
   n.grp <- length(snpMatLst)
   list.spec <- matrix(integer(),ncol=n.grp,nrow=2)
   for (cc in 1:n.grp)
   {
     snpMat <- get(paste(load(paste(snpMatLst[[cc]]))))
     list.spec[,cc] <- dim(snpMat)
   }
 } else {
   list.spec <- sapply(snpMatLst,dim)
 }
 return(list.spec)
}

convert.smp.to.chr22 <- function(snpMatLst,snpInfo,HD=F,loc="")
{
 ## convert snp.matrix list separated into different sample sets
 ## with all chromosomes, into 22 chromosome lists with all samples
 ## HD=T saves chromosome list elements to disk (loc) instead of memory
 ## in a situation where data is very large or memory is limited
 ## In HD=T the output will be a list of disk locations rather than an array
 HD <- switch(snp.mat.list.type(snpMatLst,fail=T),memory=F,disk=T,error=NULL)
 cat("\nConverting from 'n' group list of all markers to 22 chr list of all samples\n")
 snpChrLst <- vector("list", 22)
 #num.subs <- length(subIDs.actual)
 if (is(snpInfo)[1]=="RangedData")
 {
   must.use.package("IRanges",bioC=T)
   rang <- T
 } else {
   rang <- F
 }
 must.use.package("snpStats",bioC=T)
 list.spec <- get.snpMat.spec(snpMatLst)
 rownames(list.spec) <- c("Samples","SNPs")
 cat("\n"); print(list.spec); cat("\n")
 #e.g, sanger    t1d    uva
 #[1,]   4537   6808   5461   --> samples
 #[2,] 196524 196524 196524   --> snps
 if(any(diff(range(list.spec[2,]))!=0)) { cat("Warning! SNP files have different numbers of markers\n") }
 st <- 1+c(0,cumsum(list.spec[1,])[1:(ncol(list.spec)-1)])
 en <- cumsum(list.spec[1,])
 num.subs <- sum(list.spec[1,],na.rm=T)
 for (cc in 1:22)
 {
   cat(" processing chromosome ",cc,"...",sep="")
   if (rang) {
     # snpInfo is a RangedData object (IRanges)
     next.chr <- rownames(snpInfo[cc]) 
   } else {
     # snpInfo is a table or matrix, col aa = chr, col bb = snp label
     aa <- 1; bb <- 2
     cat("Warning: snpInfo object is not IRanges, assuming col",aa,"is chr and col",
         bb,"is snp label. If this is wrong please modify aa/bb in function 'convert.smp.to.chr22'\n")
     cat("\nFile preview of snp.info source file:\n")
     print(head(snpInfo))
     next.chr <- snpInfo[as.numeric(snpInfo[,aa])==cc,bb]
   }
   next.mat <- matrix(raw(),nrow=num.subs,ncol=length(next.chr))
   colnames(next.mat) <- next.chr
   if(HD) { rownames(next.mat) <- paste(1:num.subs)
   } else {
     rownames(next.mat) <- unlist(c(sapply(snpMatLst,rownames)))
   }
   next.snpmat <- new("SnpMatrix", next.mat)
   for (dd in 1:length(snpMatLst))
   {
     cat(dd,"..",sep="")
     if(HD){
       snpMat <- get(paste(load(paste(snpMatLst[[dd]]))))
       col.select <- next.chr %in% colnames(snpMat)
       next.snpmat[st[dd]:en[dd],col.select] <- snpMat[,next.chr[col.select]]
       rownames(next.snpmat)[st[dd]:en[dd]] <- rownames(snpMat)
       snpMat <- NULL
     } else {
       col.select <- next.chr %in% colnames(snpMatLst[[dd]])
       next.snpmat[st[dd]:en[dd],col.select] <- snpMatLst[[dd]][,next.chr[col.select]]
     }
   }
   if(HD) {
     fnm <- paste(loc,"chrmat_",cc,sep="")
     save(next.snpmat,file=fnm)
     snpChrLst[[cc]] <- fnm
   } else {
     snpChrLst[[cc]] <- next.snpmat
   }
   cat("..done\n")
 }
 return(snpChrLst)
}


convert.chr22.to.smp <- function(snpChrLst,snpInfo,subIDs.actual,HD=F,group.nums,loc="")
{
 ## convert snp.matrix list separated into 22 chromosome lists with all samples
 ## into different sample sets with all chromosomes. 
 ## HD=T saves chromosome list elements to disk (loc) instead of memory
 ## in a situation where data is very large or memory is limited
 ## In HD=T the output will be a list of disk locations rather than an array
 HD <- switch(snp.mat.list.type(snpChrLst,fail=T),memory=F,disk=T,error=NULL)
 cat("\nConverting from 22 chr list of all samples to 'n' group list of all markers\n")
 num.grps <- length(unique(group.nums))
 snpMatLst <- vector("list", num.grps)
 num.subs <- length(subIDs.actual)
 grp.sizes <- table(group.nums)
 if (is(snpInfo)[1]=="RangedData")
 {
   must.use.package("IRanges",bioC=T)
   rang <- T
 } else {
   rang <- F
 }
 must.use.package("snpStats",bioC=T)
 list.spec <- get.snpMat.spec(snpChrLst)
 #e.g, chr1  chr2  chr3  chr4  chr5  chr6  chr7snp  , ...   ... 
 #[1,] 16806 16806 16806 16806 16806 16806 16806 <-- samples
 #[2,] 23314 21355 12307  6444 13768 22316  7751  <-- markers
 if((any(diff(range(list.spec[1,]))!=0))) { cat("Warning! SNP files have different numbers of samples\n") }
 st <- 1+c(0,cumsum(list.spec[2,])[1:(ncol(list.spec)-1)])
 en <- cumsum(list.spec[2,])
 num.snps <- rowSums(list.spec)[2]
 for (dd in 1:num.grps) {
   cat(paste(" writing group",dd,"\n"))
   next.mat <- matrix(raw(),nrow=grp.sizes[dd],ncol=num.snps)
   colnames(next.mat) <- rownames(snpInfo)
   next.grp <- subIDs.actual[group.nums==dd]
   rownames(next.mat) <- next.grp
   next.snpmat <- new("SnpMatrix", next.mat)
   for (cc in 1:22)
   {
     if(HD){
       snpMat <- get(paste(load(paste(snpChrLst[[cc]]))))
       row.select <- next.grp %in% rownames(snpMat)
       next.snpmat[row.select,st[cc]:en[cc]] <- snpMat[next.grp[row.select],]
       snpMat <- NULL
     } else {
       cat(" processing chromosome ",cc,"...",sep="")
       row.select <- next.grp %in% rownames(snpMatLst[[cc]])
       next.snpmat[row.select,st[cc]:en[cc]] <- snpChrLst[[cc]][next.grp[row.select],]
       cat("..done\n")
     }
   }
   if(HD) {
     fnm <- paste(loc,"snpmat_",dd,sep="")
     save(next.snpmat,file=fnm)
     snpMatLst[[dd]] <- fnm
   } else {
     snpMatLst[[dd]] <- next.snpmat
   }
   cat(" complete!\n")
 }
 return(snpMatLst)
}

get.gc.markers <- function(snp.fn, mode="vcf", wndw=(10^6/2), hgV=18, 
                          ret=c("chr.info","gc","bio")[1], dir="", ...)
{
 # for a list of SNPs, get the average GC% in the 1MB surrounding window
 # load annotation package from bioconductor
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 hgFn <- paste("BSgenome.Hsapiens.UCSC.hg",hgV,sep="")
 must.use.package(c("genoset","BiocGenerics",hgFn),bioC=T)
 cat("\nGenerate GC percentage for SNPs\n")
 # recalculate: chrLens,chrStarts for bioconductor annotation (could be different to other annotation)
 # these are just the lengths and genome starting positions for each autosome
 chrLs <- as.numeric(seqlengths(Hsapiens)[1:22])
 chrSts <- c(0,cumsum((chrLs)))[1:22] # genome position of chromosome starts

 CHR.INFO <- calc.chr.ind(dir, snp.fl=NULL, mode=mode, 
                          vcf.fn=snp.fn, sav=F , ...)

 CHR.INFO <- update.chr.obj(CHR.INFO,chrLens=chrLs,chrStarts=chrSts,long.chrIDs=unlist(CHR.INFO$chrIDs))
 CHR.INFO <- remove.boundary.snps.window(CHR.INFO,window=wndw)
 lData <- convert.chr.obj.IRanges(CHR.INFO)

 cat(" loading GC data for SNPs...")
 gc.dat <- loadGC(lData, expand = wndw, bsgenome = Hsapiens)
 cat("complete\n")

 CHR.INFO[["long.GC"]] <- numeric(length(CHR.INFO$long.chrIDs))
 CHR.INFO$long.GC[match(rownames(gc.dat),CHR.INFO$long.chrIDs)] <- round(gc.dat$gc,2) # update VCF file with GC info just calculated

 CHR.INFO$long.GC[CHR.INFO$long.GC==0] <- NA # set any zero values to NA (as zero = missing)
 cat(" GC done!\n")
 # choose from 3 possible formats to return GC data
 return(switch(ret,chr.info=CHR.INFO,gc=CHR.INFO$long.GC,bio=gc.dat))
}

get.gc.human <- function(windowsize=10^6,hgV=18,ret=c("bio","gc")[1]) {
 # get the mean GC% for every 1MB window of the human genome
 # can set build ('hgV') to 18 or 19
 hgFn <- paste("BSgenome.Hsapiens.UCSC.hg",hgV,sep="")
 must.use.package(c("genoset","BiocGenerics",hgFn),bioC=T)
 # get chromosome lengths from annotation (possibly slightly different to other sources)
 chrLens <- as.numeric(seqlengths(Hsapiens)[1:22])
	chrsz <- chrLens %/% windowsize
	windowsts <- list()
	for (cc in 1:22) {	
		windowsts[[cc]] <- (((1:(chrsz[cc]-1))-1)*windowsize)+1 }
	  
	leno <- sapply(windowsts,length)
 coco <- rep(c(1:length(leno)),leno)

	stz <- as.integer(as.numeric(unlist(windowsts)))

 cat("\nExtracting GC% for each",windowsize,"window of the human genome...")
	lData <- RangedData(ranges=IRanges(start=(stz+(windowsize/2)),width=1),
                  space=coco)

	gc.dat <- loadGC(lData, expand = windowsize/2, bsgenome = Hsapiens)
 gc.dat$ranges <- (flank(gc.dat$ranges,windowsize/2,both=T)) # accurately reflect GC ranges (expand)

 cat("complete\n")
 cat("NB: please ignore warnings about previous imports and 'number of items to replace' - this is normal\n")
 det.nm <- paste("package:",hgFn,sep="")
 detach(det.nm,character.only = TRUE) # un-load large annotation package.
	return(switch(ret,bio=gc.dat,gc=gc.dat$gc))
}



remove.samp.from.plink <- function(samps,plink,why="TMcnvs")
{
	# remove specific sample from a plink file - eg because had too many cnvs
 ext <- c("cnv","fam","cnv.map")
 what <- c("cnvs","samples","")
 for (dd in 1:3) {
   fnn <- paste(plink,ext[dd],sep=".")
   inp <- readLines(fnn)
   to.rmv <- NULL
   for (cc in 1:length(samps))
   {
     to.rmv <- c(to.rmv,grep(samps[cc],inp))
   }
   to.rmv <- unique(to.rmv)
   if (length(to.rmv)>0)
   {
     out <- inp[-to.rmv]
     cat(paste("",length(to.rmv),what[dd],"removed from",fnn,"\n"))
   } else {
     out <- inp
     cat(paste(" no samples removed from",fnn,"\n"))
   }
   writeLines(out,con=paste(plink,why,ext[dd],sep="."))
 }
}



## CALCULATE LRR STATS ##
lrr.stats.tab <- function(stats.table)
{
  ## based on a table of sample-wise LRR stats (eg, mean, stdev, dlrs, etc)
  # calculate distribution indices for each LRR-sample-statistic for the whole cohort
	statz <- colnames(stats.table)
  stats.table <- as.data.frame(stats.table)  #ensure type is dataframe (eg. not matrix)
	TableRowLabels <- c("Mean","SD","Min","Q1","Median","Q3","Max","LB","UB","-2SD","+2SD","-3SD","+3SD","Low1%","Hi1%")
	nrowz <- length(TableRowLabels)
	ncolz <- length(statz)
	s.tab <- matrix(numeric(),ncol=ncolz,nrow=nrowz)
	
	for (cc in 1:ncolz)
	{
	 nxt.type <- statz[cc]
	 next.stat <- stats.table[[paste(nxt.type)]]
	 sl <- summary(next.stat)
	 IQR <- (sl[5]-sl[2])
	 pcts <- get.pctile(next.stat)
	 StDev <- sd(next.stat)
	 s.tab[,cc] <- c(sl[4],StDev,sl[c(1:3,5:6)],sl[2]-1.5*IQR,sl[5]+1.5*IQR,sl[4]+c(-2,2,-3,3)*StDev,pcts)
	}	
	rownames(s.tab) <- TableRowLabels
	colnames(s.tab) <- statz
	return(s.tab)
}


plot.extreme.samples <- function(rez,stat.table,bigMat2,CHR.INFO,dir,pref="",scl=10^6,ap="")
{
  # plot the most extreme sample for each combination of QC failure types
  dir <- validate.dir.for(dir,c("ind","big"),warn=F)
  venn.lists <- get.fail.type.list(rez,stat.table[,"Mean"]) 
  ex.id <- get.extreme.examples(stat.table,rez,venn.lists)
  cat("\nPlotting extreme samples for combinations of: Mean vs GC vs DLRS\n")
	loopz <- length(ex.id)
  warnz <- NULL
	ofn <- character(loopz)
	for (cc in 1:loopz)
	{
		if(!is.na(ex.id[cc])) {
		 loc.col <- match(ex.id[cc],colnames(bigMat2)) 
		 #LRR.dat <- sub.big.matrix(bigMat2, firstCol = loc.col, lastCol = loc.col, backingpath=dir$big)
     LRR.mat <- bigMat2[1:nrow(bigMat2),loc.col]
		 titl <- c(paste("Sample",ex.id[cc]),paste("Extreme LRR",names(ex.id)[cc]))
		 ofn[cc] <- paste(dir$ind,names(ex.id)[cc],ap,"Sample",ex.id[cc],pref,".pdf",sep="")
		 pdf(ofn[cc])
		   par(mar=c(5, 4, 10, 2))
		   plot((CHR.INFO$long.gnmIndx/scl),LRR.mat,col=CHR.INFO$long.chrCols,pch=".",cex=1.5,
		        xlab = "Genome Position (Megabases)", ylab="Log R Ratio",
		        main = titl)
		   axis(side=3,at=((CHR.INFO$chrStarts+(CHR.INFO$chrLens[1:22]/2))/scl), labels=paste(1:22))
		   mtext ("Chromosome number",side=3,line=2.5,cex=.9)
		 dev.off()
     loop.tracker(cc,loopz)
		} else { 
      #no subjects in this category, skipped..
      warnz <- c(warnz,paste(names(ex.id)[cc]))
      loop.tracker(cc,loopz)
    }
	}
  if(!is.null(warnz)) {
    cat("Warning: no subjects found for extreme values of:",paste(warnz,collapse=","),"\n")
  }
	cat(paste("produced",length(ofn),"files:\n"))
  cat(get.dir(ofn[rev(order(nchar(ofn)))][1]),";\n",sep="")
	cat(paste("  ",rmv.spc(rmv.dir(ofn)),"\n",sep="")); cat("\n")
}	

get.all.samp.fails.old <- function(dir)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl"),warn=F)
  ###src <- dir$excl;   src <- paste(dir$ano,"SAMPLE_EXCLUDE/",sep="")
  fnz <- list.files(dir$excl)
  excl.samps <- character()
  for (cc in 1:length(fnz)) {
    excl.samps <- c(excl.samps,readLines(paste(dir$excl,fnz[cc],sep="")))
  }
  excl.samps <- unique(excl.samps)
  return(excl.samps)
}


get.all.snp.fails <- function(dir,verb=F)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl2"),warn=F)
  snps.to.cut <- NULL
  all.fls <- list.files(dir$excl2)
  for (cc in 1:length(all.fls))
  {
    pr.len <- length(snps.to.cut)
    next.fl <- readLines(dir.fn.cmb(dir$excl2,all.fls[cc]))
    snps.to.cut <- unique(c(snps.to.cut,next.fl))
    if(verb) {  cat("",(length(snps.to.cut)-pr.len),
                    "Snps to remove, found in",all.fls[cc],"\n") }
  }
  return(snps.to.cut)
}


get.all.samp.fails <- function(dir,verb=F)
{
  # get all sample ids in fail lists in dir.ano/SAMPLE_EXCLUDE/
  dir <- validate.dir.for(dir,c("excl"),warn=F)
  samps.to.cut <- NULL
  all.fls <- list.files(dir$excl)
  for (cc in 1:length(all.fls))
  {
     pr.len <- length(samps.to.cut)
     next.fl <- readLines(dir.fn.cmb(dir$excl,all.fls[cc]))
     samps.to.cut <- unique(c(samps.to.cut,next.fl))
     if(verb) {  cat("",(length(samps.to.cut)-pr.len),
                     "samples to remove (across all groups), found in",all.fls[cc],"\n") }
  }
  return(samps.to.cut)
}


get.chr.ab.sample.info <- function(dir, failerMode="NOTLRR", n.bad.per.samp, max.bad.to.plot=22, chrWarns) 
{
 ### SCRIPT TAKES LIST OF SAMPLES WITH CHR ABERRATIONS (using contrasts)
 ##  AND SELECTS A SUBSET TO FLAG/PLOT WHICH CAN BE BASED ON NOT FAILING OTHER CRITERIA
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 # get previous LRR exclusion lists from files
 fnfn <- paste(dir$ano,"SAMPLE_EXCLUDE/",c("Mean","DLRS","GCWave"),".txt",sep="")  #"StDev",
 excluded.prev <- character()
 for (cc in 1:length(fnfn)) {
   excluded.prev <- c(excluded.prev,readLines(fnfn[cc]))
 }
 excluded.prev <- unique(excluded.prev)

 ## which samples to use with respect to those excluded already on other stats
 if(failerMode==("NOTLRR")) {
   checklist <- n.bad.per.samp[which(!names(n.bad.per.samp) %in% excluded.prev)]
 }
 if(failerMode==("ONLYLRR")) {
   checklist <- n.bad.per.samp[which(names(n.bad.per.samp) %in% excluded.prev)]
 }
 if(failerMode==("ALL")) { checklist <- n.bad.per.samp }

 # select sample subset with only up to a max number of bad chromosomes
 # eg, might only be interested in looking at those with 1 bad chromosome
 badcheckz <- names(checklist)[checklist<=max.bad.to.plot]

 # decode for each sample which chromosomes are bad and create labels for graphs
 chrLab <- character(length(badcheckz))
 chrN <- list()
 for (mm in 1:length(badcheckz)){
   failstr <- names(unlist(chrWarns))[which(unlist(chrWarns) %in% badcheckz[mm])]
   extr.str <- sapply(strsplit(failstr,"_",fixed=T),"[",1)
   chrN[[mm]] <- gsub("chr","",extr.str)
   if(length(failstr)==1)
   { chrLab[mm] <- paste("Chromosome",chrN[[mm]]) } else {
     chrLab[mm] <- paste("Chromosomes",paste(chrN[[mm]],collapse=","))
   }
 }
 names(chrN) <- badcheckz
 out.list <- list(badcheckz,chrN,chrLab)
 names(out.list) <- c("badcheckz","chrN","chrLab")
 return(out.list)
}

insert.fake.median <- function(pass.stats)
{
	# a hack for now to avoid fixing code so Median can be left out of tables
	cnmz <- colnames(pass.stats)
	if(!"Median" %in% cnmz)
	{
		cnmz <- cnmz[c(1,2,1,3)] ;cnmz[3] <- "Median"
		pass.stats <- pass.stats[,c(1,2,1,3)]
		colnames(pass.stats) <- cnmz	
	}
	return(pass.stats)
}


add.chr.top.axis <- function(select, badChrN, x.co, nC=22, sub=T) 
{
 # on chromosome graphs there is an option in plumbCNV to have the
 # chromosomes on the top axis, this script manages this allowing for
 # different levels of zoom which may cause overlaps. manages it better
 # than the default. also can change a target chr label to red colour.
 badChrz <- ((1:nC)[select] %in% badChrN )
 ## labels get messed as these too close together- so hack to fix
 if(any((c("18","19","21","22")) %in% badChrN ))
 {
   if("19" %in% badChrN )
   { badChrz[20] <- T } else {
     if("21" %in% badChrN )
     { badChrz[c(20,22)] <- NA
     } else { 
       if("22" %in% badChrN )
       { badChrz[21] <- T } else { badChrz[19] <- T }
     }
   }
 }
 axis(side=3,at=x.co[select][!badChrz], labels=paste(1:nC)[select][!badChrz],
      col.axis="black")
 if(length(paste(1:nC)[select][badChrz])>0) {
   axis(side=3,at=x.co[select][badChrz], labels=paste(1:nC)[select][badChrz],
      col.axis="red") }
 if (sub) { mtext ("Chromosome number",side=3,line=2.5,cex=.9) }
 return(list(select,badChrz,badChrN))
}



# Black and White Type Plot
bw.plot.lrr <- function (samp.chr.means, SX=NULL, pos.INFO=NULL, samp.chr.sds=NULL, num.labels=T,
                        CI.Z=1.96, grpAv=NA, grpLoCI=NA, grpHiCI=NA, scl=10^6, nC=22,...) 
{
 # Alternate graph without the full LRR data
 # This plots the means of each chromosome and their confidence intervals
 # also overlays the mean and std deviation bands of the reference group
 # of all good samples
 # CI.Z = confidence interval z score
 if(is.null(SX)) {  SX <- ((pos.INFO$chrStarts+(pos.INFO$chrLens[1:nC]/2))/scl) }
 YY <- samp.chr.means
 plot(SX, grpAv, cex=1.5,
      xlab = "Genome Position (Megabases)", ylab="Log R Ratio", 
      ,type="l", lty="dashed", col="black", ...)
 if (!is.na(grpLoCI[1]) & !is.na(grpHiCI[1])) {
   lines(SX,grpLoCI,lty="dotted",col="black")
   lines(SX,grpHiCI,lty="dotted",col="black")
 }
 lines(SX,YY) # subject means plot (main data)
 # confidence intervals based on specified confidence level
 if(!is.null(samp.chr.sds) & !is.null(pos.INFO)) {
   SS <- CI.Z*(samp.chr.sds/sqrt(pos.INFO$nECF))
   arrows(SX,YY,SX,YY+SS,angle=90,length=0.05)
   arrows(SX,YY,SX,YY-SS,angle=90,length=0.05)
 }
 if(num.labels) {
   text(x=SX,y=YY+(.1),labels=paste(round(YY,3)),srt=90,cex=.75) }
 return(SX)
}

# Color Type Plot

col.plot.lrr <- function (ID, bigMat, centre.chr=1:22, CHR.INFO=list(),
                         plotAdj=F, samp.chr.means=NULL, set.chr.ref=T, 
                         c.pre.n=2, c.post.n=2, m.smooth=F, scl=10^6, nC=22, ratio=10,...) 
{
 # colour type plot of LRR data for a sample over a given range, optional smoothing, overlay
 # of means, etc.
 #targ.chr.ref sets the x axis scale start point (e.g., relative to start of what chromosome)
 if(set.chr.ref) {
	targ.chr.ref <- as.numeric(centre.chr[1]) } else {
	targ.chr.ref <- NULL	
	}
 if(!is.chr.obj(CHR.INFO)) { stop("Error: not a valid CHR.INFO object") }
 XX <- (CHR.INFO$long.gnmIndx/scl)
 chr.x.labs <- ((CHR.INFO$chrStarts+(CHR.INFO$chrLens[1:nC]/2))/scl)

 # extract samples' LRR data from big matrix
 loc.col <- match(ID,colnames(bigMat))
 LRR.mat <- bigMat[1:nrow(bigMat),loc.col] #CHR.INFO$long.gnmIndx
 #LRR.dat <- sub.big.matrix(bigMat, firstCol = loc.col, lastCol = loc.col, backingpath=b.dir)
 #LRR.mat <- LRR.dat[1:nrow(LRR.dat)]
 # set indices if only plotting a limited range of genome, eg only 2 adjacent chrs, etc
 if (plotAdj) {
   ii <- range(as.integer(centre.chr))
   chr.select <- max(1,(ii[1]-c.pre.n)):min((ii[2]+c.post.n+1),nC)
   rng <- CHR.INFO$chrStarts[range(chr.select)]/scl
   rng[1] <- floor(rng[1]); rng[2] <- ceiling(rng[2])
   select <- c( head(which(XX>rng[1]),1): tail(which(XX<rng[2]),1) ) 
 } else {
   select <- c(1:length(LRR.mat)) 
   chr.select <- 1:nC
 }
 # apply smoothing (mean/median) to LRR data
 if (m.smooth) { 
   yy <- MedianScale(LRR.mat[select],ratio) 
   xx <- MedianScale(XX[select],ratio)
   ccc <- (CHR.INFO$long.chrCols[select])[seq(1,length(LRR.mat[select]),length.out=length(xx))]
 } else {
   xx <- XX[select]; yy <- LRR.mat[select] 
   ccc <- CHR.INFO$long.chrCols[select] 
 }
 # Do the large colour plot
 if (!is.null(targ.chr.ref)) 
 {    
   # ie adjust scale if zooming into a single chromosome
   offset <- (CHR.INFO$chrStarts[targ.chr.ref]/scl)
   xx <- xx-offset
   chr.x.labs <- chr.x.labs-offset
 } else { offset <- 0 }
 plot(xx,yy,col=ccc,pch=".",type="p",cex=1.5,...) 
 # add lines between the means of each chromosome
 if (!is.null(samp.chr.means)) {
  lines(chr.x.labs[chr.select],samp.chr.means[chr.select]) }
 out <- list(chr.x.labs,chr.select,select,offset)
 names(out) <- c("x.coords","chr.select","select","offset")
 return(out)
}


search.tree2 <- function(full.line,matchlist,rng)
{
 # attempt at fast way to search through long list (not used???)
 shorter <- substr(full.line,1,rng[2])
 spli <- strsplit(shorter," ",fixed=T)[[1]][1]
 return(spli %in% matchlist)
}




create.buffer <- function(matchlist) {
  ## cool search but takes twice as long as conventional method
 char.r <- range(nchar(matchlist))
 n.snp <- length(matchlist)
 find.tree <- matrix(character(),ncol=(diff(char.r)+1),nrow=n.snp)
 dd <- char.r[1]
 find.tree[,1] <- substr(matchlist,1,dd)
 for (dd in (char.r[1]+1):char.r[2])
 {
   find.tree[,dd-char.r[1]+1] <- substr(matchlist,dd,dd)
 }
 rownames(find.tree) <- paste(1:nrow(find.tree))
 return(find.tree)
}

search.tree <- function(full.line,ftree)
{
 # attempt at fast way to search through long list (not used???)
 nf <- T; coln <- 1
 len1 <- nchar(ftree[1,1])
 next.search <- substr(full.line,1,len1)
 nC <- ncol(ftree)
 while(nf & (coln<=nC))
 {
   less.rows <- which(ftree[,coln]==next.search)
   ftree <- ftree[less.rows,]
   if(length(less.rows)==1) { 
     nf <- F
   } else { 
     next.search <- substr(full.line,len1+coln,len1+coln) 
     coln <- coln + 1
   }
 }
 if(!nf) {
   str.fnd <- (paste(ftree,collapse=""))
   nc <- nchar(str.fnd)
   if (substr(full.line,1,nc)==str.fnd)
   { return(T) } else { (return(F))  }
 } else {
   return(F)
 }
}


rg2 <- function(pos,ws,chrs)
{
 # not really sure what this little function does?
 targ.pos <- chrs+pos   #-1
 reg.num <- targ.pos %/% ws + 1
 return(reg.num)
}


PC.fn <- function(next.row,nPCs,col.sel)
{
 # apply PC correction for a single SNP, allowing for missing data.
 bad1 <- which(is.na(next.row))
 if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
 next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals
 return(next.row)
}

parse.fn <- function(text.list,n=1)
{  
  # simplified way to access element #n of each element of a list
  return(as.numeric(sapply(text.list,"[",n)))  
}

smart.match <- function(loc,big,es=10^5,hint=NA,dont=F,ss=5,tol=10^4) 
{
 #fast matching in a big vector when extra info available (not used???)
 if(is.na(hint)) {
   pred.spot <- ((loc-big[1]) %/% ss) } else { es <- 10^3 ; pred.spot <- (hint+es-10) }
 st <- max(1,(pred.spot - es)); en <- (pred.spot+es)
 subt <- abs(big[st:en]-loc)
 mm <- min(subt)
 if (mm < tol) {
   rel <- match(mm,subt)[1]
   if(is.na(rel) & !dont) { 
     out <- smart.match(loc,big,es=es*100,dont=T)
   } else {
     out <-(st+rel-1)
   }
 } else {
   cat(" sorry - something weird going on, trying the slow match now\n")
   out <- match(loc,big)
   if(!is.na(out)) { cat(paste(loc,"succeeded\n")) } else { cat(paste("match to",loc,"failed\n")) }
 }
 return(out)
}

rg <- function(rng,gc.locs)
{
 # range get where position value is different to
 # array index, eg for chr position in array of markers, etc (not used???)
 tt <- range(rng)
 o1 <- smart.match(tt[1],gc.locs)
 o2 <- smart.match(tt[2],gc.locs,hint=o1)
 return(o1:o2)
}


check.hdr <- function(hdr,chr.list) 
{
 # for parsing headers in sanger GCBASE5 GC content UCSC files
 #if first line contains 'random' skip whole file
 ifrandom <- grep("random",hdr)
 ifchrom <- grep("chrom",hdr)
 next.chr <- 0
 if(length(ifrandom)==0) {
   if(length(ifchrom)>0) {
     hdrF <- gsub("variableStep chrom=chr","",hdr)
     hdrF <- gsub(" span=5","",hdrF)
     hdrF <- substr(hdrF,1,2)
     hdrF <- gsub("_","",hdrF)
     if (hdrF %in% paste(chr.list))
     { 
       next.chr <- as.numeric(hdrF)
     }
   }
 }
 return(next.chr)
}


reader <- function(fn,dir="",pref="",suf="",want.type=NULL,def="\t",qt=F,treatas=NULL,override=F,...)
{
 # try to read in data from many types of datafile
 types <- c("RData","csv","txt","tab")
 # add "/" to directory if not present
 if(dir!="" & substr(dir,nchar(dir),nchar(dir))!="/")
 { dir <- paste(dir,"/",sep="") }
 # split to get file extension
 segs <- strsplit(fn,".",fixed=T)[[1]]
 fsuf <- tail(segs,1)
 # imply main file name
 file.main <- substr(fn,1,nchar(fn)-nchar(fsuf)-1)
 # construct the full path
 full.path <- paste(dir,pref,file.main,suf,".",fsuf,sep="")
 if(!is.null(treatas) & (treatas[1] %in% types))
 {
   cat(paste(" will assume file is",treatas,"\n"))
   typ <- treatas[1]
 } else {
   # select file type to open as
   typ <- types[types %in% fsuf][1]
 }
 if(length(typ)==0)
 {
   # unknown file type
   cat(paste("Warning: file type suffix:",fsuf,"unknown\n"))
   cat("needs to be one of:\n"); print(types)
   cat("will attempt to read as lines of text\n")
   file.out <- readLines(fn)
 } else {
   if(typ==types[1])
   {
     ## .RData binary file
     file.out <- list()
     fns <- load(full.path)
     if(!qt) { cat(paste(" found",paste(fns,collapse=","),"in",full.path,"\n")) }
     if (!is.null(want.type)) {
       for (ii in 1:length(fns)) {
         if(want.type %in% is(get(fns[ii])))
         { fns <- fns[ii] ; break }
       } }
     for (cc in 1:length(fns))
     {  file.out[[cc]] <- get(fns[cc]) }
     names(file.out) <- fns
     if(length(file.out)==1) { file.out <- file.out[[1]] }
   }
   if(typ==types[2])
   {
     # csv file
     file.out <- read.csv(full.path, ...)
     file.out <- first.col.to.rownames(file.out,override)
   }    
   if(typ==types[3] | typ==types[4])
   {
     # other text .txt file
     first.2 <- readLines(full.path,n=2)
     splitto <- sapply(strsplit(first.2,def,fixed=T),length)
     if(all(splitto==c(1,1))) {
       # probably just a vector file
       file.out <- readLines(full.path)
     } else {
       # some kind of delimited file
       if (all(splitto[1]==splitto[2]))
       {
         file.out <- read.delim(full.path, ..., sep=def)
         file.out <- first.col.to.rownames(file.out,override)
       } else {
         cat(paste("Warning: *.txt file not delimited by",def,"\n"))
         cat(" will just read as text\n")
         file.out <- readLines(full.path)
       }
     }
   }
 }    
 return(file.out)
}

first.col.to.rownames <- function(dataf,override=F)
{
 # if a dataframe's first column is IDs, convert the dataframe so these become rownames
 if(is.data.frame(dataf))
 {
   typz <- sapply(sapply(dataf,is),"[[",1)
   rn <- paste(dataf[,1])
   dataout <- dataf[,-1]
   if (typz[1]=="character" & all(typz[-1] %in% c("numeric","integer")))
   {
     numerify <- T
   } else {
     suppressWarnings(tst <- length(which(is.na(as.numeric(rn))))/nrow(dataout) )
     if (tst>=.75 | override)
     { 
       suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dataout[,1])))))/nrow(dataout) )
       if(tst<=.75) {
         numerify <- T
       } else {
         numerify <- F }
     } else {
       cat("Warning: proposed rownames were mostly numbers\n")
       return(dataf)
     }
   }    
   if(numerify) { for (dd in ncol(dataout))
   { dataout[,dd] <- as.numeric(as.character(dataout[,dd])) }
   } 
   rownames(dataout) <- rn
   return(dataout)
 } else {
   if(length(dim(dataf))==2)
   {
     if("character" %in% is(dataf[,1]))
     {
       if(length(unique(dataf[,1])==nrow(dataf)))
       {
         if(ncol(dataf)>1) {
           sup <- assess.dat.type(dataf)
           if(sup<2) { 
             cat("Warning: not sure if rownames should be added from col 1\n")
             return(dataf)
           } else {
             rn <- dataf[,1]
             dataout <- dataf[,-1]
             rownames(dataout) <- paste(rn)
             if(sup>=10) { for (dd in ncol(dataout))
               { dataout[,dd] <- as.numeric(as.character(dataout[,dd])) }
             }    
             return(dataout)
           }            
         } else {
           cat("Warning: only a single column\n")
           return(dataf)
         }
       } else {
         cat("Warning: duplicate row names\n")
         return(dataf)
       }
     }
   } else {
     cat("Warning: must be formatted in rows and columns\n")
     return(dataf)
   }
   cat("Warning: couldn't change first col, not a dataframe!\n")
   return(dataf)
 }
}

assess.dat.type <- function(dat)
{
 # try to work out what format a file is in, whether IDs in column 1
 nchar.mn <- function(vec) { mean(nchar(paste(vec)),na.rm=T) }

 if(length(unique(dat[,1])==nrow(dat)))
 { support <- support+1 }
 suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dat[,1])))))/nrow(dat) )
 if (tst>=.75)
 { 
   suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dat[,2])))))/nrow(dat) )
   if(tst<=.75) {
     support <- support+10
   } else {
     support <- support+1 }
 }
 if (ncol(dat)>2 & support <2) {
   # slow and potentially inaccurate so only run if in doubt
   rc <- dim(dat)
   if (nrow(dat)>100) { r.sel <- 1:100 } else { r.sel <- 1:rc[1] } 
   if (ncol(dat)>50) { c.sel <- 1:50 } else { c.sel <- 1:rc[2] } 
   col.chr <- apply(dat[r.sel,c.sel],2,nchar.mn)
   dif.to.1 <- mean(abs(col.chr[-1]-col.char[1]))
   dif.to.o <- mean(abs(col.chr[-1]-rev(col.char[-1])))
   if (dif.to.0 < dif.to.1)
   { support <- support+1 } else { support <- support - 1 }
 }
 return(support)
}

dlrs <- function(X,na.rm=T)
{
  # calculate dlrs
  out <- (sd(diff(X),na.rm=na.rm))
  return(out/sqrt(2))
}

stdev <- function(x,na.rm=T) { sd(x,na.rm=na.rm) }

circle <- function(x,y,r) {
 # draw a circle
 theta <- seq(0, 2 * pi, length=(10000))
 xo <- x + r * cos(theta)
 yo <- y + r * sin(theta)
 return(cbind(xo,yo))
}


getBigMat <- function(fn,dir)
{
 # loads a big.matrix either using an big.matrix description object
 # , or this object in a binary file or text file, or points to a bigmatrix or matrix
 dir <- validate.dir.for(dir,c("big"),warn=F)
 if(is(fn)[1]=="big.matrix.descriptor")
 {
   bigMat2 <- attach.big.matrix(fn,path=dir$big)
 } else {
   if(is(fn)[1]=="big.matrix" | is(fn)[1]=="matrix")
   {
     if(is(fn)[1]=="matrix") {
       bigMat2 <- as.big.matrix(fn,descriptorfile="TEMPBIG",backingpath=dir$big)
     } else { bigMat2 <- fn }
   } else {
     lastchar <- substr(dir$big,nchar(dir$big),nchar(dir$big))
     if (length(grep(".RData",fn))==0) {
       fn <- rmv.dir(fn)
       if(!fn %in% list.files(dir$big)) { 
         stop(paste("Error: big.matrix file '",fn,"' not in 'dir$big'",sep=""))
       }
       cat(" loading big matrix using text description\n")
       if(lastchar=="/") { dir$big <- substr(dir$big,1,nchar(dir$big)-1) }
       bigMat2 <- attach.big.matrix(fn,path=dir$big)
     } else {
       cat(" loading big matrix using RData description\n")
       if(lastchar!="/") { dir$big <- paste(dir$big,"/",sep="") }
       filenm <- dir.fn.cmb(dir$big,fn,must.exist=T)
       dscnm <- paste(load(filenm))
       big.fn <- NULL
       for (ii in 1:length(dscnm)) {
         if("big.matrix.descriptor" %in% is(get(dscnm[ii])))
         { big.fn <- dscnm[ii] } 
       }
       if(!is.null(big.fn)) {
         descr <- get(big.fn) 
       } else {
         stop(paste("Error: didn't find bigmatrix descriptor in file",fn))
       }
       bigMat2 <- attach.big.matrix(descr,path=dir$big) 
     }
   }
 }
 return(bigMat2)
}

choose.best.sort.file <- function(search.dir,ref.list,zero.replace=T)
{
 # hopefully sort directory only contains 1 valid file
 # if there are multiple files, no files, or only invalid files with no matches to ref.list
 # (i.e, to the current big.matrix) this script will try to extract/create the best possible list
 file.list <- paste(search.dir,list.files(search.dir),sep="")
 if(length(grep("SNP_SORT",search.dir,fixed=T))>length(grep("SAMPLE_SORT",search.dir,fixed=T)))
 { guess.mode.is.snp <- T } else { guess.mode.is.snp <- F }
 goodm <- F
 if(length(file.list)>0)
 {
   if(length(file.list)>1) {
     cat(" multiple files found in sort annotation directory - should only be 1\n")
     all.snp <- list(); lens <- integer()
     for (cc in 1:length(file.list))
     {
       all.snp[[cc]] <- (readLines(file.list[cc]))
       lens[cc] <- length(which(ref.list %in% all.snp[[cc]]))
     }
     best.matches <- max(lens,na.rm=T)
     cat("",round(100*(best.matches/length(ref.list)),1),"% of big.matrix indices in sort.list\n")
     if(best.matches>0) {
       out.order <- all.snp[[which(lens==best.matches)[1]]]; goodm <- T
     } 
   } else {
     file1 <- readLines(file.list[1])
     val.len <- length(which(ref.list %in% file1))
     if (val.len>0) {
       cat("",round(100*(val.len/length(ref.list)),1),"% of big.matrix indices in sort.list\n")
       cat(" using",file.list[1],"to sort",c("samples","snps")[1+guess.mode.is.snp],"\n")
       out.order <- file1
       goodm <- T
     } else {
       cat(" Warning: no matches to reference list found in sort file, assuming invalid\n")
     }
   }
 }
 if(!goodm) {
   cat(" no valid files found in sort annotation directory - creating new from bigmatrix names\n")
   if (guess.mode.is.snp) { fnm <- "snpsort.txt" } else { fnm <- "sampsort.txt" }
   ofn <- paste(search.dir,fnm,sep="")
   writeLines(paste(ref.list),con=ofn)
   cat("wrote file:",ofn,"\n")
   out.order <- ref.list
 }
 return(out.order)
}

check.file.and.subset <- function(fnm, mat, by.row=T, min.size=1000, stop.if.fail=T)
{
 # check that file dir/fnm is a genuine subset of the row/col names (by.row) of matrix 'mat'
 # return subset or if not valid, then return all row/col names of the original 'mat'
 if(by.row)  { ref <- rownames(mat)  } else { ref <- colnames(mat) }
 if(stop.if.fail) { actn <- "aborting script" } else { actn <- "selecting all" }
 if(file.exists(fnm)) {
   out.L <- readLines(fnm)
   #remove headers if required (assumption header unlikely with same number of chars as ids)
   h1 <- out.L[1]
   if(!nchar(h1) %in% nchar(out.L[2:99])) { out.L <- out.L[-1]  ;
                                            cat("\nauto removed header:",h1) }
   if(length(which(out.L %in% ref))<min.size)
   { 
     wrn <- paste(c("column","row")[1+by.row],"subset file",fnm,
                  "did not have at least",min.size,"matches to the reference matrix --> ",actn)
     if(stop.if.fail) { stop(wrn) } else { cat(wrn,"\n"); out.L <- ref }
   }
 } else { 
   wrn <- paste(c("column","row")[1+by.row],"subset file",fnm,"was not found --> ",actn)
   if(stop.if.fail) { stop(wrn) } else { cat(wrn,"\n"); out.L <- ref } 
 }
 return(out.L)
}

sort.exclude.from.annot <- function(bigMat,dir="",dosnp=T,dosamp=T,ordR=T,ordC=T, verb=T)
{
 # sort big matrix according to annotation and exclude those in exclusion files
 dir <- validate.dir.for(dir,c("ano","sort","excl","sort2","excl2"),warn=F)
 cat(" calculating exclusions and order\n")
 samp.names <- colnames(bigMat); snp.names <- rownames(bigMat)
 # order files should only be one file - fix otherwise
 sample.order <- choose.best.sort.file(dir$sort,ref.list=samp.names,T)
 snp.order   <-  choose.best.sort.file(dir$sort2,ref.list=snp.names,T)
 ## Get SNP and SAMP exclusion lists ##
 snps.to.cut <- samps.to.cut <- NULL
 if(dosamp) { samps.to.cut <- get.all.samp.fails(dir,verb) }  
 if(dosnp) { snps.to.cut <- get.all.snp.fails(dir,verb) }  
 # use sort/exclusion lists to get reordering vectors
 if(ordR) { to.order.r <- match(snp.order,rownames(bigMat)) } else { to.order.r <- c(1:nrow(bigMat)) }
 if(ordC) { to.order.c <- match(sample.order,colnames(bigMat)) } else { to.order.c <- c(1:ncol(bigMat)) }
 to.order.r <- narm(to.order.r); to.order.c <- narm(to.order.c)
 to.remove.r <- narm(match(snps.to.cut,rownames(bigMat)[to.order.r]))
 to.remove.c <- narm(match(samps.to.cut,colnames(bigMat)[to.order.c]))
 if (length(to.remove.r)>0) { to.order.r <- to.order.r[-to.remove.r] }
 if (length(to.remove.c)>0) { to.order.c <- to.order.c[-to.remove.c] }
 out.list <- list(to.order.r,to.order.c,to.remove.r,to.remove.c,samps.to.cut,snps.to.cut)
 names(out.list) <- c("to.order.r","to.order.c","to.remove.r","to.remove.c","sample.excl","snp.excl")
 return(out.list)
}


select.samp.snp.custom <- function(bigMat,snp,samp)
{
 # based on files/vectors of snp-ids and sample-ids create selection
 # vectors to select only the ids in these lists for a matrix
 cat(" calculating selections for snps\n")
 # try to detect whether a vector of IDs, or file names
 snp.ref <- rownames(bigMat)  ; samp.ref <- colnames(bigMat) 
 
 if (length(snp)==1 & length(samp)==1 & is.character(snp) & is.character(samp))
 {
   cat(" [assuming 'samp' and 'snp' are file names containing sample and snp ids]")
   if(file.exists(snp)) {
     snp.sel <- readLines(snp)
   } else {
     if(snp=="") {
       cat(c(" snp subset file was empty, selecting all\n"))
       snp.sel <- snp.ref
     } else {
       stop("Error: argument 'snp' should be a vector of SNPs length>1 or a filename with a list of snps (no header)")
     }
   }
   if(file.exists(samp)) {
     sample.sel <- readLines(samp)
   } else {
     if(samp=="") {
       cat(c(" sample subset file was empty, selecting all\n"))
       sample.sel <- samp.ref
     } else {
       stop("Error: argument 'samp' should be a vector of Samples length>1 or a filename with a list of snps (no header)")
     }
   }
 } else { 
   #cat("[assuming 'samp' and 'snp' are vectors of sample and snp ids]")
   # if blank then assign all ids
   if(all(snp=="")) {
     snp.sel <- snp.ref
   } else {
     snp.sel <- snp
   }
   if(all(samp=="")) {
     sample.sel <- samp.ref
   } else { 
     sample.sel <- samp  
   }
 }
 # use sort/exclusion lists to get reordering vectors
 row.sel <- snp.sel ; col.sel <- sample.sel

 #print(head(row.sel));print(head(col.sel))
 to.order.r <- narm(match(row.sel,rownames(bigMat)))
 to.order.c <- narm(match(col.sel,colnames(bigMat)))
 if (!(length(to.order.r[!is.na(to.order.r)])>0 & length(to.order.c[!is.na(to.order.c)])>0))
 { cat("Warning: selection of SNPs and/or Samples has resulted in an empty dataset!\n") 
   cat("Check rownames, column names and selection lists for errors\n")}

 out.list <- list(to.order.r,to.order.c,snp.sel,sample.sel)
 names(out.list) <- c("to.order.r","to.order.c","snp.list","sample.list")
 return(out.list)
}

big.matrix.operation <- function(bigM1,bigM2,operation="-",bck,des,dir,split.to=40,low.ram=F)
{
 # apply simple arithmetic functions to a pair of big matrices
 # splits the operation into 'split.to' sets of rows. only need to increase this if
 # the datafile is more than 20 times the available RAM memory.
 dir <- validate.dir.for(dir,c("big"),warn=F)
 if(require(bigalgebra)) {
   cat(" using bigalgebra method to reduce RAM requirement for DLRS matrix\n")
   options(bigalgebra.mixed_airthmetic_returns_R_matrix=FALSE)
   use.big.al <- T
 } else {
   use.big.al <- F
   cat(" bigalgebra not installed, may increase RAM requirement for DLRS matrix\n")
 }
 if (!operation %in% c("+","-","/","*"))
 { cat(paste("Warning: operation",operation,"may not be supported\n")) }

 if(all(dim(bigM1)==dim(bigM2)))
 {
   nR <- nrow(bigM1); nC <- ncol(bigM1)
   cat(" applying function to matrices. expect this to take some time\n")
   #bigR <- eval(call("-",bigM1[1:nR,1:nC],bigM2[1:nR,1:nC]))
   bigO <- big.matrix(nR,nC, backingfile=bck,
                         backingpath=dir$big, descriptorfile=des)
   stepz <- round(seq(from=1,to=nR+1,length.out=round((split.to+1))))
   for (cc in 1:split.to)
   {
   	x1 <- stepz[cc]; x2 <- ((stepz[cc+1])-1)
    if(use.big.al & operation=="-") {
      bigO[x1:x2,1:nC] <- bigM2[x1:x2,1:nC] - bigM1[x1:x2,1:nC]
    } else {
     	bigO[x1:x2,1:nC] <- eval(call(operation,bigM1[x1:x2,1:nC],bigM2[x1:x2,1:nC]))
    }
   	loop.tracker(cc,split.to,freq=1)
    if(cc %% 4 == 0) {
      # reset memory after every 4 chunks to prevent skyrocketing RAM use #
      fl.suc <- flush(bigO) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()  
      if(low.ram) {
        RR <- describe(bigO); rm(bigO); bigO <- attach.big.matrix(RR,path=dir$big)
      }
    }
   }
   cat(" function complete, converting result to big matrix\n")
   return(bigO)
 } else {
   cat(" Warning: matrix dimensions don't match\n")
   return(NULL)
 }
}


load.data.to.bigmat <- function(dat.fn,inputType="TXT",bck,des,dir,Hopt=F,RNopt=T)
{
 # get data stored in a text or binary file and get it into a bigmatrix object format
 dir <- validate.dir.for(dir,c("ano","big"),warn=F)
 if(inputType=="RDATA")
 {
   # attach datafile bigmemory object
   cat("\nLoading RData file...")
   this.mat <- load(dat.fn)
   if (this.mat!="m.by.s.matrix") { 
     m.by.s.matrix <- get(this.mat)
     rm(this.mat) 
   }
   cat("complete\n")
   if(is.null(colnames(m.by.s.matrix)))
   {
     gen.fl <- paste(dir$ano,"snpNames.txt",sep="")
     gene.list <- readLines(gen.fl)
     if (ncol(m.by.s.matrix)==length(gene.list))
     {
       colnames(m.by.s.matrix) <- gene.list
       cat(paste("Warning: added column names to matrix from",gen.fl,"\n"))
       cat("if replacement with these names is undesired please check/modify the .R source code\n")
     } else {
       cat("Error: matrix has no column names and default file doesn't\n")
       cat("match number of columns, please check/modify the .R source code\n")
       stop()
     }
   }
   cat(" saving datafile as big matrix\n")
   bigMat <- as.big.matrix(m.by.s.matrix, backingfile=bck,
                           backingpath=dir$big, descriptorfile=des)
 } else {
   # assume inputType=="TXT"
   cat("\nLoading TAB file...(this will probably be slow)...")
   read.big.matrix(dat.fn, sep = '\t', header = Hopt,
                   has.row.names=RNopt, ignore.row.names=FALSE,
                   backingfile = bck, backingpath = dir$big,
                   descriptorfile = des, extraCols = NULL)
   cat("complete\n")
 }
}

deident <- function(sample.id.list,dir="",lookup="deidentityCodes.RData")
{
	# de identify sample IDs in a consistent manner for tables, figs
	# lookup file of the form: sample.id.list fake.id.list
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  look.tab <- reader(lookup,dir$ano,qt=T)
	fake.id.list <- look.tab[,2][match(sample.id.list,look.tab[,1])]
	return(fake.id.list)
}

reident <- function(fake.id.list,dir="",lookup="deidentityCodes.RData")
{
	# de identify sample IDs in a consistent manner for tables, figs
	# lookup file of the form: sample.id.list fake.id.list
  dir <- validate.dir.for(dir,c("ano"),warn=F)
	look.tab <- reader(lookup,dir$ano,qt=T)
	sample.id.list <- look.tab[,1][match(fake.id.list,look.tab[,2])]
	return(sample.id.list)	
}

create.deid.lookup <- function(all.ids,dir="",fn="deidentityCodes.RData",replace=F)
{
  # create lookup table of deidentified IDs
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  ful.fn <- paste(dir$ano,fn,sep="")
	nids <- length(all.ids)
	fake.ids <- paste("ID",pad.zeros(1:nids),sep="")
	rand.order <- runif(nids)
	look.tab <- cbind(all.ids[order(rand.order)],fake.ids)
	colnames(look.tab) <- c("sample.id.list","fake.id.list")
 if(replace | (!fn %in% list.files(dir$ano)))
	{
		if(get.ext(fn)=="RData") { 
			save(look.tab,file=ful.fn) 
		} else {
			write.table(look.tab,file=ful.fn,col.names=F,row.names=F,sep="\t")	
		}
   cat(" wrote new lookup deidentification table:",ful.fn,"\n")
	} else {
		cat(" ID de-identifcation table previously generated (not replaced)\n")
	}
}

pad.zeros <- function(n,numdigits=NA,char="0")
{
	## function takes a vector and makes all numbers into a string
	## with each having the same number of digits, using zeros as
	## padding for lower numbers. default is to pad up to the longest
	## number in the vector, but can manually set the length with numdigits
	if (is.na(numdigits))
	{ numdigits <- nchar(as.character(max(n,na.rm=TRUE))) }
	padn <- character(length(n))
	padn[is.na(n)] <- NA
	all0 <- paste((rep(char,times=numdigits)),collapse="")
	padn[n==0] <- all0
	nz <- nchar(as.character(n[!is.na(n)]))
	if (max(nz,na.rm=TRUE) > numdigits)
	{ cat(paste(" padding set too low for some values - will be values bigger than",
		numdigits,"digits\n")) }
	for (cc in 1:(length(n[!is.na(n)])))
	{
		padn[!is.na(n)][cc] <- paste(substr(all0,1,(numdigits-nz[cc])) , 
			as.character(n[!is.na(n)][cc]),sep="")
	}
	return(padn)
}


rmv.dir.old <- function(X) {
  # return file name (without directory if present)
  file.segs <- strsplit(X,"/",fixed=T)[[1]]
  lss <- length(file.segs)
  if (lss>1) { out <- paste(file.segs[lss]) } else { out <- X }
  return(out)
}

rmv.dir <- function(X) {
  # return file name (without directory if present)
  if(is.null(X)) { X <- "" }
  file.segs <- strsplit(X,"/",fixed=T)
  out <- sapply(file.segs,tail,1)
  return(out)
}

get.dir <- function(X) {
  # return directory without filename (if present)
  file.segs <- strsplit(X,"/",fixed=T)[[1]]
  lss <- length(file.segs)
  if (lss>1) { out <- paste(file.segs[-lss],collapse="/") } else { out <- "" }
  return(out)
}


get.ext <- function(X) {
 # get file extension from a filename character string
 file.segs <- strsplit(X,".",fixed=T)[[1]]
 lss <- length(file.segs)
 if (lss>1) { out <- paste(file.segs[lss]) } else { out <- "" }
 return(out)
}


rmv.ext <- function(X) {
 # remove file extension from a filename character string
 file.segs <- strsplit(paste(X),".",fixed=T)[[1]]
 lss <- length(file.segs)
 if (lss>1) { out <- paste(file.segs[-lss],collapse=".") } else { out <- X }
 return(out)
}


dir.fn.cmb <- function(dir="",fn,pref="",suf="",ext="",must.exist=F,force=T) 
{
  # combine a directory and a filename, taking care not to double up.
  # if 'fn' already has a dir, that will override the 'dir' argument
  # making sure the slashes are right, allowing insertion of prefix, suffix or
  # extensions, and optionally stop() if the resulting file doesn't currently exist
  
  if(length(dir)>1) { dir <- dir[1]; cat("only first dir was used\n") }
  if(length(ext)>1) { ext <- ext[1]; cat("only first extension was used\n") }
  if(length(grep("/",fn))>0) {
    dir <- get.dir(fn)  #split into dir and fn if fn has /'s
    fn <- rmv.dir(fn)
  }
  dir <- dir.force.slash(dir)
  if(ext!="") {
    #make sure ext includes the dot
    if(substr(ext,1,1)!=".")   { ext <- paste(".",ext,sep="") }
    #if ext is already built into suffix or filename, remove it from there
    fn <- rmv.ext(fn)
    suf <- rmv.ext(suf)
  }
  location <- paste(dir,pref,fn,suf,ext,sep="")
  if(any(!file.exists(location)) & must.exist) {
    warn <- paste("Warning: required file",location,"not found!")
    stop(warn)
  }
  return(location)
}


chr.lab.to.num <- function(chr.col)
{
 # some annotation files will have some chromosome information entered
 # as text (e.g mitochondrial = MT, X/Y, etc) ; converts these to numbers
 chromogrps <- list("1","2","3","4","5","6","7","8","9","10","11","12"
                    ,"13","14","15","16","17","18","19","20","21","22",
                    c("X","Y","MT","0","23","24","25","26","27","28","XY"))
 chrgrp <- numeric(length(chr.col))
 for (cc in 1:length(chromogrps))
 {
   chrgrp[paste(chr.col) %in% chromogrps[[cc]]] <- cc
 }
 if (0 %in% chrgrp)
 {
   uniq.bad <- (unique(chr.col[chrgrp==0]))
   cat(paste("Error: unknown chromosome label in file:",uniq.bad,"\n"))
   cat("change file or add value to 'chromogrps' in 'chr.lab.to.num'\n")
   stop(print(chr.col[which(!paste(chr.col) %in% unlist(chromogrps))]))
 }
 return(chrgrp)
}

parse.args <- function(arg.list,coms=c("M"),def=c(0))
{
 # parse arguments entered running R from the command line
 # defaults
 if(length(coms)>1 & length(def)==1) { def <- rep(def,length(coms)) }
 outframe <- data.frame(val=paste(def),stringsAsFactors=F)
 rownames(outframe) <- coms
 assign.cmds <- grep("=",arg.list,fixed=T)
 if(length(assign.cmds)>0)
 {
   vars.lst <- strsplit(arg.list[assign.cmds],"=",fixed=T)
   vars <- sapply(vars.lst,"[",1)
   vals <- sapply(vars.lst,tail,1)
   vals <- paste(vals)
   #vals <- as.integer(vals)

   if(any(toupper(vars) %in% coms))
   {
     which.coms <- match(toupper(vars),coms)
     for (cc in 1:length(which.coms))
     {
       if(!is.na(vals[cc]) & !is.na(which.coms[cc]))
       {
         print(paste("set",coms[which.coms[cc]],"=",vals[cc]))
         outframe[coms[which.coms[cc]],1] <- paste(vals[cc])
       } else {
         cat(paste(" skipping invalid variable",vars[cc],"or invalid value",vals[cc],"\n"))
       }
     }
   } else {
     cat("Warning: command line arguments entered but none are valid\n")
   }
 } else {
 	outframe <- NULL
 } 
 return (outframe)
}



all.counts <- function(tab,grp.range=2,char="+",mode="integer")
{
 # create table of combinations for a dataframe, with 'grp.range' vars at a time
 nn <- ncol(tab)
 combs <- all.combs.n.binary(nn,colnames(tab),colchar=char)
 combs <- matrix(as(as.matrix(combs),mode),ncol=ncol(combs),
                 dimnames=list(rownames(combs),colnames(combs)))
 use.combs <- combs[,colSums(combs) %in% grp.range]
# combs <- apply(combs,2,as.logical)
 return(combs)
}

all.combs.n.binary <- function(n,name=paste(c(1:n)),colchar="+")
{
 # return all possible logical combinations of a set of binary variables
 nc <- 2^n
 store <- character(nc)
 for (cc in ((1:nc)-1)) 
 { 
   store[cc+1] <- (paste(rev(as.integer(intToBits(cc))), collapse="")) 
 }
 maxlen <- nchar(store[1]) #32 usually?
 firstNonZero <- min(which(strsplit(tail(store,1),"")[[1]]==1))
 store <- substr(store,firstNonZero,maxlen)
 store.fr <- as.data.frame(strsplit(store,""),stringsAsFactors=F)
 for (dd in 2:nc) {
   colnames(store.fr)[dd] <- paste(name[as.logical(as.integer(store.fr[,dd]))],collapse=colchar) }
 colnames(store.fr)[1] <- "none"
 return(store.fr)
}


###LRR stats function defintions###

get.excl.list.from.stats <- function(stat.dat,stat.frame, dir="",
 lo.stat.nm=c("LB",NA,"LB"),hi.stat.nm=c("UB","UB","UB"),writeFiles=T,append=F)
{
 # calculate which samples in stat.dat should be excluded based on statistical thresholds
 # on the LRR data in stat.frame
 ### CALC EXCL SAMPLES ###
 headers <- colnames(stat.dat)
 excl.list <- list()
 sample.List <- row.names(stat.frame)
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 for (dd in 1:length(headers))
 {
   excl.list[[dd]] <- character()

   if (!is.na(lo.stat.nm[dd])) {
     lo.stat <- stat.dat[lo.stat.nm[dd],headers[dd]] 
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[stat.frame[,headers[dd]]<lo.stat]) }

   if (!is.na(hi.stat.nm[dd])) {
     hi.stat <- stat.dat[hi.stat.nm[dd],headers[dd]] 
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[stat.frame[,headers[dd]]>hi.stat]) }

   if(writeFiles) {
     fnfn <- paste(dir$ano,"SAMPLE_EXCLUDE/",headers[dd],".txt",sep="")
     if(!append) {
       writeLines(excl.list[[dd]],con=fnfn) 
       cat(paste("wrote sample list:",fnfn,"\n"))
     } else {
       if(file.exists(fnfn)) {
         exist.lns <- readLines(fnfn) ####
       } else {
         exist.lns <- NULL
       }
       to.wrt <- unique(c(exist.lns,excl.list[[dd]]))
       num.added <- length(to.wrt)-length(exist.lns)
       writeLines(to.wrt,con=fnfn) 
       cat(paste(" appended",num.added,"to sample list:",fnfn,"\n"))     
     }
   }
 }
 fail.table <- matrix(logical(),nrow=nrow(stat.frame),ncol=length(lo.stat.nm))
 for (cc in 1:ncol(fail.table)) {
   fail.table[,cc] <- sample.List %in% excl.list[[cc]] }

 rownames(fail.table) <- sample.List; colnames(fail.table) <- headers
 return(fail.table)
}


get.borders.from.stats <- function(stat.dat,stat.frame,lo.stat.nm=c("LB",NA,NA,NA),
  hi.stat.nm=c("UB","UB","UB","UB"),lo.stat.nm2=c("-2SD",NA,NA,NA),
  hi.stat.nm2=c("+2SD","+2SD","+2SD","+2SD"),writeFiles=T)
{
 # CALC SAMPLES falling between possible cutoff thresholds 
 headers <- colnames(stat.dat)
 excl.list <- list()
 sample.List <- row.names(stat.frame)
 for (dd in 1:length(headers))
 {
   excl.list[[dd]] <- character()
   x <- stat.frame[,headers[dd]]
   if (!is.na(lo.stat.nm[dd])) {
     lo.stat <- stat.dat[lo.stat.nm[dd],headers[dd]]
     lo.stat2 <- stat.dat[lo.stat.nm2[dd],headers[dd]]
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[x < lo.stat2 & x > lo.stat]) }

   if (!is.na(hi.stat.nm[dd])) {
     hi.stat <- stat.dat[hi.stat.nm[dd],headers[dd]] 
     hi.stat2 <- stat.dat[hi.stat.nm2[dd],headers[dd]] 
     excl.list[[dd]] <- c(excl.list[[dd]],sample.List[x > hi.stat2 & x < hi.stat]) }
 }
 minn <- 4  # as regardless of how many stats, for later functions, rez needs just length 3
 fail.table <- matrix(logical(),nrow=nrow(stat.frame),ncol=min(minn,length(headers)))
 for (cc in 1:ncol(fail.table)) {
   fail.table[,cc] <- sample.List %in% excl.list[[cc]] }

 rownames(fail.table) <- sample.List; colnames(fail.table) <- headers
 return(fail.table)
}


get.fail.type.list <- function(rez,chk,labels=c("Mean","DLRS","GCWave"))
{
 ## GET LISTS OF SUBJECTS FAILING ON ALL POSSIBLE SETS OF CRITERIA
 ## USE FOR INDIVIDUAL GENOME vs LRR Plot Examples
 venn.lists <- list()
 venn.lists[[1]] <- rownames(rez)[rez[,1] & !rez[,2] & !rez[,3] & chk<mean(chk,na.rm=T)]
 venn.lists[[2]] <- rownames(rez)[rez[,1] & !rez[,2] & !rez[,3] & chk>mean(chk,na.rm=T)]
 venn.lists[[3]] <- rownames(rez)[!rez[,1] & rez[,2] & !rez[,3]]
 venn.lists[[4]] <- rownames(rez)[!rez[,1] & !rez[,2] & rez[,3]]
 venn.lists[[5]] <- rownames(rez)[rez[,1] & rez[,2] & !rez[,3]]
 venn.lists[[6]] <- rownames(rez)[rez[,1] & !rez[,2] & rez[,3]]
 venn.lists[[7]] <- rownames(rez)[!rez[,1] & rez[,2] & rez[,3]]
 venn.lists[[8]] <- rownames(rez)[rez[,1] & rez[,2] & rez[,3]] 
 venn.lists[[9]] <- rownames(rez)[!rez[,1] & !rez[,2] & !rez[,3]]

 names(venn.lists)[1:2] <- paste(labels[1],c("Lo","Hi"),sep="")
 names(venn.lists)[3:4] <- labels[2:3]
 names(venn.lists)[5:7] <- c(paste(labels[1],labels[2],sep="+"),
                             paste(labels[1],labels[3],sep="+"),
                             paste(labels[2],labels[3],sep="+"))
 names(venn.lists)[8:9] <- c(paste(labels,collapse="+",sep=""),"None")
 return(venn.lists)
}

get.distinct.cols <- function(n=22)
{
 # get set of distinct colours (e.g, n>10; most builtin packages do a max of 12 distinct)
 # hand picked set of colour numbers to be distinguishable
 distinct.cols <- c(38,6,23,11,12,26,30,31,94,134,100,47,139,53,58,68,116,128,172,142,367,656,77)
 # reorder so similar colours aren't adjacent.
 distinct.cols <- distinct.cols[c(16,10,5,1,21,8,13,18,7,11,3,20,22,14,2,6,19,4,17,12,9,15)]
 colz <- colors()[distinct.cols[1:min(n,length(distinct.cols))]]
 if(n>length(distinct.cols)) { cat(paste("warning: only",length(distinct.cols),"were available\n")) }
 return(colz)
}


remove.boundary.snps.window <- function(CHR.INFO,window=0,chr.n=22)
{
 ## remove snps which when adding a window exceed the boundaries of their chromosome
 ## e.g, for analyses with windows calculations, e.g 1MB either side of each snp, etc
 if(is.chr.obj(CHR.INFO) & length(CHR.INFO$chrLens)>=chr.n) {
   for (cc in 1:chr.n) {  
     too.big <- which(CHR.INFO$chrIndx[[cc]]>(CHR.INFO$chrLens[cc]-window))
     # [1] 16210 16211 16212 16213 16214 16215 16216 16217 16218 16219 16220
     too.small <- which(CHR.INFO$chrIndx[[cc]]<(window+1))
     to.remove <- c(too.big,too.small)
     if(length(to.remove)>0) {
       CHR.INFO$chrIndx[[cc]] <- CHR.INFO$chrIndx[[cc]][-to.remove]
       CHR.INFO$chrIDs[[cc]] <- CHR.INFO$chrIDs[[cc]][-to.remove] 
     }
   }
 } else {
     cat("Warning: not a CHR.INFO object, function failed to remove boundary snps - GC code may crash\n")
 }
 return(CHR.INFO)
}


update.chr.obj <- function(CHR.INFO.OLD=NULL,long.chrIndx=NA,long.chrCols=NA,long.chrIDs=NA,long.gnmIndx=NA,
                          long.GC=NA,chrLens=NA,nECF=NA,chrStarts=NA,chrStrtsF=NA,
                        chrEndsF=NA,chrIDs=NA,chrIndx=NA,chrCols=NA,gnmIndx=NA)
{
 # update existing CHR.INFO object with some fields
 CHR.INFO.NEW <- list(long.chrIndx,long.chrCols,long.chrIDs,long.gnmIndx, long.GC,chrLens,nECF,chrStarts,
                  chrStrtsF,chrEndsF,chrIDs,chrIndx,chrCols,gnmIndx)
 all.names <- c("long.chrIndx", "long.chrCols", "long.chrIDs", "long.gnmIndx", "long.GC","chrLens",
                "nECF", "chrStarts", "chrStrtsF", "chrEndsF", "chrIDs","chrIndx", "chrCols", "gnmIndx")
 names(CHR.INFO.NEW) <- all.names
 for (cc in 1:length(CHR.INFO.OLD))
 {
   if (!(all(is.na(CHR.INFO.OLD[[all.names[cc]]]))) & all(is.na(CHR.INFO.NEW[[all.names[cc]]])) )
   {
     CHR.INFO.NEW[[all.names[cc]]] <- CHR.INFO.OLD[[all.names[cc]]]
   }
 }
 return(CHR.INFO.NEW)
}

make.chr.obj <- function(long.chrIndx=NA,long.chrCols=NA,long.chrIDs=NA,long.gnmIndx=NA,long.GC=NA,chrLens=NA,
                        nECF=NA,chrStarts=NA,chrStrtsF=NA,
                        chrEndsF=NA,chrIDs=NA,chrIndx=NA,chrCols=NA,gnmIndx=NA)
{
 # Make a CHR.INFO object using whatever data is available
 # needs objects from >> load(paste(dir.ano,"filteredChrInfo.RData",sep=""))
 # see also functions 'get.chr.filt.info' & 'calc.chr.ind'
 #long.chrIndx,  #genome positions of each snp [using chr position + chr start]
 #long.chrCols,  #chr based colour for each snp [all together, not list by chr]
 #long.chrIDs,   #ID for each snp [all together, not list by chr]
 #long.gnmIndx,  #genome positions of each snp [using chr position + chr start]
 #long.GC,       #local GC percentage for each snp 
 #chrLens,       #length of each chromosome
 #nECF,          #number of markers for each chromosome
 #chrStarts,     #genome wide start location (bp) for each chr
 #chrStrtsF,     #genome wide start location (#snp) for each chr
 #chrEndsF,      #genome wide end location (#snp) for each chr
 #chrIDs,        #list, 1 sublist for each chr with ID for each snp
 #chrIndx,       #list, 1 sublist for each chr with chr position of each snp
 #chrCols,       #list, 1 sublist for each chr with chr col for each snp
 #gnmIndx        #list, 1 sublist for each chr with full genome position for each snp

 CHR.INFO <- list(long.chrIndx,long.chrCols,long.chrIDs,long.gnmIndx,long.GC,chrLens,nECF,chrStarts,
                  chrStrtsF,chrEndsF,chrIDs,chrIndx,chrCols,gnmIndx)
 names(CHR.INFO) <- c("long.chrIndx", "long.chrCols", "long.chrIDs","long.gnmIndx","long.GC","chrLens",
                      "nECF", "chrStarts", "chrStrtsF", "chrEndsF", "chrIDs","chrIndx", "chrCols","gnmIndx")
 return(CHR.INFO)
}

print.chr.obj <- function(CHR.INFO,n=6,c=6)
{
 # nice printout of a CHR.INFO object without overloading the terminal
 if(is.chr.obj(CHR.INFO)){
   cat("\nPrintout of Object 'CHR.INFO':\n")
   num.long <- ((which(sapply(sapply(CHR.INFO,is),"[",1)=="list"))[1]-1)
   num.chr <- length(CHR.INFO)-num.long
   cat("\nVectors:\n")
   for (j in 1:num.long) {
     print(paste(paste(names(CHR.INFO)[j],":"),paste(CHR.INFO[[j]][1:n],collapse=" "),collapse=" "))
   }
   cat("\nLists:")
   for (j in (num.long+c(1:num.chr))) {
     cat("\n")
     print(paste(names(CHR.INFO)[j],":"))
     for (k in 1:c) {
       if (length(CHR.INFO[[j]])>=k) {
         cat(k,". ",sep="")
         cat(CHR.INFO[[j]][[k]][1:n],"\n") 
       } 
     }
   }
 } else {
   cat("Warning: not a CHR.INFO object!\n")
 }
}


sync.snpmat.with.info <- function(snpMatLst,infoRange=NULL,sample.info=NULL)
{
 # auto detect whether snpMatLst is a list of SnpMatrix objects or a list of RData
 # file locations and act accordingly. autodetect whether snp and/or sample info inputted.
 if(!is.null(infoRange)) {
   if (is(infoRange)[1]=="RangedData" & is(snpMatLst[[1]])[1]=="SnpMatrix")
   {
     for (cc in 1:length(snpMatLst)) {
       cat(" reordering/selecting markers from SnpMatrix",cc,"...")
       to.keep <- match(rownames(infoRange),colnames(snpMatLst[[cc]]))
       to.keep <- to.keep[!is.na(to.keep)]
       snpMatLst[[cc]] <- snpMatLst[[cc]][,to.keep]
       cat("done\n")
     }
     return(snpMatLst)
   } else {
     if(is(snpMatLst[[1]])[1]=="character") {
       if(all(file.exists(unlist(snpMatLst)))) {
         for (cc in 1:length(snpMatLst)) {
           cat(" reordering/selecting markers from SnpMatrix",cc,"...")
           snpMat <- get(paste(load(paste(snpMatLst[[cc]]))))
           to.keep <- match(rownames(infoRange),colnames(snpMat))
           to.keep <- to.keep[!is.na(to.keep)]
           snpMat <- snpMat[,to.keep]
           save(snpMat,file=paste(snpMatLst[[cc]]))
           cat("done\n")
         }
         return(snpMatLst)
       } else {
         cat(paste("Warning: invalid input parameters, could not find listed file:",unlist(snpMatLst),"\n"))
       }
     } else {
       cat("Warning: invalid input parameters, need 'RangedData' or list of file names\n")
     }
   }
 }
 if(!is.null(sample.info)) {
   if (is(sample.info)[1]=="data.frame" & is(snpMatLst[[1]])[1]=="SnpMatrix")
   {
     for (cc in 1:length(snpMatLst)) {
       cat(" reordering/selecting samples from sample.info",cc,"...")
       to.keep <- match(rownames(sample.info),rownames(snpMatLst[[cc]]))
       to.keep <- to.keep[!is.na(to.keep)]
       snpMatLst[[cc]] <- snpMatLst[[cc]][to.keep,]
       cat("done\n")
     }
     return(snpMatLst)
   } else {
     if(is(snpMatLst[[1]])[1]=="character") {
       if(all(file.exists(unlist(snpMatLst)))) {
         for (cc in 1:length(snpMatLst)) {
           cat(" reordering/selecting samples from sample.info",cc,"...")
           snpMat <- get(paste(load(paste(snpMatLst[[cc]]))))
           to.keep <- match(rownames(sample.info),rownames(snpMat))
           to.keep <- to.keep[!is.na(to.keep)]
           snpMat <- snpMat[to.keep,]
           save(snpMat,file=snpMatLst[[cc]])
           cat("done\n")
         }
         return(snpMatLst)
       } else {
         cat("Warning: invalid input parameters, need 'RangedData' and list of SnpMatrix objects\n")
       }
     } else {
       cat("Warning: invalid input parameters, need 'RangedData' and list of SnpMatrix objects\n")
     }
   }
 }
}


rmv.chr.23 <- function(CHR.INFO,silent=F)
{
 # remove chromosome 23 from a CHR.INFO object
 if(is.chr.obj(CHR.INFO)){
   lenz <- sapply(CHR.INFO,length)
   indxs.to.convert <- which(lenz==23)
   twenty3 <- names(CHR.INFO)[indxs.to.convert]
   if(length(indxs.to.convert)>0) {
     for (cc in 1:length(indxs.to.convert))
     {
       if(!silent) { cat(paste(" removing chr23 from",twenty3[cc],"\n")) }
       CHR.INFO[[indxs.to.convert[cc]]][[23]] <- NULL
       long.name <- paste("long.",twenty3[cc],sep="")
       if(!silent) { cat(paste(" re-writing",long.name,"\n")) }
       CHR.INFO[[long.name]] <- unlist(CHR.INFO[[indxs.to.convert[cc]]])
     }
   } else {
     #cat("No parameters of CHR.INFO had length = 23. CHR23 not altered in object.\n")
   }
 } else {
   stop("Error: not a CHR.INFO object!")
 }
 return(CHR.INFO)
}

is.chr.obj <- function(CHR.INFO)
{
 # is this object a 'CHR.INFO' object?
 if (is.null(CHR.INFO))
 {
   ITEST <- F
 } else {
   if(all(is.na(CHR.INFO))) {
     ITEST <- F
   } else {
     ITEST <- all(names(CHR.INFO) %in% c("long.chrIndx", "long.chrCols", "long.chrIDs", "long.gnmIndx", "long.GC", "chrLens",
                      "nECF", "chrStarts", "chrStrtsF", "chrEndsF", "chrIDs", "chrIndx", "chrCols", "gnmIndx"))
   }
 }
 return(ITEST)
}

get.chr.lens <- function(dir,len.fn="humanChrLens.txt",mode=c("hg18","hg19")[1],n=c(1:22))
{
 # retrieve chromosome lengths from local annotation file, else download from UCSC
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 chrlens.f <- paste(dir$ano,len.fn,sep="") # existing or future lengths file
 if(len.fn %in% list.files(dir$ano))
 {
   # file seems to be in annotation directory already
   chrLens <- readLines(chrlens.f)
   if (length(as.numeric(chrLens))!=length(n))
   {
     cat(paste("Length of chromosome file doesn't match expected",length(n),"\n"))
     notGot <- T
   } else { notGot <- F}
 } else { notGot <- T }
 if (notGot) {
   #download from UCSC
   cat("attempting to download chromosome lengths from UCSC\n")
   urL <- switch(mode,
             hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
             hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",)
   download.file(urL, chrlens.f)
   chrL.f <- readLines(chrlens.f)
   len.lst <- strsplit(chrL.f,"\t")
   nmz <- sapply(len.lst,"[",1)
   lnz <- sapply(len.lst,"[",2)
   want.chr.names <- match(paste("chr",n,sep=""),nmz)
   want.chr.names <- want.chr.names[!is.na(want.chr.names)]
   chrLens <- lnz[want.chr.names]
   names(chrLens) <- want.chr.names
   writeLines(chrLens,con=chrlens.f) # save file for future use
 }
 return(as.numeric(chrLens))
}


get.chr.filt.info <- function(bignames, dir, recalcF=T, recalcO=F, ret=F, pCHR.INFO=NULL,
                             vcf.file="load", nC=22, snp.info.fn="",
                             sav.fn="filteredChrInfo.RData", pre.fn="chr.snp.lists.RData", ...)
{
 # calculate chromosome indexes, in particular for the filtered version of snplist rownames(bigMat2)
 # basically does nothing if recalcF == FALSE
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 recalcO <- recalcO | ( !(pre.fn %in% list.files(dir$ano)) & is.null(pCHR.INFO) )
 recalcF <- recalcF | !(sav.fn %in% list.files(dir$ano))
 if(recalcO & !is.chr.obj(pCHR.INFO)) 
   { pCHR.INFO <- calc.chr.ind(dir,vcf.fn=snp.info.fn,sav.fn=pre.fn,sav=T,...) 
     cat(" [recalculated against latest SNP list]\n") }
 # get chr indexes: load num.chr, chr.snp.lists, chr.pos.lists
 #  [ID and position for each SNP within chromosomes]
 if(recalcF) {
   cat("\nCalculating indexes for filtered SNPs per chromosome...")
   if (!is.chr.obj(pCHR.INFO)) {
     cat(" loaded from file\n")
     CHR.INFO <- get(paste(load(file=paste(dir$ano,pre.fn,sep=""))))
   } else {
     CHR.INFO <- pCHR.INFO
   }
   chrLens <- get.chr.lens(dir,mode="hg18")
   chrStarts <- c(0,cumsum((chrLens)))[1:22] # genome position of chromosome starts
   numEachChr <- sapply(CHR.INFO$chrIndx,length) # num SNPs in array per chr.
   coloz <- get.distinct.cols(22)
   # get genome position, and assign a colour to each SNP within each chromosome
   # remove SNPs from each sublist that are no longer in dataframe (filtered out)
   chrIDs <- CHR.INFO$chrIDs
   chrIndx <- chrCols <- gnmIndx <- CHR.INFO$chrIndx
   chrIDs[[23]] <- chrIndx[[23]] <- gnmIndx[[23]] <- chrCols[[23]] <- NULL
   cat("",length(bignames),"ids entered for filtered snp set\n")
   for (dd in 1:22) 

   { 
     gnmIndx[[dd]] <- chrIndx[[dd]] + chrStarts[dd]
     chrCols[[dd]] <- rep(coloz[dd],numEachChr[dd])
     # remove snps which have been filtered in current dataset
     remvs <- which(!(paste(chrIDs[[dd]]) %in% bignames))
     if  (length(remvs)>0) {
       chrIndx[[dd]] <- chrIndx[[dd]][-remvs]
       chrCols[[dd]] <- chrCols[[dd]][-remvs] 
       chrIDs[[dd]] <- chrIDs[[dd]][-remvs]
       gnmIndx[[dd]] <- gnmIndx[[dd]][-remvs]
     }
   }
   numEachChrFilt <- sapply(chrIndx,length)
   nECF <- numEachChrFilt; names(nECF) <- NULL # auto-names cause index confusion
   # get start and end array indexs of chromosomes for a global array (filtered)
   chrStrtsF <- (1+c(0,cumsum(nECF)))[1:nC]
   chrEndsF <- cumsum(nECF)
   CHR.INFO <- make.chr.obj(chrLens=chrLens,nECF=nECF,chrStarts=chrStarts,chrStrtsF=chrStrtsF,
                            chrEndsF=chrEndsF,chrIDs=chrIDs,chrIndx=chrIndx,chrCols=chrCols,gnmIndx=gnmIndx)
   CHR.INFO <- update.longs(CHR.INFO)  # make unlisted versions
   cat(" snp.info object produced with",length(CHR.INFO$long.chrIndx),"snps\n")
   # save objects to file in annotation directory for future use
   save(CHR.INFO,file=paste(dir$ano,sav.fn,sep=""))
 } else { 
   load(paste(dir$ano,sav.fn,sep=""))
   cat(" filtered chromosome info not recalculated [as recalcF==F]\n")
 }
 if(ret) 
 {
   return(CHR.INFO)
 }
}

update.longs <- function(CHR.INFO,refresh=T) 
{
 # extend separate chromosome lists into long form if present/needed
 if(is.chr.obj(CHR.INFO)){
   ones.to.do <- c("chrIDs","chrIndx","chrCols","gnmIndx")
   long.ones <- paste("long.",ones.to.do,sep="")
   for (cc in 1:length(ones.to.do)) {
     if(refresh | (all(is.na(CHR.INFO[[long.ones[cc]]])) | is.null(CHR.INFO[[long.ones[cc]]]))) 
     {
       if(!all(is.na(CHR.INFO[[ones.to.do[cc]]][[1]])) & !is.null(CHR.INFO[[ones.to.do[cc]]][[1]])) {
         CHR.INFO[[long.ones[cc]]] <- unlist(CHR.INFO[[ones.to.do[cc]]])  
       }
     }
   }
 } else {
   stop("Error: not a CHR.INFO object")
 } 
 return(CHR.INFO)
}

convert.chr.obj.IRanges <- function(CHR.INFO,full=F)
{
 # convert CHR.INFO object to IRanges object
 must.use.package(c("IRanges","BiocGenerics"),T)
 if (is.chr.obj(CHR.INFO)) {
   CHR.INFO <- rmv.chr.23(CHR.INFO)
   CHR.INFO <- update.longs(CHR.INFO)

   # Generate chromosome numbers for each list element
   leno <- sapply(CHR.INFO$chrIndx,length)
   #for (yy in 1:22) { coco <- c(coco,rep(yy,leno[yy]))}
   coco <- rep(c(1:length(leno)),leno)

   # IRanges is sensitive to type of the position index, so ensure it is correct:
   stz <- (1*CHR.INFO$long.chrIndx) ;  names(stz) <- NULL
   stz <- as.integer(as.numeric(stz))
   nmz <- paste(CHR.INFO$long.chrIDs)

   if(length(nmz)==length(stz) & length(nmz)==length(coco))
   {
     cat(" converting CHR.INFO to IRanges, RangedData object...")
     if(full){
       # add all available data to the RangedData object
       colorz <- paste(CHR.INFO$long.chrCols)
       gnmindz <- as.numeric(CHR.INFO$long.gnmIndx)
       gcz <- as.numeric(CHR.INFO$long.GC)
       maintwo <- (length(nmz)==length(colorz) & length(nmz)==length(gnmindz))
       must.use.package("IRanges",bioC=T)
       if(length(gcz)==length(nmz)) {
         if(maintwo){
           lData <- RangedData(ranges=IRanges(start=stz,end=stz,width=1,names=nmz),space=coco,gindx=gnmindz,color=colorz,gc=gcz)
         } else {
           lData <- RangedData(ranges=IRanges(start=stz,end=stz,width=1,names=nmz),space=coco,gc=gcz)
         }
       } else {
         if(maintwo){
           lData <- RangedData(ranges=IRanges(start=stz,end=stz,width=1,names=nmz),space=coco,gindx=gnmindz,color=colorz)
         } else {
           lData <- RangedData(ranges=IRanges(start=stz,end=stz,width=1,names=nmz),space=coco)
         }
       }
     } else {
       # just insert components necessary for GC calculations
       lData <- RangedData(ranges=IRanges(start=stz,end=stz,width=1,names=nmz),space=coco)
     }
     cat("complete\n")
   } else {
     stop("Error: lengths of snp IDs, positions and chromosomes do not match, review file")
   }
 } else {
   cat("Warning: not a CHR.INFO object, null IRanges object produced\n")
   lData <- NULL
 }
 return(lData)
}


derive.file.delim <- function(fn,n=10,comment.char="#",delims=c("\t"," ",","),large=10)  
{
 # test top 'n' lines to determine what delimeter the file uses
 comment.char <- "#"
 test.bit <- readLines(fn,n)
 cmnt <- which(substr(test.bit,1,1)==comment.char)
 if(length(cmnt)>0) { test.bit <- test.bit[-cmnt] }
 num.del <- list()
 for (cc in 1:length(delims)) {
   num.del[[cc]] <- sapply(strsplit(test.bit,delims[[cc]],T),length)
 }
 #print(num.del)
 need.0 <- sapply(num.del,function(X) { sum(diff(X)) })
 num.del <- sapply(num.del,"[",1)
 if(any(!need.0)) {
   #rng <- range(num.del)
   candidates <- which(num.del>1 & num.del<=large & !need.0)
   if(length(candidates)>0) { out <- candidates[1] 
   } else {
     candidates <- which(num.del>large & !need.0)
     if(length(candidates)>0) { out <- candidates[1]
     } else {
       candidates <- which(num.del==1 & !need.0)
       if(length(candidates)>0) { out <- candidates[1]
       } else {
         cat("Warning: no delimiters tried were able to produce a valid file spec\n")
         out <- NULL
       }
     }
   }
 } else {
   cat("Warning: no delimiters tried were able to produce a valid file spec\n")
   out <- NULL
 }
 return(delims[out])
}

calc.chr.ind <- function(dir,snp.fl=NULL, mode=c("bim","vcf","map","map3")[1], vcf.file="load",
                        vcf.fn="", sav=F,
                        sav.fn="chr.snp.lists.RData", snp.col=NA, pos.col=NA, chr.col=NA )
{
 # calculate chromosome-wise position and ID for each SNP in list()
 # dir should be contain the annotation directory assumed to contain the snp list (dir.ano)
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 # read in gene annotation (vcf file/bim file) if not passed as fn. arg.
 # sort info by chromosome and position
 cat("\nRetrieving SNP information (chr,pos,id) from",mode,"file\n")
 loadvcf <- F
 if(is.null(dim(vcf.file))) { if(vcf.file[1]=="load") { loadvcf <- T  } }
 if(loadvcf) {
   if(file.exists(vcf.fn)) {
     del <- derive.file.delim(vcf.fn)
     vcf.file <- switch(mode,
          bim=read.delim(vcf.fn,comment.char="#",header=F,stringsAsFactors=F,sep=del),
          vcf=read.delim(vcf.fn,comment.char="#",header=T,stringsAsFactors=F,sep=del),
          map=read.delim(vcf.fn,header=F,stringsAsFactors=F,sep=del),
          map3=read.delim(vcf.fn,header=F,stringsAsFactors=F,sep=del),           
          )
     cat(" file preview:\n"); print(head(vcf.file,4)); cat("\n")
   } else {
     linez <- character()
     linez[1] <- (paste("Error: expecting file: '",vcf.fn," from function parameter 'vcf.fn'.",sep=""))
     linez[2] <- ("This file should be a vcf, bim, map, or map3 file (see plink website")
     linez[3] <- ("for description of map and map3 formats). 'Map3' is the simplest with")
     linez[4] <- ("chromosome number in column 1, snp id in column 2, snp position in column 3.")
     linez[5] <- ("Can be space, tab or comma delimited.")
     cat(linez,"\n"); stop()
   }
 }
 # if no user value entered, set column with SNP labels to default for given mode
 if(is.na(snp.col)) { snp.col <- switch(mode,bim=2,vcf=2,map=2,map3=2) }
 # if no user value entered, set column with SNP positions to default for given mode
 if(is.na(pos.col)) { pos.col <- switch(mode,bim=4,vcf=4,map=4,map3=3) }
 if(is.na(chr.col)) { chr.col <- 1 } # (these main file types all have chr in col 1)
 cat(" assuming columns are:  SNP-id:",snp.col,"; Chr:",chr.col,"; Pos:",pos.col,"\n")
 cat(" if this does not match file preview above, please stop and change file 'mode', or set")
 cat(" the values of:\n snp.col, pos.col, chr.col, explicitly in functions passing args to 'calc.chr.ind'.\n")
 if(is.null(snp.fl) | length(snp.fl)>10)
 {
   if(is.null(snp.fl))
   {
     cat(" no subset of snps was selected, default is to use every snp in map/bim/vcf file\n")
     snp.list <- vcf.file[,snp.col]
     match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
   } else {
     cat(" using vector of snp IDs 'snp.fl' as subset\n")
     snp.list <- snp.fl
     match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
     pc.missing <- (length(which(is.na(match.list.to.vcf)))/length(snp.list))
     cat("",paste(round(100*pc.missing,1),"% snps missing from annotation file\n"))
   }
 } else {
   if(!snp.fl %in% list.files(dir$ano)) 
   {
     linez <- character()
     linez[1] <- (paste("Error: expecting file: '",snp.fl,"' from parameter 'snp.fl' in ",dir$ano,sep=""))
     linez[2] <- ("This file should be a list of snp ids (one per line) to include in the current process.")
     linez[3] <- ("This could be a list of all snps, or any subset of snps in the map/bim/vcf file.")
     linez[4] <- ("Alternatively pass in a character() list of ids, or NULL to include all in the map/bim/vcf file.")
     cat(linez,"\n") ; stop()
   } else {
     snp.list <- readLines(paste(dir$ano,snp.fl,sep=""))
   }
   # match snp list to vcf file
   match.list.to.vcf <- match(snp.list,vcf.file[,snp.col])
   pc.missing <- (length(which(is.na(match.list.to.vcf)))/length(snp.list))
   cat("",paste(round(100*pc.missing,1),"% snps missing from annotation file\n"))
   if(pc.missing>.5) { stop("Error: too many missing. comment out [#] line in function 'calc.chr.ind' if intentional")}
 }
 chrnums <- chr.lab.to.num(vcf.file[,chr.col])
 chrnum <- chrnums[match.list.to.vcf]
 num.chr <- length(unique(chrnum))
 chrpos <- vcf.file[,pos.col][match.list.to.vcf]
 # reorder by genome position
 pos.order <- order(chrpos)
 chrnum <- chrnum[pos.order]
 chrpos <- chrpos[pos.order]
 snps.ordered <- snp.list[pos.order]
 chr.snp.lists <- tapply(snps.ordered,factor(chrnum),list)
 chr.pos.lists <- tapply(chrpos,factor(chrnum),list)

 CHR.INFO <- make.chr.obj(chrIDs=chr.snp.lists,chrIndx=chr.pos.lists)
 if(sav) {
   ofn <- paste(dir$ano,sav.fn,sep="")
   save(CHR.INFO, file=ofn)
   cat(paste(" chromosome positions file saved to:\n",ofn,"\n"))
   return(NULL)
 } else {
   cat(" chromosome positions return as list-object (not saved to file)\n")
   return(CHR.INFO)
 }
}


sample.bad.count.table <- function(dir,sub.list,refresh=T,type=1,extralist=NA,addExtraTo="Mean")
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason for each plate
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 nsamp <- length(sub.list)
 if(type==1) {
   new.fr <- data.frame(CallRate=integer(nsamp),
                        Mean=integer(nsamp),DLRS=integer(nsamp),ChrAb=integer(nsamp),
                        GCWave=integer(nsamp),TOTAL=integer(nsamp))  #StDev=integer(nsamp),
 } else {
   new.fr <- data.frame(Mean=integer(nsamp),
                        DLRS=integer(nsamp),ChrAb=integer(nsamp),
                        GCWave=integer(nsamp),TOTAL=integer(nsamp)) #StDev=integer(nsamp),
 }
 rownames(new.fr) <- sub.list

 if (!is.na(extralist[1]) & addExtraTo %in% colnames(new.fr))
 { extraCol <- match(addExtraTo,colnames(new.fr)) } else { extraCol <- 0 }
 excl.samp.dir <- paste(dir$ano,"SAMPLE_EXCLUDE/",sep="")
 dirfilz <- list.files(excl.samp.dir)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- agrep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   if (length(nxt.col)>0) {
     next.bad[[cc]] <- readLines(paste(excl.samp.dir,dirfilz[cc],sep=""))
     if (extraCol==nxt.col) { next.bad[[cc]] <- c(next.bad[[cc]],extralist) }
     bad.cnts <- (table(next.bad[[cc]],exclude=NA,useNA="no"))
     new.fr[names(bad.cnts),nxt.col] <- as.vector(bad.cnts)
   }
 }
 TOTAL.bad <- unique(unlist(next.bad))
 bad.cnts <- (table(TOTAL.bad,exclude=NA,useNA="no"))
 new.fr[names(bad.cnts),ncol(new.fr)] <- as.vector(bad.cnts)

 return(new.fr)
}



list.qc.fail.types <- function(dir,filt.list)
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 ns <- 1
 new.fr <- data.frame(CallRate=integer(ns),
                        Mean=integer(ns),DLRS=integer(ns),ChrAb=integer(ns),
                        GCWave=integer(ns),TOTAL=integer(ns))  #StDev=integer(ns),
 excl.samp.dir <- paste(dir$ano,"SAMPLE_EXCLUDE/",sep="")
 dirfilz <- list.files(excl.samp.dir)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- agrep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   next.bad[[cc]] <- readLines(paste(excl.samp.dir,dirfilz[cc],sep=""))
   if (length(which(!next.bad[[cc]] %in% filt.list))>0) {
     next.bad[[cc]] <- next.bad[[cc]][-which(!next.bad[[cc]] %in% filt.list)]
   } 
   bad.cnts <- length(next.bad[[cc]])
   new.fr[1,nxt.col] <- bad.cnts
 }
 TOTAL.bad <- unique(unlist(next.bad))
 bad.cnts <- length(TOTAL.bad)
 new.fr[1,ncol(new.fr)] <- bad.cnts

 return(new.fr)
}


update.plate.bad.count.table <- function(dir,plate.index=NULL,plate.list="plate.list.txt",refresh=T,filt.list=NULL)
{
 # look at all sample exclusion files and make table
 # of how many excluded for each reason for each plate
 # plate.index should contain sample ids in column 1, plate ids in column 2
 # would be easy to run for wells by switching cols 2&3
 dir <- validate.dir.for(dir,c("ano"),warn=F)
 if (!is.null(filt.list)) {
   fl <- find.id.col(plate.index,filt.list,ret="index")
   plate.index.sub <- plate.index[fl$index,]
 } else {
   plate.index.sub <- plate.index
 }
 if (!refresh & ("plate.stats.txt" %in% list.files(dir$ano))) {
   new.fr <- read.delim(file=paste(dir$ano,"plate.stats.txt",sep=""))
   plts <- paste(new.fr[,1])
 } else { 
   if(length(plate.list)==1) {
     if(plate.list %in% list.files(dir$ano))
     {
       plat.fn <- paste(dir$ano,plate.list,sep="")
       plts  <- readLines(plat.fn)
     } else {
       plts <- paste(unique(plate.index[,2])) # plate ids should be in col 2
     }
   } else {
     plts <- plate.list # assume this is the list of plates
   }
   npl <- length(plts)
   new.fr <- data.frame(ID=plts,CallRate=integer(npl),
                        Mean=integer(npl),DLRS=integer(npl),ChrAb=integer(npl),
                         GCWave=integer(npl),TOTAL=integer(npl)) 
   if(!is.null(filt.list) & !is.null(plate.index))
   {
   	per.plate.cnt <- table(plate.index.sub[,2])
   	new.fr[["SIZE"]] <- per.plate.cnt[paste(plts)]
   }
 }
 dirfilz <- list.files(dir$excl)
 to.search <- gsub(".txt","",dirfilz,fixed=T)
 next.bad <- list()
 for (cc in 1:length(to.search))
 {
   nxt.col <- grep(toupper(to.search[cc]),toupper(colnames(new.fr)))
   if(length(nxt.col)==1) {
	   next.bad[[cc]] <- readLines(paste(dir$excl,dirfilz[cc],sep=""))
	   plt.bad <- (plate.index[(plate.index[,1] %in% next.bad[[cc]]),2])
	   plt.bad2 <- (plate.index.sub[(plate.index.sub[,1] %in% next.bad[[cc]]),2])
   	 bad.cnts <- (table(plt.bad))
     if(length(bad.cnts)!=0) {
   	   new.fr[match(names(bad.cnts),new.fr[,1]),nxt.col] <- as.vector(bad.cnts)
     }
   	 cat(dirfilz[cc]," ",sum(table(plt.bad2)),"/",length(next.bad[[cc]])," in this cohort\n",sep="")
	} else {
		if(length(nxt.col)>1) { 
			stop("Error: more than 1 file matched header. check files") 
		} else {
			cat(paste("Warning: file",dirfilz[cc],"skipped in counts because not in table header\n"))
			next.bad[[cc]] <- ""
		}
	}
   if (!is.null(filt.list) & nchar(next.bad[[cc]][1])>1) {
     excl <- which(!next.bad[[cc]] %in% filt.list)
     if(length(excl)>1) {
       next.bad[[cc]] <- next.bad[[cc]][-excl]  }
   }
 }
 TOTAL.bad <- unique(unlist(next.bad))
 plt.bad <- (plate.index[(plate.index[,1] %in% TOTAL.bad),2])
 plt.bad <- plt.bad[plt.bad!=""]
 bad.cnts <- (table(plt.bad))
 new.fr[match(names(bad.cnts),new.fr[,1]),"TOTAL"] <- as.vector(bad.cnts)

 return(new.fr)
}


plate.box.plots <- function(plot.stats,pass.tab=NULL,plate.lookup,pref,dir="")
{
  # make box plots across a cohort showing the distribution between plates
  dir <- validate.dir.for(dir,c("qc.pl"),warn=F)
  if(!(all(c("id","plate") %in% colnames(plate.lookup))))
    { stop("plate.lookup must contain 'plate' and 'id' columns")}
  if(is.matrix(plate.lookup)) { plate.lookup <- as.data.frame(plate.lookup) }
  plts <- paste(unique(plate.lookup$plate))
  pltnums <- 1:length(plts) #match(substr(rownames(plot.stats),1,6),plts)
  samp.plts <- match(plate.lookup$plate,plts)
  plt.samp <- samp.plts[match(rownames(plot.stats),plate.lookup$id)]
  statz <- colnames(plot.stats) # e.g, c("Mean","DLRS","GCWave")
  ylz <- list( c(-.2,.1), c(0,.6) , c(-.2,.1) )
  spotz <- c("topleft","topright")
  xllz <- c(-.2,.1,-.2); xlhz <- c(.1,.5,1)
  
  ofn <- paste(dir$qc.pl,"PlateBoxPlotsMDG",pref,".pdf",sep="")
  pdf(ofn)
  
  for (tt in 1:length(statz))
  { 
    stat <- paste(statz[tt])
    yl <- ylz[[tt]]
    spott <- spotz[tt]
    boxplot(plot.stats[[stat]]~factor(plt.samp),xlab="plate number",
            ylab=paste("LRR",stat), lwd=1, pch=".", bty="l",
            main=paste("LRR",stat,"distribution"),ylim=yl,bty="l")
    # use segments here instead of abline to avoid overshadowing some boxes
    if(!is.null(pass.tab)) {
      abline(h=pass.tab[c("LB","UB"),stat],col="lightblue") #,lty="dotted")
      legend("top",legend=c("Lower and upper bounds"),
           lty="solid",col="lightblue",bty="n")
    }
  }
  dev.off()
  
  print(paste("produced file:",ofn))
}


three.way.comparison <- function(clean.stats,sample.info,raw.fn="StatsPerSample",dir,
                                 batch.comps=c("plate","grp"),bxplot=T,pref="auto",...) 
{
  ## function assuming data is PC-clean and raw stats per sample summary files are
  # in the default location with default names. if so, makes nice visualisation
  dir <- validate.dir.for(dir,c("qc.pc"),warn=F)
  grpnumz <- as.numeric(unique(sample.info$grp))
  raw.stat.table <- combine.raw.stats(grpnumz,dir,base.fn=raw.fn,pref=pref)
  sampqc.stat.table <- remove.qc.failers(raw.stat.table,sample.info)
  three.part.comp <- list(raw.stat.table,sampqc.stat.table,clean.stats)
  names(three.part.comp) <- c("Raw","Sample.QC","PC.Corrected")
  vars <- colnames(clean.stats)
  vars <- vars[vars %in% colnames(raw.stat.table)]
  if(length(vars)<2) { return(three.part.comp) } # don't bother w/ plot if not at least 2 stats 
  raw.stat.table <- raw.stat.table[,match(vars,colnames(raw.stat.table))]
  sampqc.stat.table <- sampqc.stat.table[,match(vars,colnames(sampqc.stat.table))]
  
  if(bxplot) {
    for (dd in 1:length(batch.comps)) {
      ofn <- dir.fn.cmb(dir$qc.pc,"ThreeStageComparison",suf=batch.comps[dd],ext=".pdf")
      pdf(ofn,...)
      par(mfcol=c(ncol(clean.stats),length(three.part.comp)))
      # some batch effects diagnostic plots to see how well the correction has worked
      for (cc in 1:length(three.part.comp)) {
        batch.box.plots(three.part.comp[[cc]],NULL,sample.info,batch=batch.comps[dd],
                        pref=paste(batch.comps[dd],names(three.part.comp)[cc],sep="_"),
                        subtitle=paste("[",names(three.part.comp)[cc],"]"),dir=dir$qc.pc,to.file=F)
      }
      dev.off()
      cat("Wrote plots to:",ofn,"\n")
    }
  }
  return(three.part.comp)
}


combine.raw.stats <- function(grpnumz,dir,base.fn="StatsPerSample",pref="auto") 
{
  # combine set of sample statistics from multiple files (cohorts) into 1 big file
  # or in the case of length grpnumz = 1, simply return the only file.
  file.info <- get.file.specs(dir)
  if(pref=="auto") {
    ## try to automatically deduce the prefix for the raw stats file(s)
    if(length(unique(file.info$GRP))==1) {
      pref <- "LRR"
    } else {
      pref <- paste("LRR",grpnumz,sep="")
    }
  } 
  dir <- validate.dir.for(dir,c("qc.lrr"),warn=F)
  if(length(grpnumz)==1) {
    ifn <- dir.fn.cmb(dir$qc.lrr,base.fn,pref=pref,ext=".tab")
    ifn2 <- dir.fn.cmb(dir$qc.lrr,base.fn,pref=paste(pref,"1",sep=""),ext=".tab")
    if(!file.exists(ifn)) { 
      if(file.exists(ifn2)) {
        ifn <- ifn2
      } else {
        in.dir <- list.files(dir$qc.lrr)
        tryanything <- grep(base.fn,in.dir,fixed=T)
        if(length(tryanything)>0) {
          cat("Warning: Didn't find expected sample-wise statistics text file in",dir$qc.lrr,"\n")
          cat("but did find:",in.dir[tryanything][1],"\nattempting to use as stats file. \n")
          cat("stats file with 1 grp should have been named:\n",ifn,"\n")
          ifn <- dir.fn.cmb(dir$qc.lrr,in.dir[tryanything[1]])
        } else {
          Emsg=paste("Error: could not find any sample-wise stats file named *",base.fn,"* in:\n",dir$qc.lrr,"\n",sep="")
          stop(Emsg)
        }
      }
    }
    stat.table <- read.table(ifn)
  } else {
    stat.table <- list()
    for (kk in grpnumz) {
     ifn <- dir.fn.cmb(dir$qc.lrr,base.fn,pref=pref[kk],ext=".tab")
      if(file.exists(ifn)) {
        stat.table[[kk]] <- read.table(ifn)
      } else {
        cat("Warning: Didn't find expected samplewise stats file:",ifn,"\n")
      }
    }
    stat.table <- do.call("rbind",stat.table)
  }
  return(stat.table)
}


batch.box.plots <- function(plot.stats,pass.tab=NULL,lookup,batch="plate",pref="All",subtitle="",dir="",to.file=T)
{
  ## plot box plot of batch effect, using stats summary and looking up batch categories from table
  batch <- paste(batch[1])
  index <- find.id.col(lookup,ids=rownames(plot.stats))$index
  if(!(all(c(batch) %in% colnames(lookup))))
  { stop(paste("'lookup' expected to contain '",batch,"' and 'id' named columns"),sep="")}
  if(is.matrix(lookup)) { lookup <- as.data.frame(lookup) }
  plts <- paste(unique(lookup[,batch]))
  pltnums <- 1:length(plts) #match(substr(rownames(plot.stats),1,6),plts)
  samp.plts <- match(lookup[,batch],plts)
  plt.samp <- samp.plts[index]  #match(rownames(plot.stats),lookup[,"id"])]
  statz <- colnames(plot.stats) # e.g, c("Mean","DLRS","GCWave")
  ylz <- list( c(-.2,.1), c(0,.6) , c(-.2,.1) )
  spotz <- c("topleft","topright")
  xllz <- c(-.2,.1,-.2); xlhz <- c(.1,.5,1)
  
  if(to.file) {
    dir <- validate.dir.for(dir,c("qc.pl"),warn=F)
    ofn <- paste(dir$qc.pl,batch,"BoxPlotsMDG",pref,".pdf",sep="")
    pdf(ofn)
  }
  
  for (tt in 1:length(statz))
  { 
    stat <- paste(statz[tt])
    yl <- ylz[[tt]]
    spott <- spotz[tt]
    boxplot(plot.stats[[stat]]~factor(plt.samp),xlab=paste(batch,"number"),
            ylab=paste("LRR",stat), lwd=1, pch=".", bty="l",
            main=paste("LRR",stat,"distribution"),sub=subtitle,ylim=yl,bty="l")
    # use segments here instead of abline to avoid overshadowing some boxes
    if(!is.null(pass.tab)) {
      abline(h=pass.tab[c("LB","UB"),stat],col="lightblue") #,lty="dotted")
      legend("top",legend=c("Lower and upper bounds"),
             lty="solid",col="lightblue",bty="n")
    }
  }
  if(to.file) {  dev.off() ; cat(paste("produced file:",ofn,"\n")) }
}


get.extreme.examples <- function(stat.table,failures,venn.lists)
{
 # find most extreme subject for each part of venn list
 # [NB: or for 'none' group - the least extreme ]
 last.n <- length(venn.lists)
 z.scores <- apply(stat.table,2,StandardizeX)
 sample.dev.scores <- rowSums(abs(z.scores)*(failures)) # scores to get MOST extreme
 none.grp <- match(venn.lists[[last.n]],names(sample.dev.scores)) # group with no outliers
 sample.dev.scores[none.grp] <- rowSums(abs(z.scores))[none.grp] # scores to get LEAST extreme

 ex.id <- integer(last.n)
 # find most extreme subject for each part of venn list

 for (cc in 1:(last.n))
 {
   samp.nos <- match(venn.lists[[cc]],names(sample.dev.scores)) 
   # find max for 1-8, min for 9 [none group]
   if (length(samp.nos)>0) {
     if(cc<length(venn.lists))
     {
       mx <- which(sample.dev.scores[samp.nos]==max(sample.dev.scores[samp.nos],na.rm=T))
     } else {
       mx <- which(sample.dev.scores[samp.nos]==min(sample.dev.scores[samp.nos],na.rm=T))
     }
     ex.id[cc] <- names(sample.dev.scores)[samp.nos[mx[1]]]
   } else { ex.id[cc] <- NA }
 }
 names(ex.id) <- names(venn.lists)   
 return(ex.id)
}



rescale.dens <- function(dnsty,V,W,X=NULL,Y=NULL,ret=F,lim=T)
{
 # modify the units of a density object so the plot can be put into a
 # plot area defined by parameters V,W,X,Y
 if(is.null(X) ) { X <- range(dnsty$x) }
 if(is.null(Y) ) { Y <- range(dnsty$y) }
 dnsty$x <- ((V[2]-V[1])/(X[2]-X[1]))*(dnsty$x - X[1]) + V[1]
 dnsty$y <- ((W[2]-W[1])/(Y[2]-Y[1]))*(dnsty$y - Y[1]) + W[1]
 if (lim)
 { 
   dnsty$x[dnsty$x < V[1]] <- V[1]
   dnsty$x[dnsty$x > V[2]] <- V[2]
   dnsty$y[dnsty$y < W[1]] <- W[1]
   dnsty$y[dnsty$y > W[2]] <- W[2]
 }
 if (ret)
 { return(list(dnsty,X,Y)) } else { return(dnsty) }
}

col.extremes <- function(lat,selBlue="extr",selRed="extr",sdc=2)
{
 #if ="extr" then colour extreme values in a (pre) latex object
 #otherwise colour red/blue according to logical vars selBlue, selRed
 lat.num <- as.numeric(lat); dim(lat.num) <- dim(lat)
 colsds <- rep(apply(lat.num,2,sd,na.rm=T),each=nrow(lat))
 colmns <- rep(apply(lat.num,2,mean,na.rm=T),each=nrow(lat))
 dim(colsds) <- dim(colmns) <- dim(lat)
 Ff <- !is.na(lat.num)
 if (!is.logical(selBlue) ) {
   selBlue <- which(lat.num[Ff]<(colmns[Ff]-(sdc*colsds[Ff]))) }
 if (!is.logical(selRed)) {
   selRed <- which(lat.num[Ff]>=(colmns[Ff]+(sdc*colsds[Ff]))) }
 lat[Ff][selBlue] <- paste("\\textcolor{blue}{",lat[selBlue],"}",sep="")
 lat[Ff][selRed] <-  paste("\\textcolor{red}{",lat[selRed],"}",sep="")
 return(lat)
}

StandardizeX <- function(x)
{
 # standardize some data according to itself
 u <- mean(x,na.rm=T)
 s <- sd(x,na.rm=T)
 (x-u)/s
}


get.pctile <- function(dat,pc=0.01)
{
 # get top/bottom percentile 'pc' subset of data 'dat'.
 rr <- rank(dat,na.last=NA)
 tpbt <- round(c(pc,1-pc)*max(rr,na.rm=T))
 ord <- dat[!is.na(dat)][order(dat[!is.na(dat)])]
 if(tpbt[1]==0) { tpbt[1] <- 1 }
 pcts <- ord[tpbt]
 return(pcts)
}


MedianScale <- function(X, n=100, meth="median",ratio=10)
{
 ## smooth data series by median / mean
 rem <- length(X) %% n
 YorN <- (rem != 0) 
 n.chunks <- (length(X) %/% n) + as.integer(YorN)
 if (YorN) {  X <- c(X,rep(NA,times=(n-rem))) }
 Xmat <- matrix(X,nrow=n)
 if(meth=="median") {
   meds.out <- apply(Xmat,2,median,na.rm=T) 
 } else {
   meds.out <- apply(Xmat,2,mean,na.rm=T)
 }
 meds.out <- rep(meds.out,each=n)[1:length(X)]
 meds.out <- ((ratio*meds.out)+X)/(ratio+1)
 return(meds.out)
}


initialise.excl.files <- function(dir,reset.files=c("ChrAb.txt","DLRS.txt","GCWave.txt","Mean.txt","BadPlates.txt"))
{
  # delete existing sample exclusion file if present
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  excl.dir <- paste(dir$ano,"SAMPLE_EXCLUDE/",sep="")
  cur.files <- list.files(excl.dir)
  for (cc in 1:length(reset.files))
  {
    if(reset.files[cc] %in% cur.files)
    {
      cat(" deleting file:",reset.files[cc],"for initialisation purposes\n")
      unlink(paste(excl.dir,reset.files[cc],sep=""))
    }
  }
}

big.exclude.sort <- function(des.fn="LRRdescrFile", dir="", deepC=T, tranMode=2, pref="LRRFilt",
                             f.snp="", f.samp="", verbose=T)
{
  # sort and exclude snps/samples from a big.matrix
  must.use.package("bigmemory")
  dir <- validate.dir.for(dir,c("big","ano"),warn=F)
  if (tranMode==1) {  cat("\nFiltering samples/snps with SNP-QC failures:\n")
  } else { cat("\nExcluding samples/snps with SNP-QC failures, sorting SNPs by chr, pos:\n") }
  
  # bigmatrix file names for re-ordered filtered matrix (which is the final output of this script)
  bck.fn.o <- paste(pref,"Sort","bckfile",sep="")
  des.fn.o <- paste(pref,"Sort","descrFile",sep="")
  ofn <- dir.fn.cmb(dir$big,des.fn.o,ext=".RData")
  
  bigMat <- getBigMat(des.fn,dir)
  cat(paste(" attached matrix with dims:",paste(dim(bigMat),collapse=","),"\n"))
  # get list of deleting/reordering vectors using annotation files
  if (tranMode==1) {
    #if(f.snp=="") { sif1 <- F } else { sif1 <- T }
    #if(f.samp=="") { sif2 <- F } else { sif2 <- T }
    #snp.L <- check.file.and.subset(fnm=paste(dir$ano,f.snp,sep=""), bigMat, by.row=(rowsAre=="SNP"), stop.if.fail=sif1)
    #samp.L <- check.file.and.subset(fnm=paste(dir$ano,f.samp,sep=""), bigMat, by.row=(rowsAre=="SAMP"), stop.if.fail=sif2)
    trans.list <- select.samp.snp.custom(bigMat,snp=f.snp,samp=f.samp)
    cat(paste(" selected",length(trans.list[[4]]),"listed samples and",length(trans.list[[3]]),"Snps\n"))
  } else {  
    cat(" excluding and ordering listed samples based on annotation directory files; SNP/SAMP\n")
    trans.list <- sort.exclude.from.annot(bigMat,dir=dir,dosnp=T,dosamp=T,ordR=T,ordC=T, verb=verbose)
    ###sort.exclude.from.annot.old(bigMat,snp=T,samp=F,dir=dir$ano,crOff=F)  - original version!
  }
  wrn <- "warning: trans.list was already attached, detaching now..\n"
  while("trans.list" %in% search()) { detach(trans.list); cat(wrn) }
  attach(trans.list)
  if(verbose) {
    cat("\nReordering SNPs and Samples...\n")
    cat("\nINDEXES SUMMARY\n")
    cat(paste(length(to.order.r),"row indexes range is from",min(to.order.r),"to",max(to.order.r),"\n"))
    cat("-->",head(to.order.r),sep=", "); cat ("\n")
    cat(paste(length(to.order.c),"col indexes range is from",min(to.order.c),"to",max(to.order.c),"\n"))
    cat("-->",head(to.order.c),sep=", ")
    cat("\n\n raw big.matrix summary before ordering and exclusion based on SNP-QC:\n\n")
    printBigSummary(bigMat,"bigMat")
  }
  if(!deepC)
  {
    # this is fast with available RAM (like 20 secs)
    cat(" running reorder in system memory\n")
    system.time(bigMat1 <- bigMat[to.order.r,to.order.c])
    if(colnames(bigMat1)[1]!=colnames(bigMat)[to.order.c[1]])
    { cat(" adding colnames\n") ; colnames(bigMat1) <- colnames(bigMat)[to.order.c] }
    if(rownames(bigMat1)[1]!=rownames(bigMat)[to.order.r[1]])
    { cat(" adding rownames\n") ; rownames(bigMat1) <- rownames(bigMat)[to.order.r] }
    cat(" converting matrix to big.matrix\n")
    bigMat2 <- as.big.matrix(bigMat1, backingfile=bck.fn.o,
                             backingpath=dir$big, descriptorfile=des.fn.o)
    cat(paste(" matrix descr saved as standard description file:",des.fn.o,"\n"))
    descr <- describe(bigMat2)
  } else {
    #this is slow but creates backing file and will speed up ops later
    cat(" starting deep copy...")
    bigMat2 <- deepcopy(bigMat, cols = to.order.c, rows = to.order.r,
                        backingfile=bck.fn.o,backingpath=dir$big, descriptorfile=des.fn.o )
    cat("done\n")
    cat("\nAdding names\n")
    options(bigmemory.allow.dimnames=TRUE)
    colnames(bigMat2) <- colnames(bigMat)[to.order.c]
    cat(" added colnames\n")
    rownames(bigMat2) <- rownames(bigMat)[to.order.r]  
    cat(" added rownames\n")
    descr <- describe(bigMat2)
    flush(bigMat2) # hopefully this will ensure the row/colnames are added to the file backing
    cat(paste(" due to use of deep copy option, recommend only to use descr saved as rbinary description file\n"))
  }
  cat(paste(" created big.matrix description file:",des.fn.o,"\n"))
  cat(paste(" created big.matrix backing file:",bck.fn.o,"\n"))
  
  save(descr,file=ofn)
  cat(paste(" created big.matrix binary description file:",rmv.dir(ofn),"\n"))
  
  while("trans.list" %in% search()) { detach(trans.list) }
  return(descr) 
}


head.list <- function(alist,n=6)
{
  ## uses recursion if necessary to 
  # sensibly apply the 'head' function
  # to a list which may have sublists
  for(cc in 1:length(alist))
  {
    if(!is.null(names(alist))) 
    { cat("$",names(alist)[[cc]],":\n",sep="") } else { (cat("[[",cc,"]]\n",sep="")) }
    if (is(alist[[cc]],"list"))
    {
      head.list(alist[[cc]],n)
    } else {
      print(head(alist[[cc]],n))
    }
  }
}


plain.vec.to.big.matrix <- function(dir, grp=NULL, snp.fn="snpNames.txt", sample.fn=NULL, input.fn=NULL, pref="LRR", input.is.vec=T, delete.existing=T)
{
  # import from a text (hopefully long format) datafile to a big.matrix
  max.mem <- 4096 # stupidly large
  dat.file.suf <- "LRR.dat" ## note this must match initial bash script that extracts data!
  ## Define data types for big.matrix
  dat.type <- double(1)
  dat.typeL <- "double"  # label
  dir <- validate.dir.for(dir,c("ano","big","col"),warn=F)
  #
  if(all(!is.null(sample.fn))) {
    sample.fn <- dir.fn.cmb(dir$ids,sample.fn,must.exist=T)  # sample id file name
    cat("Reading sample and snp lists from text files...\n")
    cat(paste("Reading samples from",sample.fn,"\n"))
    ID.list <- lapply(sample.fn,readLines)  #readLines(sample.fn)
  } else {
    ID.list <- get.subIDs(dir,"group")
  }
  if(!is.null(grp)) { 
    ID.list <- ID.list[[grp]]
  } 
  cmb.ID.list <- paste(do.call("c",ID.list))
  if(all(is.null(input.fn))) {
    input.fn <- get.lrr.files(dir,grp,suffix=dat.file.suf)    
    ifn <- dir.fn.cmb(dir$col,input.fn,must.exist=T)
  }
  snp.fn <- dir.fn.cmb(dir$ano,snp.fn,must.exist=T)  # snp annotation file name
  numfls <- length(ID.list)
  #print(head(cmb.ID.list))
  cat(paste(" reading snps from",snp.fn,"\n"))
  gene.list <- readLines(snp.fn)
  print(head(gene.list))
  if(length(input.fn)>1) {
    if(length(input.fn)==numfls) {
      cat("Warning: reading a single cohort from",numfls,"source files. Edit file.spec.txt if this is unexpected\n")
    } else {
      stop("Error: when reading a single cohort from multiple source files, need same number of id files")
    }
  }
  num.sub <- length(cmb.ID.list) #ID.list)
  smp.szs <- sapply(ID.list,length)
  fil.ofs <- c(0,cumsum(smp.szs)) #note last element is the end of the last file
  num.gen <- length(gene.list)
  cat(paste("found",num.sub,"samples and",num.gen,"markers\n"))

  cls <- num.sub; rws <- num.gen
  cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
  memory.estimate <- as.double((as.double(rws)*as.double(cls))/cells.per.gb)
  if (memory.estimate > max.mem) {
    cat("Error: Insufficient disk space availability expected for this import method\n")
    cat("Please free up some disk space and try again\n")
    stop()
  }
  # use 'pref' as the name of the big.matrix backing files for this cohort
  bck.fn <- paste(pref,"bckfile",sep="")
  des.fn <- paste(pref,"descrFile",sep="")
  if ((!des.fn %in% list.files(dir$big)) | delete.existing )
  {
    if(delete.existing & (des.fn %in% list.files(dir$big)))
    {
      dfn <- paste(dir$big,des.fn,sep="")
      cat("deleting",dfn,"\n")
      unlink(dfn)
    } else {
      #all clear, no files already exist with same name
    }
  } else {
    cat(paste("Error: Big matrix description file",des.fn,"already exists in",dir$big,"\n"))
    cat("Please delete, rename or move this file, or use option DELETE=1, before re-running this script\n")
    stop()
  }
  cat("\nCreating big matrix object...")
  cat("\n predicted disk use: ",round(memory.estimate,1),"GB\n")
  bigVar <- big.matrix(nrow=rws,ncol=cls, backingfile=bck.fn, dimnames=list(gene.list,cmb.ID.list),
                       type=dat.typeL, backingpath=dir$big, descriptorfile=des.fn)
  for(ff in 1:numfls) {
    ##ifn <- dir.fn.cmb(dir$col,input.fn[ff],must.exist=T)
    dat.file <- file(ifn[ff])
    open(con=dat.file,open="r")
    cat(paste(" opening connection to ",c("matrix","long")[1+input.is.vec],
              " format datafile (",ff,"/",numfls,"):",rmv.dir(ifn[ff]),"\n",sep=""))
    cat("\nLoading text data into big matrix object:\n")
    nxt.rng <- (fil.ofs[ff]+1):(fil.ofs[ff+1])
    cls <- length(nxt.rng)
    if(!input.is.vec)
    {
      ## read from matrix format tab file
      twty.pc <- round(rws/5)
      for (cc in 1:rws) {
        if ((cc %% twty.pc)==0)  { fl.suc <- flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
        loop.tracker(cc,rws)
        next.line <- readLines(dat.file,n=1)
        next.row <- strsplit(next.line,"\t",fixed=T)[[1]]
        if (cc==1) { if (length(next.row)!=cls) {
          stop("dimensions of import file do not match gene and snp list spec, exiting")
          break; break; } }
        bigVar[cc,nxt.rng] <- next.row
      }
    } else {
      ## read from (long) vector format tab file
      twty.pc <- round(cls/5)
      for (cc in 1:cls) {
        loop.tracker(cc,cls)
        if ((cc %% twty.pc)==0)  { fl.suc <- flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
        bigVar[,(cc+fil.ofs[ff])] <- as(readLines(dat.file,n=rws),dat.typeL)
      }
    }
    close(dat.file)
  }
  cat("\n")
  cat(paste("created big.matrix description file:",des.fn,"\n"))
  cat(paste("created big.matrix backing file:",bck.fn,"\n"))
  return(describe(bigVar))
  cat("...complete!\n")
}


dir.force.slash <- function(dir) {
  # make sure 'dir' directory specification ends in a / character
  the.test <- (dir!="" & substr(dir,nchar(dir),nchar(dir))!="/")
  dir[the.test] <- paste(dir[the.test],"/",sep="")
  return(dir)
}


def.dirs.fn <- function(dir.raw="LRRQC_Raw_Files",
                        dir.sup="LRRQC_Support_Files",
                        dir.base="plumbCNV_LRRQC")
{
  #create a directory structure for plumbCNV to work upon starting with a rawdata and base directory
  #raw data:  dir.raw <- "/ipswich/data/Immunochip/FinalReports/"
  #location of support files directory dir.sup <- "/chiswick/data/store/metabochip/PLINK/"
  #location to write summaries of call rate: dir.base <- "/chiswick/data/ncooper/ImmunochipReplication/"
  #make sure each ends in a slash
  dir.raw <- dir.force.slash(dir.raw)
  dir.sup <- dir.force.slash(dir.sup)
  dir.base <- dir.force.slash(dir.base)
  
  dir.cr <- paste(dir.base,"CALLRATES/",sep="")
  dir.cr.plk <- paste(dir.cr,"Plink/",sep="")
  #location of big matrix main data files (contains LRR data for all samples/snps)
  dir.big <- paste(dir.base,"BigMatrixFiles/",sep="")
  #location of annotation files directory
  dir.ano <- paste(dir.base,"ANNOTATION/",sep="")
  
  #qc directory
  dir.qc <- paste(dir.base,"LRRQC/",sep="")
  dir.qc.lrr <- paste(dir.qc,"Mean_DLRS_GC/",sep="")
  dir.qc.cab <- paste(dir.qc,"ChrAb/",sep="")
  dir.qc.gc <- paste(dir.qc,"GCWave/",sep="")
  dir.qc.pl <- paste(dir.qc,"Plates/",sep="")
  dir.qc.pc <- paste(dir.qc,"PCComparison/",sep="")
  dir.qc.qc <- paste(dir.qc,"PassQC/",sep="")
  dir.ind <- paste(dir.qc,"Examples/",sep="")
  
  #scripts directories
  dir.scr <- paste(dir.base,"Scripts/",sep="")
  dir.ss <- paste(dir.scr,"LRRSubScripts/",sep="")
  dir.sa <- paste(dir.scr,"ChrAbSubScripts/",sep="")
  # sample exclude directory
  dir.excl <- paste(dir.ano,"SAMPLE_EXCLUDE/",sep="")
  dir.sort <- paste(dir.ano,"SAMPLE_SORT/",sep="")
  dir.excl2 <- paste(dir.ano,"SNP_EXCLUDE/",sep="")
  dir.sort2 <- paste(dir.ano,"SNP_SORT/",sep="")
  # principle components directories
  dir.pc <- paste(dir.base,"PC/",sep="")
  #penn cnv files
  dir.cnv <- paste(dir.base,"PENNCNV/",sep="")
  dir.lrr.dat <- paste(dir.base,"LRRDATA/",sep="")
  dir.baf <- paste(dir.base,"BAFDATA/",sep="")
  dir.col <- paste(dir.lrr.dat,"ColumnData/",sep="")
  dir.baf.col <- paste(dir.baf,"ColumnData/",sep="")
  dir.ids <- paste(dir.lrr.dat,"ColumnIds/",sep="")
  
  if ("dir" %in% ls()) { rm(dir) }
  all.locs <- ls()[grep("dir.",ls())]
  nloc <- length(all.locs)
  dir <- vector("list", nloc)
  for (cc in 1:nloc) {  dir[[cc]] <- get(all.locs[cc]) }
  names(dir) <- substr(all.locs,5,nchar(all.locs))
  return(dir)
}


init.dirs.fn <- function(dir,overwrite=F,ignore=c("raw","sup"),
                         silent=F,info.dir=NULL,file.spec="file.spec.txt",
                      plate.info="plate.lookup",pheno="pheno.lookup")
{
  # create the directory structure implied by the object 'dir'
  # provides options to delete existing directories if overwrite=T
  if(!is.list(dir)) { failnow <- T } else { failnow <- F }
  if(!all(sapply(lapply(dir,is),"[",1) %in% "character")) { failnow <- T }
  if(failnow) { stop("object 'dir' should be a list of directory locations") }
  ## don't try to create 'ignore' directories, ie., might be read-only, etc.
  if(!is.null(ignore)) {
    ignore <- ignore[ignore %in% names(dir)]
    for (dd in 1:length(ignore)) {
      dir[[paste(ignore[dd])]] <- NULL
    }
  }
  ## sort list by number of 'slash's, as this will mean that the
  # higher level directories will already exist when creating sub dirs
  countslash <- function(text) { length(text[text=="/"]) } 
  out1 <- sapply(dir,strsplit,split="")
  slashcounts <- sapply(out1,countslash)
  dirs.ordered <- unlist(dir)[order(slashcounts[match(names(dir),names(slashcounts))])]
  ldd <- length(dirs.ordered)
  passfail <- logical(ldd)
  ## iterate through creating each directory, taking note of whether it already exists
  # option to keep or erase contents if it does (eg, to create a fresh setup)
  for (cc in 1:ldd) {
    next.dir <- dirs.ordered[cc]
    if(!file.exists(next.dir)) {
      passfail[cc] <- dir.create(next.dir)
    } else {
      if(!silent) { cat("directory ",getwd(),"/",next.dir,"already existed, ",sep="") }
      if(overwrite) {
        cat("contains",length(list.files(next.dir)),"files and/or sub-directories\n")
        choicez <- c("I am sure - DELETE","DO NOT delete")
        choice <- select.list(choicez,preselect=choicez[2],title=paste("delete contents of",next.dir,"?"))
        if(choice==choicez[1]) {
          cat("attempt to overwrite...")
          unlink(next.dir,recursive=T)
          passfail[cc] <- dir.create(next.dir)
          cat(c("succeeded","failed")[2-as.numeric(passfail[cc])],"\n")
        } else {
          cat("Directory deletion cancelled by user\n")
        }
      } else {
        if(!silent) { cat("was left unmodified\n") }
      }
    }
  }
  if(!silent) {
    cat(length(which(passfail)),"of",length(passfail),"directories created successfully\n") }
  if(!is.null(info.dir)) {
    if(file.spec!="") {
      o_spcf <- dir.fn.cmb(info.dir,file.spec)
      n_spcf <- dir.fn.cmb(dir$lrr.dat,file.spec)
      if(file.exists(o_spcf)) {
        if(!file.exists(n_spcf)) {
          file.copy(from=o_spcf, to=n_spcf, overwrite = F, recursive = F, copy.mode = T)
          cat(" copied file",rmv.dir(o_spcf),"\ninto:",dir$lrr.dat,"\nfrom:",get.dir(o_spcf),"\n")
        } }
    }
    if(plate.info!="") {
      o_spcf <- dir.fn.cmb(info.dir,plate.info)
      n_spcf <- dir.fn.cmb(dir$ano,plate.info)
      if(file.exists(o_spcf)) {
        if(!file.exists(n_spcf)) {
          file.copy(from=o_spcf, to=n_spcf, overwrite = F, recursive = F, copy.mode = T)
          cat(" copied file",rmv.dir(o_spcf),"\ninto:",dir$lrr.dat,"\nfrom:",get.dir(o_spcf),"\n")
        } }
    }
    if(pheno!=""){
      o_spcf <- dir.fn.cmb(info.dir,pheno)
      n_spcf <- dir.fn.cmb(dir$ano,pheno)
      if(file.exists(o_spcf)) {
        if(!file.exists(n_spcf)) {
          file.copy(from=o_spcf, to=n_spcf, overwrite = F, recursive = F, copy.mode = T)
          cat(" copied file",rmv.dir(o_spcf),"\ninto:",dir$lrr.dat,"\nfrom:",get.dir(o_spcf),"\n")
        } }
    }
  }
  return(dirs.ordered)
}


## compiles certain functions if compiler package available
if(cmp.there){
  smart.match <- cmpfun(smart.match)
  rg <- cmpfun(rg)
  parse.fn <- cmpfun(parse.fn)
  #PC.fn <- cmpfun(PC.fn) #compiled version actually seems slower!
}

index.fn.file <- function(fn.in,fn.out="fn.index.txt") 
{
  # makes html index of each function in a large functions file
  grp <- function(what,ins) { grep(what,ins,fixed=T) }
  if(file.exists(fn.in))  {
    fl <- readLines(fn.in)
    fl <- rmv.spc(fl)
    fn.lines <- unique(c(grp("<- function",fl),grp("<-function",fl)))
    nfn <- length(fn.lines)
    fn.list <- vector("list",nfn)
    for (kk in 1:nfn) {
      first.ln <- fl[fn.lines[kk]]
      n <- 1; while(substr(first.ln,n,n)!="<" & substr(first.ln,n,n)!=" ") { n <- n+1 }
      fn.nm <- substr(first.ln,1,n-1)
      names(fn.list)[kk] <- paste("<p></p><b>",fn.nm,"</b>",sep=""); descr <- c()
      lnn <- fn.lines[kk]; while(length(grp("{",fl[lnn]))==0) { lnn <- lnn+1 }
      #print(fl[lnn])
      lnn <- lnn+1 ; 
      while(length(grp("#",fl[lnn]))>0) { 
        descr <- c(descr,gsub("#","",fl[lnn],fixed=T))
        lnn <- lnn+1 
      }
      fn.list[[kk]] <- rmv.spc(paste(descr))
    }
  } else {
    cat("Warning: could not find function file to index\n")
  }
  fil.cont <- sapply(fn.list,paste,collapse="\n")
  write.table(fil.cont,file=fn.out,quote=F,col.names=F)
  return(fn.list)
}


dummylst <- index.fn.file("/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsCNVAnalysis.R",
              "/chiswick/data/ncooper/ImmunochipReplication/Scripts/FunctionsIndex.htm")

rm(dummylst)
