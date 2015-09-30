####### INTERNAL FUNCTIONS
# Create network matrix from adjacency matrix and nodes list
netfromadj = function(A, names) {
  nedges = sum(A)
  ngenes = length(names)
  out = matrix(nrow=nedges, ncol=2)
  k=1
  for (i in 1:(ngenes) ) {
    for (j in (1:ngenes)[-i] ) {
      if (A[i,j]==1) {
        out[k,1]=names[i]
        out[k,2]=names[j]
        k=k+1
      }
    }
  }
  return(out)
}

# compute pvalue
pvalue = function(naf, nout, oa, i_f) {
  oac = nout - oa
  p1 = phyper(naf, m = oa, n = oac, k = i_f, lower.tail = FALSE) # p > naf
  p2 = phyper(naf-1, m = oa, n = oac, k = i_f, lower.tail = TRUE) # p < naf
  p3 = 1-p1-p2
  pval = 2*min(p1,p2)+p3
  overenr = (p1<p2)
  return(list('p'=pval,'overenr'=overenr))
}


########## pnea CLASS
pneac = function(x,...) UseMethod('pnea', x)

pnea = function(agslist, fgslist = NULL, network, nodes, alpha = NULL, agsnames = NULL, fgsnames = NULL) {
  if (class(network) == "matrix") {
    if ( ncol(network)>2 ) net = netfromadj(network, nodes)
    else if ( ncol(network)==2 ) net=network
  }
  else if (class(network) == 'igraph') {
    requireNamespace("igraph", quietly = TRUE)
    adjacency = as.matrix(igraph::get.adjacency(network))
    net = netfromadj(adjacency, nodes)
  }
  # NB: from this point on, 'net' is the network matrix to be used!!!
  if (is.null(alpha) == FALSE) {if ( alpha<=0 | alpha>=1 ) stop('alpha must be in (0,1)')}
  if (is.factor(nodes) == T) {nodes = as.character(nodes)}
  oa = numeric()
  i_f = numeric()
  agslogic = vector("list", length(agslist)) 
  netred = vector("list", length(agslist)) # rows where the first column entry is in A
  fgslogic = vector("list", length(fgslist))
  naf = numeric()
  p = naf
  expect = naf
  o = vector()
  concl = o
  from = o
  to = o
  if ( is.null( names(agslist) ) ) names(agslist) = agsnames
  for (i in 1:length(agslist)) {
    if (is.factor(agslist[[i]]) == T) {agslist[[i]] = as.character(agslist[[i]])}
    agslogic[[i]] = (net[,1] %in% agslist[[i]])
    oa[i] = sum(agslogic[[i]])
    netred[[i]] = net[agslogic[[i]],]
  }
  # first case: each ags vs each fgs (fgslist provided!)
  if (is.null(fgslist) == FALSE) {
    if ( is.null( names(fgslist) ) ) names(fgslist) = fgsnames
    for (i in 1:length(fgslist)) {
      if (is.factor(fgslist[[i]]) == T) {fgslist[[i]] = as.character(fgslist[[i]])}
      fgslogic[[i]] = (net[,2] %in% fgslist[[i]])
      i_f[i] = sum(fgslogic[[i]])
    }
    nout = dim(net)[1]
    k=1
    for (i in 1:length(agslist)) {
      for (j in 1:length(fgslist)) {
        from[k] = names(agslist)[i]
        to[k] = names(fgslist)[j]
        naf[k] = sum( netred[[i]][,2] %in%  fgslist[[j]])
        temp = pvalue(naf[k], nout = nout, oa = oa[i], i_f[j])
        p[k] = temp$p
        expect[k] = round(i_f[j] * oa[i] / nout, 4)
        o[k] = temp$overenr
        if (is.null(alpha) == FALSE) {
          if (p[k] > alpha) concl[k] = 'No enrichment'
          else if (p[k] <= alpha & o[k]==T) concl[k] = 'Overenrichment'
          else if (p[k] <= alpha & o[k]==F) concl[k] = 'Underenrichment'
        }
        k = k+1
      }
    }
  }
  # second case: only agslist provided
  if (is.null(fgslist) == TRUE) {
    for (i in 1:length(agslist)) {
      fgslogic[[i]] = (net[,2] %in% agslist[[i]])
      i_f[i] = sum(fgslogic[[i]])
    }
    nout = dim(net)[1]
    k=1
    for (i in 1:(length(agslist)-1)) {
      for (j in (i+1):length(agslist)) {
        from[k] = names(agslist)[i]
        to[k] = names(agslist)[j]
        naf[k] = sum( netred[[i]][,2] %in%  agslist[[j]])
        temp = pvalue(naf[k], nout = nout, oa = oa[i], i_f[j])
        p[k] = temp$p
        expect[k] = round(i_f[j] * oa[i] / nout, 4)
        o[k] = temp$overenr
        if (is.null(alpha) == FALSE) {
          if (p[k] > alpha) concl[k] = 'No enrichment'
          else if (p[k] <= alpha & o[k]==T) concl[k] = 'Overenrichment'
          else if (p[k] <= alpha & o[k]==F) concl[k] = 'Underenrichment'
        }
        k = k+1
      }
    }
  }
  # final common code
  if (is.null(alpha) == FALSE) {
    results = data.frame(from, to, naf, expect, p, concl)
    names(results) = c('AGS', 'FGS', 'naf', 'expected_naf', 'pvalue', 'conclusion')
  }
  else {
    results = data.frame(from, to, naf, expect, p)
    names(results) = c('AGS', 'FGS', 'naf', 'expected_naf', 'pvalue')
  }
  class(results)=c('pnea','data.frame')
  results
}

summary.pnea = function(object, ...) {
  cat("Number of comparisons:",dim(object)[1],"\n")
  cat("Enrichments at 1% level:",sum(object$pvalue<0.01),"\n")
  cat("Enrichments at 5% level:",sum(object$pvalue<0.05),"\n")
  cat("Kolmogorov-Smirnov test for uniformity of p-values:",round(ks.test(x=object$pvalue, y='punif')$p.value,4),"\n")
}

plot.pnea = function(x, nbreaks=10, ...) {
  par(mfrow=c(1,2))
  hist(x$pvalue,breaks=nbreaks, xlim=c(0,1), xlab='p-value', main='Histogram of p-values')
  plot(ecdf(x$p),main='P-p plot of p-values',xlab='x',ylab ='F(x)',xlim=c(0.037,0.963),ylim=c(0.037,0.963))
  abline(v=1);abline(a=0,b=1,col='red')
  legend(0.6,0.2,c('Distribution of','p-values','Uniform','distribution'),col=c('black','white','red','white'),lty=1,cex=0.6)
  }

print.pnea = function(x,nrows=10,...) {
  class(x) = 'data.frame'
  if (nrows == 'ALL') {
    nrows = dim(x)[1]
  }
  head(x,nrows)
}


