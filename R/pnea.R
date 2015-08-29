####### INTERNAL FUNCTIONS
pneatest = function(ags, fgs, net, nodes) {
  # ags, fgs, nodes: vectors, net: matrix with 2 columns
  Apres = ags[ags %in% nodes] #AGS
  Fpres = fgs[fgs %in% nodes] #FGS
  AC = nodes[ (nodes %in% ags)==F ] #complement of A
  a = length(Apres)
  f = length(Fpres)
  if (a==0) stop("None of the genes of A is present in the network")
  if (f==0) stop("None of the genes of F is present in the network")
  i_f = sum(net[,2] %in% Fpres)
  naf = sum( (net[,1] %in% Apres) ==T & (net[,2] %in% Fpres) ==T )
  oa = sum( (net[,1] %in% Apres) ==T )
  oac = sum( (net[,1] %in% AC) ==T )
  p1 = phyper(naf, m = oa, n = oac, k = i_f, lower.tail = FALSE) # p > naf
  p2 = phyper(naf-1, m = oa, n = oac, k = i_f, lower.tail = TRUE) # p < naf
  p3 = 1-p1-p2
  pval = 2*min(p1,p2)+p3
  overenr = (p1<p2)
  result = list('naf'=naf,'pval'=pval,'overenr'=overenr) # create result object
  result # give result object as last command
}

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

########## pnea CLASS
pneac = function(x,...) UseMethod('pneac', x)

pnea = function(agslist, fgslist, network, nodes, alpha=NULL, agsnames, fgsnames) {
  if (class(network) == "matrix") {
    if (ncol(network)>2) net = netfromadj(network, nodes)
    else if (ncol(network)==2) net=network
    }
  else if (class(network) == 'igraph') {
    requireNamespace("igraph", quietly = TRUE)
    adjacency = as.matrix(igraph::get.adjacency(network))
    net = netfromadj(adjacency, nodes)
    }
  if ( is.null( names(agslist) ) ) names(agslist) = agsnames
  if ( is.null( names(fgslist) ) ) names(fgslist) = agsnames
  from = character()
  to = character()
  nafvec = numeric()
  pvalvec = numeric()
  concl = character()
  k=1
  for (i in 1:length(agslist)) {
    for (j in 1:length(fgslist)) {
      from[k] = names(agslist)[i]
      to[k] = names(fgslist)[j]
      pnea = pneatest(agslist[[i]], fgslist[[j]], net, nodes)
      nafvec[k] = pnea$naf
      pvalvec[k] = round(pnea$pval,5)
      if (is.null(alpha) == FALSE) {
        if ( alpha<=0 | alpha>=1 ) stop('alpha must be in (0,1)')
        if (pnea$pval > alpha) concl[k] = 'No enrichment'
        else if (pnea$pval<= alpha & pnea$overenr==T) concl[k] = 'Overenrichment'
        else if (pnea$pval<= alpha & pnea$overenr==F) concl[k] = 'Underenrichment'
      }
      k=k+1
    }
  }
  if (is.null(alpha) == FALSE) {
    results = data.frame(from, to, nafvec, pvalvec, concl)
    names(results) = c('AGS', 'FGS', 'naf', 'pvalue', 'conclusion')
  }
  else {
    results = data.frame(from, to, nafvec, pvalvec)
    names(results) = c('AGS', 'FGS', 'naf', 'pvalue')
  }
  class(results)=c('pneac','data.frame')
  results
}

summary.pneac = function(object, ...) {
  cat("Number of AGSs tested:",length(unique(object$AGS)),"\n")
  cat("Number of FGSs tested:",length(unique(object$FGS)),"\n")
  cat("Number of comparisons:",dim(object)[1],"\n")
  cat("Enrichments at 1% level:",sum(object$pvalue<0.01),"\n")
  cat("Enrichments at 5% level:",sum(object$pvalue<0.05),"\n")
  cat("Kolmogorov-Smirnov test for uniformity of p-values:",round(ks.test(x=object$pvalue, y='punif')$p.value,4),"\n")
}

plot.pneac = function(x, nbreaks=10, ...) {
  par(mfrow=c(1,2))
  hist(x$pvalue,breaks=nbreaks,xlab='p-value', main='Histogram of p-values')
  plot(ecdf(x$p),main='P-p plot of p-values',xlab='x',ylab ='F(x)',xlim=c(0.037,0.963),ylim=c(0.037,0.963))
  abline(v=1);abline(a=0,b=1,col='red')
  legend(0.6,0.2,c('Distribution of','p-values','Uniform','distribution'),col=c('black','white','red','white'),lty=1,cex=0.6)
  }

print.pneac = function(x,nrows=10,...) {
  class(x) = 'data.frame'
  head(x,nrows)
}

