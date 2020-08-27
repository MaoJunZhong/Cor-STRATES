

################################

# R code from Cooney & Thomas 'Heterogeneous relationships between rates of speciation and 
# body size evolution across vertebrate clades'
# Nature Ecology & Evolution

################################


### Branch-rate metrics


DR_stat <- function(x, return.mean = FALSE){ 
  rootnode <- length(x$tip.label) + 1
  sprates <- numeric(length(x$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]    
      qx <- qx + el* (1 / 2^(index-1))
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  if (return.mean){
    return(mean(sprates))   
  }else{
    names(sprates) <- x$tip.label
    return(sprates)
  }
}


TR_stat <- function(x, return.mean = FALSE){ 
  # Calculate metric based on weighted averages of branch rates (TRes)
  rootnode <- length(x$tip.label) + 1
  sprates <- numeric(length(x$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    rx <- c()
    wx <- c()
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]    
      rx <- c(rx, el) ### modified
      wx <- c(wx,  1/2^(index-1))
      index <- index + 1
    }
    sprates[i] <- weighted.mean(rx, wx) # take weighted mean of bls to avoid node height bias associated with summing values
  }
  if (return.mean){
    return(mean(sprates))   
  }else{
    names(sprates) <- x$tip.label
    return(sprates)
  }
}


TB_stat <- function(phy) {
  # length of terminal branches of rate scaled tree
  tb <- setNames(phy$edge.length[sapply(1:length(phy$tip.label), function(x,y) which(y==x), y=phy$edge[,2])], phy$tip.label)
  return(tb)
}


# --------------------------- #


##¢ Simulation test functions


get.null.rates <- function (j, phy, null.trait.rates, type) {
  frates <- c()
  ftree <- phy
  ftree$edge.length <- null.trait.rates[j,]
  if (type == "tree") {
    frates <- ftree$edge.length
  }
  if (type == "tips") {
    frates <- TB_stat(ftree)
  }
  if (type == "es") {
    frates <- TR_stat(ftree)
  }
  return(frates)
}


sim.test <- function (spp.rates, trait.rates, null.trait.rates, cor.method = "spearman") {
  obs.cor <- cor(log(trait.rates), log(spp.rates), method = cor.method)
  sims <- nrow(null.trait.rates)
  sim.cor <- c()
  for (k in 1:sims) {
    frates <- null.trait.rates[k,]
    sim.cor <- c(sim.cor, cor(log(frates), log(spp.rates), method = cor.method))
  }
  upper <- (length(sim.cor[sim.cor >= obs.cor])+1)/(sims+1)
  lower <- (length(sim.cor[sim.cor <= obs.cor])+1)/(sims+1)
  pval <- 2*min(c(upper,lower)) # Calculate the two-tailed p value (remove "2" for one-tailed)
  if (pval == 0) { pval <- 2*(1/sims) }
  cis <- quantile(sim.cor, probs = c(0.025, 0.975))
  ses <- (obs.cor - mean(sim.cor)) / sd(sim.cor)
  out <- c(obs.cor, mean(sim.cor), cis[1], cis[2], pval, ses)
  names(out) <- c("obs.cor", "sim.cor.mean", "sim.cor.lci", "sim.cor.uci", "pval", "ses")
  return(out)
}


# --------------------------- #

