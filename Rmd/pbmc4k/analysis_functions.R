liblist <-
  c(
    "ggplot2",
    "MASS",
    "plyr",
    "reshape2",
    "ggpubr",
    "gridExtra",
    "ggpmisc",
    "cowplot",
    "corrplot",
    "ggrepel",
    "optparse",
    "rsvd",
    "scran",
    "igraph",
    "rtracklayer",
    "knitr",
    "mclust",
    "Rtsne",
    "factoextra",
    "edgeR",
    "gmodels",
    "RANN",
    "NMF"
  )
for (libname in liblist) {
  suppressWarnings(suppressMessages(library(libname, character.only = TRUE)))
}
addTheme <- function(p) {
  p <- p + theme(legend.position = "top")
  p <-
    p + theme(plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = 'bold'
    ))
  p <-
    p + theme(
      axis.title.x = element_text(
        color = 'black',
        size = 10,
        face = 'bold'
      ),
      axis.title.y = element_text(
        size = 10,
        color = 'black',
        face = 'bold'
      )
    )
  p <-
    p + theme(
      axis.text.y = element_text(
        color = 'black',
        size = 10,
        face = 'bold'
      ),
      axis.text.x = element_text(
        color = 'black',
        size = 10,
        face = 'bold'
      )
    )
  p <-
    p + theme(
      legend.text = element_text(size = 10, face = 'bold'),
      legend.title = element_text(size = 10, face = 'bold')
    )
  return(p)
}
SummaryQuantificationByType <- function(data, genes, match_name,nthreshold) {
  biotypes <- genes[match(rownames(data), genes[,match_name]), 'gene_type']
  headers <- colnames(data)
  summ <- data.frame()
  for (i in 1:ncol(data)) {
    x <- data[, i]
    y <-
      data.frame(aggregate(x == nthreshold, by = list(biotypes), FUN = sum))
    colnames(y) <- c('biotype', headers[i])
    if (i == 1) {
      summ <- y
    } else{
      summ <- merge(summ, y, by = 'biotype')
    }
  }
  return(summ)
}

ParseGene <- function(gtf_file) {
  gtf_gencode <-
    readGFF(
      gtf_file,
      version = 2L,
      tags = c("gene_name", "gene_id", "transcript_id", "gene_type")
    )
  genes <- subset(gtf_gencode, gtf_gencode$type == "gene")
  return(genes)
}

rpm<-function(x,log=FALSE){
  tot.size<-apply(x,2,sum)
  rpm<-apply(x,1,function(x){x*10e6/tot.size})
  if(log){
    rpm<-log(rpm+1,base=2)
  }
  return(list('rpm'=t(rpm),'tot.size'=tot.size))
}
matrixZscore<-function(x,dim,alpha){
  if(dim==2){
    x<-scale(x)
  }else{
    x<-t(scale(t(x)))
  }
  qval<-abs(qnorm(alpha/2))
  x<-ifelse(x>qval,1,0)
  return(x)
}
compareNormExpr<-function(x,names,genes,log=TRUE){
  n<-length(names)
  delta.genes<-list()
  for (i in 1:(n-1)){
    x1<-x[[names[i]]]
    for(j in (i+1):n){
      x2<-x[[names[j]]]
      lab<-paste(names[i],names[j],sep=':')
      mlist<-match(rownames(x1),rownames(x2))
      nlist<-match(colnames(x1),colnames(x2))
      x2<-x2[mlist,nlist]
      if(log==FALSE){
        delta<-x1-x2
      }else{
        x1<-log(x1+1,base=2)
        x2<-log(x2+1,base=2)
        delta<-x1-x2
      }
      z.delta.cell<-matrixZscore(delta,2,0.05)
      diff.norm<-SummaryQuantificationByType(z.delta.cell, genes, 'gene_name',0)
      delta.genes[[lab]]<-diff.norm
    }
  }
  return(delta.genes)
}
plotCompNormExpr <- function(expr) {
  names <- names(expr)
  for (n in names) {
    X <- expr[[n]]
    rownames(X)<-X$biotype
    pseudo.genes <- grep(pattern = "pseudogene", x = rownames(X), value = TRUE)
    pseudos<-apply(X[pseudo.genes,-1],2,sum)
    X<-subset(X,!(rownames(X) %in% pseudo.genes))
    X<-rbind(c('pseudogenes',pseudos),X)
    umi.biotype.frac.m <- melt(X,id.vars='biotype')
    umi.biotype.frac.m$value<-as.numeric(umi.biotype.frac.m$value)
    p <- ggboxplot(
      umi.biotype.frac.m,
      x = "biotype",
      y = "value",
      orientation = "horizontal"
    ) + theme_bw() +
      scale_y_log10() +
      ggtitle(stringr::str_wrap(
        paste(
          n,
          ": # genes with significant foldchanges between two normalization"
        ),
        40
      ))
    p <- addTheme(p)
    print(p)
  }
}
plotTSNE <- function(data, genes) {
  names <- names(data)
  plots <- list()
  Y<-list()
  
  for (n in names) {
    x <- data[[n]]
    tsne <-
      Rtsne(
        t(x),
        check_duplicates = FALSE,
        dims = 2,
        perplexity = 30,
        max_iter = 500,
        verbose = FALSE
      )
    mlist <- match(genes, rownames(x))
    x <- data.frame(tsne$Y, t(x[mlist[!is.na(mlist)], ]))
    mx <- melt(x, id.vars = c('X1', 'X2'))
    p <-
      ggscatter(
        data = mx,
        x = 'X1',
        y = 'X2',
        color = 'value',
        facet.by = 'variable'
      ) + gradient_color(c("blue", "white", "red")) + 
      ggtitle(paste(n, " normalization", sep =''))
    plots[[n]] <- p
    Y[[n]]<-x
  }
  out<-list('plots'=plots,'Y'=Y)
  return(out)
}
doPCA<-function(data,pcafunction){
  names<-names(data)
  pcs<-list()
  for(n in names){
    if(pcafunction=="rpca"){
      mat.pcs <- rpca(data[[n]], scale = T,center=T)
    }else if (pcafunction == "fast"){
      mat.pcs<-fast.prcomp(data[[n]],center=T, scale=T)
    }else if (pcafunction == "regular"){
      mat.pcs<-prcomp(data[[n]],center=T, scale=T)
    }
    pcs[[n]] <- mat.pcs
  }
  return(pcs)
}
doScreePlots<-function(pcas,npcs){
  names<-names(pcas)
  fracs<-c()
  for(n in names){
    p<-pcas[[n]]
    f<-p$sdev^2/sum(p$sdev^2)
    subf<-f[1:npcs]
    names(subf)<-paste('PC',1:npcs,sep='')
    if(length(fracs)==0){
      fracs<-subf
    }else{
      fracs<-cbind(fracs,subf)
    }
  }
  colnames(fracs)<-names
  m<-melt(fracs)
  p <-
    ggplot(
      data = m,
      aes(Var1,value),
      facet_grid = Var2,
      xlab="PCs",
      ylab="% of variance explained by PCs "
    ) +geom_bar()+
    ggtitle(paste(" Screeplot of PCA", sep =''))
  return(p)
}

doEdgesCluster<-function(X,nn=30,do.jaccard=TRUE,method="Louvain") {
  nearest<-nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim <- 1*(nearest$nn.dists >= 0 )
  edges <- melt(t(nearest$nn.idx))
  colnames(edges) <- c("B", "A", "C")
  edges <- edges[,c("A","B","C")]
  edges$B <- edges$C
  edges$C<-1
  
  #Remove repetitions
  edges <- unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  if (do.jaccard){
    NN = nearest$nn.idx
    jaccard_dist <- apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    edges$C <- jaccard_dist
    edges <- subset(edges, C != 0) ## skip no-overlapping
    edges$C <- edges$C/max(edges$C) ## convert to ratio
  }
  Adj <- matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) <- rownames(X)
  colnames(Adj) < rownames(X)
  Adj[cbind(edges$A,edges$B)] <- edges$C ## up
  Adj[cbind(edges$B,edges$A)] <- edges$C ## down
  g<-graph.adjacency(Adj, mode = "undirected", weighted=TRUE)
  if (method=="Louvain") graph.out <- cluster_louvain(g)
  if (method=="Infomap") graph.out <- cluster_infomap(g)
  return(graph.out)
}
doSNNCluster<-function(X){
  # Run SNN graph
  grps <- buildSNNGraph(t(X),rand.seed = 1000)
  # extract memebership
  clusters <- cluster_fast_greedy(grps)
  memb <-  membership(clusters)
  return(memb)
}
plotTSneClustering<-function(pcas, tsne, methodname) {
  names<-names(pcas)
  ps<-list()
  mem<-data.frame()
  for(n in names) {
    X <- pcas[[n]]$rotation[, c(1:50)]
    X1<-tsne$Y[[n]]$X1
    X2<-tsne$Y[[n]]$X2
    e <- doEdgesCluster(X,
                        nn = 30,
                        do.jaccard = TRUE,
                        method = methodname)
    if(length(mem)==0){
      mem<-e$membership
    }else{
    mem<-cbind(mem,e$membership)
    }
    d <-
      data.frame(
        'X1' = X1,
        'X2' = X2,
        'membership' = as.factor(e$membership)
      )
    md <- melt(d, id.vars = c("X1", "X2"))
    ps[[n]]<-ggscatter(
      data = md,
      x = 'X1',
      y = 'X2',
      color = 'value',
      font.label = 8
    ) + theme_bw() + ggtitle(paste(n,' ',methodname,' Clustering'))
    
  }
  colnames(mem)<-names
  return(list('plots'=ps,'mem'=mem))
}