# Rank Rank Hypergeometric Overlap (Plaisier et al., Nucleic Acids Research, 2010)
# Compares 2 lists of differentially expressed genes, each from independent gene expression experiments
# Facilitates generation of multiple RRHO maps, each having the same scale


# Main function to call
generate.rrho<-function(pval.data,logfc.data,list,outdir,BY=FALSE){
  
  # pval.data: A data frame, col1 has gene ID; subsequent cols = p-values in each expt
  # logfc.data: A data frame, col1 has gene ID; subsequent cols = logfc in each expt
  # list: A matrix, 1st col has name of expt 1, 2nd col has name of expt 2
  # outdir: name of output directory
  # BY: Benjamini-Yekutieli multiple hypothesis correction
  
  max.scale<-list()  
  
  for(i in 1:nrow(list)) # for every expt pair, calculates rrho.max(list1, list2)
  {
    
    list1<-cbind(rownames(pval.data),-1*log10(pval.data[,as.character(list[i,1])])*sign(logfc.data[,as.character(list[i,1])]))
    list2<-cbind(rownames(pval.data),-1*log10(pval.data[,as.character(list[i,2])])*sign(logfc.data[,as.character(list[i,2])])) 
    
    # list 1 = data frame from expt 1 with two cols, col 1 has Gene ID 
    # and col 2 has signed ranking value (e.g. signed -log10(p-value))
    # list 2 = data.frame from expt 2 with two cols, col 1 has Gene ID 
    # and col 2 has signed ranking value (e.g. signed -log10(p-value))
    
    max.scale<-append(max.scale,rrho.max(list1,list2)) 
  }
  
  # Creates RRHO map using max(max.scale) for each expt pair
  for(j in 1:nrow(list)){ 
    
    list1<-cbind(rownames(pval.data),-1*log10(pval.data[,as.character(list[j,1])])*sign(logfc.data[,as.character(list[j,1])]))
    list2<-cbind(rownames(pval.data),-1*log10(pval.data[,as.character(list[j,2])])*sign(logfc.data[,as.character(list[j,2])])) 
    
    rrho.mod(list1,list2,maximum=max(unlist(max.scale)),labels=c(as.character(list[j,1]),as.character(list[j,2])),outputdir=outdir,BY=BY)
  }
}


###########RRHO Dependency functions###############################

# Calculates step size (no. of genes to increase by in each iteration of algorithm)
defaultftepSize <-function(list1, list2){
  n1<- dim(list1)[1]
  n2<- dim(list2)[1]
  result <- ceiling(min(sqrt(c(n1,n2))))  

}  

# Computes counts and p-values of RRHO Matrix
numericListOverlap<- function(sample1, sample2, stepsize){
  # sample1: genes in ordered list1
  # sample2: genes in ordered list2
  
  n<- length(sample1)

  # returns no. of common genes and hyper geometric -log10(p-value) for subset of genes
  overlap<- function(a,b) { 
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b]))) # checks if genes of a subset in list1 are in subset of list2, counts no. of such instances
    log.pval<- -phyper(q=count-1, m=a, n=n-a, k=b, lower.tail=FALSE, log.p=TRUE) # performs hypergeometric test for subset of genes in list1 and list2
    return(c(counts=count, log.pval=log.pval))    
  } 
  
  indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize)) # creates matrix where each cell represents i and j genes
  overlaps<- apply(indexes, 1, function(x) overlap(x['i'],x['j'])) # calculates values of counts and -log10(p-value) for each cell of matrix
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  # matrix having no. of common genes in each subset
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows) # matrix having -log10(p-values) of each subset 
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals))  
}
                   

#################gets max value to scale subsequent rrho maps #########################
rrho.max<-function (list1, list2, stepsize = defaultftepSize(list1, list2))
{
  if (length(list1[, 1]) != length(unique(list1[, 1])))  
    stop("Non-unique gene identifier found in list1")
  if (length(list2[, 1]) != length(unique(list2[, 1]))) 
    stop("Non-unique gene identifier found in list2")
  
  
  result <- list(hypermat = NA, hypermat.counts = NA, n.items = nrow(list1), 
                 stepsize = stepsize, hypermat.by = NA, call = match.call())

  list1 <- list1[order(as.numeric(list1[, 2]), decreasing = TRUE), ] # rearrange rows acc to decreasing order of signed -log10(p-value)
  list2 <- list2[order(as.numeric(list2[, 2]), decreasing = TRUE), ] # rearrange rows acc to decreasing order of signed -log10(p-value)
  
  
  nlist1 <- length(list1[, 1]) 
  nlist2 <- length(list2[, 1]) 
  
  
  N <- max(nlist1, nlist2) 
  
  .hypermat <- numericListOverlap(list1[, 1], list2[, 1], stepsize)  
  
  result$hypermat <- hypermat <- .hypermat$log.pval 
  result$hypermat.counts <- .hypermat$counts
  
  return(max(hypermat, na.rm = TRUE)) 
  
}


#################Actual RRHO function (draws rrho map)#################################
rrho.mod<-function (list1, list2, stepsize = defaultftepSize(list1, list2), maximum,
                    labels, plots = TRUE, outputdir , BY = FALSE) 
{
  if (length(list1[, 1]) != length(unique(list1[, 1])))
    stop("Non-unique gene identifier found in list1")
  if (length(list2[, 1]) != length(unique(list2[, 1]))) 
    stop("Non-unique gene identifier found in list2")
  if (plots && (missing(outputdir) || missing(labels))) 
    stop("When plots=TRUE, outputdir and labels are required.")
  
  result <- list(hypermat = NA, hypermat.counts = NA, n.items = nrow(list1), 
                 stepsize = stepsize, hypermat.by = NA, call = match.call())
  
  list1 <- list1[order(as.numeric(list1[, 2]), decreasing = TRUE), ] # reorder list 1 based on signed -log10(p-values) (decreasing)
  list2 <- list2[order(as.numeric(list2[, 2]), decreasing = TRUE), ] # reorder list 2 based on signed -log10(p-values) (decreasing)
  nlist1 <- length(list1[, 1]) 
  nlist2 <- length(list2[, 1]) 
  N <- max(nlist1, nlist2) # max of lengths
  
  .hypermat <- numericListOverlap(list1[, 1], list2[, 1], stepsize) # stores -log10(p-values) and count of overlap
  result$hypermat <- hypermat <- .hypermat$log.pval 
  result$hypermat.counts <- .hypermat$counts 
  
  if (BY) {
    hypermatvec <- matrix(hypermat, nrow = nrow(hypermat) * 
                            ncol(hypermat), ncol = 1) 
    hypermat.byvec <- p.adjust(exp(-hypermatvec), method = "BY") 
    hypermat.by <- matrix(-log(hypermat.byvec), nrow = nrow(hypermat), 
                          ncol = ncol(hypermat))
    result$hypermat.by <- hypermat.by
  }
  
  if (plots) { 
    require(VennDiagram)
    require(grid)
    color.bar <- function(lut, min, max = -min, nticks = 11, 
                          ticks = seq(min, max, len = nticks), title = "") {
      scale <- (length(lut) - 1)/(max - min) 
      plot(c(0, 10), c(min, max), type = "n", bty = "n", 
           xaxt = "n", xlab = "", yaxt = "n", ylab = "")
      mtext(title, 2, 2.3, cex = 0.8)
      axis(2, round(ticks, 0), las = 1, cex.lab = 0.8)
      for (i in 1:(length(lut) - 1)) {
        y <- (i - 1)/scale + min
        rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
      }
    }
    .filename <- paste("RRHOMap", labels[1], "_VS_", labels[2], 
                       ".jpg", sep = "")
    jpeg(paste(outputdir, .filename, sep = "/"), width = 8, 
         height = 8, units = "in", quality = 100, res = 150)
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                     "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
    #image(hypermat, xlab = "", ylab = "", col = jet.colors(100), 
    #axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map",zlim=c(min(hypermat, na.rm = TRUE),maximum));
    
    #test minimum
    image(hypermat, xlab = "", ylab = "", col = jet.colors(100), 
          axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map",zlim=c(0,maximum));
 
    
    mtext(labels[2], 2, 0.5)
    mtext(labels[1], 1, 0.5)
    #     color.bar(jet.colors(100), min = min(hypermat, na.rm = TRUE), 
    #               max = maximum, nticks = 6, title = "-log(P-value)")
    
    
    #test minimum
    color.bar(jet.colors(100), min =0, 
              max = maximum, nticks = 6, title = "-log(P-value)")
    
    dev.off()
    
    # scatter plot
    list2ind <- match(list1[, 1], list2[, 1])
    list1ind <- 1:nlist1
    corval <- cor(list1ind, list2ind, method = "spearman")
    .filename <- paste("RankScatter", labels[1], "_VS_", 
                       labels[2], ".jpg", sep = "")
    jpeg(paste(outputdir, .filename, sep = "/"), width = 8, 
         height = 8, units = "in", quality = 100, res = 150)
    plot(list1ind, list2ind, xlab = paste(labels[1], "(Rank)"), 
         ylab = paste(labels[2], "(Rank)"), pch = 20, main = paste("Rank-Rank Scatter (rho = ", 
                                                                   signif(corval, digits = 3), ")", sep = ""), cex = 0.5)
    model <- lm(list2ind ~ list1ind)
    lines(predict(model), col = "red", lwd = 3)
    dev.off()
    maxind.ur <- which(max(hypermat[ceiling(nrow(hypermat)/2):nrow(hypermat), 
                                    ceiling(ncol(hypermat)/2):ncol(hypermat)], na.rm = TRUE) == 
                         hypermat, arr.ind = TRUE)
    indlist1.ur <- seq(1, nlist1, stepsize)[maxind.ur[1]]
    indlist2.ur <- seq(1, nlist2, stepsize)[maxind.ur[2]]
    genelist.ur <- intersect(list1[indlist1.ur:nlist1, 1], 
                             list2[indlist2.ur:nlist2, 1])
    maxind.lr <- which(max(hypermat[1:(ceiling(nrow(hypermat)/2) - 
                                         1), 1:(ceiling(ncol(hypermat)/2) - 1)], na.rm = TRUE) == 
                         hypermat, arr.ind = TRUE)
    indlist1.lr <- seq(1, nlist1, stepsize)[maxind.lr[1]]
    indlist2.lr <- seq(1, nlist2, stepsize)[maxind.lr[2]]
    genelist.lr <- intersect(list1[1:indlist1.lr, 1], list2[1:indlist2.lr, 
                                                            1])
    .filename <- paste(outputdir, "/RRHO_GO_MostDownregulated", 
                       labels[1], "_VS_", labels[2], ".csv", sep = "")
    write.table(genelist.ur, .filename, row.names = F, quote = F, 
                col.names = F)
    .filename <- paste(outputdir, "/RRHO_GO_MostUpregulated", 
                       labels[1], "_VS_", labels[2], ".csv", sep = "")
    write.table(genelist.lr, .filename, row.names = F, quote = F, 
                col.names = F)
    .filename <- paste(outputdir, "/RRHO_VennMost", labels[1], 
                       "__VS__", labels[2], ".jpg", sep = "")
    jpeg(.filename, width = 8.5, height = 5, units = "in", 
         quality = 100, res = 150)
    vp1 <- viewport(x = 0.25, y = 0.5, width = 0.5, height = 0.9)
    vp2 <- viewport(x = 0.75, y = 0.5, width = 0.5, height = 0.9)
    pushViewport(vp1)
    h1 <- draw.pairwise.venn(length(indlist1.ur:nlist1), 
                             length(indlist2.ur:nlist2), length(genelist.ur), 
                             category = c(labels[1], labels[2]), scaled = TRUE, 
                             lwd = c(0, 0), fill = c("cornflowerblue", "darkorchid1"), 
                             cex = 1, cat.cex = 1.2, cat.pos = c(0, 0), ext.text = FALSE, 
                             ind = FALSE, cat.dist = 0.01)
    grid.draw(h1)
    grid.text("Down Regulated", y = 1)
    upViewport()
    pushViewport(vp2)
    h2 <- draw.pairwise.venn(length(1:indlist1.lr), length(1:indlist2.lr), 
                             length(genelist.lr), category = c(labels[1], labels[2]), 
                             scaled = TRUE, lwd = c(0, 0), fill = c("cornflowerblue", 
                                                                    "darkorchid1"), cex = 1, cat.cex = 1.2, cat.pos = c(0, 
                                                                                                                        0), ext.text = FALSE, main = "Negative", ind = FALSE, 
                             cat.dist = 0.01)
    grid.draw(h2)
    grid.text("Up Regulated", y = 1)
    dev.off()
  }
  return(result)
  
}
