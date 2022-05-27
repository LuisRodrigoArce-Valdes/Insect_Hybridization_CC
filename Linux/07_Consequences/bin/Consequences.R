# With this script we will do some plots for the consecuencues of hybridization
rm(list = ls())
#library(fmsb)
library(tidyr)
# Radarchart code ####
radarchart <- function(df, axistype=0, seg=4, pty=16, pcol=1:8, plty=1:6, plwd=1,
                       pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                       cglcol="navy", axislabcol="blue", title="", maxmin=TRUE,
                       na.itp=TRUE, centerzero=FALSE, vlabels=NULL, vlcex=NULL,
                       caxislabels=NULL, calcex=NULL,
                       paxislabels=NULL, palcex=NULL, ...) {
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }
  if (maxmin==FALSE) { # when the dataframe does not include max and min as the top 2 rows.
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1, ...) # define x-y coordinates without any plot
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+CGap)/(seg+CGap), yy*(i+CGap)/(seg+CGap), lty=cglty, lwd=cglwd, border=cglcol)
    if (axistype==1|axistype==3) CAXISLABELS <- paste(i/seg*100,"(%)")
    if (axistype==4|axistype==5) CAXISLABELS <- sprintf("%3.2f",i/seg)
    if (!is.null(caxislabels)&(i<length(caxislabels))) CAXISLABELS <- caxislabels[i+1]
    if (axistype==1|axistype==3|axistype==4|axistype==5) {
      if (is.null(calcex)) text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol) else
        text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol, cex=calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  else {
    arrows(xx/(seg+CGap), yy/(seg+CGap), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  PAXISLABELS <- df[1,1:n]
  if (!is.null(paxislabels)) PAXISLABELS <- paxislabels
  if (axistype==2|axistype==3|axistype==5) {
    if (is.null(palcex)) text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol) else
      text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol, cex=palcex)
  }
  VLABELS <- colnames(df)
  
  if (!is.null(vlabels)) VLABELS <- vlabels
  
  
  
  
  
  ##--------------------------------------------------
  ## Modified by Killbill-(Me)
  ##--------------------------------------------------
  # Main code:
  # if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
  #   text(xx*1.2, yy*1.2, VLABELS, cex=vlcex, adj=adjVec)
  
  
  # Modified code:
  # Create a variable that round 'xx' value to 0 and 1 for non zero and 0.5 for 0 values.
  adjVec <- ifelse(round(xx) < 0, 1, ifelse(round(xx) > 0, 0, 0.5))
  
  #apply 'adjVec' variable to "adj" parameters of text.
  
  for (i in seq_along(xx)){
    if (is.null(vlcex)) text(xx[i]*1.1, yy[i]*1.1, VLABELS[i], adj=adjVec[i]) else
      text(xx[i]*1.1, yy[i]*1.1, VLABELS[i], cex=vlcex, adj=adjVec[i])
  }
  
  ##-------------------------------------------------
  ## End
  ##-------------------------------------------------
  
  
  
  
  series <- length(df[[1]])
  SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg+CGap)+(df[i,]-df[2,])/(df[1,]-df[2,])*seg/(seg+CGap)
    if (sum(!is.na(df[i,]))<3) { cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n",i,df[i,])) # for too many NA's (1.2.2012)
    } else {
      for (j in 1:n) {
        if (is.na(df[i, j])) { # how to treat NA
          if (na.itp) { # treat NA using interpolation
            left <- ifelse(j>1, j-1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left>1, left-1, n)
            }
            right <- ifelse(j<n, j+1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right<n, right+1, 1)
            }
            xxleft <- xx[left]*CGap/(seg+CGap)+xx[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            yyleft <- yy[left]*CGap/(seg+CGap)+yy[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            xxright <- xx[right]*CGap/(seg+CGap)+xx[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            yyright <- yy[right]*CGap/(seg+CGap)+yy[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft; yytmp <- yyleft;
              xxleft <- xxright; yyleft <- yyright;
              xxright <- xxtmp; yyright <- yytmp;
            }
            xxs[j] <- xx[j]*(yyleft*xxright-yyright*xxleft)/(yy[j]*(xxright-xxleft)-xx[j]*(yyright-yyleft))
            yys[j] <- (yy[j]/xx[j])*xxs[j]
          } else { # treat NA as zero (origin)
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j]*CGap/(seg+CGap)+xx[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
          yys[j] <- yy[j]*CGap/(seg+CGap)+yy[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], col=pfcols[i-2])
      } else {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], 
                density=pdensities[i-2], angle=pangles[i-2], col=pfcols[i-2])
      }
      points(xx*scale, yy*scale, pch=ptys[i-2], col=pcols[i-2])
    }
  }
}
# Beginning script ####

# Reading data
consequences <- read.table("../../../Windows/Busquedas/Consequences/Unified_Consequences.txt", sep = "\t", header = T)
consequences <- consequences[,-5]

# In odonates, removing from general consequences those in climatic change
# Sorting
consequences <- consequences[order(consequences$Group, consequences$Order, consequences$Sp1, consequences$Sp2),]

# Removing duplicates
consequences <- consequences[!duplicated(consequences[,-1]),]

# Adding ID
consequences$ID <- 1:nrow(consequences)

# Removing other repeated species [SCRIPT SPECIFIC FOR THIS DATASET] {mannualy because lazyness}
consequences <- consequences[!consequences$ID %in% c(31,44,70,71,72,73,82,92,95,100),]

# Removing extra rows
consequences <- consequences[consequences$Include=="Y",]

# Tableing
abs <- as.data.frame(table(consequences$Unified.Consequence, consequences$Order, consequences$Group))
rel <- as.data.frame(prop.table(table(consequences$Unified.Consequence, consequences$Order, consequences$Group), margin = 2:3))

# Tableing for all insects
abs.all <- as.data.frame(table(consequences$Unified.Consequence, consequences$Group))
abs.all$Group <- "Insects"
abs.all <- abs.all[,c(1,4,2,3)]
colnames(abs.all) <- colnames(abs)
abs <- rbind(abs, abs.all)

rel.all <- as.data.frame(prop.table(table(consequences$Unified.Consequence, consequences$Group), margin = 2))
rel.all$Group <- "Insects"
rel.all <- rel.all[,c(1,4,2,3)]
colnames(rel.all) <- colnames(rel)
rel <- rbind(rel, rel.all)


# Merging
consequences <- cbind(abs, rel=rel$Freq)
colnames(consequences)[1:3] <- c("Consequence","Order","Group") 
rm(abs, rel, abs.all, rel.all)

# Subsetting
insects <- list()
for (i in c(unique(consequences$Order))) {
  insects[[i]] <- consequences[consequences$Order==i,]
  insects[[i]] <- insects[[i]][,-c(2,4)]
  insects[[i]] <- insects[[i]][order(insects[[i]]$Consequence),]
  insects[[i]] <- spread(insects[[i]], Consequence, rel)
  row.names(insects[[i]]) <- insects[[i]][,1]
  insects[[i]] <- insects[[i]][,-1]
  # Adding maximum and minimum
  insects[[i]] <- rbind(rep(1,ncol(insects[[i]])), rep(0,ncol(insects[[i]])), insects[[i]])
  # Ordering columns
  insects[[i]] <- insects[[i]][,c(1,8,7,6,5,4,3,2)]
}

# Color vector
colors_border=c(rgb(0.8,0.2,0.5,0.9), rgb(0.2,0.5,0.5,0.9))
colors_in=c(rgb(0.8,0.2,0.5,0.4), rgb(0.2,0.5,0.5,0.4))

# Plotting
dir.create("../figures/", showWarnings = F)

# Spider plot
png("../figures/01_SpiderPlot.png", width = 32, height = 26, units = "cm", res = 300)
par(mfrow=c(3,2), family="serif")
for(i in unique(consequences$Order)) {
  radarchart(insects[[i]], axistype = 1,
             # custom polygon
             pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
             #custom labels
             vlcex=2.0, title = i)
}
dev.off()
