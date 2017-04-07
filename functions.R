genes.plot.pseudo <- function (object,genes.use, pseudo="pseudo.time",pch.use = 16, cex.use = 1.5,
                               spline.span = 0.75,name.x=NULL,col.use=NULL,do.logit=FALSE,
                               inset=c(0,0),ylim=NULL,do.label=T,lwd=1,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }

  cell.ids = set.ifnull(cell.ids, object@cell.names)
  xlab=paste(name.x,"Pseudotime")
  data.use = data.frame(t(FetchData(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.scaled = FALSE)))
  pseudo.order=order(data.use[pseudo,cell.ids])
  pseudo.data=unlist(data.use[pseudo,cell.ids][pseudo.order])
  object@scale.data=object@scale.data[,pseudo.order]
  ylim=set.ifnull(ylim,range(data.use[genes.use,cell.ids]))
  
  plot(0,0,xlim = range(pseudo.data),ylim = ylim,type = "l",xlab=xlab,ylab="log(FPKM)",
       cex = cex.use, pch = pch.use, font.lab=2,cex.lab=1.2,bty="l",...)
  col.use=set.ifnull(col.use,rainbow(length(genes.use)))
  
  for (i in 1:length(genes.use)){
    g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
    loess.fit = loess(g2 ~ pseudo.data, span = spline.span)
    lines(pseudo.data, loess.fit$fitted,type="l",col=col.use[i],lwd=lwd)
  }
  if(do.logit){
    gene.norm=apply(object@data[genes.use[i],cell.ids],1, FUN = function(X) (X - min(X))/diff(range(X)))
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm[pseudo.order]~pseudo.data,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn)]
    turn.mat=c(turn.mat,turn.point)
    abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
  }
  if(do.label){
    legend("topright",legend = genes.use, 
           lty=1, col=col.use, cex=.75,xpd=T,inset = inset)
  }
}

genePlot.pseudo <- function (object,pseudo="pseudo.time", gene, cell.ids = NULL, col.use = NULL, 
                             pch.use = 16, cex.use = 1.5, do.ident = FALSE, conf.int=FALSE,
                             do.spline = FALSE,spline.span = 0.75,name.x=NULL, do.logit=FALSE,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = use.scale)))
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  name.x = set.ifnull(name.x, pseudo)
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = as.numeric(data.use[gene, cell.ids])
  ident.use = as.factor(object@ident[cell.ids])
  pseudo.data=as.matrix(FetchData(object,pseudo))
  
  if(do.logit){
    gene.norm=apply(object@data[gene,colnames(object@data)],1, FUN = function(X) (X - min(X))/diff(range(X)))
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm~pseudo.data,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn),1]
  }
  
  if (length(col.use) > 1) {
    col.use = col.use[as.numeric(ident.use)]
  }
  else {
    col.use = set.ifnull(col.use, as.numeric(ident.use))
  }

  plot(g1, g2, xlab = name.x, ylab = gene, col = col.use, cex = cex.use, 
       main = "", pch = pch.use, font.lab=2,cex.lab=1.2, ...)
  
  if (do.logit){
    abline(v=turn.point,lwd=4,col="purple",lty="longdash")
  }
  
  if (do.spline) {
    loess.fit = loess(g2 ~ g1, span = spline.span)
    lines(g1[order(pseudo.data)], loess.fit$fitted[order(pseudo.data)],col = "black",lwd=3)
    if(conf.int){
      prid = predict(loess.fit,se=T)
      lines(g1[order(g1)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(g1)],col = "darkblue",lwd=1.2,lty=2)
      lines(g1[order(g1)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(g1)],col = "darkblue",lwd=1.2,lty=2)
    }
  }
}
