# Shiny Ciona Server
library('shiny')
library('Seurat')
library('RColorBrewer')
library('gplots')
source("./functions.R")
load("./temporal.Robj")
load("./spatial.Robj")
heartGene=rownames(subset(read.table("./panHP_20.markers.txt"), power>0.5 & avg_diff>1))
asmGene=rownames(subset(read.table("./ASM_20.markers.txt"), power>0.5 & avg_diff>1))
shpSpecific=rownames(subset(read.table("./SHPspecific_20.markers.txt"),power>0.5&avg_diff>1))
fhpSpecific=rownames(subset(read.table("./FHPspecific_20.markers.txt"),power>0.5&avg_diff>1))

shinyServer(function(input, output) {
  
  genes <- reactive({
    unlist(strsplit(input$genes,"\\, |\\,| "))
  })
  
  geneNames <- reactive({
    unlist(strsplit(input$geneHeatmap,"\\, |\\,| "))
  })
  
  trajectObj <- reactive({
    as.character(input$trajectory)
  })
  
  hpfgeneList <- reactive({
    switch(input$hpfgeneList,
           All = c(heartGene[1:15], fhpSpecific[1:15], shpSpecific,asmGene[1:15]),
           ASM = asmGene[1:15],
           Heart = heartGene[1:15],
           FHP_Specific = fhpSpecific[1:15],
           SHP_Specific = shpSpecific)
  })
  
  geneList <- reactive({
    switch(input$geneList,
           ASM = asmGene[1:15],
           Heart = heartGene[1:15],
           FHP_Specific = fhpSpecific[1:15],
           SHP_Specific = shpSpecific)
  })
    
  trajectory <- reactive({
    switch(input$trajectory, 
           ASM = "asm.test", 
           FHP = "fhp.test", 
           SHP = "shp.test") 
    })
  
  datasetInput <-reactive({
    switch(input$dataset,
           ASM = "ASM_20.markers.txt",
           Heart = "panHP_20.markers.txt",
           FHP = "FHPspecific_20.markers.txt",
           SHP = "SHPspecific_20.markers.txt",
           KH2013_UniqueNAME = "KH2013_UniqueNAME.txt"
           )
  })
  
  output$tsne <- renderPlot({
    TSNEPlot(hpf20,do.label = T,colors.use = c("red","orange","lightblue","blue"),pt.size = 1.5)
  })
  
  output$featurePlot <- renderPlot({
    if(length(genes()) != 0){
      if(length(genes()) == 1){
        FeaturePlot(hpf20,features.plot = genes(), nCol = 1,pt.size = 1.5)
      } else{
        FeaturePlot(hpf20,features.plot = genes(),pt.size = 1.5)
      }
    }
  })
  
  output$violinPlot <- renderPlot({
    if(length(genes()) != 0){
      VlnPlot(hpf20, genes(),cols.use = c("red","orange","lightblue","blue"),size.x.use = 0.8)
      }
    })
  
  output$hpfheatmap <- renderPlot({
    if(length(genes())==0){
      if(input$hpfgeneList=="All"){
        DoHeatmap(hpf20,genes.use = hpfgeneList(),order.by.ident = T,
                  slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                  key.xlab = "", RowSideColors = c(rep("pink",15),rep("red",15),rep("orange",4),rep("blue",15)), key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                  sepwidth = c(1,0.2),srtCol=0,cex.col = 1,cexRow=0.7,keysize=1,density.info=c("none"))
        
      } else {
        DoHeatmap(hpf20,genes.use = hpfgeneList(),order.by.ident = T,
                  slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                  key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                  sepwidth = c(1,0.2),srtCol=0,cex.col = 1,cexRow=0.7,keysize=1,density.info=c("none"))
      }
    } else {
      if(length(genes()) == 1){
        heatmap.2(rbind(hpf20@scale.data[genes(),order(hpf20@ident)],
                        hpf20@scale.data[genes(),order(hpf20@ident)]),
                  Rowv = NA, Colv = NA, trace = "none",
                  key = T,labRow = genes(),col=col,
                  key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                  colsep = cumsum(table(hpf20@ident)),sepcolor = "black",sepwidth = c(2,1),
                  labCol=c(rep("",50),"20FHP",rep("",80),"20SHP",rep("",100),"20ASM1",
                           rep("",50),"20ASM2"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
      } else {
        DoHeatmap(hpf20,genes.use = genes(),order.by.ident = T,
                  slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                  key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                  sepwidth = c(1,1),srtCol=0,cex.col = 1,cexRow=1,keysize=1,density.info=c("none"))
      }
    }
    })

  output$heatmap <- renderPlot({
    if(length(geneNames()) == 0){
      if(input$trajectory == "SHP"){
        if(input$smooth){
          heatmap.2(as.matrix(shp.test@raw.data[geneList(),]),
                    Rowv = NA, Colv = NA, trace = "none",
                    key = T,labRow = geneList(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                    key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                    colsep = cumsum(table(shp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                    labCol=c(rep("",20),"TVC",rep("",50),"STVC",rep("",55),"SHP1",
                             rep("",70),"SHP2"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
        } else {
          DoHeatmap(shp.test,genes.use = geneList(),order.by.ident = F,
                    slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                    key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                    sepwidth = c(1,1),srtCol=0,cex.col =1,cexRow=1,keysize=1,density.info=c("none"))
        }
      }
      if(input$trajectory == "FHP"){
        if(input$smooth){
          heatmap.2(as.matrix(fhp.test@raw.data[geneList(),]),
                    Rowv = NA, Colv = NA, trace = "none",
                    key = T,labRow = geneList(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                    key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                    colsep = cumsum(table(fhp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                    labCol=c(rep("",35),"TVC",rep("",90),"FHP1",rep("",120),"FHP2",
                             rep("",100),"FHP3"),srtCol=0,cexCol=1,cexRow=1,density.info=c("none"))
        } else {
          DoHeatmap(fhp.test,genes.use = geneList(),order.by.ident = F,
                    slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                    key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                    sepwidth = c(1,1),srtCol=0,cex.col=1,cexRow=1,keysize=1,density.info=c("none"))
        }
      }
      if(input$trajectory == "ASM"){
        if(input$smooth){
          heatmap.2(as.matrix(asm.test@raw.data[geneList(),]),
                    Rowv = NA, Colv = NA, trace = "none",
                    key = T,labRow = geneList(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                    key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                    colsep = cumsum(table(asm.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                    labCol=c(rep("",30),"TVC",rep("",80),"STVC",rep("",90),"ASM1",
                             rep("",80),"ASM2",rep("",55),"ASM3"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
        } else {
          DoHeatmap(asm.test,genes.use = geneList(),order.by.ident = F,
                    slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                    key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                    sepwidth = c(1,1),srtCol=0,cex.col=1,cexRow=1,keysize=1,density.info=c("none"))
        }
      } 
    } else {
      if(length(geneNames()) == 1){
        if(input$trajectory == "SHP"){
          if(input$smooth){
            heatmap.2(as.matrix(rbind(shp.test@raw.data[geneNames(),],
                                      shp.test@raw.data[geneNames(),])),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(shp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",20),"TVC",rep("",50),"STVC",rep("",55),"SHP1",
                               rep("",70),"SHP2"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          } else {
            heatmap.2(rbind(shp.test@scale.data[geneNames(),],
                            shp.test@scale.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col,
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(shp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",20),"TVC",rep("",50),"STVC",rep("",55),"SHP1",
                               rep("",70),"SHP2"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
        if(input$trajectory == "FHP"){
          if(input$smooth){
            heatmap.2(as.matrix(rbind(fhp.test@raw.data[geneNames(),],
                                      fhp.test@raw.data[geneNames(),])),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(fhp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",35),"TVC",rep("",90),"FHP1",rep("",120),"FHP2",
                               rep("",100),"FHP3"),srtCol=0,cexCol=1,cexRow=1,density.info=c("none"))
          } else {
            heatmap.2(rbind(fhp.test@scale.data[geneNames(),],
                            fhp.test@scale.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col,colsep = cumsum(table(fhp.test@ident)),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      sepcolor="black",sepwidth = c(2,1),labCol=c(rep("",35),"TVC",
                                                                  rep("",90),"FHP1",rep("",120),"FHP2",rep("",100),"FHP3"),
                      srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
        if(input$trajectory == "ASM"){
          if(input$smooth){
            heatmap.2(as.matrix(rbind(asm.test@raw.data[geneNames(),],
                                      asm.test@raw.data[geneNames(),])),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(asm.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",30),"TVC",rep("",80),"STVC",rep("",90),"ASM1",
                               rep("",80),"ASM2",rep("",55),"ASM3"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          } else {
            heatmap.2(rbind(asm.test@scale.data[geneNames(),],
                            asm.test@scale.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col,sepcolor = "black",colsep = cumsum(table(asm.test@ident)),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),sepwidth = c(2,1),
                      labCol=c(rep("",30),"TVC",rep("",80),"STVC",rep("",90),"ASM1",
                               rep("",80),"ASM2",rep("",55),"ASM3"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
      } else {
        if(input$trajectory == "SHP"){
          if(input$smooth){
            heatmap.2(as.matrix(shp.test@raw.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(shp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",20),"TVC",rep("",50),"STVC",rep("",55),"SHP1",
                               rep("",70),"SHP2"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          } else {
            DoHeatmap(shp.test,genes.use = geneNames(),order.by.ident = F,
                      slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                      key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                      sepwidth = c(1,1),labCol=c(rep("",20),"TVC",rep("",50),"STVC",rep("",55),"SHP1",
                                                 rep("",70),"SHP2"),srtCol=0,cex.col=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
        if(input$trajectory == "FHP"){
          if(input$smooth){
            heatmap.2(as.matrix(fhp.test@raw.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(fhp.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",35),"TVC",rep("",90),"FHP1",rep("",120),"FHP2",
                               rep("",100),"FHP3"),srtCol=0,cexCol=1,cexRow=1,density.info=c("none"))
          } else {
            DoHeatmap(fhp.test,genes.use = geneNames(),order.by.ident = F,
                      slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                      key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                      sepwidth = c(1,1),srtCol=0,cex.col=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
        if(input$trajectory == "ASM"){
          if(input$smooth){
            heatmap.2(as.matrix(asm.test@raw.data[geneNames(),]),
                      Rowv = NA, Colv = NA, trace = "none",
                      key = T,labRow = geneNames(),col=col, breaks = unique(c(seq(0,2,length=200),seq(2,6,length=200))),
                      key.title = "Expression Level Scale", key.xlab = "", mar=c(8,8),
                      colsep = cumsum(table(asm.test@ident)),sepcolor = "black",sepwidth = c(2,1),
                      labCol=c(rep("",30),"TVC",rep("",80),"STVC",rep("",90),"ASM1",
                               rep("",80),"ASM2",rep("",55),"ASM3"),srtCol=0,cexCol=1,cexRow=1,keysize=1,density.info=c("none"))
          } else {
            DoHeatmap(asm.test,genes.use = geneNames(),order.by.ident = F,
                      slim.col.label = T,draw.line = T,key.title = "Expression Level Scale", 
                      key.xlab = "", key.ylab = "", mar=c(8,8),col=col,sepcolor="black",
                      sepwidth = c(1,1),srtCol=0,cex.col=1,cexRow=1,keysize=1,density.info=c("none"))
          }
        }
      }
    }
  })
  
  output$genePlot <- renderPlot({
    if(length(geneNames()) == 0){
      if(trajectObj() == "ASM"){
        genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneList(), lwd=1.5, name.x = "ASM",inset=c(0,-0.1))
        abline(v=ps.asm,lwd=2,lty=2)
      }
      if(trajectObj() == "FHP"){
        genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneList(), lwd=1.5, name.x = "FHP",inset=c(0,-0.1))
        abline(v=ps.fhp,lwd=2,lty=2)
      }
      if(trajectObj() == "SHP"){
        genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneList(), lwd=1.5, name.x = "SHP",inset=c(0,-0.1))
        abline(v=ps.shp,lwd=2,lty=2)
      }
    } else {
      if(length(geneNames()) == 1){
        if(trajectObj() == "ASM"){
          genePlot.pseudo(eval(as.name(trajectory())),gene=geneNames(),col.use = c("green","yellow","lightblue","blue1","blue4"),do.spline = T,name.x = "ASM Trajectory",
                          cex.use = 0.8,cex.lab=1,do.logit = T)
          abline(v=ps.asm,lwd=2,lty=2)
        }
        if(trajectObj() == "FHP"){
          genePlot.pseudo(eval(as.name(trajectory())),gene=geneNames(),col.use = c("green","pink1","red1","red4"),do.spline = T,name.x = "FHP Trajectory",
                          cex.use = 0.8,cex.lab=1,do.logit = T)
          abline(v=ps.fhp,lwd=2,lty=2)
        }
        if(trajectObj() == "SHP"){
          genePlot.pseudo(eval(as.name(trajectory())),gene=geneNames(),col.use = c("green","yellow","orange1","orange4"),do.spline = T,name.x = "SHP Trajectory",
                          cex.use = 0.8,cex.lab=1,font.lab=2,do.logit = T)
          abline(v=ps.shp,lwd=2,lty=2)
        }
      } else {
        if(trajectObj() == "ASM"){
          genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneNames(), lwd=1.5, name.x = "ASM",inset=c(0,-0.1))
          abline(v=ps.asm,lwd=2,lty=2)
        }
        if(trajectObj() == "FHP"){
          genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneNames(), lwd=1.5, name.x = "FHP",inset=c(0,-0.1))
          abline(v=ps.fhp,lwd=2,lty=2)
        }
        if(trajectObj() == "SHP"){
          genes.plot.pseudo(eval(as.name(trajectory())),genes.use = geneNames(), lwd=1.5, name.x = "SHP",inset=c(0,-0.1))
          abline(v=ps.shp,lwd=2,lty=2)
        }
      }
    }
  })
  
  output$vlnPlot <- renderPlot({
    if(length(geneNames()) != 0){
      if(trajectObj() == "ASM"){
        VlnPlot(eval(as.name(trajectory())), geneNames(),cols.use = c("green","yellow","lightblue","blue1","blue4"),size.x.use = 8)
      }
      if(trajectObj() == "FHP"){
        VlnPlot(eval(as.name(trajectory())), geneNames(),cols.use = c("green","pink","red1","red4"),size.x.use = 8)
      }
      if(trajectObj() == "SHP"){
        VlnPlot(eval(as.name(trajectory())), geneNames(),cols.use = c("green","yellow","orange1","orange4"),size.x.use = 8)
      }
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function(){paste(input$dataset,"txt",sep = ".")},
    content = function(file){file.copy(datasetInput(),file)},
    contentType = "text/plain"
  )
})
