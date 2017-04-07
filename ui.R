# Shiny Ciona UI
library(shiny)
shinyUI(fluidPage(
  titlePanel(h1("Cardiopharyngeal Specification 
             scRNA-seq Data",align = "center"),windowTitle = "Ciona scRNA-seq"),
  sidebarLayout(
    sidebarPanel("",h5("Email",span("xn262@nyu.edu", style = "color:blue"),"if you
        have any questions regarding of this website.")),
    mainPanel(
      h4("This website hosts the Spatial-Temporal Cardiopharyngeal Specification Single-Cell RNA-seq data 
        in", HTML("<I> Ciona Intestinalis.</I>"), "The data was generated in",a("Satija Lab",href="http://satijalab.org",target="_blank"), 
         "and", a("Christiaen Lab",href="http://biology.as.nyu.edu/object/LionelChristiaen.html",target="_blank"),"at", a("NYGC",href="http://www.nygenome.org",target="_blank"), 
         "/",a("NYU",href="http://biology.as.nyu.edu",target="_blank"), "and was analysed using",a("Seurat.",href="http://www.satijalab.org/seurat",target="_blank")),
      h4("Warning: For better performance, it is recommended to run this Shiny app localy.",a("",href="http://biology.as.nyu.edu",target="_blank"),"If you encountered a crash of this page, wait 3-5 mins and retry.")
      )
  ),
  
  sidebarLayout(
    sidebarPanel(h3("Spatial Dynamics (using hpf20 data)"),
                 textInput("genes", value = NULL, label = h3("Enter Gene Names")),
                 selectInput("hpfgeneList", label = h3("Select Gene List"),
                             choices = list("All","ASM", "Heart", "FHP_Specific", "SHP_Specific"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(h3("tsne Plot"), plotOutput("tsne",width = "80%")),
        tabPanel(h3("Feature Plot"), plotOutput("featurePlot",width = "80%")),
        tabPanel(h3("Violin Plot"), plotOutput("violinPlot",width="80%")),
        tabPanel(h3("Heatmap"), plotOutput("hpfheatmap",width="80%")
      )))),
    
  sidebarLayout(
    sidebarPanel(h3("Temporal Dynamics (with regulatory states)"),
                 selectInput("trajectory", label = h3("Select Trajectory"),
                                      choices = list("ASM", "FHP", "SHP")),
                 selectInput("geneList", label = h3("Select Gene List"),
                                      choices = list("ASM", "Heart", "FHP_Specific", "SHP_Specific")),
                 textInput("geneHeatmap", value = NULL, label = h3("Enter Gene Names")),
                 checkboxInput("smooth","Smoothed Data"),
                 selectInput("dataset", label = h3("Download Gene List"), choices=c("ASM","Heart","FHP","SHP","KH2013_UniqueNAME")),
                 downloadButton("downloadData", "Download")
                 ),

  mainPanel(
    tabsetPanel(
      tabPanel(h3("Heatmap"),plotOutput("heatmap",width = "80%")),
      tabPanel(h3("Gene Plot"), plotOutput("genePlot",width = "80%")),
      tabPanel(h3("Violin Plot"), plotOutput("vlnPlot",width="80%"))
    )))
))