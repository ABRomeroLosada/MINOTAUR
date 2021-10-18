#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(shiny)
#library(clusterProfiler)
#library(pathview)
#library(ChIPseeker)
#library(ChIPpeakAnno)
#library(rtracklayer)
#library(seqinr)
library(shinythemes)
library(shinyjs)
## Load microalgae annotation packages
#library(org.Otauri.eg.db) 
#install.packages(pkgs = "./org.Otauri.eg.db/",repos = NULL,type="source")

## Load microalgae genome annotation packages
# library(TxDb.Otauri.JGI)

## Auxiliary functions
## Auxiliary function to compute enrichments
compute.enrichments <- function(gene.ratios, bg.ratios)
{
    gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
    bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
    enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
    enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
    
    return(enrichments.text)  
}

## Auxiliary function to split a string using commas
split.commas <- function(annotation.str)
{
    return(strsplit(annotation.str,split=",")[[1]])
}

## Ostreococcus tauri gene link to ORCAE
## https://bioinformatics.psb.ugent.be/orcae/annotation/OsttaV2/current/ostta15g02520
ostta.gene.link <- function(gene.name)
{
    orcae.link <- paste0("https://bioinformatics.psb.ugent.be/orcae/annotation/OsttaV2/current/",gene.name)
    gene.link <- paste(c("<a href=\"",
                         orcae.link,
                         "\" target=\"_blank\">",
                         gene.name, "</a>"),
                       collapse="")
    return(gene.link)
}


## Gene Ontology term link
# http://amigo.geneontology.org/amigo/term/GO:0015979
go.link <- function(go.term)
{
    link <- paste0("http://amigo.geneontology.org/amigo/term/", go.term)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             go.term, "</a>"),
                           collapse = "")
    return(complete.link)
}

## KO link
#https://www.genome.jp/dbget-bin/www_bget?ko:K00276
ko.link <- function(ko.term)
{
    link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", ko.term)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             ko.term, "</a>"),
                           collapse = "")
    return(complete.link)
}

## KOG link
## https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720
kog.link <- function(kog.term)
{
    link <- paste0("https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720", kog.term)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             kog.term, "</a>"),
                           collapse = "")
    return(complete.link)
}

## ENZYME link
## https://www.brenda-enzymes.org/enzyme.php?ecno=1.4.3.21
enzyme.link <- function(ec.term)
{
    link <- paste0("https://www.brenda-enzymes.org/enzyme.php?ecno=", ec.term)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             ec.term, "</a>"),
                           collapse = "")
    return(complete.link)
}


## KEGG pathway link
## https://www.genome.jp/kegg-bin/show_pathway?cre04136
kegg.pathway.link <- function(kegg.pathway)
{
    link <- paste0("https://www.genome.jp/kegg-bin/show_pathway?",kegg.pathway)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             kegg.pathway, "</a>"),
                           collapse = "")
    return(complete.link)
}

## KEGG module link
## https://www.genome.jp/kegg-bin/show_module?cre04136
kegg.module.link <- function(kegg.module)
{
    link <- paste0("https://www.genome.jp/kegg-bin/show_module?",kegg.module)
    complete.link <- paste(c("<a href=\"",
                             link,
                             "\" target=\"_blank\">",
                             kegg.module, "</a>"),
                           collapse = "")
    return(complete.link)
}


####SD plot
plot.sd <- function(gene.id, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                             paste(sd.zt,2,sep="_"),
                                                             paste(sd.zt,3,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd)
  max.expression <- max(current.gene.expression.sd)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(as.numeric(current.gene.expression.sd[1,]),type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=18),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),3),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
}

plot.sd.ll <- function(gene.id, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                             paste(sd.zt,2,sep="_"),
                                                             paste(sd.zt,3,sep="_"),
                                                             paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(as.numeric(current.gene.expression.sd.ll[1,]),type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
}

plot.sd.dd <- function(gene.id, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                             paste(sd.zt,2,sep="_"),
                                                             paste(sd.zt,3,sep="_"),
                                                             paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(as.numeric(current.gene.expression.sd.ll),type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
}

####LD plots
plot.ld.ll <- function(gene.id, gene.expression)
{
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                             paste(ld.zt,2,sep="_"),
                                                             paste(ld.zt,3,sep="_"),
                                                             paste(ld.zt,4,sep="_"),
                                                             paste(ld.zt,5,sep="_"))]
  
  min.expression <- min(current.gene.expression.ld.ll)
  max.expression <- max(current.gene.expression.ld.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(as.numeric(current.gene.expression.ld.ll),type="o",lwd=3,col="blue",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
  polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
}

# gene.id<-selected.gene
# gene.expression <- total.gene.expression
plot.ld.dd <- function(gene.id, gene.expression)
{
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld.dd <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                             paste(ld.zt,2,sep="_"),
                                                             paste(ld.zt,3,sep="_"),
                                                             paste(ld.zt,6,sep="_"),
                                                             paste(ld.zt,7,sep="_"))]
  
  min.expression <- min(current.gene.expression.ld.dd)
  max.expression <- max(current.gene.expression.ld.dd)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(as.numeric(current.gene.expression.ld.ll),type="o",lwd=3,col="blue",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
  polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
}

#### SD and LD plots
#gene.id<-selected.gene



# Define UI
ui <- shinyUI(fluidPage(#theme= "bootstrap.css",
    theme = shinytheme("flatly"),
    
    fluidRow(
        column(
            width = 2,
            img(src='logo_1.png', align = "center", width=200),
            tags$br(),
            radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                         label="",
                         choices=c(
                             "Home" = "home",
                             "Cluster Functional Analysis" = "clusters",
                             "Individual Gene Exploration" = "individual",
                            "Tutorials" = "tutorials",
                             "GitHub repository" = "github",
                             "Citation and Contact" = "citation"
                         ))),
        column(
            width = 8,
            tags$div(align = "center", 
                     tags$h1(tags$b("OsttaCIRC"), tags$br()),
                     tags$h2("microALGAE FUNctional enrichment tool for MicroAlgae RnA-seq and Chip-seq AnalysiS")),
            tags$br(),tags$br(),
            conditionalPanel(condition = "input.navigation_bar == 'home'",
                             tags$div(align = "justify", "Welcome to", tags$b("ALGAEFUN")," with ", tags$b("MARACAS"), "a microalgae web based tool for the analysis of ", 
                                      tags$b("RNA-seq"), "and ", tags$b("ChIP-seq"), "data and the", tags$b("functional annotation"), "of the resulting gene sets and genomic loci. ",
                                      tags$b("ALGAEFUN"), "with ", tags$b("MARACAS"), "supports the analysis for a wide collection 
               of microalgae that includes", tags$i("Chlamydomonas reinhardtii, Ostreococcus tauri, Phaeodactylum tricornutum"), "and ", 
                                      tags$i("Nannochlorpsis gaditana."), "Please select from the navigation bar on the left the type of analysis you want to perform. You can also 
               see our", tags$b("video tutorials"), "on how to analyse RNA-seq and
               ChIP-seq data as well as on how to functionally annotate gene sets and genomic loci. Our
               code is freely available at", tags$b("Github."), "Please cite our work if you find it useful in your research."),
                             
                             tags$div(align = "justify", "Below you can find the phylogenetic relationship between the different microalgae species supported in ALGAEFUN with MARACAS: "),
                             tags$br(),tags$br(),
                             # 
                             tags$div(align ="center",img(src='phylogeny.png', align = "center", width=600))
            ),
            
            conditionalPanel(condition = "input.navigation_bar == 'clusters'",
                             tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("functional annotation"), 
                                      "over gene sets.", tags$b("Gene Ontology (GO) enrichment"), "analysis as well as", tags$b("KEGG (Kyoto Encyclopedia
                       of Genes and Genomes) pathway enrichment"), "analysis are supported. The gene set of interest can be obtained, for example,
                       as the result of a differential expression analysis carried out using", tags$b("MARACAS."), " See our", tags$b("video tutorial"),
                                      "for details or follow the next steps to perform your analysis:",
                                      tags$ol(
                                          tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the type of enrichment analysis 
                                          to perform and the", tags$b("p-value threshold.")),
                                          tags$li("Insert your ", tags$b("gene set"), " in the text box or load it from a file using the",
                                                  tags$b("Browse …"), " button. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set."),
                                          tags$li("Users can choose between the default", tags$b("background"), " gene provided by AlgaeFUN of a custom one 
                                          that can be specified."),
                                          tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified functional enrichment analysis. The
                                          results will be shown in the different tabs below.")
                                      )
                             )),
            
            conditionalPanel(condition = "input.navigation_bar == 'individual'",
                             tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("annotation analysis 
                                of genomic loci or regions."), "These are typically generated from", tags$b("ChIP-seq"), "studies 
                                of the genome-wide distribution of", tags$b("epigenetic marks or transcription factor binding sites."),
                                      "Our tool", tags$b("MARACAS"), "can be used to perform this type of analysis. The set of marked genes 
                                can be obtained as well as the distribution of the genomic loci overlapping specific genes parts. Also, individual marked genes 
                                and the average signal level around the TSS (Transcription Start Site) and TES (Transcription
                                End Site) over the complete set of marked genes can be visualized.", " See our", tags$b("video tutorial"),
                                      "for details or follow the next steps to perform your analysis:",
                                      tags$ol(
                                          tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the gene", 
                                                  tags$b("promoter length"), "and the",  tags$b("gene parts"), "that will be considered
                                          when determining the marked genes."),
                                          tags$li("Insert in the text box your ", tags$b("set of genomic regions"), " as a table consisting 
                                          of three tab-separated columns representing the chromosome, the start and end position of 
                                          the regions. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set. Alternatively, using the",
                                                  tags$b("Browse..."), "button, the genomic regions can be uploaded from a file in BED format as 
                                          described previously containing as least three columns. This file can be obained using our tool", 
                                                  tags$b("MARACAS.")),
                                          tags$li("Optionally, users can upload the genome wide signal level of a epigenetic mark or transcription 
                                          factor binding in a BigWig file. This file can be obained using our tool", tags$b("MARACAS.")),
                                          tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified analysis. The
                                          results will be shown in the different tabs below.")
                                      )
                             )),
            

            conditionalPanel(condition = "input.navigation_bar == 'github'",
                             tags$div(align = "justify", tags$b("AlgaeFUN,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                      tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                      tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN/issues","issue")), 
                                      "in GitHub and we will address it."),
                             tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN",
                                                                           "AlgaeFUN at GitHub")))),
                             tags$br(),
                             tags$br(),
                             
                             tags$div(align = "justify", tags$b("MARACAS,"), "is developed using bash scripting and several 
        bioconductor R packages. The source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                      tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                      tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS/issues","issue")), 
                                      "in GitHub and we will address it."),
                             tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS",
                                                                           "MARACAS at GitHub")))),
            ),
            
            conditionalPanel(condition = "input.navigation_bar == 'citation'",
                             tags$div(align = "justify", "We are strongly committed to", tags$b("open access software"), 
                                      "and", tags$b("open science."),"Following our philosophy we have deposited our GitHub code 
                       into", tags$a(href="https://zenodo.org/record/4754516#.YJxLPSaxUws", target="_blank",tags$b("Zenodo")), ", a
                       general-purpose open-access repository developed under the", 
                                      tags$a(href="https://www.openaire.eu/", target="_blank", tags$b("European OpenAIRE program.")), "Meanwhile we publish 
                       our work in a journal if you find", tags$b("AlgaeFUN with MARACAS"), "useful in your research we would be most grateful if you cite 
                       our GitHub repository with a,", tags$b("DOI"),  "as follows:",
                                      tags$br(),
                                      tags$br(),
                                      tags$div(tags$b("Romero-Losada, A.B., Arvanitidou, C., de los Reyes, P., 
                                García-González, M., Romero-Campero, F.J. (2021) AlgaeFUN with MARACAS, microAlgae FUNctional 
                                enrichment tool for MicroAlgae RnA-seq and Chip-seq AnalysiS v1.0, Zenodo, doi:10.5381/zenodo.4754516 doi:10.5381/zenodo.4752818"))),
                             
                             tags$br(),
                             tags$br(),
                             #tags$div(align="center", img(src='smiley.png', align = "center", width=200,hight=200)),
                             tags$br()
                             
            ),
            
            conditionalPanel(condition = "input.navigation_bar == 'tutorials'",
                             tags$div(align="center",uiOutput("video_tutorial")),
                             tags$div(align = "justify", 
                                      tags$br(),
                                      tags$br(),
                                      tags$div(tags$h4(tags$b("Above you can find a video tutorial on how to use the different tools implemented 
                                in AlgaeFUN with MARACAS."))))
                             
            ),
            
        ),
        column(
            width = 2,
            img(src='logo_ibvf.jpg', align = "center", width=100),
            img(src='logo_us.png', align = "center", width=100),
            tags$br(),tags$br(),tags$br(),
            img(src='logo_csic.jpg', align = "center", width=100),
            tags$br(),tags$br(),
            tags$div(align="center",width=60,
                     HTML("<script type=\"text/javascript\" src=\"//rf.revolvermaps.com/0/0/8.js?i=5jamj0c2y0z&amp;m=7&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33\" async=\"async\"></script>"))
        )
    ),
    
    
    tags$br(),tags$br(),
    
    #Interface where the user can choose his/her preferencies, separated by columns
    fluidRow(
        column(width = 2),
        column(width = 4,
                  #Choose your favourite photoperiod
                 conditionalPanel(condition = "input.navigation_bar == 'clusters' || input.navigation_bar == 'individual' ",
                                 selectInput(inputId = "season", label="Choose your favourite photoperiod", 
                                            
                                            choices=c("Long days (Summer)" = "LD", 
                                                      "Short days (Winter)" = "SD",
                                                      "Both" = "cicle_comparison")
                                            )),
               
               #Choose your favourite time of the day
               conditionalPanel(condition = "input.navigation_bar == 'clusters'",
                                selectInput(inputId = "zt", label="Choose your favourite time of the day", 
                                            
                                            choices=c("Sunrise(ZT0)" = "0", 
                                                      "ZT4" = "4",
                                                      "ZT8" = "8",
                                                      "ZT12" = "12",
                                                      "ZT16" = "16", 
                                                      "ZT20" = "20"
                                                       ))),
               #Choose your favourite gene
               conditionalPanel(condition = "input.navigation_bar == 'individual'",
                                textInput(inputId= "gene",
                                          label= "Choose your favourite gene",
                                          value = "",
                                          placeholder = "ostta07g03440"))
               ),
        column(width = 4,
               
               #Transcriptome or proteome?
               conditionalPanel(condition = "input.navigation_bar == 'individual' ",    
                                #Choose the kind of analysis that you are interested in
                                radioButtons(inputId = "omics",
                                             label="Are you interested in proteomics or transcriptomics?",
                                             choices=c("Transcriptomics rules!" = "rna",
                                                       "Proteomics nerd" = "prot",
                                                       "Multi-omics integration for the win!" = "integration"
                                             ))),
               conditionalPanel(condition = "input.navigation_bar == 'individual' && input.omics == 'rna' ",    
                                #Choose the kind of analysis that you are interested in
                                radioButtons(inputId = "continuo",
                                             label="Would you like to combine your favourite photoperiod with 
                                             continuous light or darkness conditions?",
                                             choices=c("Continuous light" = "LL",
                                                       "Continuous darkness" = "DD",
                                                       "My favourite photoperiod is enough" = "3days"
                                             ))),
              #choose the analysis to execute
               conditionalPanel(condition = "input.navigation_bar == 'clusters'",    
                                #Choose the kind of analysis that you want us to execute 
                                radioButtons(inputId = "analysis",
                                             label="Choose your desirable analysis",
                                             choices=c("GO terms enrichment" = "go",
                                                       "KEGG pathways enrichment analysis" = "kegg",
                                                       "Both" = "both"
                                             ))),
               
               conditionalPanel(condition= "(input.analysis == 'go' || input.analysis == 'both') &&
                                   input.navigation_bar == 'clusters'",
                                radioButtons(inputId = "ontology",
                                             label="Choose gene ontology:",
                                             choices = c("Biological process" = "BP",
                                                         "Cellular Component" = "CC",
                                                         "Mollecular Function" = "MF"))),
               #Choose a p-value
               conditionalPanel(condition = "input.navigation_bar == 'clusters'",
                                
                                numericInput(inputId = "pvalue", 
                                             label= "Which will be your chosen p-value?", 
                                             value= 0.05)),
                #Button for functional enrichment over the chosen cluster of genes.
                conditionalPanel(condition = "input.navigation_bar == 'clusters'",
                                 actionButton(inputId = "go.button",label = "Have fun!", icon("send") )                       
                                ),
              conditionalPanel(condition = "input.navigation_bar == 'individual'",
                               actionButton(inputId = "circ.button",label = "Have fun!", icon("send") )                       
              )
               )
              ),
              
        
    tags$br(), tags$br(),
    mainPanel(width = 13,
              
              #Main panel containing the results organized in different tabs: GO map, Go terms data table, and 
              #KEGG pathway maps for gene set enrichment analysis
              
              conditionalPanel(condition = "input.navigation_bar == 'clusters'",
                               tabsetPanel(type ="tabs",
                                           tabPanel(tags$b("GO ENRICHMENT"),
                                                    shinyjs::useShinyjs(),
                                                    hidden(div(id='loading.enrichment.go',h3('Please be patient, computing GO enrichment ...'))), 
                                                    hidden(div(id='ready.enrichment.go',h3('Your GO enrichment is ready!'))), 
                                                    htmlOutput(outputId = "gene_sanity_go"),
                                                    htmlOutput(outputId = "wrong_genes_go"),
                                                    tags$br(),
                                                    htmlOutput(outputId = "intro_go"),
                                                    tags$br(),
                                                    tags$br(),
                                                    tabsetPanel(type = "tabs",
                                                                tabPanel(tags$b("GO Enrichment Table"),
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "textGOTable"),
                                                                         tags$br(), tags$br(),
                                                                         dataTableOutput(outputId = "output_go_table"),
                                                                         uiOutput(outputId = "download_ui_for_go_table")#,
                                                                         #htmlOutput(outputId = "revigo")
                                                                ),
                                                                tabPanel(tags$b("GO Map"),
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "go_graph"),
                                                                         tags$br(), tags$br(),
                                                                         div(style= "overflow:scroll; height:500px; text-align: center;", 
                                                                             plotOutput(outputId = "go.plot", inline = T)),
                                                                         tags$br(),
                                                                         tags$br(), tags$br()),
                                                                tabPanel(tags$b("GO Barplot"),
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "barplot_text"),
                                                                         tags$br(),
                                                                         div(style= "text-align: center;",
                                                                             plotOutput(outputId = "bar.plot",inline=TRUE)),
                                                                         tags$br(),
                                                                         tags$br(), tags$br()),
                                                                tabPanel(tags$b("GO Dotplot"),
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "dotplot_text"),
                                                                         tags$br(),
                                                                         div(style= "text-align: center;",
                                                                             plotOutput(outputId = "dot.plot",inline=TRUE)),
                                                                         tags$br(),
                                                                         tags$br(), tags$br()),
                                                                tabPanel(tags$b("GO Concept Map"),
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "cnetplot_text"),
                                                                         tags$br(),
                                                                         div(style= "text-align: center;",
                                                                             plotOutput(outputId = "cnet.plot",inline=TRUE)),
                                                                         tags$br(),
                                                                         tags$br(),tags$br())
                                                    ) # close tabsetPanel for go result
                                           ), # close tabPanel for KEGG ENRICHMENT
                                           tabPanel(tags$b("KEGG PATHWAY ENRICHMENT"),
                                                    shinyjs::useShinyjs(),
                                                    hidden(div(id='loading.enrichment.kegg',h3('Please be patient, computing KEGG pathway enrichment ...'))), 
                                                    hidden(div(id='ready.enrichment.kegg',h3('Your KEGG enrichment is ready!'))), 
                                                    htmlOutput(outputId = "gene_sanity_kegg"),
                                                    htmlOutput(outputId = "wrong_genes_kegg"),
                                                    tags$br(),
                                                    htmlOutput(outputId = "intro_kegg"),
                                                    tags$br(), tags$br(),
                                                    tabsetPanel(type = "tabs",
                                                                tabPanel(tags$b("KEGG Pathway Enrichment Table"), 
                                                                         tags$br(),
                                                                         htmlOutput(outputId = "textKEGGTable"),
                                                                         htmlOutput(outputId = "kegg_pathway_table_text"),
                                                                         tags$br(),
                                                                         dataTableOutput(outputId = "output_pathway_table"),
                                                                         br(), br(), br()),
                                                                tabPanel(tags$b("KEGG Pathway Visualization"),
                                                                         tags$br(), tags$br(),
                                                                         htmlOutput(outputId = "textKEGGImage"),
                                                                         br(), br(),
                                                                         uiOutput(outputId = "kegg_selectize"),
                                                                         div(style= "text-align: center;",
                                                                             imageOutput("kegg_image", inline = T)),
                                                                         tags$br(), tags$br(),
                                                                         br(), br()),
                                                                tabPanel(tags$b("KEGG Module Enrichment Table"),
                                                                         tags$br(), tags$br(),
                                                                         htmlOutput(outputId = "text_module_kegg"),
                                                                         br(), br(),
                                                                         dataTableOutput(outputId = "output_module_table"),
                                                                         tags$br(), tags$br(),
                                                                         tags$br(), tags$br(),
                                                                         uiOutput(outputId = "kegg_module_selectize"),
                                                                         imageOutput("kegg_module_image"),
                                                                         tags$br(), tags$br())
                                                    )
                                                    
                                           )
                               )
              ),
              
              ## Main panel containing the results for graphical representation of expression/amount of protein profiles 
              ##and statistical analysis of their rhythmicity
              
              conditionalPanel(condition = "input.navigation_bar == 'individual'",
                               hidden(div(id='loading.circ',h3('Please be patient, we are working on it ...'))), 
                               hidden(div(id='ready.circ',h3('Here are the results!!'))),
                               tabsetPanel(type = "tabs",
                                           tabPanel(tags$b("Graphical representation"),
                                                    tags$br(),
                                                    div(style= "text-align: center;",
                                                        plotOutput(outputId = "circadian.plot",inline=TRUE))
                                                    ),
                                           tabPanel(tags$b("Statistical analysis"),
                                                    tags$br(), tags$br(),
                                                    dataTableOutput(outputId = "output_statistical_table"),
                                                    uiOutput(outputId = "download_ui_for_statistical_table"),
                                                    tags$br(),
                                                    div(style= "text-align: center;",
                                                        plotOutput(outputId = "circacompare",inline=TRUE)),
                                                    tags$br(), tags$br(),
                                                    br(), br())
                                           ))
             
                                                    
)))

server <- shinyServer(function(input, output, session) {
    
    ## video tutorial
    output$video_tutorial <- renderUI({
        HTML("<iframe width=\"560\" height=\"315\" src=\"https://www.youtube.com/embed/ZCWrqOxrdJM\" title=\"YouTube video player\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture\" allowfullscreen></iframe>")
    })  
  
 
    ## Actions to perform after click the go button
    observeEvent(input$go.button , {
        shinyjs::showElement(id = 'loading.enrichment.go')
        shinyjs::hideElement(id = 'ready.enrichment.go')
        
        # Remove previous results
        output$intro_go <- renderText(expr = "")
        output$gene_sanity_go <- renderText(expr = "")
        output$wrong_genes_go  <- renderText(expr = "")
        output$textGOTable <- renderText(expr = "")
        output$output_go_table <- renderDataTable("")
        output$download_ui_for_go_table<- renderUI(expr = NULL)
        #output$revigo<- renderUI(expr = NULL)
        output$go_graph <- renderText(expr = "")
        output$go.plot <- renderPlot(expr = NULL)
        output$barplot_text <- renderText("")
        output$bar.plot <- renderPlot(expr = NULL)
        output$dotplot_text <- renderText("")
        output$dot.plot <- renderPlot(expr = NULL)
        #output$emapplot_text <- renderText("")
        #output$emap.plot <- renderPlot(expr = NULL)
        output$cnetplot_text <- renderText("")
        output$cnet.plot <- renderPlot(expr = NULL)
        
        output$intro_kegg <- renderText(expr = "")
        output$output_pathway_table <- renderDataTable(expr = NULL)
        output$textKEGGImage <- renderText(expr = "")
        output$kegg_selectize <- renderUI(expr = NULL)
        output$kegg_image <- renderImage(expr = NULL,deleteFile = T)
        output$text_module_kegg <- renderText(expr = "")
        output$output_module_table <- renderDataTable(expr = NULL)
        
        # Load libraries
        library(clusterProfiler)
        library(pathview)
        library(org.Otauri.eg.db)
            
        org.db <- org.Otauri.eg.db
        microalgae.genes <- read.table(file = "otauri_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
        gene.link.function <- ostta.gene.link
        
        ## Extract genes from cluster text files
        file.name <- paste(c("cluster_","peak_",as.character(input$zt), ".txt"), collapse="")
        path <- paste(c("clusters_",input$season), collapse="")
        complete_path<- paste(c(path,file.name),collapse="/")
        target.genes <- read.table(file=complete_path, header=F,as.is = T, comment.char="")
      
        ## GO term enirchment analysis
        if(input$analysis == "kegg")
        {
            output$intro_go <- renderText(expr = "<p style=\"color:blue\"><b> You have chosen only to perform a KEGG enrichment analysis.
                                    Please check the content of the KEGG ENRICHMENT tab.</b></p>")
        } else if((input$analysis == "go" || input$analysis == "both") && (length(target.genes$V1) == 0))
        {
            output$intro_go <- renderText(expr = "<p style=\"color:red\"><b> You forgot to select your favourite time of the day. Please, select one. p</b></p>")
        } else if((input$analysis == "go" || input$analysis == "both") && (length(target.genes$V1) > 0))
        {
            ## Intro text for GO enrichment
            go.intro.text <- paste(c("The tabs below present the results from the <b>GO enrichment analysis</b> 
                                      performed over the genes that peak at (<i>ZT", input$zt,
                                     "</i> ...) from the microalgae <b> <i> Ostreococcus tauri </i> </b>"),collapse="") 
            
            output$intro_go <- renderText(expr = go.intro.text)
            
            ## Perform GO enrichment
            enrich.go <- enrichGO(gene          = target.genes$V1,
                                  universe      = microalgae.genes,
                                  OrgDb         = org.db,
                                  ont           = input$ontology,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = input$pvalue,
                                  readable      = TRUE,
                                  keyType = "GID")
            
            
            ## Generate ouput table
            enrich.go.result <- as.data.frame(enrich.go)
            
            if(nrow(enrich.go.result) > 0)
            {
                ## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
                go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                                           bg.ratios = enrich.go.result$BgRatio)
                
                go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                                              enrich.go.result$pvalue, enrich.go.result$qvalue,
                                              go.term.enrichments, 
                                              gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                                              stringsAsFactors = FALSE)
                
                colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                                               "Enrichment (Target Ratio; BG Ration)","Genes")
                
                go.result.table.with.links <- go.result.table
                ## Add links to the genes
                genes.in.go.enrichment <- go.result.table$Genes
                
                ## Add link to genes
                for(i in 1:length(genes.in.go.enrichment))
                {
                    go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
                }
                
                ## Add links to GO ids
                go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
                
                ## Introductory text for GO enrichment table
                go.table.text <- "The table below summarizes the result of the GO term
enrichment analysis. Each row represents a GO term significantly enriched in the target
gene set with respect to the selected gene universe. The first column represents the GO term
identifier. The second column contains a human readable description. For more details on the 
corresponding GO term, click on the identifier in the first column. The third and fourth 
column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
n is the number of genes with annotation from the target set, N is the number of genes with
annotation from the gene universe, m is the number of genes from the target set annotated with the
corresponding GO term and M is the number of genes from the gene universe annotated with
the GO term associated with the corresponding row. The enrichment is then computed as
E = (m/n) / (M/N). Finally, the last column, contains the genes from the target set
annotated with the GO term represented in the corresponding row."
                output$textGOTable <- renderText(expr = go.table.text)
                
                ## Output table with GO enrichment result
                output$output_go_table <- renderDataTable({
                    go.result.table.with.links #go.result.table
                },escape=FALSE,options =list(pageLength = 5))
                
                ## Generate UI to download go enrichment table and creating a downlodable table
                output$download_ui_for_go_table<- renderUI(
                    tagList(downloadButton(outputId= "downloadGOTable", "Download GO Enrichment Table"),tags$br(),tags$br())
                )
                
                ## Download result
                output$downloadGOTable<- downloadHandler(
                    filename= function() {
                        paste("go_enrichment_table_zt",input$zt ,"_",input$season,".tsv", sep="")
                    },
                    content= function(file) {
                        write.table(x = go.result.table,quote = F,sep = "\t",
                                    file=file,row.names=FALSE,col.names=TRUE)
                    })
                
                
            go.graph.text <- "The following acyclic graph represents the GO term enrichment
        in the target gene set. Each node stands for a GO term. The color of each node
        indicates the level of significance from grey, non-significant, to intense red,
        highly significant. An arrow is drawn from GO term A to GO term B when A is a more
        general GO term than B or B is more specific than A. Right click on the image to download it."
                
                output$go_graph <- renderText(expr = go.graph.text)
                
                ## GO plot
                output$go.plot <- renderPlot(
                    width     = 940,
                    height    = 900,
                    res       = 120,
                    expr = {
                        goplot(enrich.go,showCategory = 10)
                    })
                
                output$barplot_text <- renderText("In the following barplot each bar represents a significantly enriched 
        GO term. The length of the bar corresponds to the number of genes in the
        target set annotated with the given GO term. The bar color captures the level
        of significance from blue, less significant, to red, more significant. Right click on the image to download it.")
                
                ## Barplot
                output$bar.plot <- renderPlot(
                    width     = 870,
                    height    = 600,
                    res       = 120,
                    expr = {
                        barplot(enrich.go,drop=TRUE,showCategory = 10)
                    })
                
                output$dotplot_text <- renderText("In the following dotplot each dot represents a significantly enriched 
        GO term. The x-position of the dot corresponds to the ratio between the number of genes annotated with the
corresponding GO term and the total number of annotated genes in the target set. The dot color captures the level
        of significance from blue, less significant, to red, more significant. Right click on the image to download it.")
                
                ## Dotplot
                output$dot.plot <- renderPlot(
                    width     = 870,
                    height    = 600,
                    res       = 120,
                    expr = {
                        dotplot(enrich.go)
                    })
              
                output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term. Right click on the image to download it.")
                
                
                ##CNET plot
                output$cnet.plot <- renderPlot(
                    width     = 870,
                    height    = 600,
                    res       = 120,
                    expr = {
                        cnetplot(enrich.go)
                    })
            } else
            {
                output$textGOTable <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
                output$go_graph <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
                output$barplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
                output$dotplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")  
                output$cnetplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")        
            }
            
            shinyjs::hideElement(id = 'loading.enrichment.go')
            shinyjs::showElement(id = 'ready.enrichment.go')
        }
        
        ## KEGG pathways enrichment analysis
        if(input$analysis == "go")
        {
            output$intro_kegg <- renderText(expr = "<p style=\"color:blue\"><b> You have chosen only to perform a GO enrichment analysis.
                                    Please check the content of the GO ENRICHMENT tab.</b></p>")
        } else if((input$analysis == "kegg" || input$analysis == "both") && (length(target.genes$V1) == 0))
        {
            output$intro_kegg <- renderText(expr = "<p style=\"color:red\"><b> You forgot to select your favourite time of the day. Please, select one. p</b></p>")
        } else if( (input$analysis == "kegg"  || input$analysis == "both") && (length(target.genes$V1) > 0))
        {
            shinyjs::showElement(id = 'loading.enrichment.kegg')
            shinyjs::hideElement(id = 'ready.enrichment.kegg')
            
            ## Update target genes and universe depending on the microalgae
                target.genes <- paste0("OT_",target.genes$V1)
                gene.universe <- paste0("OT_",microalgae.genes)
                organism.id <- "ota"
            
            ## Compute KEGG pathway enrichment
               pathway.enrichment <- enrichKEGG(gene = target.genes, organism = organism.id, keyType = "kegg",
                                                 universe = gene.universe,qvalueCutoff = input$pvalue)
            }
            shinyjs::showElement(id = 'ready.enrichment.kegg')
            shinyjs::hideElement(id = 'loading.enrichment.kegg')
            
            pathway.enrichment.result <- as.data.frame(pathway.enrichment)
            if(nrow(pathway.enrichment.result) > 0)
            {
                kegg.intro.text <- paste(c("The tabs below present the results from the <b>KEGG pathways/modules enrichment analysis</b> 
                                      performed over the genes that peak at (<i>ZT", input$zt,
                                           "</i> ...) from the microalgae <b> <i> Ostreococcus tauri </i> </b>"),collapse="") 
                output$intro_kegg <- renderText(expr = kegg.intro.text)
                
                pathways.enrichment <- compute.enrichments(gene.ratios = pathway.enrichment.result$GeneRatio,
                                                           bg.ratios = pathway.enrichment.result$BgRatio)
                
                ## Generate gene lists for enriched pathways
                 kegg.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = pathway.enrichment.result$geneID))
               
                pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                                    pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                                    pathways.enrichment, 
                                                    kegg.enriched.genes,
                                                    stringsAsFactors = FALSE)
                
                colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                                     "Enrichment (Target Ratio; BG Ration)","Genes")
                
                kegg.result.table.with.links <- pathways.result.table
                
                for(i in 1:length(kegg.enriched.genes))
                {
                    kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
                }
                
                ## Add links to kegg pathways
                kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)
                
                
                ## Introductory text for GO enrichment table
                kegg.table.text <- "The table below summarizes the result of the KEGG pathway
enrichment analysis. Each row represents a pathway significantly enriched in the target
gene set with respect to the selected gene universe. The first column represents the KEGG pathway 
identifier. The second column contains a human readable description. For more details on the 
corresponding pathway, click on the identifier in the first column. The third and fourth 
column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
n and N are the number of genes associated to any pathway in the target set and in the gene universe
respectively; m and M are the number of genes associated to the pathway represented in the corresponding
row in the target gene set and in the gene universe respectively. The enrichment is then computed as
E = (m/n) / (M/N). Finally, the last column, contains the genes from the target gene set
assocated to the enriched pathway represented in the corresponding row." 
                
                output$textKEGGTable <- renderText(expr = kegg.table.text)
                
                output$output_pathway_table <- renderDataTable({
                    kegg.result.table.with.links
                },escape=FALSE,options =list(pageLength = 5))
                
                
                ## Figures for KEGG pathway enrichment analysis
                
                ## Prepare gene set for representation
                   genes.pathway <- rep(0, length(microalgae.genes))
                   names(genes.pathway) <- microalgae.genes
                    
                  genes.pathway[target.genes$V1] <- 1
               
                pathways.for.select <- paste(pathways.result.table[["KEGG ID"]], pathways.result.table[["Description"]], sep=" - ")
                
                kegg.image.text <- "<b> The enriched pathways detected above can be visualized using the dropdown menu below. 
      Genes in the target set associated to the corresponding pathway will be highlighted as red rectangles: </b>"
                
                output$textKEGGImage <- renderText(expr = kegg.image.text)
                
                output$kegg_selectize <- renderUI({
                    selectInput(inputId = "kegg_pathway", label="Choose Pathway for Representation",multiple = FALSE,selected = pathways.for.select[1],
                                choices=pathways.for.select)
                })  
            } else
            {
                output$textKEGGTable <- renderText(expr = "<b>No significant KEGG 
                                           pathway enrichment detected in 
                                           the input gene set.")
            }
            
            #KEGG modules enrichment  
            modules.enrichment <- enrichMKEGG(gene = target.genes$V1, universe = microalgae.genes, organism = organism.id, keyType = "kegg",minGSSize = 4)
            
          
            modules.enrichment.result <- as.data.frame(modules.enrichment)
            if(nrow(modules.enrichment.result) > 0)
            {
                modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                          bg.ratios = modules.enrichment.result$BgRatio)
                modules.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = modules.enrichment.result$geneID))
                
                modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                                   modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                                   modules.enrichment, 
                                                   modules.enriched.genes,
                                                   stringsAsFactors = FALSE)
                
                colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                                    "Enrichment (Target Ratio; BG Ration)","Genes")
                
                modules.result.table.with.links <- modules.result.table
                
                for(i in 1:length(modules.enriched.genes))
                {
                    modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
                }
                
                ## Add links to kegg pathways
                modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
                
                
                ## Introductory text for GO enrichment table
                module.table.text <- "The table below summarizes the result of the KEGG module
        enrichment analysis. Modules are sometimes easier to interpret than pathways. Each row represents a module significantly enriched in the target
        gene set with respect to the selected gene universe. The first column represents the KEGG module 
        identifier. The second column contains a human readable description. For more details on the 
        corresponding module, click on the identifier in the first column. The third and fourth 
        column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
        of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
        n and N are the number of genes associated to any module in the target set and in the gene universe
        respectively; m and M are the number of genes associated to the module represented in the corresponding
        row in the target gene set and in the gene universe respectively. The enrichment is then computed as
        E = (m/n) / (M/N). Finally, the last column, contains the genes from the target gene set
        associated to the enriched module represented in the corresponding row." 
                
                output$text_module_kegg <- renderText(expr = module.table.text)
                
                output$output_module_table <- renderDataTable({
                    modules.result.table.with.links
                },escape=FALSE,options =list(pageLength = 5)) 
            } else 
            {
                output$text_module_kegg <- renderText(expr = "<b>No significant KEGG 
                                              module enrichment detected in the 
                                              input gene set.")
            }
         
        
        enriched.pathway.id <- reactive({ 
            if(is.null(input$kegg_pathway))
            {
                return()
            } else
            {
                return(strsplit(input$kegg_pathway,split=" - ")[[1]][1] )
            }
        })
        
        observeEvent(enriched.pathway.id(), {
            
          
             organism.id <- "ota"
             output$kegg_image <- renderImage({
                pathview(gene.data = sort(genes.pathway,decreasing = TRUE),kegg.dir = "pathways",
                         pathway.id = enriched.pathway.id(),
                         species = organism.id,
                         limit = list(gene=max(abs(genes.pathway)), cpd=1),
                         gene.idtype ="kegg")
                
                list(src = paste(c(enriched.pathway.id(),"pathview","png"), collapse="."),
                     contentType="image/png",width=1200,height=900)
            },deleteFile = T)
        })
        
    })
   
    
    observeEvent(input$circ.button, {
  shinyjs::showElement(id = 'loading.circ')
  shinyjs::hideElement(id = 'ready.circ')
  
  #Remove previous results
  #output$intro_go <- renderText(expr = "")
  output$output_statistical_table <- renderDataTable("")
  output$download_ui_for_statistical_table<- renderUI(expr = NULL)
  output$circadian.plot <- renderPlot(expr = NULL)
  output$circacompare <- renderPlot(expr = NULL)
  
  if(input$omics == "rna" || input$omics == "integration")
  {
    #extract gene expression levels of the target gene.
  selected.gene <- as.character(input$gene)
  #selected.gene <-"ostta01g00060"
  total.gene.expression <- read.table(file = "gene_expression.tsv", header =T)
  gene.expression<- total.gene.expression[selected.gene,]
  
  }else if (input$omics == "prot" || input$omics == "integration")
  {
    #Extract protein abundance levels of the corresponding target gene
    target.prot <- as.character(input$gene)
    #target.prot <-"ostta01g00060"
    swath.normalized.data.SD <- read.table(file = "sd_swath_processed_data.tsv",header=T,sep="\t")
    rownames(swath.normalized.data.SD)<- swath.normalized.data.SD$X
    swath.normalized.data.SD$X <- NULL
    
    swath.normalized.data.LD <- read.table(file = "swath_processed_data.tsv",header=T,sep="\t")
    swath.normalized.data.LD[,"zt20_2"] <- swath.normalized.data.LD[,"zt16_2"]
    swath.normalized.data.LD[,"zt16_2"] <- swath.normalized.data.LD[,"zt20_2"]
    
    zts.omics <- c(paste(paste(paste("zt", seq(from=0, to=20, by=4), sep = ""), sep = "_"), 1, sep = "_"),
                   paste(paste(paste("zt", seq(from=0, to=20, by=4), sep = ""), sep = "_"), 2, sep = "_"),
                   paste(paste(paste("zt", seq(from=0, to=20, by=4), sep = ""), sep = "_"), 3, sep = "_"))
    proteins.carotenoids.SD <- na.omit(swath.normalized.data.SD[target.prot,])
    proteins.carotenoids.SD <- proteins.carotenoids.SD[,zts.omics]
    proteins.carotenoids.LD <- na.omit(swath.normalized.data.LD[target.prot,])
    proteins.carotenoids.LD <- proteins.carotenoids.LD[,zts.omics]
  }
  
  if(input$omics == "rna")
  {
    if(input$season == "SD" && input$continuo == "3days")
    {
      gene.expression.SD <- gene.expression[,43:60]
      output$circadian.plot<- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.sd(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
          
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("sd_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.SD[,new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      
      rain.results <- matrix(ncol=3, nrow=1)
      rownames(rain.results) <- c("SD")
      colnames(rain.results) <- c("", "Period 24h", "Period 12h")
      rain.results[,1] <-"Rhythmicity under short day (SD) conditions"
      rain.results["SD","Period 24h"] <- rain24.sd$pVal
      rain.results["SD","Period 12h"] <- rain12.sd$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results 
      },escape=FALSE,options =list(pageLength = 5))
      
    }else if (input$season == "SD" && input$continuo == "LL")
    {
      gene.expression.SD.LL <- gene.expression[,43:72]
      output$circadian.plot<- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.sd.ll(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
          
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("sd_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.SD.LL[,new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      ###rain for LL conditions
      new.time.points.order <- c(paste("sd_", "zt00_", seq(from=2, to=5), sep=""),
                                 paste("sd_", "zt04_", seq(from=2, to=5), sep=""),
                                 paste("sd_", "zt08_", seq(from=2, to=5), sep=""),
                                 paste("sd_", "zt12_", seq(from=2, to=5), sep=""),
                                 paste("sd_", "zt16_", seq(from=2, to=5), sep=""),
                                 paste("sd_", "zt20_", seq(from=2, to=5), sep=""))
      gene.expression.rain <- gene.expression.SD.LL[,new.time.points.order]
     
      new.rain.order <-c(paste( "zt0_", seq(from=2, to=5), sep=""),
                         paste( "zt4_", seq(from=2, to=5), sep=""),
                         paste( "zt8_", seq(from=2, to=5), sep=""),
                         paste( "zt12_", seq(from=2, to=5), sep=""),
                         paste( "zt16_", seq(from=2, to=5), sep=""),
                         paste( "zt20_", seq(from=2, to=5), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      
      library(rain)
      rain24.sd.ll<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.sd.ll<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      rain.results <- matrix(ncol=3, nrow=2)
      rownames(rain.results) <- c("SD", "SD+LL")
      colnames(rain.results) <- c("","Period 24h", "Period 12h")
      rain.results[,1] <- c("Rhythmicity under short day (SD) conditions", 
                            "Rhythmicity under constant light conditions (after SD)")
      rain.results["SD","Period 24h"] <- rain24.sd$pVal
      rain.results["SD","Period 12h"] <- rain12.sd$pVal
      rain.results["SD+LL","Period 24h"] <- rain24.sd.ll$pVal
      rain.results["SD+LL","Period 12h"] <- rain12.sd.ll$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results #go.result.table
      },escape=FALSE,options =list(pageLength = 5))
      
    }else if (input$season == "SD" && input$continuo == "DD")
    {
      gene.expression.SD.DD <- c(as.numeric(gene.expression[,43:60]), as.numeric(gene.expression[,73:84]))
      names <-c(names(gene.expression[,43:60]), names(gene.expression[,73:84]))
      names(gene.expression.SD.DD) <- names
      output$circadian.plot<- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.sd.dd(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
          
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("sd_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("sd_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.SD.DD[new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      names(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.sd<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      ###rain for DD conditions
      new.time.points.order <- c(paste("sd_", "zt00_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt00_", seq(from=6, to=7), sep=""),
                                 paste("sd_", "zt04_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt04_", seq(from=6, to=7), sep=""),
                                 paste("sd_", "zt08_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt08_", seq(from=6, to=7), sep=""),
                                 paste("sd_", "zt12_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt12_", seq(from=6, to=7), sep=""),
                                 paste("sd_", "zt16_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt16_", seq(from=6, to=7), sep=""),
                                 paste("sd_", "zt20_", seq(from=2, to=3), sep=""),
                                 paste("sd_", "zt20_", seq(from=6, to=7), sep=""))
      gene.expression.rain <- gene.expression.SD.DD[new.time.points.order]
      
      new.rain.order <-c(paste( "zt0_", seq(from=2, to=5), sep=""),
                         paste( "zt4_", seq(from=2, to=5), sep=""),
                         paste( "zt8_", seq(from=2, to=5), sep=""),
                         paste( "zt12_", seq(from=2, to=5), sep=""),
                         paste( "zt16_", seq(from=2, to=5), sep=""),
                         paste( "zt20_", seq(from=2, to=5), sep=""))
      names(gene.expression.rain)<- new.rain.order
      
      library(rain)
      rain24.sd.dd<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.sd.dd<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      rain.results <- matrix(ncol=3, nrow=2)
      rownames(rain.results) <- c("SD", "SD+DD")
      colnames(rain.results) <- c("","Period 24h", "Period 12h")
      rain.results[,1] <- c("Rhythmicity under short day (SD) conditions", 
                            "Rhythmicity under constant darkness conditions (after SD)")
      rain.results["SD","Period 24h"] <- rain24.sd$pVal
      rain.results["SD","Period 12h"] <- rain12.sd$pVal
      rain.results["SD+DD","Period 24h"] <- rain24.sd.dd$pVal
      rain.results["SD+DD","Period 12h"] <- rain12.sd.dd$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results #go.result.table
      },escape=FALSE,options =list(pageLength = 5))
      
      
    }else if (input$season == "LD" && input$continuo == "3days")
    { 
      gene.expression.LD <- gene.expression[,1:18]
      output$circadian.plot <- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.ld.ll(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("ld_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.LD[,new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      rain.results <- matrix(ncol=3, nrow=1)
      rownames(rain.results) <- c("LD")
      colnames(rain.results) <- c("","Period 24h", "Period 12h")
      rain.results[,1] <- "Rhythmicity under long day conditions"
      rain.results["LD","Period 24h"] <- rain24.ld$pVal
      rain.results["LD","Period 12h"] <- rain12.ld$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results 
      },escape=FALSE,options =list(pageLength = 5))
      
      
       }else if (input$season== "LD" && input$continuo == "LL")
    {
      gene.expression.LD.LL <- gene.expression[,1:30]
      output$circadian.plot <- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.ld.ll(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("ld_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.LD.LL[,new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      ###rain for DD conditions
      new.time.points.order <- c(paste("ld_", "zt00_", seq(from=2, to=5), sep=""),
                                 paste("ld_", "zt04_", seq(from=2, to=5), sep=""),
                                 paste("ld_", "zt08_", seq(from=2, to=5), sep=""),
                                 paste("ld_", "zt12_", seq(from=2, to=5), sep=""),
                                 paste("ld_", "zt16_", seq(from=2, to=5), sep=""),
                                 paste("ld_", "zt20_", seq(from=2, to=5), sep=""))
      gene.expression.rain <- gene.expression.LD.LL[,new.time.points.order]
      
      new.rain.order <-c(paste( "zt0_", seq(from=2, to=5), sep=""),
                         paste( "zt4_", seq(from=2, to=5), sep=""),
                         paste( "zt8_", seq(from=2, to=5), sep=""),
                         paste( "zt12_", seq(from=2, to=5), sep=""),
                         paste( "zt16_", seq(from=2, to=5), sep=""),
                         paste( "zt20_", seq(from=2, to=5), sep=""))
      colnames(gene.expression.rain) <- new.rain.order
      
      library(rain)
      rain24.ld.ll<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.ld.ll<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      rain.results <- matrix(ncol=3, nrow=2)
      rownames(rain.results) <- c("LD", "LD+LL")
      colnames(rain.results) <- c("","Period 24h", "Period 12h")
      rain.results[,1] <- c("Rhythmicity under long day (LD) conditions", 
                            "Rhythmicity under constant light conditions (after LD)")
      rain.results["LD","Period 24h"] <- rain24.ld$pVal
      rain.results["LD","Period 12h"] <- rain12.ld$pVal
      rain.results["LD+LL","Period 24h"] <- rain24.ld.ll$pVal
      rain.results["LD+LL","Period 12h"] <- rain12.ld.ll$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results 
      },escape=FALSE,options =list(pageLength = 5))
      
          
      
    }else if (input$season== "LD" && input$continuo == "DD")
    {
      
      gene.expression.LD.DD <- c(as.numeric(gene.expression[,1:18]), as.numeric(gene.expression[,31:42]))
      names <-c(names(gene.expression[,1:18]), names(gene.expression[,31:42]))
      names(gene.expression.LD.DD) <- names
      output$circadian.plot <- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          plot.ld.dd(gene.id=selected.gene, 
                     gene.expression=total.gene.expression)
        })
      #Prepare data for rain analysis
      new.time.points.order <- c(paste("ld_", "zt00_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt04_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt08_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt12_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt16_", seq(from=1, to=3), sep=""),
                                 paste("ld_", "zt20_", seq(from=1, to=3), sep=""))
      gene.expression.rain <- gene.expression.LD.DD[new.time.points.order]
      head(gene.expression.rain)
      
      
      new.rain.order <-c(paste( "zt0_", seq(from=1, to=3), sep=""),
                         paste( "zt4_", seq(from=1, to=3), sep=""),
                         paste( "zt8_", seq(from=1, to=3), sep=""),
                         paste( "zt12_", seq(from=1, to=3), sep=""),
                         paste( "zt16_", seq(from=1, to=3), sep=""),
                         paste( "zt20_", seq(from=1, to=3), sep=""))
      names(gene.expression.rain) <- new.rain.order
      library(rain)
      rain24.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.ld<- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      ###rain for DD conditions
      new.time.points.order <- c(paste("ld_", "zt00_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt00_", seq(from=6, to=7), sep=""),
                                 paste("ld_", "zt04_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt04_", seq(from=6, to=7), sep=""),
                                 paste("ld_", "zt08_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt08_", seq(from=6, to=7), sep=""),
                                 paste("ld_", "zt12_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt12_", seq(from=6, to=7), sep=""),
                                 paste("ld_", "zt16_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt16_", seq(from=6, to=7), sep=""),
                                 paste("ld_", "zt20_", seq(from=2, to=3), sep=""),
                                 paste("ld_", "zt20_", seq(from=6, to=7), sep=""))
      gene.expression.rain <- gene.expression.LD.DD[new.time.points.order]
      
      new.rain.order <-c(paste( "zt0_", seq(from=2, to=5), sep=""),
                         paste( "zt4_", seq(from=2, to=5), sep=""),
                         paste( "zt8_", seq(from=2, to=5), sep=""),
                         paste( "zt12_", seq(from=2, to=5), sep=""),
                         paste( "zt16_", seq(from=2, to=5), sep=""),
                         paste( "zt20_", seq(from=2, to=5), sep=""))
      names(gene.expression.rain) <- new.rain.order
      
      library(rain)
      rain24.ld.dd <- rain(as.numeric(gene.expression.rain), deltat=4, period=24, verbose=T, nr.series=3)
      rain12.ld.dd <- rain(as.numeric(gene.expression.rain), deltat=4, period=12, verbose=T, nr.series=3)
      
      rain.results <- matrix(ncol=3, nrow=2)
      rownames(rain.results) <- c("LD", "LD+DD")
      colnames(rain.results) <- c("","Period 24h", "Period 12h")
      rain.results[,1] <-c("Rhythmicity under long day (LD) conditions",
                           "Rhythmicity under constant darkness conditions (after LD)")
      rain.results["LD","Period 24h"] <- rain24.ld$pVal
      rain.results["LD","Period 12h"] <- rain12.ld$pVal
      rain.results["LD+DD","Period 24h"] <- rain24.ld.dd$pVal
      rain.results["LD+DD","Period 12h"] <- rain12.ld.dd$pVal
      
      output$output_statistical_table <- renderDataTable({
        rain.results 
      },escape=FALSE,options =list(pageLength = 5))
      
      
    }else if (input$season== "cicle_comparison" && input$continuo == "3days")
    {
      gene.expression.SD <- gene.expression[,43:60]
      gene.expression.LD <- gene.expression[,1:18]
      
      # output$circadian.plot<- renderPlot(
      #   width     = 870,
      #   height    = 600,
      #   res       = 120,
      #   expr = {
      #     plot.sd.ll(gene.id=selected.gene, 
      #                gene.expression=total.gene.expression)
      #     
        # })
   #####Circacompare analysis
      library(circacompare)
      time.points <- seq(from=0,by=4,length.out = 18)
      circacompare.SD.LD <- matrix(nrow=15,ncol=2)
      current.gene <- selected.gene
      circacomp.data <- data.frame(time=c(time.points,time.points),
                                   measure=c(t(gene.expression.LD/max(gene.expression.LD)),
                                             t(gene.expression.SD/max(gene.expression.SD))),
                                   group=c(rep("Selected gene under LD conditions",18),rep("selected gene under SD conditions",18)))
      
      result.i<- circacompare(x = circacomp.data, 
                              col_time = "time", 
                              col_group = "group", 
                              col_outcome = "measure",
                              alpha_threshold = 1)
      circacompare.SD.LD[,2] <- result.i[[2]][,2]
      colnames(circacompare.SD.LD) <- c("","")
      rownames(circacompare.SD.LD) <- result.i[[2]][,1]
      circacompare.SD.LD[,1] <- result.i[[2]][,1]
      
      output$output_statistical_table <- renderDataTable({
        circacompare.SD.LD 
      },escape=FALSE,options =list(pageLength = 15))
      
      output$circacompare<- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            result.i$plot
        })
      
    }else if (input$season== "cicle_comparison" && input$continuo == "LL")
    {ji<-1}else if (input$season== "cicle_comparison" && input$continuo == "DD")
    {ji<-1}  
  }
  shinyjs::showElement(id = 'ready.circ')
  shinyjs::hideElement(id = 'loading.circ')
})
     })
    
# Run the application 
shinyApp(ui = ui, server = server)