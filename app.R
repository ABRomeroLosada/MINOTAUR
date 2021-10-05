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
# library(org.Otauri.eg.db) 
##install.packages(pkgs = "./packages/annotation_packages/org.Otauri.eg.db/",repos = NULL,type="source")

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
                     tags$h1(tags$b("ALGAEFUN with MARACAS"), tags$br()),
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
                                                      "Short days (Winter)" = "SD"
                                            ))),
               
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
               
               conditionalPanel(condition = "input.navigation_bar == 'individual' ",    
                                #Choose the kind of analysis that you are interested in
                                radioButtons(inputId = "continuo",
                                             label="Would you like to combine your favourite photoperiod with 
                                             continuous light or darkness conditions?",
                                             choices=c("Continuous light" = "LL",
                                                       "Continuous darkness" = "DD",
                                                       "My favourite photoperiod is enough" = "3days"
                                             ))),
               #Transcriptome or proteome?
               conditionalPanel(condition = "input.navigation_bar == 'individual' ",    
                                #Choose the kind of analysis that you are interested in
                                radioButtons(inputId = "omics",
                                             label="Are you interested in proteomics or transcriptomics?",
                                             choices=c("Transcriptomics rules!" = "rna",
                                                       "Proteomics nerd" = "prot",
                                                       "Multi-omics integration for the win!" = "integration"
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
              
              conditionalPanel(condition = "(input.navigation_bar == 'clusters')",
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
                                                                # tabPanel(tags$b("GO Emap"),
                                                                #          tags$br(),
                                                                #          htmlOutput(outputId = "emapplot_text"),
                                                                #          tags$br(),
                                                                #          div(style= "text-align: center;",
                                                                #              plotOutput(outputId = "emap.plot",inline=TRUE)),
                                                                #          tags$br(),
                                                                #          tags$br(), tags$br()),
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
        file.name <- paste(c("cluster_","peak_",input$zt, ".txt"), collapse="")
        path <- paste(c("clusters_",input$season), collapse="")
        complete_path<- paste(c(path,filename),collapse="/")
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
               pathway.enrichment <- enrichKEGG(gene = target.genes$V1, organism = organism.id, keyType = "kegg",
                                                 universe = microalgae.genes,qvalueCutoff = input$pvalue)
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
                if (input$microalgae == "otauri")
                {
                    kegg.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = pathway.enrichment.result$geneID))
                } else if (input$microalgae == "creinhardtii")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(cre.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if (input$microalgae == "vcarteri")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(vocar.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if (input$microalgae == "ptricornutum")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(phatri.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if (input$microalgae == "ngaditana")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(naga.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if (input$microalgae == "knitens")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
                    }
                } else if (input$microalgae == "hlacustris")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
                    }
                } else if (input$microalgae == "czofingiensis" | 
                           input$microalgae == "mpusilla" |
                           input$microalgae == "bprasinos" |
                           input$microalgae == "mendlicherianum" | 
                           input$microalgae == "smuscicola" |
                           input$microalgae == "dsalina" |
                           input$microalgae == "csubellipsoidea")
                {
                    kegg.enriched.genes <- pathway.enrichment.result$geneID
                    for(i in 1:length(kegg.enriched.genes))
                    {
                        kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
                    }
                }
                
                pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                                    pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                                    pathways.enrichment, 
                                                    kegg.enriched.genes,
                                                    stringsAsFactors = FALSE)
                
                colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                                     "Enrichment (Target Ratio; BG Ration)","Genes")
                
                kegg.result.table.with.links <- pathways.result.table
                ## Add links to genes
                # if(input$microalgae == "otauri")
                # {
                #   gene.link.function <- ostta.gene.link
                # } else if(input$microalgae == "creinhardtii" | input$microalgae == "vcarteri" | 
                #               input$microalgae == "csubellipsoidea" | input$microalgae == "dsalina" |
                #               input$microalgae == "mpusilla")
                # {
                #   gene.link.function <- phytozome.gene.link
                # } else if(input$microalgae == "bprasinos")
                # {
                #   gene.link.function <- bathy.gene.link
                # } else if(input$microalgae == "ptricornutum")
                # {
                #   gene.link.function <- phaeodactylum.gene.link
                # } else if(input$microalgae == "ngaditana")
                # {
                #   gene.link.function <- ngaditana.gene.link
                # } else if(input$microalgae == "knitens")
                # {
                #   gene.link.function <- knitens.gene.link
                # }else if(input$microalgae == "hlacustris")
                # {
                #   gene.link.function <- ncbi.gene.link
                # }else if(input$microalgae == "czofi")
                # {
                #   gene.link.function <- zofi.gene.link
                # }
                
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
                if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
                   input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                   input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
                   input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea")
                {
                    genes.pathway <- rep(0, length(ko.universe))
                    names(genes.pathway) <- ko.universe
                    
                    genes.pathway[target.ko] <- 1
                } else 
                {
                    genes.pathway <- rep(0, length(gene.universe))
                    names(genes.pathway) <- gene.universe
                    
                    genes.pathway[target.genes] <- 1
                }
                
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
            
            if( input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
                input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
                input$microalgae == "dsalina"| input$microalgae == "csubellipsoidea")
            {
                modules.enrichment <- enrichMKEGG(gene = target.ko, universe = ko.universe, organism = "ko", keyType = "kegg",minGSSize = 4)
            } else
            {
                modules.enrichment <- enrichMKEGG(gene = target.genes, universe = gene.universe, organism = organism.id, keyType = "kegg",minGSSize = 4)
            }
            
            modules.enrichment.result <- as.data.frame(modules.enrichment)
            if(nrow(modules.enrichment.result) > 0)
            {
                modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                          bg.ratios = modules.enrichment.result$BgRatio)
                
                if (input$microalgae == "otauri")
                {
                    modules.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = modules.enrichment.result$geneID))
                } else if (input$microalgae == "bprasinos")
                {
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
                    }
                }
                else if (input$microalgae == "creinhardtii")
                {
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        modules.enriched.genes[i] <- paste(cre.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if(input$microalgae == "vcarteri")
                {
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        modules.enriched.genes[i] <- paste(vocar.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if(input$microalgae == "ptricornutum")
                {
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        modules.enriched.genes[i] <- paste(phatri.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if(input$microalgae == "ngaditana")
                {
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        modules.enriched.genes[i] <- paste(naga.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
                    }
                } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
                          input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                          input$microalgae == "smuscicola"| input$microalgae == "mendlicherianum" |
                          input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea")
                {
                    microalga.ko <- AnnotationDbi::select(org.db,columns = c("KO"),keys=keys(org.db,keytype = "GID"))
                    modules.enriched.genes <- modules.enrichment.result$geneID
                    for(i in 1:length(modules.enriched.genes))
                    {
                        current.Ks <- strsplit(modules.enriched.genes[i],split="/")[[1]]
                        current.genes <- c()
                        for(j in 1:length(current.Ks))
                        {
                            current.genes <- c(current.genes,subset(microalga.ko, KO == current.Ks[j])$GID)
                        }
                        modules.enriched.genes[i] <- paste(current.genes,collapse=" ")
                    }
                }
                
                modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                                   modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                                   modules.enrichment, 
                                                   modules.enriched.genes,
                                                   stringsAsFactors = FALSE)
                
                colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                                    "Enrichment (Target Ratio; BG Ration)","Genes")
                
                modules.result.table.with.links <- modules.result.table
                
                ## Add links to genes
                # if(input$microalgae == "otauri")
                # {
                #   gene.link.function <- ostta.gene.link
                # } else if(input$microalgae == "creinhardtii" | input$microalgae == "vcarteri"  | 
                #               input$microalgae == "cocsu" | input$microalgae == "dsalina" |
                #               input$microalgae == "mpusilla")
                # {
                #   gene.link.function <- phytozome.gene.link
                # } else if(input$microalgae == "ptricornutum")
                # {
                #   gene.link.function <- phaeodactylum.gene.link
                # } else if(input$microalgae == "ngaditana")
                # {
                #   gene.link.function <- ngaditana.gene.link
                # } else if(input$microalgae == "knitens")
                # {
                #   gene.link.function <- knitens.gene.link
                # }else if(input$microalgae == "hlacustris")
                # {
                #   gene.link.function <- ncbi.gene.link
                # }else if(input$microalgae == "zofi")
                # {
                #   gene.link.function <- zofi.gene.link
                # }
                
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
            
            if(input$microalgae == "creinhardtii")
            {
                organism.id <- "cre"
            } else if(input$microalgae == "otauri")
            {
                organism.id <- "ota"
            } else if(input$microalgae == "vcarteri")
            {
                organism.id <- "vcn"
            } else if(input$microalgae == "ptricornutum")
            {
                organism.id <- "pti"
            } else if(input$microalgae == "ngaditana")
            {
                organism.id <- "ngd"
            } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris" |
                      input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                      input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
                      input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea"
            )
            {
                organism.id <- "ko"
            }
            
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
    
    ## Actions to perform after click on the genomic button
    observeEvent(input$genomic_button,{
        
        ## Remove previous results
        output$textTableAnnotatedGenes <- renderText(expr = "")
        output$output_gene_chip_table <- renderDataTable(expr = NULL)
        output$download_gene_chip_table <- renderUI(expr = NULL)
        output$piechart.text <- renderText(expr = "")
        output$annotation.pie.chart <- renderPlot(expr = NULL)
        output$tss.distance.text <- renderText(expr = "")
        output$distance.to.tss <- renderPlot(expr = NULL)
        output$gene.signal.text <- renderText(expr = "")
        output$annotated_genes <- renderUI(expr = NULL)
        output$tss.signal.text <- renderText(expr = "")
        output$tss_signal <- renderPlot(expr = NULL)
        shinyjs::hideElement(id = 'loading.tss.signal')
        shinyjs::hideElement(id = 'ready.tss.signal')
        shinyjs::hideElement(id = 'loading.chip')
        shinyjs::hideElement(id = 'ready.chip')
        
        shinyjs::showElement(id = 'loading.chip')
        shinyjs::hideElement(id = 'ready.chip')
        
        ## Load libraries
        library(ChIPseeker)
        library(ChIPpeakAnno)
        library(rtracklayer)
        library(Biostrings)
        library(seqinr)
        
        ## Select txdb 
        if(input$microalgae == "otauri")
        {
            gene.link.function <- ostta.gene.link
            library("org.Otauri.eg.db")
            library("TxDb.Otauri.JGI")
            txdb <- TxDb.Otauri.JGI
            org.db <- org.Otauri.eg.db
            microalgae.annotation <- read.table(file = "annotations/otauri_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/otauri_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "mpusilla")
        {
            gene.link.function <- mpusilla.gene.link
            library("org.MpusillaCCMP1545.eg.db")
            library("TxDb.MpusillaCCMP1545.Phytozome")
            org.db <- org.MpusillaCCMP1545.eg.db
            txdb <- TxDb.MpusillaCCMP1545.Phytozome
            microalgae.annotation <- read.table(file = "annotations/mpusilla_gene_annotation.tsv",sep="\t",header = T,as.is=T,comment.char = "", quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/mpusilla_gene_annotation_links.tsv",sep="\t",header = T,as.is=T,comment.char = "", quote="\"")
        } else if (input$microalgae == "bprasinos")
        {
            gene.link.function <- bathy.gene.link
            library("org.Bprasinos.eg.db")
            library("TxDb.Bprasinos.ORCAE")
            org.db <- org.Bprasinos.eg.db
            txdb <- TxDb.Bprasinos.ORCAE
            microalgae.annotation <- read.table(file = "annotations/bprasinos_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/bprasinos_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "csubellipsoidea")
        {
            gene.link.function <- csubellipsoidea.gene.link
            library("org.Csubellipsoidea.eg.db")
            library("TxDb.Csubellipsoidea.Phytozome")
            org.db <- org.Csubellipsoidea.eg.db
            txdb <- TxDb.Csubellipsoidea.Phytozome
            microalgae.annotation <- read.table(file = "annotations/csubellipsoidea_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/csubellipsoidea_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "creinhardtii")
        {
            gene.link.function <- chlamy.gene.link
            library("org.Creinhardtii.eg.db")
            library("TxDb.Creinhardtii.Phytozome")
            txdb <- TxDb.Creinhardtii.Phytozome
            org.db <- org.Creinhardtii.eg.db
            microalgae.annotation <- read.table(file = "annotations/creinhardtii_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/creinhardtii_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "vcarteri")
        {
            gene.link.function <- vcarteri.gene.link
            library("org.Vcarteri.eg.db")
            library("TxDb.Vcarteri.Phytozome")
            txdb <- TxDb.Vcarteri.Phytozome
            org.db <- org.Vcarteri.eg.db
            microalgae.annotation <- read.table(file = "annotations/vcarteri_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/vcarteri_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "dsalina")
        {
            gene.link.function <- dsalina.gene.link
            library("org.Dsalina.eg.db")
            library("TxDb.Dsalina.Phytozome")
            txdb <- TxDb.Dsalina.Phytozome
            org.db <- org.Dsalina.eg.db
            microalgae.annotation <- read.table(file = "annotations/dsalina_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/dsalina_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "hlacustris")
        {
            gene.link.function <- ncbi.gene.link
            library("org.Hlacustris.eg.db")
            library("TxDb.Hlacustris.NCBI")
            org.db <- org.Hlacustris.eg.db
            txdb <- TxDb.Hlacustris.NCBI
            microalgae.annotation <- read.table(file = "annotations/hlacustris_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/hlacustris_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "czofingiensis")
        {
            gene.link.function <- czofingiensis.gene.link
            library("org.Czofingiensis.eg.db")
            library("TxDb.Czofingiensis.Phytozome")
            org.db <- org.Czofingiensis.eg.db
            txdb <- TxDb.Czofingiensis.Phytozome
            microalgae.annotation <- read.table(file = "annotations/czofingiensis_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/czofingiensis_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "knitens")
        {
            gene.link.function <- no.gene.link
            library("org.Knitens.eg.db")
            library("TxDb.Knitens.Phycocosm")
            org.db <- org.Knitens.eg.db
            txdb <- TxDb.Knitens.Phycocosm
            microalgae.annotation <- read.table(file = "annotations/knitens_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/knitens_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "mendlicherianum")
        {
            gene.link.function <- no.gene.link
            library("org.Mendlicherianum.eg.db")
            library("TxDb.Mendlicherianum.pub")
            org.db <- org.Mendlicherianum.eg.db
            txdb <- TxDb.Mendlicherianum.pub
            microalgae.annotation <- read.table(file = "annotations/mendlicherianum_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/mendlicherianum_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "smuscicola")
        {
            gene.link.function <- no.gene.link
            library("org.Smuscicola.eg.db")
            library("TxDb.Smuscicola.pub")
            org.db <- org.Smuscicola.eg.db
            txdb <- TxDb.Smuscicola.pub
            microalgae.annotation <- read.table(file = "annotations/smuscicola_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/smuscicola_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "ptricornutum")
        {
            gene.link.function <- phaeodactylum.gene.link
            library("org.Ptricornutum.eg.db")
            library("TxDb.Ptricornutum.Ensembl.Protists")
            txdb <- TxDb.Ptricornutum.Ensembl.Protists
            org.db <- org.Ptricornutum.eg.db
            microalgae.annotation <- read.table(file = "annotations/ptricornutum_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/ptricornutum_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } else if (input$microalgae == "ngaditana")
        {
            gene.link.function <- ngaditana.gene.link
            library("org.Ngaditana.eg.db")
            library("TxDb.Ngaditana.JGI")
            txdb <- TxDb.Ngaditana.JGI
            org.db <- org.Ngaditana.eg.db
            microalgae.annotation <- read.table(file = "annotations/ngaditana_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
            microalgae.annotation.links <- read.table(file = "annotations/ngaditana_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
        } 
        
        
        ## Load reference genome for the chosen microalgae
        microalgae.genome.data <- read.fasta(file = paste(c("genomes/",input$microalgae,".fa"),collapse=""),
                                             seqtype = "DNA")
        microalgae.genome <- getSequence(microalgae.genome.data)
        names(microalgae.genome) <- getName(microalgae.genome.data)
        
        ## Extract genomic regions from text box or uploaded file
        if(is.null(input$genomic_regions_file))
        {
            genomic.regions <- as.vector(unlist(strsplit(input$genomic_regions, split="\n",
                                                         fixed = TRUE)[1]))
            
            chrs <- vector(mode = "character", length=length(genomic.regions))
            start.points <- vector(mode = "character", length=length(genomic.regions))
            end.points <- vector(mode = "character", length=length(genomic.regions))
            
            for(i in 1:length(genomic.regions))
            {
                current.splitted.row <- strsplit(genomic.regions[i],split="\\s+")[[1]]
                current.splitted.row <- current.splitted.row[current.splitted.row != ""]
                chrs[i] <- current.splitted.row[1]
                start.points[i] <- current.splitted.row[2]
                end.points[i] <- current.splitted.row[3]
            }
            
            genomic.df <- data.frame(chr=chrs,start=start.points,end=end.points)
            genomic.df <- genomic.df[complete.cases(genomic.df),]
            print(head(genomic.df))
            genomic.regions <- makeGRangesFromDataFrame(df = genomic.df, 
                                                        seqnames.field = "chr",
                                                        start.field = "start",
                                                        end.field = "end")
        } else
        {
            genomic.regions <- readPeakFile(peakfile = input$genomic_regions_file$datapath,header=FALSE)
        }
        
        ## Define promoter region around TSS
        promoter <- getPromoters(TxDb=txdb, 
                                 upstream=input$promoter_length, 
                                 downstream=input$promoter_length)
        
        ## Annotate genomic loci
        peakAnno <- annotatePeak(peak = genomic.regions, 
                                 tssRegion=c(-input$promoter_length, input$promoter_length),
                                 TxDb=txdb)
        
        ## Introductory text for pie chart
        output$piechart.text <- renderText(expr = "<b>The following piechart represents the distribution of the 
    input genomic loci overlapping different gene features. For example, if a section appears corresponding to 
    Promoter (90.75%), this means that 90.75% of the input genomic loci overlap a gene promoter:</b>")
        
        ## Plot pie chart with annotation
        output$annotation.pie.chart <- renderPlot(width = 940, height = 900, res = 120, {
            plotAnnoPie(peakAnno)
        })
        
        ## Introductory text for tss distance distribution
        output$tss.distance.text <- renderText(expr = "<b>The following graph represents the distance distribution 
    either upstream or downstream from genes TSS of the input genomic loci:</b>")
        
        ## Plot distance to tss
        output$distance.to.tss <- renderPlot(width = 940, height = 450, res = 120, {
            plotDistToTSS(peakAnno,
                          title="Distribution of genomic loci relative to TSS",
                          ylab = "Genomic Loci (%) (5' -> 3')")
        })
        
        ## Extract genes 
        peak.annotation <- as.data.frame(peakAnno)
        
        simple.annotation <- sapply(X = as.vector(peak.annotation$annotation),FUN = extract.annotation)
        names(simple.annotation) <- NULL
        
        genes.promoter <- peak.annotation$geneId[simple.annotation == "Promoter"]
        genes.5utr <- peak.annotation$geneId[simple.annotation == "5' UTR"]
        genes.3utr <- peak.annotation$geneId[simple.annotation == "3' UTR"]
        genes.exon <- peak.annotation$geneId[simple.annotation == "Exon"]
        genes.intron <- peak.annotation$geneId[simple.annotation == "Intron"]
        
        ## Select final gene set
        genes <- c()
        if( "Promoter" %in% input$selected_genomic_features )
        {
            genes <- c(genes,genes.promoter)
        }
        if( "5' UTR" %in% input$selected_genomic_features )
        {
            genes <- c(genes,genes.5utr)
        }
        if( "3' UTR" %in% input$selected_genomic_features )
        {
            genes <- c(genes,genes.3utr)
        }
        if( "Exon" %in% input$selected_genomic_features )
        {
            genes <- c(genes,genes.exon)
        }
        if( "Intron" %in% input$selected_genomic_features )
        {
            genes <- c(genes,genes.intron)
        }
        
        genes <- unique(genes)
        
        ## Output table with gene annotation
        genes.annotation.download <- subset(microalgae.annotation, Gene.ID %in% genes) 
        genes.annotation.links <- microalgae.annotation.links[genes,]  
        rownames(genes.annotation.links) <- NULL
        genes.annotation.links <- genes.annotation.links[!is.na(genes.annotation.links$Gene.ID),]
        
        ## Introductory text for target genes 
        annotated.genes.table.text <- "<b>The table below enumerates the potential gene targets associated with the input
    genomic loci. A gene is associated as a target of a genomic locus when it overlaps at least one of the selected
    gene features (i.e. promoter, exon, etc.). Each row represents a gene with its available annotation. Click on the 
    gene id to access information from the corresponding data base. Click on the annotation ids to access a detailed 
    description. You can select the number of genes displayed in each page of the table with the <i>Show</i> 
    dropdown menu below. You can use the <i>Search</i> box to identify a specific gene or annotation. The search
    boxes at the end of the table can be used to search in a specific column. Finally, this table can be downloaded
    by clicking on the button below.</b>"
        
        output$textTableAnnotatedGenes <- renderText(expr = annotated.genes.table.text)
        
        ## Output table with annotated genes
        output$output_gene_chip_table <- renderDataTable({
            genes.annotation.links
        },escape=FALSE,options =list(pageLength = 10)) 
        
        ## Generate UI to download gene table from genomic annotation and creating a downlodable table
        output$download_gene_chip_table<- renderUI(
            tagList(downloadButton(outputId= "downloadGeneTargets", "Download Potential Gene Targets"),tags$br(),tags$br())
        )
        
        ## Download result
        output$downloadGeneTargets <- downloadHandler(
            filename= function() {
                paste("annotated_genes_",microalgae.names[input$microalgae] , ".tsv", sep="")
            },
            content= function(file) {
                write.table(x = genes.annotation.download, quote=F, sep="\t",#genes.annotation.download,quote = F,sep = "\t",
                            file=file,row.names=FALSE,col.names=TRUE)
            })
        
        ## Extraction of the genomic features of the specified genes.
        genes.data <- subset(genes(txdb), gene_id %in% genes)
        
        ## Extract info from genes for profile representations
        genes.data.df <- as.data.frame(genes.data)
        exons.data <- as.data.frame(exons(txdb))
        cds.data <- as.data.frame(cds(txdb))
        
        output$annotated_genes <- renderUI({
            fluidRow(
                column (width = 3,  
                        selectInput(inputId = "selected_annotated_gene", 
                                    label="Choose a gene to inspect:",
                                    multiple = FALSE,selected = genes[1],
                                    choices=genes), 
                        ## Selectize to choose target gene to represent
                        selectizeInput(inputId = "selected.motifs",
                                       label = "Select Motifs",
                                       choices = motif.names,
                                       multiple = TRUE),
                        
                        ## Checkbox to select all available motifs
                        checkboxInput(inputId = "all.motifs",
                                      label = "Select All Motifs:",
                                      value = FALSE),
                        
                        ## Numeric input for PWM score
                        numericInput(inputId = "min.score.pwm", 
                                     label = "Minimum Score for Motif Identification:",
                                     value = 100, 
                                     min = 80,
                                     max = 100,
                                     step = 5),
                        actionButton(inputId = "individual_gene_mark",label = "Go")
                ),
                column(width = 9,
                       plotOutput(outputId = "individual_gene_profile")
                )
            )
        })
        
        selected.annotated.gene.id <- reactive({ 
            if(is.null(input$selected_annotated_gene))
            {
                return()
            } else
            {
                return(input$selected_annotated_gene)
            }
        })
        
        observeEvent(input$individual_gene_mark, { 
            gene.name <- selected.annotated.gene.id()
            
            target.gene.body <- genes.data.df[gene.name,]
            target.gene.chr <- as.character(target.gene.body$seqnames)
            target.gene.start <- target.gene.body$start
            target.gene.end <- target.gene.body$end
            
            target.gene.strand <- as.character(target.gene.body$strand)
            
            ## Extract cds annotation
            cds.data.target.gene <- subset(cds.data, 
                                           seqnames == target.gene.chr & 
                                               (start >= target.gene.start & 
                                                    end <= target.gene.end))
            
            ## Extract exons annotation
            exons.data.target.gene <- subset(exons.data, 
                                             seqnames == target.gene.chr & 
                                                 (start >= target.gene.start & 
                                                      end <= target.gene.end))
            
            ## Determine the genome range to plot including promoter, gene body and 5' UTR
            ## This depends on whether the gene is on the forward or reverse strand
            range.to.plot <- target.gene.body
            
            if(target.gene.strand == "+")
            {
                range.to.plot$start <- range.to.plot$start - input$promoter_length
                range.to.plot$end <- range.to.plot$end + input$promoter_length
            } else if (target.gene.strand == "-")
            {
                range.to.plot$end <- range.to.plot$end + input$promoter_length
                range.to.plot$start <- range.to.plot$start - input$promoter_length
            }
            
            ## Compute the length of the genome range to represent
            current.length <- range.to.plot$end - range.to.plot$start
            
            ## Compute profile in gene
            if(!is.null(input$bw_file))
            {
                selected.bigwig.files <- input$bw_file$data
                #selected.bed.files <- input$genomic_regions_file$data
                
                ## Since ChIPpeakAnno needs more than one region to plot our region
                ## is duplicated 
                regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
                
                ## Import signal from the bigwig files
                cvglists <- sapply(selected.bigwig.files, import, 
                                   format="BigWig", 
                                   which=regions.plot, 
                                   as="RleList")
                
                ## Compute signal in the region to plot
                chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                                    upstream=ceiling(current.length/2), 
                                                    downstream=ceiling(current.length/2),
                                                    n.tile=current.length) 
                
                ## Compute mean signal 
                if(target.gene.strand == "+")
                {
                    chip.signal.mean <- colMeans(chip.signal[[1]],na.rm = TRUE)
                } else if (target.gene.strand == "-")
                {
                    chip.signal.mean <- rev(colMeans(chip.signal[[1]],na.rm = TRUE))
                }
                
                ## Normalization
                chip.signal.mean <- 20 * chip.signal.mean / max(chip.signal.mean)
                
                ## Colors to draw signal
                line.colors <- "blue"
                area.colors <- "lightblue"
            } #else
            # {
            #   chip.signal.mean <- rep(20,length(cord.x)) 
            #   ## Colors to draw signal
            #   line.colors <- "white"
            #   area.colors <- "white"
            # }
            
            ## Determine upper limit of the graph
            upper.lim <- 10#21
            
            ## Height to draw DNA strand
            gene.height <- -25
            cord.x <- 1:current.length
            
            ## Exon width to plot
            exon.width <- 1#2
            
            ## Cds width to plot
            cds.width <- 1.5#3
            
            ## Width of the rectangule representing the peak reagion
            peak.width <- 0.5#1
            
            ## Extract exons for target gene
            exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
            
            ## Transform exon coordinates to current range
            min.pos <- min(exons.data.target.gene$start)
            
            exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter_length
            exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter_length
            
            ## Extract cds for target gene
            cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
            
            cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter_length
            cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter_length
            
            output$individual_gene_profile <- renderPlot({#width = 940, height = 700, res = 120, {
                plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
                     cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
                     ylim=c(-35,upper.lim),xlim=c(-3000,max(cord.x)))
                
                ## Represent exons
                for(i in 1:nrow(exons.data.target.gene))
                {
                    # Determine start/end for each exon
                    current.exon.start <- exons.data.target.gene$start[i]
                    current.exon.end <- exons.data.target.gene$end[i]
                    
                    ## Determine coordinates for each exon polygon and represent it
                    exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
                    exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
                    
                    polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
                }
                
                if(nrow(cds.data.target.gene) > 0)
                {
                    for(i in 1:nrow(cds.data.target.gene))
                    {
                        # Determine current cds start/end
                        current.cds.start <- cds.data.target.gene$start[i]
                        current.cds.end <- cds.data.target.gene$end[i]
                        
                        # Determine curret cds coordinates for the polygon and represent it
                        cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
                        cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
                        polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
                    }
                }
                
                ## Draw arrow to represent transcription direction 
                if(target.gene.strand == "+")
                {
                    lines(c(input$promoter_length,input$promoter_length,input$promoter_length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
                    lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+6,gene.height+5),lwd=3)
                    lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+4,gene.height+5),lwd=3)
                } else if (target.gene.strand == "-")
                {
                    lines(c(current.length - input$promoter_length, current.length - input$promoter_length, current.length - input$promoter_length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
                    lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
                    lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
                }
                
                ## Draw promoter range
                if(target.gene.strand == "+")
                {
                    axis(side = 1,labels = c(- input$promoter_length, - input$promoter_length / 2,"TSS"),at = c(1,input$promoter_length/2,input$promoter_length),lwd=2,cex=1.5,las=2,cex=2)
                } else if(target.gene.strand == "-")
                {
                    axis(side = 1,labels = c("TSS",- input$promoter_length / 2,- input$promoter_length),at = c(current.length-input$promoter_length,current.length-input$promoter_length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
                }
                
                ## Draw gene name
                text(x = current.length / 2, y = -33 , 
                     labels = bquote(italic(.(gene.name))),cex = 1.7,font = 3)
                
                ## Extract bed file name 1 and read it
                current.peaks <- as.data.frame(genomic.regions)#read.table(file=input$genomic_regions_file$data,header = F, as.is = T)
                peak.coordinates <- subset(current.peaks, seqnames == as.character(range.to.plot$seqnames) & start >= range.to.plot$start & end <= range.to.plot$end) 
                current.peaks.to.plot <- peak.coordinates[,2:3]
                
                ## Transform coordinates 
                current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
                
                ## Check if there are peaks for the target gene
                if(nrow(current.peaks.to.plot) > 0)
                {
                    #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
                    for(j in 1:nrow(current.peaks.to.plot))
                    {
                        ## Extract start and end point of each peak region
                        current.peak.start <- current.peaks.to.plot[j,1]
                        current.peak.end <- current.peaks.to.plot[j,2]
                        
                        ## Computer coordinates for polygon and draw it
                        peak.x <- c(current.peak.start,current.peak.end,
                                    current.peak.end,current.peak.start)
                        peak.y <- c(peak.width - 12,   peak.width - 12, 
                                    - peak.width - 12, - peak.width - 12)  
                        
                        polygon(x = peak.x, y = peak.y, col = "bisque", border = "darkred",lwd=2)
                    }
                }
                
                ## Draw profiles if bw file is provided
                if(!is.null(input$bw_file))
                {
                    ## Compute base line for current TF
                    current.base.line <- - 10
                    
                    ## Represent signal from the current TF
                    lines(chip.signal.mean+current.base.line,type="l",col=line.colors,lwd=3)
                    ## Determine polygon coordinates and represent it
                    cord.y <- c(current.base.line,chip.signal.mean+current.base.line,current.base.line)
                    cord.x <- 1:length(cord.y)
                    polygon(cord.x,cord.y,col=area.colors)
                }
                
                ## Determine TFBS motifs to search for
                if(input$all.motifs)
                {
                    selected.motifs.pwm <- motifs.pwm
                } else if(!is.null(input$selected.motifs))
                {
                    selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
                }
                
                if( input$all.motifs | !is.null(input$selected.motifs))
                {
                    selected.motif.names <- names(selected.motifs.pwm)
                    selected.motif.ids <- motif.ids[selected.motif.names]
                    
                    ## Initialize data frame containing TF binding sequences in the peak regions
                    df.hits <- data.frame(0,"","","")
                    colnames(df.hits) <- c("position","id","name","seq")
                    
                    ## Identify TF binding DNA motifs 
                    if(nrow(current.peaks.to.plot) > 0)
                    {
                        #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
                        for(j in 1:nrow(current.peaks.to.plot))
                        {
                            ## Genomic coordinates of the current peak
                            peak.chr <- peak.coordinates[j, 1]
                            peak.start <- peak.coordinates[j, 2]
                            peak.end <- peak.coordinates[j, 3]
                            
                            ## Extract start and end point of each peak region in our plot
                            current.peak.start <- current.peaks.to.plot[j,1]
                            current.peak.end <- current.peaks.to.plot[j,2]
                            
                            ## Extract peak sequence
                            peak.sequence <- c2s(microalgae.genome[[peak.chr]][peak.start:peak.end])
                            peak.rev.comp.sequence <- reverse.complement(peak.sequence)
                            
                            for(k in 1:length(selected.motifs.pwm))
                            {
                                motif.pwm <- selected.motifs.pwm[[k]]
                                
                                hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                                    min.score = paste0(input$min_score_pwm,"%"))
                                hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
                                hits.fw <- as(hits.fw, "IRanges")
                                hits.fw.start <- start(hits.fw)
                                hits.fw.end <- end(hits.fw)
                                
                                if(length(hits.fw.start) > 0)
                                {
                                    df.hits.fw <- data.frame(((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                                             rep(selected.motif.ids[k],length(hits.fw.start)),
                                                             rep(selected.motif.names[k],length(hits.fw.start)),
                                                             hits.fw.seqs)
                                    colnames(df.hits.fw)  <- c("position","id","name","seq")
                                    df.hits <- rbind(df.hits,df.hits.fw)
                                }
                                
                                hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                                                     min.score = paste0(input$min.score.pwm,"%"))
                                hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
                                hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
                                names(hits.rev.seqs) <- NULL
                                
                                hits.rev <- as(hits.rev, "IRanges")
                                hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
                                hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
                                
                                if(length(hits.rev.start) > 0)
                                {
                                    df.hits.rev <- data.frame(((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                                              rep(selected.motif.ids[k],length(hits.rev.start)),
                                                              rep(selected.motif.names[k],length(hits.rev.start)),
                                                              hits.rev.seqs)
                                    colnames(df.hits.rev)  <- c("position","id","name","seq")
                                    df.hits <- rbind(df.hits,df.hits.rev)
                                }
                            }
                        }
                    }
                    
                    ## Remove first line of the data frame added just for technical reason
                    df.hits <- df.hits[-1,]
                    nrow(df.hits)
                    
                    ## Draw TF binding sites
                    detected.tfbs <- unique(as.vector(df.hits$name))
                    
                    # TF binding sites colors and symbol shapes
                    symbol.shapes <- c(17, 18, 19, 15)
                    symbol.color <- c("blue", "red", "darkgreen", "magenta")
                    
                    number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
                    necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
                    necessary.colors <- rep(symbol.color,number.of.shapes)
                    
                    if(length(detected.tfbs) > 0)
                    {
                        for(i in 1:length(detected.tfbs))
                        {
                            current.tfbs <- detected.tfbs[i]
                            current.shape <- necessary.shapes[i]
                            current.color <- necessary.colors[i]
                            
                            positions <- subset(df.hits, name == current.tfbs)
                            
                            for(j in 1:nrow(positions))
                            {
                                pos.to.draw <- positions$position[j]
                                
                                points(x = pos.to.draw, y = -15,
                                       pch = current.shape, col = current.color, cex = 1)
                            }
                        }
                        
                        ## Add legend for TFBS
                        legend.step <- 5
                        for(i in 1:length(detected.tfbs))
                        {
                            points(x = -3000, y = upper.lim - (i-1)*legend.step, 
                                   pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
                            
                            
                            current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
                            current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
                            
                            text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
                                 adj = 0,cex = 0.7)
                        }
                    }
                }
            }) # close output$individual_gene_profile
        }) # observeEvent input$individual_gene_mark
        
        
        individual.gene.signal.text <- "<b> Here you can visualize an individual marked
    gene following these steps:
    <ol>
      <li>Choose or type your gene of interest from the left dropdown menu where
      autocompletion can be used.</li>
      <li>Optionally, choose or type the name of a DNA motif or transcription factor consensus
      binding sequence to be searched for in the genomic loci associated with your
      gene of interest. Alternatively tick the \"Select All Motifs\" box when you
      want to search for DNA motifs in our database. Here autocompletion is also activated.</li>
      <li>Choose a a minimun score for the selected DNA motif identification.</li>
      <li>Click on \"Go\" to visualize the selected genes with its associated input
      genomic loci, the identification of the selected DNA motifs according to the
      specified score and the signal profile over the selected gene when a BigWig
      file is provided.</li>
    </ol>
    The visualization is available on the right. The DNA sequence is represented by a
    black line. The widest blue rectangles represent exons. Untranslated regions at 5' and
    3' end of the transcript are depicted as thinner blue rectangles. Black lines also 
    represent introns. An arrow is used to mark the TSS (Transcriptional Start Site) 
    and the transcription sense. The red boxes located on top to the gene representation
    correspond to the specific input genomic loci associated with the selected gene. When
    a BigWig file is provided representing, for example, the signal of a epigenetic mark
    the signal profile over the selected gene is reprenting in lighblue. 
    </b>
    "
        
        ## Plot signal around tss and tes
        if(!is.null(input$bw_file))
        {
            # shinyjs::showElement(id = 'loading.tss.signal')
            # shinyjs::hideElement(id = 'ready.tss.signal')
            
            ## Extraction of the TSS 
            genes.tss <- resize(genes.data, width=1, fix='start')
            
            ## Centering around TSS with promoter length
            around.genes.tss <- genes.tss
            start(around.genes.tss) <- start(genes.tss) - input$promoter_length
            end(around.genes.tss) <- end(genes.tss) + input$promoter_length
            
            ## Importing bigWig file
            cvglists <- sapply(input$bw_file$data, import,
                               format="BigWig",
                               which=around.genes.tss,
                               as="RleList")
            
            ## Extracting the signal around TSS with promoter length
            number.tiles <- 2*input$promoter_length/20
            tss.sig <- featureAlignedSignal(cvglists, around.genes.tss,
                                            upstream=input$promoter_length,
                                            downstream=input$promoter_length,
                                            n.tile=number.tiles)
            
            ## Extraction of the TES
            genes.tes <- resize(genes.data, width=1, fix='end')
            
            ## Centering around TES 
            around.genes.tes <- genes.tes
            start(around.genes.tes) <- start(genes.tes) - input$promoter_length
            end(around.genes.tes) <- end(genes.tes) + input$promoter_length
            
            print("around TES")
            print(around.genes.tes)
            
            
            ## Extracting the signal 2kb from the center with 50 tile (This may be slow)
            number.tiles.tes <- 2 * input$promoter_length /20
            tes.sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                                            upstream=input$promoter_length, 
                                            downstream=input$promoter_length,
                                            n.tile=number.tiles.tes) 
            
            
            print("TSS signal:")
            print(tss.sig)
            output$tss_signal <- renderPlot({
                
                par(mfrow=c(1,2))
                ## Plotting the results around TSS (you may have to change the names of the conditions)
                profile.around.tss <- colMeans(tss.sig[[1]],na.rm = TRUE)
                profile.around.tes <- colMeans(tes.sig[[1]],na.rm = TRUE)
                max.y <- max(c(profile.around.tss, profile.around.tes))
                plot(profile.around.tss,type="l",col="blue",lwd=3,ylab="",
                     cex.lab=2,axes=FALSE,xlab="",main="Signal around TSS",cex.main=1.5,ylim=c(0,max.y))
                polygon(c(1,1:length(profile.around.tss),length(profile.around.tss)),
                        c(0,profile.around.tss,0),col="lightblue")
                
                axis(side = 1,
                     labels = c(-input$promoter_length,-input$promoter_length/2,
                                "TSS",
                                input$promoter_length/2,input$promoter_length),
                     at = c(1,number.tiles/4,number.tiles/2,3*number.tiles/4,number.tiles),lwd=2,cex=1.5,las=2,cex=2)
                
                plot(profile.around.tes,type="l",col="red4",lwd=3,ylab="",
                     cex.lab=2,axes=FALSE,xlab="",main="Signal around TES",cex.main=1.5,ylim=c(0,max.y))#,ylim=c(0,830))
                polygon(c(1,1:length(profile.around.tes),length(profile.around.tes)),
                        c(0,profile.around.tes,0),col="bisque")
                
                axis(side = 1,
                     labels = c(-input$promoter_length,-input$promoter_length/2,
                                "TES",
                                input$promoter_length/2,input$promoter_length),
                     at = c(1,number.tiles.tes/4,number.tiles.tes/2,3*number.tiles.tes/4,number.tiles.tes),lwd=2,cex=1.5,las=2,cex=2)
                
            })
            
            output$tss.signal.text <- renderText(expr = "<b>The graph below represents
      the average signal around the TSS (Transcription Start Site) and TES (Transcription
      End Site) of marked genes with your input genomic loci.</b>")
            output$gene.signal.text <- renderText(expr = individual.gene.signal.text)
        } else
        {
            output$tss.signal.text <- renderText(expr = "<b> Average signal levels around
      genes TSS (Transcription Start Site) and TES (Trascription End Site) could not 
      be computed since NO BigWig file was selected.</b>")
            output$gene.signal.text <- renderText(expr = paste(individual.gene.signal.text,"<br><b> Average signal levels for 
      specific marked genes cannot be computed since NO BigWig file was selected.</b>"))
        }
        shinyjs::hideElement(id = 'loading.chip')
        shinyjs::showElement(id = 'ready.chip')
    })
})

# Run the application 
shinyApp(ui = ui, server = server)