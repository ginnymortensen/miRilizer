####################################################################
## Title: app
## Date edited: 08/15/23
## Author: Genevieve Mortensen
## Purpose: web-hosted tool to visualize microRNA-mRNA interactome
####################################################################

# Source library paths ----
# assign(".lib.loc", "/geode2/home/u030/gamorten/Carbonate/R/x86_64-pc-linux-gnu-library/4.1", 
#         envir = environment(.libPaths))
# library(BiocManager)
# options(repos = BiocManager::repositories())
#source correct library path
#source("./.Rprofile.site")

# Import libraries ----
library(shiny)
library(maps)
library(mapproj)
library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(gridExtra)
library(rsconnect)
library(periscope)
library(IRanges)
library(cowplot)
library(grid)

# Set working directory ----
# setwd("/N/project/mirilizer/tool-development/app/")

# Source helper functions ----
source("scripts/merge_n_plot_helper.R")

# Increase size upload limit ----
options(shiny.maxRequestSize=100*1024^2)

# Define UI ----
ui <- fluidPage(
  #This is the title of the web browser for the app
  title = "miRilizer",
  #This is the title shown on the upper left.
  # the title is a png made in Adobe illustrator - it has to be in the www folder.
  titlePanel(img(src = "mirilizer-logo.png", height = 120, width = 400)),
  hr(),
  #this defines the sidebar
  # sidebarLayout(position = "right",
  #   #These are the actual elements of the panel in the sidebar
  #   # functions for text style are just HTML
  #   sidebarPanel(h3("Visualizer Options"),
  #                helpText("Choose elements of each variable to display 
  #                         in the visualization"),
                 #These are the checkbox input specifying each thing the
                 # which could be included or excluded from the plot.
                 # fluidRow(
                 #   column(4,
                 #        checkboxGroupInput("miRrun", 
                 #                    label = "miRNA read abundance", 
                 #                    choices = list("Run 1" = 1, 
                 #                                   "Run 2" = 2, 
                 #                                   "Run 3" = 3,
                 #                                   "Run 4" = 4,
                 #                                   "Run 5" = 5,
                 #                                   "Run 6" = 6))),
                 # column(4,
                 #        checkboxGroupInput("mRNArun", 
                 #                    label = "mRNA read abundance", 
                 #                    choices = list("Run 1" = 1, 
                 #                                   "Run 2" = 2, 
                 #                                   "Run 3" = 3,
                 #                                   "Run 4" = 4,
                 #                                   "Run 5" = 5,
                 #                                   "Run 6" = 6)))),
                 #Starts the next row of inputs
                 # fluidRow(
                 #   column(4,
                 #        checkboxGroupInput("seeds", 
                 #                    label = "seeds", 
                 #                    choices = list("5/6 binding" = 1, 
                 #                                   "6/6 binding" = 2))),
                 # column(4,
                 #        checkboxGroupInput("thermos", 
                 #                    label = "miRNA thermodynamics", 
                 #                    choices = list("binding energy" = 1, 
                 #                                   "opening energy" = 2,
                 #                                   "duplex energies" = 3)))),
                 # ),
    #This is the main panel where images are shown.
    # This also includes the file upload button.
  tabsetPanel(type = "tabs",
              tabPanel("Home",
                       hr(),
    h5("miRilizer is a simple tool allowing users to upload miRNA and/or RNA sequencing data and visualize the gene and miR interactome. Please cite the tool if it is used for your research. If miRilizer does not contain data needed to fill out required visualization arguments, please use the standalone pipeline to generate your own data. Additionally, you can submit a request to add processed files to the data repository on GitHub so that the project can grow increasingly inclusive for all researchers. Science is better when we work together!"),
    hr(),
    h4(strong("Upload files")),
    h5("A caveat of using this tool is that all files must be uploaded in the same run order, however, not all files need to be present for visualization. The only required file(s) to be uploaded are the bed file(s)."),
    hr(),
    fluidRow(column(3,
                    fileInput("file", NULL,
                        label = "Upload bed file(s)",
                        multiple = TRUE)
                       ),
                column(3,
                       fileInput("file1", NULL,
                                 label = "Upload rnaseq bdg file(s)",
                                 multiple = TRUE)
                       ),
             # column(3,
             #        fileInput("file2", NULL,
             #                  label = "Upload seed match output",
             #                  multiple = TRUE)
             #        ),
             # column(3,
             #        fileInput("file3", NULL,
             #                  label = "Upload rnaUP file(s)",
             #                  multiple = TRUE)
             #        ),
                ),
              fluidRow(
                column(3, 
                       tableOutput("files")
                       ),
                column(3,
                       tableOutput("morefiles")
                       ),
                # column(3,
                #        tableOutput("seedfiles")
                #        ),
                # column(3,
                #        tableOutput("rnaupfiles")
                #        )
                ),
    # fluidRow(column(3,
    #                 fileInput("file4", NULL,
    #                           label = "Upload rnaduplex file(s)",
    #                           multiple = TRUE)
    # ),
    # column(3,
    #        fileInput("file5", NULL,
    #                  label = "Upload genomic file(s)",
    #                  multiple = TRUE)
    # ),
    # column(3,
    #        fileInput("file6", NULL,
    #                  label = "Upload canonical gene output",
    #                  multiple = TRUE)
    # ),
    # ),
    # fluidRow(
    #   column(3, 
    #          tableOutput("4file")
    #   ),
    #   column(3,
    #          tableOutput("5file")
    #   ),
    #   column(3,
    #          tableOutput("6file")
    #   ),
    # ),
              hr(),
    ),
    tabPanel("Visualization",
             hr(),
              h5(strong("View the counts of molecules present in your data.")),
              h5("Select a gene by clicking it and then view the counts of the miRs bound to that gene in the adjacent table
                 which will load after making your selection. You can only select one gene at a time and view the miRs for that one gene.
                 When viewing the bound miRs, click on any number of miRs appearing in the table and then confirm the selection by clicking
                 the 'Make Selection' button at the top right of the miR counts table."),
              h5(em("You can download the full dataset or the currently shown page of data.")),
              fluidRow(
                column(6,
                  DT::dataTableOutput("genes")
                ),
                column(6,
                  DT::dataTableOutput("mirs")
                )
              ),
              hr(),
              h5("This text output confirms the selected gene and miRs which will be loaded into the dropdown selections below."),
              verbatimTextOutput("mirplot"),
              hr(),
              h5(strong("Visualize each gene/miR interactome")),
              h5("Use the dropdowns below to select the gene and miRs you'd like visualized.
                 Note that only one gene and microRNA can be visualized at a time."),
              h5(em("You can download the plot using the 'Save plot' button at the bottom of the screen.")),
              fluidRow(
                column(4,
                       selectInput("geneselect", 
                                   "Choose a gene", 
                                   choices = NULL)
                       ),
                column(4,
                       selectInput("mirselect", 
                                   "Choose a miR", 
                                   choices = NULL)
                       )
              ),
                  plotOutput("bigplot", 
                             width = "130%", 
                             height = "500px",
                             click = "plot_click",  # Equiv, to click=clickOpts(id="plot_click")
                             hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                             brush = brushOpts(id = "plot_brush")
                  ),
              br(),
                       plotOutput("legendplot",
                                  width = "auto",
                                  height = "100px"
                                  ),
              hr(),
              # plotOutput("tblplot"),
              downloadButton("saveplot", "Save plot"),
              hr(),
              br()
    ),
    tabPanel("Help",
             hr(),
             h4("Help is here"),
             h5("miRilizer is an R Shiny application which allows for upload of miR eclip bed files and RNA sequencing files. 
             Additionally, data is stored and incorporated into the application. 
             For the purposes of this demonstration, only CYP3A4 and hsa-miR-122-5p is visualized with an uploaded file. Once data storage support is established, 
                the individual file upload options will become available and the tool will be available for public use."),
             hr(),
             h5("miRilizer requires a number of outputs to generate nucleotide-level resolution of interactome visualization. gene-specific RNA sequencing files provide
                the gray trace above the x-axis (representing the transcript of the selected gene with nucleotide positions) in the relevant plot, and the colored lune traces are generated from the uploaded miR eclip assay bed files.
                The blue trace below the x axis shows rnaUP free energies, while the blue dots show rnaduplex binding free energies. Seed matches
                are shown with purple or green dots, and the gray trace below the x-axis shows the u19s output generated from the rnaUP tool."))
    ),
)

# Define server logic ----
server <- function(input, output, session) {
  #This is the server logic that renders the output shown in the UI
  ############################
  ## Create file upload table
  ############################
  output$files <- renderTable({
    #ensure a bed file was uploaded
    bedpath <-input$file
    ext <- tools::file_ext(bedpath$datapath)
    req(bedpath)
    validate(need(ext == "bed", "Please upload a bed file"))
    #create a data frame with the information from the files to show they were uploaded
    striped = TRUE
    hover = TRUE
    bordered = TRUE
    df1 <- data.frame(input$file)
    df2 <- data.frame(input$file2)
    df3 <- rbind(df1, df2)
    df3 <- df3[,1:3]
    df3
  })
  
  output$morefiles <- renderTable({
    #ensure a bed file was uploaded
    rnapath <-input$file1
    ext <- tools::file_ext(rnapath$datapath)
    req(rnapath)
    validate(need(ext == "bdg", "Please upload a bdg file"))
    #create a data frame with the information from the files to show they were uploaded
    striped = TRUE
    hover = TRUE
    bordered = TRUE
    df1 <- data.frame(input$file1)
    df2 <- data.frame(input$file2)
    df3 <- rbind(df1, df2)
    df3 <- df3[,1:3]
    df3
  })
  
  ###############################
  ## Create gene table tab output
  ###############################
  output$genes<- DT::renderDataTable(server = FALSE, {
    #require the file to prevent an error
      bedpath <-input$file
      ext <- tools::file_ext(bedpath$datapath)
      req(bedpath)
      validate(need(ext == "bed", "Please upload a bed file"))
      
    #Access function
      a <- getCountInfo(bedfiles = bedpath$datapath)
      
    #Create the table with download options
      DT::datatable(
        a$genefreq,
        class = c("compact", "cell-border stripe", "hover"),
        selection = "single",
        rownames = FALSE,
        extensions = c('Buttons'),
        options = list(
          dom = 'Bfrltip',
          pageLength = 5,
          lengthMenu = c(5,10,30,50),
          buttons = list(
            list(extend = "csv", 
                 filename = "all_gene_counts",
                 text = "Download all (.csv)",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            ),
            list(extend = "csv", 
                 filename = "selected_gene_counts",
                 text = "Download current (.csv)",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            )
        )
      )
    )
  })
  
  ##############################
  ## Create miR table tab output
  ##############################
  # server = FALSE argument necessary to download full dataset
  output$mirs<- DT::renderDataTable(server = FALSE, {
    #require the file to prevent an error
    bedpath <-input$file
    ext <- tools::file_ext(bedpath$datapath)
    req(bedpath)
    validate(need(ext == "bed", "Please upload a bed file"))
    
    #Establish function access
    a <- getCountInfo(bedfiles = bedpath$datapath)
    #Filter based on selection
    info1 <- input$genes_rows_selected
    mydata1 <- as.data.frame(a$genefreq)
    genestext <- mydata1$Gene
    genesselected <- c(as.character(genestext[info1]))
    #Access specific info based on gene selection
    b <- getSelectInfo(gene = genesselected, 
                       bedmaster = a$bedmaster, 
                       allgenes = a$allgenes, 
                       allmirs = a$allmirs)
    
    
    callback <- c(
      "var dt = table.table().node();",
      "$(dt).selectable({",
      "  distance : 10,",
      "  selecting: function(evt, ui){",
      "    $(this).find('tbody tr').each(function(i){",
      "      if($(this).hasClass('ui-selecting')){",
      "        table.row(':eq(' + i + ')').select();",
      "      }",
      "    });",
      "  }",
      "}).on('dblclick', function(){table.rows().deselect();});"
    )
    
    #Create the table with download options
    dtable <- DT::datatable(
      b$tblg,
      class = c("compact", "cell-border stripe", "hover"),
      callback = DT::JS(callback),
      selection = "multiple",
      rownames = FALSE,
      extensions = c('Buttons', 'Select'),
      options = list(
        dom = 'Bfrltip',
        pageLength = 5,
        lengthMenu = c(5,10,30,50),
        buttons = list(
          list(extend = "csv", 
               filename = "all_miR_counts", 
               text = "Download all (.csv)",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ),
          list(extend = "csv", 
               filename = "selected_miR_counts", 
               text = "Download current (.csv)",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(
            extend = 'collection',
            text = 'Make selection',
            action = DT::JS(
              "function() {
              var node = this[0].node;
              var value = $(node).attr('data-value') || 0;
              value ++;
              $(node).attr('data-value', value);
              Shiny.setInputValue('selectbtn', value, {priority: 'event'});
              }"
            )
          )
        )
      )
    )
    dep <- htmltools::htmlDependency("jqueryui", "1.12.1",
                                     "www/shared/jqueryui",
                                     script = "jquery-ui.min.js",
                                     package = "shiny")
    dtable$dependencies <- c(dtable$dependencies, list(dep))
    dtable
    
  })
  
  ##################################
  ## Show selected genes and miRs ##
  ##################################
  # make the button work only if a gene or miR has been selected
    infotoplot <- eventReactive(input$selectbtn, {
        #Some of this logic is repeated from above - this can be improved by wrapping the
        # function calls inside of reactive() or observe().
        invisible({
        bedpath <-input$file
        a <- getCountInfo(bedfiles = bedpath$datapath)
        b <- getSelectInfo(gene = input$genes_rows_selected, 
                           bedmaster = a$bedmaster, 
                           allgenes = a$allgenes, 
                           allmirs = a$allmirs)
        
        #Collect gene data
        info1 <- input$genes_rows_selected
        mydata1 <- as.data.frame(a$genefreq)
        genestext <- mydata1$Gene
        
        #Collect miR data
        info2 <- input$mirs_rows_selected
        mydata2 <- as.data.frame(b$tblg)
        mirstext <- mydata2$miR
        
        #Arrange the data to be viewable
        genesselected <- c(as.character(genestext[info1]))
        mirsselected <- c(as.character(mirstext[info2]))
        })
        
        #Build text confirmation
        if (length(info1) | length(info2)) {
          cat('This gene was selected:\n')
          cat(genesselected, sep = ', ')
          cat('\n\nThese miRs were selected:\n')
          cat(mirsselected, sep = ', ')
          cat('\n')
          invisible(updateSelectInput(session = session, inputId = "geneselect",
                            label = "Gene",
                            choices = genesselected))
          invisible(updateSelectInput(session = session, inputId = "mirselect",
                            label = "miRs",
                            choices = mirsselected)
                    )
        }
        })

    output$mirplot <- renderPrint({
      infotoplot()
      })
    
## Get Plot
    actualplot <- reactive({
      bedpath <-input$file
      ext <- tools::file_ext(bedpath$datapath)
      req(bedpath)
      validate(need(ext == "bed", "Please upload a bed file"))
      
      #Establish function access
      a <- getCountInfo(bedfiles = bedpath$datapath)
      
      #Get selected mir from drop-down
      req(input$geneselect)
      req(input$mirselect)
      
      #Conditional plotting
      if (length(input$file1)) {
        rnapath <- input$file1
        ext2 <- tools::file_ext(rnapath$datapath)
        req(rnapath)
        validate(need(ext2 == "bdg", "Please upload an RNA file"))
        
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath,
          rnafiles = rnapath$datapath
        )
        
      } else {
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath
        )
      }
      show(c$mainpltnoleg)
    })
    
    #Get Legend
    actuallegend <- reactive({
      bedpath <-input$file
      ext <- tools::file_ext(bedpath$datapath)
      req(bedpath)
      validate(need(ext == "bed", "Please upload a bed file"))
      
      #Establish function access
      a <- getCountInfo(bedfiles = bedpath$datapath)
      
      #Get selected mir from drop-down
      req(input$geneselect)
      req(input$mirselect)
      
      #Conditional plotting
      if (length(input$file1)) {
        rnapath <- input$file1
        ext2 <- tools::file_ext(rnapath$datapath)
        req(rnapath)
        validate(need(ext2 == "bdg", "Please upload an RNA file"))
        
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath,
          rnafiles = rnapath$datapath
        )
        
      } else {
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath
        )
      }
      show(grid.draw(c$mainpltleg))
    })
    
    # Plot the plot
    output$bigplot <- renderPlot(
      actualplot(),
      res = 130
      )
    
    # Plot the legend
    output$legendplot <- renderPlot(
      actuallegend(),
      res = 90
    )
    
    #Download plot
    actualplotr <- function(){
      bedpath <-input$file
      ext <- tools::file_ext(bedpath$datapath)
      req(bedpath)
      validate(need(ext == "bed", "Please upload a bed file"))
      
      #Establish function access
      a <- getCountInfo(bedfiles = bedpath$datapath)
      
      #Get selected mir from drop-down
      req(input$geneselect)
      req(input$mirselect)
      
      #Conditional plotting
      if (length(input$file1)) {
        rnapath <- input$file1
        ext2 <- tools::file_ext(rnapath$datapath)
        req(rnapath)
        validate(need(ext2 == "bdg", "Please upload an RNA file"))
        
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath,
          rnafiles = rnapath$datapath
        )
        
      } else {
        c <- plotter(
          gene_symbol = input$geneselect,
          mir = input$mirselect, 
          bedmaster = a$bedmaster, 
          allgenes = a$allgenes, 
          allmirs = a$allmirs,
          bedfilelist = a$bedfilelist,
          bedfiles = bedpath$datapath
        )
      }
      show(c$mainplot)
    }
    
    filenamer = reactive(paste0(input$geneselect, "_", input$mirselect, ".png"))
    output$saveplot <- downloadHandler(
      filename = filenamer, # variable with filename
      content = function(file) {
        png(file)
        print(actualplotr())
        dev.off()
      })
}


# Run the app ----
shinyApp(ui = ui, server = server)