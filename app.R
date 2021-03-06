library(shiny)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(dplyr)
library(Matrix)
library(survival)
library(survminer)
setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/PADACS/')

## load input choices
sc.dataset.table <- read.csv('Data/sc_comp-table.csv', header = TRUE)
bulk.dataset.table <- read.csv('Data/bulk_comp-table.csv', header = TRUE)
plot.choices <- c('Scatter Plot', 'Expression by Group', 'Differential Expression / Activity')
bulk.plot.choices <- c('Expression by Group', 'Differential Expression / Activity', 'Survival Analysis')

## UI
ui <- fluidPage(
  # webpage title
  titlePanel('PADCAS: Visualization of PDAC Gene Expression and Protein Activity'),
  
  # create tabs
  tabsetPanel(
    
    # about panel
    tabPanel(title = 'About',
             includeMarkdown('Data/about-tab.rmd')
    ),
    
    # bulk analysis panel
    tabPanel('Bulk Analysis', h1("Bulk RNASeq Analysis"),
             sidebarLayout(
               sidebarPanel(
                 # dataset input
                 selectInput('bulkData', 'Dataset', choices = unique(bulk.dataset.table$Dataset)),
                 # data modality
                 selectInput('bulkType', 'Select data modality:', choices = NULL),
                 # phenotype
                 selectInput('bulkPheno', 'Phenotype', choices = NULL),
                 # comparison
                 selectInput('bulkComp', 'Comparison', choices = NULL),
                 # plots to generate
                 #checkboxGroupInput('bulkPlots', "What plots would you like to view?", bulk.plot.choices),
                 # gene to plot
                 textInput('bulkGene', label = 'Choose a gene to display.', value = '', placeholder = NULL),
                 # action button to plot
                 actionButton('bulkPlot', label = 'Plot')
               ),
               mainPanel(
                 plotOutput("bulkBox", height = '3in'),
                 br(),
                 plotOutput("bulkDot", height = '3in'),
                 br(),
                 plotOutput("bulkSurvival", height = '3.5in'),
               )
             )
    ),
    
    # single-cell panel
    tabPanel('SC Analysis', h1("Single Cell RNASeq Analysis"),
             sidebarLayout(
               sidebarPanel(
                 # dataset input
                 selectInput('scData', 'Dataset', choices = unique(sc.dataset.table$Dataset)),
                 # data modality
                 selectInput('scType', 'Data Modality', choices = NULL),
                 # comparison input
                 selectInput('scComp', 'Comparison', choices = NULL),
                 # phenotype input
                 selectInput('scPheno', 'Phenotype', choices = NULL),
                 # plots to generate
                 #checkboxGroupInput('scPlots', "What plots would you like to view?", plot.choices),
                 # gene to plot
                 textInput('scGene', label = 'Choose a gene to display.', value = '', placeholder = NULL),
                 # action button to plot
                 actionButton('scPlot', label = 'Plot')
               ),
               mainPanel(
                 plotOutput("scScatter", height = '3.5in'),
                 br(),
                 plotOutput("scViolin", height = '3in'),
                 br(),
                 plotOutput("scDotPlot", height = '3in'),
               )
             )
    )
    
    # citation panel
    #tabPanel('Citation and Contact',
    #         includeMarkdown('Data/citation-contact-tab.rmd')
    #)
  )
)

## server functions
server <- function(input, output) {
  ## sc: reactive input functions
  ###############
  scTypeInput <- reactive({
    filter(sc.dataset.table, Dataset == input$scData)
  })
  observeEvent(scTypeInput(), {
    choices <- unique(scTypeInput()$Modality)
    updateSelectInput(inputId = "scType", choices = choices)
  })
  
  scPhenoInput <- reactive({
    filter(sc.dataset.table, Dataset == input$scData, 
           Modality == input$scType)
  })
  observeEvent(scPhenoInput(), {
    choices <- unique(scPhenoInput()$Phenotype)
    updateSelectInput(inputId = "scPheno", choices = choices) 
  })
  
  scCompInput <- reactive({
    filter(sc.dataset.table, Dataset == input$scData, 
           Modality == input$scType,
           Phenotype == input$scPheno)
  })
  observeEvent(scCompInput(), {
    choices <- unique(scCompInput()$Comparison)
    updateSelectInput(inputId = "scComp", choices = choices) 
  })
  ###############
  
  ## sc plots
  ###############
  scPlotData <- eventReactive(input$scPlot, {
    plot.type <- 1
    # load the appropriate data
    if (input$scData == 'Peng') {
      dge.file <- 'Data/peng/peng_'
      
      if (input$scPheno == 'Tumor') {
        # prepare expression vector
        gexp.cpm <- readRDS('Data/peng/peng_tumor_cpm.rds')
        validate(need(input$scGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        gexp.vec <- gexp.cpm[input$scGene,]
        rm(gexp.cpm)
        # load metadata
        gexp.meta <- readRDS('Data/peng/peng_tumor_metaData.rds')
        meta.vec <- gexp.meta 
        # load densemap
        gexp.dmap <- readRDS('Data/peng/peng_tumor_densmap.rds')
        # modify dge string
        dge.file <- paste(dge.file, 'tumor_', sep = '')
        
      } else if (input$scPheno == 'Normal') {
        # prepare expression vector
        gexp.cpm <- readRDS('Data/peng/peng_normal_cpm.rds')
        validate(need(input$scGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        gexp.vec <- gexp.cpm[input$scGene,]
        rm(gexp.cpm)
        # load metadata
        gexp.meta <- readRDS('Data/peng/peng_normal_metaData.rds')
        # load densemap
        gexp.dmap <- readRDS('Data/peng/peng_normal_densmap.rds')
        # modify dge string
        dge.file <- paste(dge.file, 'normal_', sep = '')
        
      } else if (input$scPheno == 'Tumor vs. Normal') {
        # prepare expression vector
        gexp.cpm <- readRDS('Data/peng/peng_normal_cpm.rds')
        validate(need(input$scGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        normal.vec <- gexp.cpm[input$scGene,]
        rm(gexp.cpm)
        gexp.cpm <- readRDS('Data/peng/peng_tumor_cpm.rds')
        tumor.vec <- gexp.cpm[input$scGene,]
        rm(gexp.cpm)
        gexp.vec <- c(normal.vec, tumor.vec)
        # load metadata
        gexp.meta <- readRDS('Data/peng/peng_tvn_metaData.rds')
        # load densemap
        gexp.dmap <- readRDS('Data/peng/peng_normal_densmap.rds')
        # modify dge string
        dge.file <- paste(dge.file, 'tvn_', sep = '')
        plot.type <- 2
        
      }
      
    } else if (input$scData == 'Lin') {
      dge.file <- 'Data/lin/lin_'
      # prepare expression vector
      gexp.cpm <- readRDS('Data/lin/lin_cpm.rds')
      validate(need(input$scGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
      gexp.vec <- gexp.cpm[input$scGene,]
      rm(gexp.cpm)
      # load metadata
      gexp.meta <- readRDS('Data/lin/lin_metadata.rds')
      # load densemap
      gexp.dmap <- readRDS('Data/lin/lin_densmap.rds')
    } else if (input$scData == 'Tuveson') {
      if (input$scType == 'Protein Activity') {
        pact.mat <- readRDS('Data/tuveson/tuveson_pact.rds')
        validate(need(input$scGene %in% rownames(pact.mat), "Selected protein is not present in this data."))
        pact.vec <- pact.mat[input$scGene,]
        rm(pact.mat)
        # load metadata 
        pact.meta <- readRDS('Data/tuveson/tuveson_metadata.rds')
        # load densemap
        pact.dmap <- readRDS('Data/tuveson/tuveson_pact_dmap.rds')
        dpa.obj <- readRDS('Data/tuveson/tuveson_Tumor_Subtype-dpa.rds')
      } else {
        dge.file <- 'Data/tuveson/tuveson_'
        # prepare expression vector
        gexp.cpm <- readRDS('Data/tuveson/tuveson_cpm.rds')
        validate(need(input$scGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        gexp.vec <- gexp.cpm[input$scGene,]
        rm(gexp.cpm)
        # load metadata
        gexp.meta <- readRDS('Data/tuveson/tuveson_metadata.rds')
        # load densemap
        gexp.dmap <- readRDS('Data/tuveson/tuveson_densmap.rds')
      }
    }
    
    if (input$scType == 'Protein Activity') {
      scDataList <- list('gexp.vec' = pact.vec,
                         'gexp.meta' = pact.meta,
                         'gexp.dmap' = pact.dmap,
                         'dge' = dpa.obj,
                         'plot.type' = plot.type)
    } else {
      # load dge
      dge.file <- paste(dge.file, gsub(' ', '_', input$scComp), '-dge.rds', sep = '')
      dge.obj <- readRDS(dge.file)
      
      # prepare data list and return
      scDataList <- list('gexp.vec' = gexp.vec,
                         'gexp.meta' = gexp.meta,
                         'gexp.dmap' = gexp.dmap,
                         'dge' = dge.obj,
                         'plot.type' = plot.type)
    }
    return(scDataList)
  })
  
  output$scScatter <- bindEvent(renderPlot({
    scData <- scPlotData()
    # process data
    comp.vec <- scData$gexp.meta[, input$scComp]
    names(comp.vec) <- rownames(scData$gexp.meta)
    use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
    # create plot.df
    plot.df <- data.frame('Expression' = scale(scData$gexp.vec[use.samps]),
                          'Meta' = as.factor(comp.vec[use.samps]),
                          'DMAP1' = scData$gexp.dmap[use.samps,1],
                          'DMAP2' = scData$gexp.dmap[use.samps,2])
    #median.exp <- median(scData$gexp.vec[use.samps][which(scData$gexp.vec[use.samps] > 0)])
    # generate plot
    comp.scatter <- ggplot(plot.df, aes(DMAP1, DMAP2)) +
      geom_point(aes(color = Meta), size = 0.5) + 
      labs(x = 'DMAP1', y = 'DMAP2', title = input$scComp) +
      scale_color_discrete(name = input$scComp)
    exp.scatter <- ggplot(plot.df, aes(DMAP1, DMAP2)) +
      geom_point(aes(color = Expression)) + 
      labs(x = 'DMAP1', y = 'DMAP2')
    if (input$scType == 'Gene Expression') {
      exp.scatter <- exp.scatter + scale_colour_gradient2(low = 'forestgreen', mid = 'darkgrey', high = 'purple', name = 'Expression') +
        ggtitle('Scaled Expression') 
    } else if (input$scType == 'Protein Activity') {
      exp.scatter <- exp.scatter + scale_colour_gradient2(low = 'blue', mid = 'darkgrey', high = 'red', name = 'Activity') +
        ggtitle('Protein Activity') 
    }
    return(ggarrange(plotlist = list(comp.scatter, exp.scatter),
                     nrow = 1, ncol = 2))
  }), scPlotData())
  
  output$scViolin <- bindEvent(renderPlot({
    scData <- scPlotData()
    if (scData$plot.type == 1) {
      # process data
      comp.vec <- scData$gexp.meta[, input$scComp]
      names(comp.vec) <- rownames(scData$gexp.meta)
      use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
      # create plot.df
      plot.df <- data.frame('Expression' = scData$gexp.vec[use.samps],
                            'Meta' = as.factor(comp.vec[use.samps]))
      # generate plot
      violin.plot <- ggplot(plot.df, aes(x = Meta, y = Expression)) +
        geom_violin(aes(fill = Meta), color = 'black') +
        xlab(input$scComp) + ylab('Log2 CPM') + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.title = element_blank())
      
      if (input$scType == 'Protein Activity') {
        violin.plot <- violin.plot + ggtitle(paste(input$scGene, ': Group Activity', sep = ''))
      } else {
        violin.plot <- violin.plot + ggtitle(paste(input$scGene, ': Group Expression', sep = ''))
      }
      
      return(violin.plot)
    } else if (scData$plot.type == 2) {
      # process data
      compartment.vec <- scData$gexp.meta[, 'Compartment']
      names(compartment.vec) <- rownames(scData$gexp.meta)
      comp.vec <- scData$gexp.meta[, input$scComp]
      names(comp.vec) <- rownames(scData$gexp.meta)
      use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
      # create plot.df
      plot.df <- data.frame('Expression' = scData$gexp.vec[use.samps],
                            'Meta' = as.factor(comp.vec[use.samps]),
                            'Compartment' = as.factor(compartment.vec[use.samps]))
      # generate plot
      violin.plot <- ggplot(plot.df, aes(x = Meta, y = Expression)) +
        geom_violin(aes(fill = Compartment), color = 'black') +
        ggtitle(paste(input$scPheno, ': ', input$scGene, ' Expression', sep = '')) + 
        xlab(input$scComp) + ylab('Log2 CPM') + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.title = element_blank())
      
      return(violin.plot)
    }
  }), scPlotData())
  
  output$scDotPlot <- bindEvent(renderPlot({
    scData <- scPlotData()
    # process data
    comp.vec <- scData$gexp.meta[, input$scComp]
    names(comp.vec) <- rownames(scData$gexp.meta)
    use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
    # create plot.df
    plot.df <- data.frame('p.val' = scData$dge$p.val[input$scGene,],
                          'rbsc' = scData$dge$rbsc[input$scGene,],
                          'level' = colnames(scData$dge$p.val),
                          'gene' = rep(input$scGene, length = ncol(scData$dge$p.val)))
    # generate plot
    dot.plot <- ggplot(plot.df, aes(level, gene)) +
      geom_point(aes(color = rbsc, size = p.val)) + 
      scale_size(trans = 'reverse') + 
      scale_color_gradient(low = 'blue', high = 'red') +
      xlab(input$scComp) + ylab(input$gene) + 
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
    return(dot.plot)
  }), scPlotData())
  ###############
  
  ### bulk: reactive input functions
  ###############
  bulkTypeInput <- reactive({
    filter(bulk.dataset.table, Dataset == input$bulkData)
  })
  observeEvent(bulkTypeInput(), {
    choices <- unique(bulkTypeInput()$Modality)
    updateSelectInput(inputId = "bulkType", choices = choices)
  })
  
  bulkPhenoInput <- reactive({
    filter(bulk.dataset.table, Dataset == input$bulkData, 
           Modality == input$bulkType)
  })
  observeEvent(bulkPhenoInput(), {
    choices <- unique(bulkPhenoInput()$Phenotype)
    updateSelectInput(inputId = "bulkPheno", choices = choices) 
  })
  
  bulkCompInput <- reactive({
    filter(bulk.dataset.table, Dataset == input$bulkData, 
           Modality == input$bulkType,
           Phenotype == input$bulkPheno)
  })
  observeEvent(bulkCompInput(), {
    choices <- unique(bulkCompInput()$Comparison)
    updateSelectInput(inputId = "bulkComp", choices = choices) 
  })
  ###############
  
  ### bulk: plots
  ###############
  bulkPlotData <- eventReactive(input$bulkPlot, {
    plot.type <- 1
    if (input$bulkData == "CUMC LCM") {
      if (input$bulkPheno == "Tumor") {
        # prepare expression vector
        gexp.cpm <- readRDS('Data/cumc_lcm/cumc-lcm_epi_tpm.rds')
        validate(need(input$bulkGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        gexp.vec <- gexp.cpm[input$bulkGene,]
        rm(gexp.cpm)
        # load metadata
        gexp.meta <- readRDS('Data/cumc_lcm/cumc-lcm_epi_meta.rds')
        # load dge
        gexp.dge <- readRDS('Data/cumc_lcm/cumc-lcm_epi_dge.rds')
        
      } else if (input$bulkPheno == "Stroma") {
        # prepare expression vector
        gexp.cpm <- readRDS('Data/cumc_lcm/cumc-lcm_stroma_tpm.rds')
        validate(need(input$bulkGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
        gexp.vec <- gexp.cpm[input$bulkGene,]
        rm(gexp.cpm)
        # load metadata
        gexp.meta <- readRDS('Data/cumc_lcm/cumc-lcm_stroma_meta.rds')
        # load dge
        gexp.dge <- readRDS('Data/cumc_lcm/cumc-lcm_stroma_dge.rds')
        
      }
    } else if (input$bulkData == "TCGA") {
      # prepare expression vector
      gexp.cpm <- readRDS('Data/tcga/tcga_tpm.rds')
      validate(need(input$bulkGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
      gexp.vec <- gexp.cpm[input$bulkGene,]
      rm(gexp.cpm)
      # load metadata
      gexp.meta <- readRDS('Data/tcga/tcga_meta.rds')
      # load dge or NaRnEA
      if (input$bulkType == 'Protein Activity') {
        gexp.dge <- readRDS('Data/tcga/tcga_narnea.rds')
        plot.type <- 2
      } else {
        gexp.dge <- readRDS('Data/tcga/tcga_dge.rds')
      }
    } else if (input$bulkData == "UNC") {
      # prepare expression vector
      gexp.cpm <- readRDS('Data/unc/unc_tpm.rds')
      validate(need(input$bulkGene %in% rownames(gexp.cpm), "Selected gene is not present in this data."))
      gexp.vec <- gexp.cpm[input$bulkGene,]
      rm(gexp.cpm)
      # load metadata
      gexp.meta <- readRDS('Data/unc/unc_meta.rds')
      # load dge
      if (input$bulkType == 'Protein Activity') {
        gexp.dge <- readRDS('Data/unc/unc_narnea.rds')
        plot.type <- 2
      } else {
        gexp.dge <- readRDS('Data/unc/unc_dge.rds')
      }
    }
    
    # prepare data list and return
    bulkDataList <- list('gexp.vec' = gexp.vec,
                         'gexp.meta' = gexp.meta,
                         'gexp.dge' = gexp.dge,
                         'plot.type' = plot.type)
    return(bulkDataList)
  })
  
  output$bulkBox <- bindEvent(renderPlot({
    bulkData <- bulkPlotData()
    
    # process data
    comp.vec <- isolate(bulkData$gexp.meta[, input$bulkComp])
    names(comp.vec) <- rownames(bulkData$gexp.meta)
    use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
    
    # create plot.df
    plot.df <- data.frame('Expression' = bulkData$gexp.vec[use.samps],
                          'Meta' = as.factor(comp.vec[use.samps]))
    # generate plot
    box.plot <- ggplot(plot.df, aes(x = Meta, y = Expression)) +
      geom_boxplot(aes(fill = Meta), color = 'black') +
      ggtitle(paste(input$bulkGene, ': Group Expression', sep = '')) + 
      xlab(input$bulkComp) + ylab('Log2 TPM') + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank())
    
    return(box.plot)
  }), bulkPlotData())
  
  output$bulkDot <- bindEvent(renderPlot({
    bulkData <- bulkPlotData()
    
    if (bulkData$plot.type == 1) {
      # extract data from dge object
      dge.gene.p <- sapply(bulkData$gexp.dge[[input$bulkComp]], function(x) {
        x[input$bulkGene, 'p.val']
      })
      dge.gene.rbs <- sapply(bulkData$gexp.dge[[input$bulkComp]], function(x) {
        x[input$bulkGene, 'rbs.cor']
      })
      # creat plot df
      plot.df <- data.frame('p.val' = dge.gene.p,
                            'rbsc' = dge.gene.rbs,
                            'level' = names(dge.gene.p),
                            'gene' = rep(input$bulkGene, length(dge.gene.p)))
      # generate plot
      dot.plot <- ggplot(plot.df, aes(level, gene)) +
        geom_point(aes(color = rbsc, size = p.val)) + 
        scale_size(trans = 'reverse') + 
        scale_color_gradient(low = 'blue', high = 'red') +
        xlab(input$bulkComp) + ylab(input$bulkGene) + 
        theme(axis.text.x = element_text(angle=45, hjust=1))
      
      return(dot.plot)
    } else {
      # extract data from dge object
      dge.gene.p <- sapply(bulkData$gexp.dge[[input$bulkComp]], function(x) {
        x[input$bulkGene, 'NES']
      })
      dge.gene.rbs <- sapply(bulkData$gexp.dge[[input$bulkComp]], function(x) {
        x[input$bulkGene, 'PES']
      })
      # creat plot df
      plot.df <- data.frame('NES' = abs(dge.gene.p),
                            'PES' = dge.gene.rbs,
                            'level' = names(dge.gene.p),
                            'gene' = rep(input$bulkGene, length(dge.gene.p)))
      # generate plot
      dot.plot <- ggplot(plot.df, aes(level, gene)) +
        geom_point(aes(color = PES, size = NES)) + 
        scale_color_gradient(low = 'blue', high = 'red') +
        xlab(input$bulkComp) + ylab(input$bulkGene) + 
        theme(axis.text.x = element_text(angle=45, hjust=1))
      
      return(dot.plot)
      
    }
  }), bulkPlotData())
  
  output$bulkSurvival <- bindEvent(renderPlot({
    bulkData <- bulkPlotData()
    
    # set label and meta data 
    comp.name <- input$bulkComp
    meta.mat <- bulkData$gexp.meta
    # process data
    comp.vec <- meta.mat[, comp.name]
    names(comp.vec) <- rownames(meta.mat)
    use.samps <- names(comp.vec)[which(!is.na(comp.vec))]
    km.df <- data.frame('time' = meta.mat[use.samps, 'Survival'],
                        'status' = as.numeric(meta.mat[use.samps, 'censor.status']),
                        'var' = comp.vec[use.samps])
    km.df <- km.df[complete.cases(km.df),]
    # fit the model
    km.fit <- survfit(formula = Surv(time, status) ~ var, data = km.df)
    km.dif <- survdiff(Surv(time) ~ var, data = km.df)
    km.p <- round(pchisq(km.dif$chisq, df = 2, lower.tail = FALSE), 3)
    km.p <- format(km.p, scientific = TRUE, digits = 3)
    # plot
    comp.title <- gsub('_', ' ', comp.name)
    legend.labs <- sapply(names(km.fit$strata), function(x) {strsplit(x, '=')[[1]][2]})
    survival.plot <- ggsurvplot(km.fit, data = km.df,
                                legend.labs = legend.labs,
                                legend.title = comp.title) +
      labs(title = paste('Kaplan-Meier: ', comp.title, sep = ''),
           subtitle = paste('Chi-Square p-value = ', km.p, sep = '')) +
      theme_survminer(font.main = c(16, "plain", "black"),
                      font.submain = c(12, "plain", "black"),
                      legend = 'right')
    return(survival.plot)
  }), bulkPlotData())
  ###############
}

shinyApp(ui, server)

