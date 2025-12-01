library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(plotly)
library(gplots)
library(dplyr)
library(WGCNA)
library(dynamicTreeCut)
# 建议加上这两个，防止报错
options(stringsAsFactors = FALSE)
enableWGCNAThreads() 

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  titlePanel("Complete WGCNA Analysis Pipeline (Go-Dissim to R-Official)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Step 1: Load & Cluster"),
      actionButton("loadBtn", "Load Dissimilarity Matrix", 
                   class = "btn-primary btn-lg",
                   style = "width: 100%; margin-bottom: 10px;"),
      actionButton("clusterBtn", "Run Clustering (hclust)",
                   class = "btn-success btn-lg",
                   style = "width: 100%; margin-bottom: 20px;"),
      
      conditionalPanel(
        condition = "output.clusteringDone",
        hr(),
        
        h4("Step 2: Detect & Merge Modules"),
        p(class = "text-info", "Aligns with WGCNA 'blockwiseModules' logic."),
        sliderInput("cutHeight", "Merge Threshold (MEDiss):",
                    min = 0, max = 1, value = 0.25, step = 0.05),
        helpText("0.25 corresponds to correlation 0.75"),
        numericInput("minModuleSize", "Min Module Size:",
                     value = 30, min = 10, max = 200, step = 10),
        actionButton("cutBtn", "Detect Modules (DynamicTree + Merge)",
                     class = "btn-info",
                     style = "width: 100%; margin-bottom: 20px;"),
        
        hr(),
        
        h4("Step 3: Deep Analysis"),
        actionButton("runAnalysis", "Calc Eigengenes & Hubs",
                     class = "btn-warning btn-lg",
                     style = "width: 100%; margin-bottom: 10px;"),
        actionButton("runEnrichment", "Run GO/KEGG Enrichment",
                     class = "btn-danger btn-lg",
                     style = "width: 100%; margin-bottom: 20px;"),
        
        hr(),
        
        h4("Export Results"),
        downloadButton("downloadModules", "Module Assignments", style = "width: 100%; margin-bottom: 5px;"),
        downloadButton("downloadHub", "Hub Genes (kME)", style = "width: 100%; margin-bottom: 5px;"),
        downloadButton("downloadEdges", "Cytoscape Edges", style = "width: 100%;")
      )
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "mainTabs",
        
        # Tab 1: Overview & Log
        tabPanel("Overview",
                 icon = icon("info-circle"),
                 br(),
                 h3("Pipeline Status Log"),
                 verbatimTextOutput("statusLog"),
                 hr(),
                 conditionalPanel(
                   condition = "output.modulesDetected",
                   fluidRow(
                     column(6, h4("Module Counts"), tableOutput("moduleCountTable")),
                     column(6, h4("Module Size Dist"), plotlyOutput("moduleSizePlot", height = "300px"))
                   )
                 )
        ),
        
        # Tab 2: Dendrogram (Official Style)
        tabPanel("Dendrogram & Colors",
                 icon = icon("tree"),
                 br(),
                 h3("Clustering Tree & Module Colors"),
                 p("This matches the 'plotDendroAndColors' from standard WGCNA."),
                 plotOutput("dendroPlot", height = "700px")
        ),
        
        # Tab 3: Eigengene Networks (Heatmap)
        tabPanel("Eigengene Network",
                 icon = icon("th"),
                 br(),
                 h3("Eigengene Adjacency Heatmap"),
                 p("Visualizes the correlation between merged modules (Meta-modules)."),
                 plotOutput("eigengeneHeatmap", height = "600px")
        ),
        
        # Tab 4: Module Table
        tabPanel("Module Genes",
                 icon = icon("list"),
                 br(),
                 h3("Gene-Module Assignments"),
                 DTOutput("moduleTable")
        ),
        
        # Tab 5: Hub Genes
        tabPanel("Hub Genes (kME)",
                 icon = icon("star"),
                 br(),
                 fluidRow(
                   column(4, selectInput("hubModule", "Select Module Color:", choices = NULL)),
                   column(4, numericInput("topN", "Show Top N:", value = 20))
                 ),
                 DTOutput("hubTable")
        ),
        
        # Tab 6: Enrichment
        tabPanel("Enrichment",
                 icon = icon("dna"),
                 br(),
                 fluidRow(
                   column(4, selectInput("enrichModule", "Select Module:", choices = NULL)),
                   column(4, radioButtons("enrichType", "Type:", choices = c("GO" = "go", "KEGG" = "kegg")))
                 ),
                 DTOutput("enrichTable")
        )
      )
    )
  )
)

# ==============================================================================
# Helper: Check Bioconductor packages
# ==============================================================================
check_bioc_packages <- function() {
  required <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
  missing <- c()
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  return(list(available = (length(missing) == 0), missing = missing))
}

# ==============================================================================
# Server
# ==============================================================================
server <- function(input, output, session) {
  
  # Reactive values to store pipeline state
  vals <- reactiveValues(
    dissim_matrix = NULL,
    gene_names = NULL,
    hclustObj = NULL,
    
    expr = NULL,        # (Samples x Genes)
    
    # Module results
    dynamicMods = NULL, # Pre-merge labels
    mergedMods = NULL,  # Post-merge labels (colors)
    MEs = NULL,         # Module Eigengenes
    
    datKME = NULL,      # kME matrix
    hub_genes = NULL,   # Calculated Hubs
    
    log = "Ready. Please Load Dissimilarity Matrix.\n",
    
    # State flags
    dataLoaded = FALSE,
    clusteringDone = FALSE,
    modulesDetected = FALSE,
    analysisComplete = FALSE,
    enrichmentComplete = FALSE
  )
  
  # Output flags for UI conditionals
  output$dataLoaded <- reactive({ vals$dataLoaded })
  output$clusteringDone <- reactive({ vals$clusteringDone })
  output$modulesDetected <- reactive({ vals$modulesDetected })
  output$analysisComplete <- reactive({ vals$analysisComplete })
  output$enrichmentComplete <- reactive({ vals$enrichmentComplete })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  outputOptions(output, "clusteringDone", suspendWhenHidden = FALSE)
  outputOptions(output, "modulesDetected", suspendWhenHidden = FALSE)
  outputOptions(output, "analysisComplete", suspendWhenHidden = FALSE)
  outputOptions(output, "enrichmentComplete", suspendWhenHidden = FALSE)
  
  # Update Log Helper
  addLog <- function(msg) {
    vals$log <- paste0(vals$log, msg, "\n")
  }
  
  output$statusLog <- renderText({ vals$log })
  
  # --------------------------------------------------------------------------
  # Step 1: Load Data
  # --------------------------------------------------------------------------
  observeEvent(input$loadBtn, {
    req(file.exists("dissimilarity_matrix.csv"))
    
    withProgress(message = 'Loading Data...', value = 0, {
      addLog("=== Loading Data ===")
      
      # 1. Load Dissimilarity (Go output)
      # Assuming CSV has headers and row names are the first column or implied
      raw_dis <- fread("dissimilarity_matrix.csv", data.table = FALSE)
      if(class(raw_dis[,1]) == "character") {
        rownames(raw_dis) <- raw_dis[,1]
        raw_dis <- raw_dis[,-1]
      }
      vals$dissim_matrix <- as.matrix(raw_dis)
      vals$gene_names <- colnames(vals$dissim_matrix)
      addLog(sprintf("Loaded Dissimilarity Matrix: %d genes", ncol(vals$dissim_matrix)))
      
      # 2. Load Expression Data (For Eigengenes calculation)
      # WGCNA requires Expr for Eigengenes even if tree is built from ext. matrix
      if(file.exists("clean_thyroid_matrix.csv")) {
        expr_raw <- fread("clean_thyroid_matrix.csv", data.table = FALSE)
        # Transpose logic matches your provided script
        # Assuming input is Rows=Genes, Cols=Samples -> Convert to WGCNA standard (Rows=Samples)
        row.names(expr_raw) <- expr_raw[,1]
        expr_raw <- expr_raw[,-1]
        
        # Transpose
        datExpr <- as.data.frame(t(expr_raw))
        
        # Ensure numeric
        datExpr[] <- lapply(datExpr, function(x) as.numeric(as.character(x)))
        
        # Match genes with dissimilarity matrix
        common_genes <- intersect(vals$gene_names, colnames(datExpr))
        vals$expr <- datExpr[, common_genes]
        vals$dissim_matrix <- vals$dissim_matrix[common_genes, common_genes]
        vals$gene_names <- common_genes
        
        addLog(sprintf("Matched Expression Matrix: %d samples, %d genes", nrow(vals$expr), ncol(vals$expr)))
        vals$dataLoaded <- TRUE
      } else {
        addLog("ERROR: clean_thyroid_matrix.csv not found.")
      }
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 2: Clustering
  # --------------------------------------------------------------------------
  observeEvent(input$clusterBtn, {
    req(vals$dataLoaded)
    withProgress(message = 'Clustering...', value = 0.5, {
      addLog("Running hclust(method='average')...")
      vals$hclustObj <- hclust(as.dist(vals$dissim_matrix), method = "average")
      vals$clusteringDone <- TRUE
      addLog("Clustering Complete.")
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 3: Module Detection (Official WGCNA Logic)
  # --------------------------------------------------------------------------
  observeEvent(input$cutBtn, {
    req(vals$clusteringDone)
    
    withProgress(message = 'Detecting & Merging Modules...', value = 0, {
      addLog("=== Module Detection ===")
      
      # 1. Dynamic Tree Cut
      incProgress(0.3, detail = "Cutting Tree...")
      dynamicMods <- cutreeDynamic(
        dendro = vals$hclustObj,
        distM = vals$dissim_matrix,
        deepSplit = 4,
        pamRespectsDendro = FALSE,
        minClusterSize = input$minModuleSize
      )
      
      dynamicColors <- labels2colors(dynamicMods)
      
      # Debug 1: 初始模块数量
      counts <- table(dynamicMods)
      counts_str <- paste(names(counts), counts, collapse = "; ")
      addLog(sprintf("Initial dynamicMods = %d modules", length(unique(dynamicMods))))
      addLog(paste("Initial module counts:", counts_str))
      
      # 2. Merge Close Modules
      incProgress(0.6, detail = "Merging close modules...")
      MEList <- moduleEigengenes(vals$expr, colors = dynamicColors)
      MEs <- MEList$eigengenes
      
      # 安全检查：是否有 NA
      if (any(is.na(MEs))) {
        addLog("ERROR: Module eigengenes contain NA. Merge aborted.")
        return(NULL)
      }
      
      MEDiss <- 1 - cor(MEs)
      METree <- hclust(as.dist(MEDiss), method = "average")
      
      # Debug 2: 画 MEs 聚类树
      dir.create("www", showWarnings = FALSE)
      png("www/ME_clustering_debug.png", width = 1000, height = 700)
      plot(METree, main = "MEs clustering before merge")
      abline(h = input$cutHeight, col = "red", lwd = 2)
      dev.off()
      addLog("Saved ME_clustering_debug.png in /www/")
      
      # 真正 merge
      merge <- mergeCloseModules(
        vals$expr, dynamicColors,
        cutHeight = input$cutHeight,
        verbose = 0
      )
      
      vals$mergedMods <- merge$colors
      vals$MEs <- merge$newMEs
      
      addLog(sprintf("Merged Modules: %d final modules",
                     length(unique(vals$mergedMods))))
      
      vals$modulesDetected <- TRUE
    })
  })
  
      
  
  # --------------------------------------------------------------------------
  # Step 4: Analysis (kME, Hubs)
  # --------------------------------------------------------------------------
  observeEvent(input$runAnalysis, {
    req(vals$modulesDetected)
    
    withProgress(message = 'Calculating kME...', value = 0.5, {
      addLog("Calculating signedKME (Module Membership)...")
      
      # Official WGCNA function: signedKME
      vals$datKME <- signedKME(vals$expr, vals$MEs)
      
      # Prepare Hub Genes table
      gene_colors <- vals$mergedMods
      hubs <- data.frame(GeneID = colnames(vals$expr), Module = gene_colors)
      
      # Attach kME for the assigned module
      kme_values <- c()
      for(i in 1:nrow(hubs)) {
        mod <- hubs$Module[i]
        colname <- paste0("kME", mod)
        if(colname %in% colnames(vals$datKME)) {
          kme_values[i] <- vals$datKME[i, colname]
        } else {
          kme_values[i] <- NA # Grey module usually
        }
      }
      hubs$kME <- kme_values
      vals$hub_genes <- hubs
      
      vals$analysisComplete <- TRUE
      addLog("kME Calculation Done.")
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 5: Enrichment (Fixing the cutoff part)
  # --------------------------------------------------------------------------
  observeEvent(input$runEnrichment, {
    req(vals$analysisComplete)
    
    checks <- check_bioc_packages()
    if(!checks$available) {
      showModal(modalDialog(
        title = "Missing Packages",
        paste("Please install:", paste(checks$missing, collapse=", ")),
        easyClose = TRUE
      ))
      return()
    }
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    
    vals$enrichmentComplete <- TRUE # Flag to show UI
    addLog("Enrichment libraries loaded. Ready to display results.")
  })
  
  # ==============================================================================
  # OUTPUTS & PLOTS
  # ==============================================================================
  
  # --- 1. Dendrogram & Colors (Official WGCNA Style) ---
  output$dendroPlot <- renderPlot({
    req(vals$hclustObj, vals$mergedMods)
    
    # Using the standard WGCNA plotting function
    # plotDendroAndColors requires standard margins setup inside
    plotDendroAndColors(vals$hclustObj, 
                        vals$mergedMods,
                        "Module Colors",
                        dendroLabels = FALSE, 
                        hang = 0.03,
                        addGuide = TRUE, 
                        guideHang = 0.05,
                        main = "Gene Dendrogram and module colors (Official Pipeline)")
  })
  
  # --- 2. Eigengene Network Heatmap (Official WGCNA Style) ---
  output$eigengeneHeatmap <- renderPlot({
    req(vals$MEs)
    
    # Plotting the relationship between modules
    plotEigengeneNetworks(vals$MEs, 
                          "Eigengene adjacency heatmap", 
                          marHeatmap = c(3,4,2,2),
                          plotDendrograms = TRUE, 
                          xLabelsAngle = 90)
  })
  
  # --- 3. Module Stats Table ---
  output$moduleCountTable <- renderTable({
    req(vals$mergedMods)
    counts <- table(vals$mergedMods)
    as.data.frame(counts)
  })
  
  # --- 4. Module Size Plot ---
  output$moduleSizePlot <- renderPlotly({
    req(vals$mergedMods)
    df <- as.data.frame(table(vals$mergedMods))
    colnames(df) <- c("Module", "Count")
    # Sort by count desc
    df <- df[order(-df$Count),]
    
    p <- ggplot(df, aes(x=reorder(Module, -Count), y=Count, fill=Module)) +
      geom_bar(stat="identity") +
      scale_fill_identity() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Module Color", y = "Gene Count")
    
    ggplotly(p)
  })
  
  # --- 5. Module Assignments DT ---
  output$moduleTable <- renderDT({
    req(vals$mergedMods)
    df <- data.frame(GeneID = vals$gene_names, 
                     Module = vals$mergedMods)
    datatable(df, options = list(pageLength = 10))
  })
  
  # --- 6. Hub Genes DT ---
  output$hubTable <- renderDT({
    req(vals$hub_genes, input$hubModule)
    
    # Filter by selected module
    sub_df <- vals$hub_genes %>% 
      filter(Module == input$hubModule) %>%
      arrange(desc(kME)) %>%
      head(input$topN)
    
    datatable(sub_df, options = list(pageLength = 10)) %>%
      formatRound(columns=c("kME"), digits=4)
  })
  
  # --- 7. Enrichment DT (Dynamic Calc) ---
  output$enrichTable <- renderDT({
    req(vals$enrichmentComplete, input$enrichModule)
    
    # Get genes in module
    mod_genes <- vals$gene_names[vals$mergedMods == input$enrichModule]
    
    withProgress(message = 'Running Enrichment...', {
      
      # Assuming GeneIDs are SYMBOLS. Convert to ENTREZID.
      # Modify 'keyType' if your input genes are ENSEMBL or other.
      gene_entrez <- bitr(mod_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      
      if(nrow(gene_entrez) == 0) return(NULL)
      
      res_df <- NULL
      
      if(input$enrichType == "go") {
        ego <- enrichGO(gene = gene_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP", # Biological Process
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE)
        if(!is.null(ego)) res_df <- as.data.frame(ego)
      } else {
        ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                            organism = 'hsa',
                            pvalueCutoff = 0.05)
        if(!is.null(ekegg)) res_df <- as.data.frame(ekegg)
      }
      
      if(is.null(res_df) || nrow(res_df) == 0) {
        return(data.frame(Result = "No significant enrichment found."))
      }
      
      # Return simple table
      res_df[, c("ID", "Description", "p.adjust", "geneID")]
    })
  }, selection = 'single')
  
  # ==============================================================================
  # Downloads
  # ==============================================================================
  output$downloadModules <- downloadHandler(
    filename = function() { "module_assignments.csv" },
    content = function(file) {
      df <- data.frame(GeneID = vals$gene_names, Module = vals$mergedMods)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$downloadHub <- downloadHandler(
    filename = function() { "hub_genes_kME.csv" },
    content = function(file) {
      write.csv(vals$hub_genes, file, row.names = FALSE)
    }
  )
  
  # Export for Cytoscape (Edges)
  # This requires recalculating TOM or using the Dissim matrix
  output$downloadEdges <- downloadHandler(
    filename = function() { paste0("Cytoscape_Edges_", input$hubModule, ".txt") },
    content = function(file) {
      # Use the dissimilarity matrix as the weight source (1 - dissim = sim)
      # Extract genes for the specific module
      mod_genes <- vals$gene_names[vals$mergedMods == input$hubModule]
      if(length(mod_genes) > 1000) {
        # Limit to top 1000 hubs to prevent crash
        mod_hubs <- vals$hub_genes %>% 
          filter(Module == input$hubModule) %>% 
          arrange(desc(kME)) %>% 
          head(1000)
        mod_genes <- mod_hubs$GeneID
      }
      
      # Subset matrix
      mod_sim <- 1 - vals$dissim_matrix[mod_genes, mod_genes]
      
      # Export using WGCNA function
      exportNetworkToCytoscape(mod_sim,
                               edgeFile = file,
                               nodeFile = NULL,
                               weighted = TRUE,
                               threshold = 0.02, # Cutoff
                               nodeNames = mod_genes)
    }
  )
}

shinyApp(ui, server)