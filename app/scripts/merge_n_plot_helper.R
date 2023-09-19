########################################################
## Title: merge_n_plot_helper
## Date: last edited 10Aug2023
## Author: Genevieve Mortensen from Nick Powell
## Purpose: Create helper functions for web-hosted app
########################################################

## Set working directory
#setwd("/N/project/mirilizer/tool-development/app/")
# 
# # ## Testing variables
# bedfiles <- c("/N/project/mirilizer/prim_hep/run1/eclipse_run1.bed",
#                "/N/project/mirilizer/prim_hep/run2/eclipse_run2.bed",
#               "/N/project/mirilizer/prim_hep/run3/prim_hep_kit_run1.bed",
#               "/N/project/mirilizer/prim_hep/run4/eclipse_run3.bed",
#               "/N/project/mirilizer/prim_hep/run5/eclipse_run3_gene_specific1.bed",
#               "/N/project/mirilizer/prim_hep/run6/eclipse_run3_gene_specific2.bed",
#               "/N/project/mirilizer/prim_hep/run7/prim_hep_kit_run2.bed",
#               "/N/project/mirilizer/prim_hep/run8/prim_hep_kit_run3.bed"
# )
# gene = "CYP3A4"
# mir = "hsa-miR-122-5p"
# rnafiles = "/N/project/mirilizer/tool-development/app/testing/rnaseq_CYP3A4.bdg"

##############################################################
# GET COUNT INFORMATION FOR ALL GENES AND MIRS IN BEDFILE(S) # 
##############################################################
getCountInfo <- function(bedfiles){
  
  #############################################
  # Create master data frame from the bed files
  #############################################
  bedfilelist <- list()
  for(i in 1:length(bedfiles)){
    bedfilelist[[i]] <- read_table(bedfiles[i], col_names = FALSE, skip = 1)
    bedfilelist[[i]][5:6] <- str_split_fixed(bedfilelist[[i]]$X4,"_",2)
    bedfilelist[[i]][7] <- paste0("run", i)
    colnames(bedfilelist[[i]]) <- c("chr","start","stop","name","mir_name","gene_name","run")
  }
  
  ## Condition for single or multiple bed files
  if (length(bedfilelist) > 1){
    bedmaster <- bedfilelist[[1]]
    for(i in 2:length(bedfilelist)){
      bedmaster <- rbind(bedmaster, bedfilelist[[i]])
    }
  } else {
    bedmaster <- bedfilelist[[1]]
  }
  
  ##########################
  # Create counts tables
  ##########################
  ## create miR count table
  mirfreq <- as.data.frame(table(bedmaster$mir_name))
  mirfreq <- mirfreq[order(mirfreq$Freq,decreasing = TRUE),]
  names(mirfreq)[names(mirfreq) == "Var1"] <- "miR"
  names(mirfreq)[names(mirfreq) == "Freq"] <- "Counts"
  
  ## create gene count table
  genefreq <- as.data.frame(table(bedmaster$gene_name))
  genefreq <- genefreq[order(genefreq$Freq,decreasing = TRUE),]
  names(genefreq)[names(genefreq) == "Var1"] <- "Gene"
  names(genefreq)[names(genefreq) == "Freq"] <- "Counts"
  
  #######################################
  ## Create all genes and all mirs tables
  #######################################
  allgenes <- as.data.frame(genefreq$Gene)
  names(allgenes)[names(allgenes)=="genefreq$Gene"] <- "X1"
  allmirs <- as.data.frame(mirfreq$miR)
  names(allmirs)[names(allmirs)=="mirfreq$miR"] <- "X1"
  
  return(list("bedmaster" = bedmaster, 
              "mirfreq" = mirfreq, 
              "genefreq" = genefreq, 
              "allgenes" = allgenes, 
              "allmirs" = allmirs, 
              "bedfilelist" = bedfilelist))
}




##########################################################
# GET INFORMATION BASED ON SELECTED GENE OR SELECTED MIR # 
##########################################################
getSelectInfo <- function(gene, 
                          mir, 
                          bedmaster, 
                          allgenes, 
                          allmirs){
  
  ###############################################
  ## Subset relevant genes if the argument exists
  ###############################################
  if(hasArg(gene)){
    
    if(is.numeric(gene)){ gene_symbol <- allgenes$X1[gene] } #If the gene begins with a number, subset from allgenes
    if(is.character(gene)){ gene_symbol <- gene } #if it begins with a character
    bed_gene <- subset(bedmaster, bedmaster$gene_name==gene_symbol) #Subset gene of interest from bedmaster
    
    tblg <- as.data.frame(table(bed_gene$mir_name)) #Subset miRs from gene-specific bed file
    tblg <- tblg[order(tblg$Freq,decreasing = TRUE),] #Sort by decreasing abundance
    names(tblg)[names(tblg) == "Var1"] <- "miR"
    names(tblg)[names(tblg) == "Freq"] <- "Counts"
  } else {
    tblg <- NULL
  }
    
  #############################################
  # Subset relevant mirs if the argument exists
  #############################################
  if(hasArg(mir)){
    
    if(is.numeric(mir)){ mir_symbol <- allmirs$X1[mir] } #If the gene begins with a number, subset from allgenes
    if(is.character(mir)){ mir_symbol <- mir } #if it begins with a character
    
    bed_mir <- subset(bedmaster, bedmaster$mir_name==mir_symbol) #subset miR of interest from bedmaster
    
    tblm <- as.data.frame(table(bed_mir$gene_name)) #Subset miRs from gene-specific bed file
    tblm <- tblm[order(tblm$Freq,decreasing = TRUE),] #Sort by decreasing abundance
    names(tblm)[names(tblm) == "Var1"] <- "Gene"
    names(tblm)[names(tblm) == "Freq"] <- "Counts"
  } else {
    tblm <- NULL
  }
  
  return(list("tblm" = tblm, 
              "tblg" = tblg, 
              "gene_symbol" = gene_symbol))
}

######################################
# COLLECT FILES, TIDY DATA, AND PLOT #
######################################

plotter <- function(gene_symbol, 
                    mir, 
                    bedmaster, 
                    allgenes, 
                    allmirs, 
                    bedfilelist, 
                    bedfiles,
                    rnafiles){
  
  ####################################################################
  ## Load rnaseq files - user required upload - contigent on existence
  ####################################################################
  if (hasArg(rnafiles)){
    #Create list from rnaseq files if multiple are uploaded
    rnafilelist <- list()
    for(i in 1:length(rnafiles)){
      rnafilelist[[i]] <- read_table(rnafiles[[i]], col_names = FALSE, skip = 1)
    }
    #Create rnaseq list
    rnaseq_list <- list() #caveat that these must be uploaded with the bedfiles addressed later
    #Handle if one or more rnaseq file exist
    if (length(rnafilelist) < 2){
      rnaseq_list[[1]] <- data.frame(read_table(rnafiles, col_names = FALSE))[2:3]
      colnames(rnaseq_list[[1]]) <- c("gen_coord",paste0("run1_depth"))
    } else {
      for(n in 1:length(rnafiles)){ #these are specific to the run and exist as one file per bedfile
        rnaseq_list[[n]] <- read_table(rnafile, col_names = FALSE)
        colnames(rnaseq_list[[n]]) <- c("gen_coord",paste0("run",n,"_depth"))
        }
      }
    } else {
      rnaseq_list <- NULL
    }
  
  ###########################################################
  # Load run bdg files - rendered from user uploaded bedfiles
  ###########################################################
  #Arrange bedfiles for runs_list
  bedpaths <- list()
  for(n in 1:length(bedfiles)){
    bedpaths[[n]] <- read_table(bedfiles[[n]], col_names = FALSE, skip = 1)
  }
  
  #set chromosome variable names
  chr <- c(paste0("chr",seq(1,23)),"chrX","chrY")
  
  # Load each bdg file for each miR and arrange by run
  runs_list <- list()
  
  for(n in 1:length(bedpaths)){ #caveat that bedfiles need to be in the correct run order
    #create subset of bedmaster 
    bedsub <- subset(bedpaths[[n]], bedpaths[[n]][["X4"]]==paste0(mir,"_",gene_symbol))
    if(nrow(bedsub)==0){next}
    names(bedsub) <- c("chrom","chromStart","chromEnd","name")
    bedsub <- subset(bedsub, bedsub$chrom %in% chr) #cleaning step to make sure only one chromosome exists in the data for that gene
    bedsub <- bedsub[order(bedsub$chrom, bedsub$chromStart),]
    rownames(bedsub) <- NULL
    
    #make bedgraph using IRanges (which will be in 1-based coords for start and end positions)
    x <- IRanges(start=bedsub$chromStart,end=bedsub$chromEnd)
    xcov <- coverage(x)
    start <- as.numeric(as.character(bedsub[1,2]))
    end <- length(xcov)
    xcov <- xcov[start:end]
    xcov <- as.data.frame(xcov)
    xcov$gen_coord <- seq(from=start, to=end)
    xcov <- subset(xcov, xcov$value>0)
    
    #Attach data as data frame to the runs list
    runs_list[[n]] <- data.frame(xcov[c(2,1)])
    colnames(runs_list[[n]]) <- c("gen_coord",paste0(mir,"_run",n)) #paste list of mirs
  }

  ###############################################################
  # Load allchimeric reads - rendered from user uploaded bedfiles
  ###############################################################
  ## Load individual run bdg files for all chimeric reads
  runs_allchim_list <- list()
  for(n in 1:length(bedfiles)){ #this exists as one file, so no list nesting
    
      bedpaths[[n]][5:6] <- str_split_fixed(bedpaths[[n]]$X4,"_",2)
      #subset bed to selected gene
      bedsub <- subset(bedpaths[[n]], bedpaths[[n]][["V2"]]==paste0(gene_symbol))
      if(nrow(bedsub)==0){next}
      names(bedsub) <- c("chrom","chromStart","chromEnd","name","mir","gene")
      bedsub <- subset(bedsub, bedsub$chrom %in% chr)
      bedsub <- bedsub[order(bedsub$chrom, bedsub$chromStart),]
      rownames(bedsub) <- NULL
      
      #make bedgraph using IRanges (which will be in 1-based coords for start and end positions)
      x <- IRanges(start=bedsub$chromStart,end=bedsub$chromEnd)
      xcov <- coverage(x)
      start <- as.numeric(as.character(bedsub[1,2]))
      end <- length(xcov)
      xcov <- xcov[start:end]
      xcov <- as.data.frame(xcov)
      xcov$gen_coord <- seq(from=start, to=end)
      xcov <- subset(xcov, xcov$value>0)
      
      runs_allchim_list[[n]] <- data.frame(xcov[c(2,1)])
      colnames(runs_allchim_list[[n]]) <- c("gen_coord",paste0("allchim_run",n))
  }
  
  ##############################################
  # Load canon and genom files - processed by us
  ##############################################
  # Subset canonical gene information
  a <- read_delim(paste0("data/genom_files/genom_", gene_symbol, ".txt"), delim = " ", col_types = cols())
  b <- read_delim(paste0("data/canon_files/canon_", gene_symbol, ".txt"), delim = " ", col_types = cols())
  colnames(b)[3] <- "canon_region"
  gen_canon <- full_join(a, b, by="gen_coord")
  
  ####################################
  # Load rnaup files - processed by us
  ####################################
  a <- paste0("data/rnaup/", gene_symbol, "/", mir, "_", gene_symbol, "_w25_u3.out")
  a <- read_table(a, comment = "#", col_names = FALSE)
  # Handle rnaup output and store results for each miR in a list
  rnaup_list <- list()
  # for(n in 1:length(mirs_to_plot)){
  rnaup_list[[1]] <- read_table(paste0("data/rnaup/",gene_symbol,"/",mir,"_",gene_symbol,"_w25_u3.out"),comment = "#",col_names = FALSE)
  rnaup_list[[1]] <- rnaup_list[[1]][c(1:nrow(b)),]
  rnaup_list[[1]] <- rnaup_list[[1]][-c(2,3,6)]
  names(rnaup_list[[1]]) <- c("canon_pos","u19s",paste0(mir,"_dG"))
  # }
  
  ####################################
  ## Load seed files - processed by us
  ####################################
  seeds <- read.table(paste0("data/seeds/seeds_",gene_symbol,".txt"),header = TRUE) #table of all seeds for the gene
  seed_list <- list()
  mir1 <- gsub("-",".",mir) #sub out dashes with dots
  # for(n in 1:length(mirs_to_plot)){ #create list of seeds subsetted by miR
  seed_list[[1]] <- seeds[colnames(seeds) %in% c("gen_coord",paste0(mir1,"_seedmatch"))]
  # }
  
  #####################################
  # Load duplex files - processed by us
  #####################################
  dup <- read.table(paste0("data/rnaduplex/duplex_",gene_symbol,".tsv"),header = TRUE) #dataframe of duplexes for gene
  duplex_list <- list()
  # for(n in 1:length(mirs_to_plot)){ #list of duplexes subsetted by miR
  duplex_list[[1]] <- dup[colnames(dup) %in% c("gen_coord",paste0(mir1,"_duplex_dG"))]
  # }
  
  #########################
  #  Tidy statistical data
  #########################
  bed_mir <- subset(bedmaster, bedmaster$mir_name==mir)
  bed_gene <- subset(bedmaster, bedmaster$gene_name==gene_symbol)

  #####################
  # Run statistics data
  #####################
  gg_tbl <- full_join(bed_mir %>% 
                        group_by(run) %>% 
                        summarize(distinct_genes=n_distinct(gene_name),
                                  avg_reads_per_gene=round(mean(table(gene_name)),
                                                           1)),
                      bed_gene %>% group_by(run) %>% 
                        summarize(distinct_mirs=n_distinct(mir_name),
                                  avg_reads_per_mir=round(mean(table(mir_name)),
                                                          1)),
                      by="run")
  # distinct_things is the number of distinct things for that run
  # avg_reads_per_thing is the rounded mean (nearest tenth) of the number of reads of the thing for that run.
  # Each thing, gene or miR, is subsetted from the selected gene or miR appearing in the data.
  #     e.g.) for gene CYP3A4 there are an average of 8.4 reads of miR hsa-miR-122-5p
  
  ############################################
  # Handle tidy data contigent on RNA seq data
  ############################################
  # Create conditional statement for handling rnafiles
  #contigency on rnafile argument existence in function
  if (hasArg(rnafiles)){
    if (length(rnafiles) < 2) {
      #join rnaseq data and canonical gene information
      rnaseq <- left_join(gen_canon[c(2,7)],rnaseq_list[[1]],by="gen_coord")
    } else {
      rnaseq <- left_join(gen_canon[c(2,7)],rnaseq_list[[1]],by="gen_coord")
      for(i in 2:length(rnafiles)){
        rnaseq <- left_join(rnaseq,rnaseq_list[[i]],by="gen_coord")
      }
    }
    # Only merge if rnaseq exists
    rnaseq <- rnaseq[-c(1)]
    rnaseq <- reshape2::melt(rnaseq,id = "pos")
    # Create plotting factor 2 with rnaseq data
    if (length(bedfiles) < 2){
      mir_reads <- left_join(gen_canon[c(2,7)],runs_list[[1]],by="gen_coord")
      allchim_reads <- left_join(gen_canon[c(2,7)],runs_allchim_list[[1]],by="gen_coord")
    } else {
      mir_reads <- left_join(gen_canon[c(2,7)],runs_list[[1]],by="gen_coord")
      allchim_reads <- left_join(gen_canon[c(2,7)],runs_allchim_list[[1]],by="gen_coord")
      for (i in 2:length(bedfiles)){
        mir_reads <- left_join(mir_reads,runs_list[[i]],by="gen_coord")
        allchim_reads <- left_join(allchim_reads,runs_allchim_list[[i]],by="gen_coord")
      }
    }
    allchim_reads <- allchim_reads[-c(1)]
    allchim_reads <- reshape2::melt(allchim_reads,id = "pos")
    mir_reads <- mir_reads[-c(1)]
    mir_reads <- reshape2::melt(mir_reads,id = "pos")
    factor2 <- max(allchim_reads$value,na.rm = TRUE)/max(rnaseq$value,na.rm = TRUE)
  } else {
    # Handle other files
    if (length(bedfiles) < 2){
      mir_reads <- left_join(gen_canon[c(2,7)],runs_list[[1]],by="gen_coord", copy =TRUE)
      allchim_reads <- left_join(gen_canon[c(2,7)],runs_allchim_list[[1]],by="gen_coord")
    } else {
      mir_reads <- left_join(gen_canon[c(2,7)],runs_list[[1]],by="gen_coord", copy = TRUE)
      allchim_reads <- left_join(gen_canon[c(2,7)],runs_allchim_list[[1]],by="gen_coord")
      for (i in 2:length(bedfiles)){
        mir_reads <- left_join(mir_reads,runs_list[[i]],by="gen_coord")
        allchim_reads <- left_join(allchim_reads,runs_allchim_list[[i]],by="gen_coord")
      }
    }
    # Now push everything together
    mir_reads <- mir_reads[-c(1)]
    mir_reads <- reshape2::melt(mir_reads,id = "pos") #tidy the reads
    allchim_reads <- allchim_reads[-c(1)]
    allchim_reads <- reshape2::melt(allchim_reads,id = "pos")
  }
  
    #################################
    ## Tidy native data for plotting
    #################################
    x <- left_join(gen_canon,duplex_list[[1]],by="gen_coord")
    colnames(x)[10] <- "duplex_dG"
    x <- left_join(x,rnaup_list[[1]],by="canon_pos")
    colnames(x)[11:12] <- c("u19s","rnaup_dG")
    x <- left_join(x,subset(seed_list[[1]], seed_list[[1]][,2]==5),by="gen_coord")
    colnames(x)[13] <- "seed_5"
    x <- left_join(x,subset(seed_list[[1]], seed_list[[1]][,2]==6),by="gen_coord")
    colnames(x)[14] <- "seed_6"
    factor <- max(allchim_reads$value,na.rm = TRUE)/45
    
    #Restructure seed columns in x data for plotting
    x$seed_5[is.na(x$seed_5)] = ""
    x$seed_6[is.na(x$seed_6)] = ""
    x$seeds <- as.numeric(paste(x$seed_5, x$seed_6))
  
  ################
  # Plotting data
  ################
    #Create main plot
  xx1 <- ggplot()+
      
      ## 5' UTR
      geom_rect(fill = "antiquewhite",
                alpha= 0.5,
                aes(
                  xmin = min(which(x$canon_region=="5'UTR")),
                  xmax = max(which(x$canon_region=="5'UTR")),
                  ymin = -Inf,
                  ymax = Inf))+
      
      ## 3' UTR
      geom_rect(fill = "antiquewhite",
                alpha= 0.5,
                aes(
                  xmin = min(which(x$canon_region=="3'UTR")),
                  xmax = max(which(x$canon_region=="3'UTR")),
                  ymin = -Inf,
                  ymax = Inf))+
      
      ## horizontal axis
      geom_hline(yintercept = 0,
                 color="black",
                 size=0.2)+
      
      ## duplex_dG - use guides to get continuously colored points in legend
      geom_point(data=x, 
                 aes(
                   x=pos,
                   y=duplex_dG*factor, 
                   color=duplex_dG))+
      scale_color_continuous(limits = c(-30,-10))+
      guides(color = guide_legend(title = "duplex dG",
                                  reverse = TRUE))
    
      if (hasArg(rnafiles)){
        xx <- xx1 +
          
          ## RNA sequencing reads
          geom_ribbon(data=rnaseq, 
                      aes(
                        x=pos, 
                        ymin=0,
                        ymax=value*factor2, 
                        fill=variable),
                      alpha=1)+
          scale_fill_grey(start = 0.7, end = 0.9, 
                          name = "RNA reads")+
          
          ## all chimeric reads
          new_scale_color()+
          geom_line(data=allchim_reads, 
                    aes(
                      x=pos,
                      y=value,
                      color=variable),
                    size=0.4,
                    linetype="dashed",
                    alpha=0.5)+
          scale_color_manual(values = c("darkgreen",
                                        "black",
                                        "goldenrod3",
                                        "#CC79A7",
                                        "#56B4E9",
                                        "#009E73",
                                        "#F0E442",
                                        "#0072B2",
                                        "#D55E00"), 
                             name = "all chimeric reads")+
          
          ## miR reads
          new_scale_color()+
          geom_line(data=mir_reads, 
                    aes(
                      x=pos,
                      y=value,
                      color=variable),
                    size=0.8,
                    alpha=0.7)+
          scale_color_manual(values = c("darkgreen",
                                        "black",
                                        "goldenrod3",
                                        "#CC79A7",
                                        "#56B4E9",
                                        "#009E73",
                                        "#F0E442",
                                        "#0072B2",
                                        "#D55E00"), 
                             name = "miR reads")+
          
          ## u19s
          new_scale_fill()+
          geom_ribbon(data=x, 
                      aes(
                        x=pos,
                        ymin=0,
                        ymax=-u19s*factor, 
                        fill = ""))+
          scale_fill_manual(values = c("grey50"), 
                            name = "u19s")+
          
          ## rnaup_dG
          new_scale_color()+
          geom_line(data=x, 
                    aes(
                      x=pos,
                      y=rnaup_dG*factor, 
                      color = ""))+
          scale_color_manual(values = c("navyblue"), 
                             name = "rnaup_dG")+
          
          ## plot theme and axis formatting
          theme_minimal()+
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(colour = "black",
                                         size = 0.2),
                text = element_text(size=8))+
          scale_x_continuous(expand = c(0,0))+
          ylab("deltaG                               microRNA read pileup values")+
          xlab("transcript position (all exons combined)")+
          
          ## Seed matches
          new_scale_color()+
          geom_point(data=x, 
                     aes(
                       x=pos,
                       y=seeds*factor/300, 
                       color=as.factor(seeds)), 
                     size = 2)+
          scale_color_manual(values = c("5" = "purple", 
                                        "6" = "green"), 
                             name = "seed match")
      } else {
          xx <- xx1 + 
            
            ## all chimeric reads
            new_scale_color()+
            geom_line(data=allchim_reads, 
                      aes(
                        x=pos,
                        y=value,
                        color=variable),
                      size=0.4,
                      linetype="dashed",
                      alpha=0.5)+
            scale_color_manual(values = c("darkgreen",
                                          "black",
                                          "goldenrod3",
                                          "#CC79A7",
                                          "#56B4E9",
                                          "#009E73",
                                          "#F0E442",
                                          "#0072B2",
                                          "#D55E00"), 
                               name = "all chimeric reads")+
            
            ## miR reads
            new_scale_color()+
            geom_line(data=mir_reads, 
                      aes(
                        x=pos,
                        y=value,
                        color=variable),
                      size=0.8,
                      alpha=0.7)+
            scale_color_manual(values = c("darkgreen",
                                          "black",
                                          "goldenrod3",
                                          "#CC79A7",
                                          "#56B4E9",
                                          "#009E73",
                                          "#F0E442",
                                          "#0072B2",
                                          "#D55E00"), 
                               name = "miR reads")+
            
            ## u19s
            new_scale_fill()+
            geom_ribbon(data=x, 
                        aes(
                          x=pos,
                          ymin=0,
                          ymax=-u19s*factor, 
                          fill = ""))+
            scale_fill_manual(values = c("grey50"), 
                              name = "u19s")+
            
            ## rnaup_dG
            new_scale_color()+
            geom_line(data=x, 
                      aes(
                        x=pos,
                        y=rnaup_dG*factor, 
                        color = ""))+
            scale_color_manual(values = c("navyblue"), 
                               name = "rnaup_dG")+
            
            ## plot theme and axis formatting
            theme_minimal()+
            theme(
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.line = element_line(colour = "black",
                                           size = 0.2),
                  text = element_text(size=8))+
            scale_x_continuous(expand = c(0,0))+
            ylab("deltaG                               microRNA read pileup values")+
            xlab("transcript position (all exons combined)")+
            
            ## Seed matches
            new_scale_color()+
            geom_point(data=x, 
                       aes(
                         x=pos,
                         y=seeds*factor/300, 
                         color=as.factor(seeds)), 
                       size = 2)+
            scale_color_manual(values = c("5" = "purple", 
                                          "6" = "green"), 
                               name = "seed match")
      }
    
    #create plot with no legend
    xx2 <- xx + theme(legend.position="none")

    #Create table plot
        xxx <- tableGrob(gg_tbl, rows = NULL, theme = ttheme_default(base_size = 8))
        
    #Merge plots for easy return
        xxxx <- ggarrange(xx, xxx, nrow = 2, heights = c(10,2))
        
    #Legend wrangling
        xxleg <- xx + theme(legend.direction = "horizontal")+
          guides(guide_legend(ncol = 2))
        legxx <- get_legend(xxleg)
        grid.draw(legxx)
        
    return(list("mainplot" = xx, "mainpltnoleg" = xx2, "tableplot" = xxx, "combinedplot" = xxxx, "mainpltleg" = legxx))
}
