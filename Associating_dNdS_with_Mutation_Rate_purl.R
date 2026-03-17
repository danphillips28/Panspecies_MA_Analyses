library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
#library(patchwork)
library(Rsamtools)
library(ggpubr)
library(clinfun)
library(segmented)
library(mgcv)
library(tidyverse)
library(stringr)
library(dplyr)
`%ni%` <- Negate(`%in%`)

load(file = "/home/ocdm0351/DPhil/R_Data/Annotated_VCFs_Raw")
load(file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_Ratios")
load(file = "/home/ocdm0351/DPhil/R_Data/Species_Genes_gtf")
load(file = "/home/ocdm0351/DPhil/R_Data/Species_Canonical_Transcripts_gtf")
load(file = "/home/ocdm0351/DPhil/R_Data/Species_Canonical_Mappings")

# Subset for Debugging Purposes
#Annotated_VCFs_Raw <- Annotated_VCFs_Raw[names(Annotated_VCFs_Raw) %ni% c(
#                                                                          "Weller_2014_MA_Ppacificus_El_paco")]
#Annotated_VCFs_Raw <- Annotated_VCFs_Raw[1:6]

Species_dNdS_Ratios <- lapply(Species_dNdS_Ratios, function(df) {
  df$query_id <- sub("^(.*)[.].*", "\\1", df$query_id) 
  df
})

Annotated_VCFs_Raw <- lapply(Annotated_VCFs_Raw, function(df) {
  df <- df[, 1:12] # SnpSift is adding annoying empty column names for some reason. Will remove upstream.
  df$Variant_ID <- paste0(df$CHROM, "-", df$POS, ":", df$REF, ">", df$ALT)
  df$MutationID <- paste0(df$REF,">", df$ALT)
  df <- df %>% 
  mutate(MutationalCategory = 
           ifelse(MutationID %in% c("C>A","G>T"), "C>A | G>T",
           ifelse(MutationID %in% c("C>G","G>C"), "C>G | G>C",
           ifelse(MutationID %in% c("C>T","G>A"), "C>T | G>A",
           ifelse(MutationID %in% c("T>A","A>T"), "T>A | A>T",
           ifelse(MutationID %in% c("T>C","A>G"), "T>C | A>G",
           ifelse(MutationID %in% c("T>G","A>C"), "T>G | A>C",
                  "non-SNP")))))))
  df$TRANSCRIPT_ID <- sub("^(.*)[.].*", "\\1", df$TRANSCRIPT_ID) 
  df
})


Annotated_VCFs_Gene_Body_Variants <- lapply(
  Annotated_VCFs_Raw,
  function(df) df %>% filter(!EFFECT %in% c("upstream_gene_variant",
                                 "downstream_gene_variant",
                                 "intergenic_region")))           

for (element_name in unique(names(Annotated_VCFs_Gene_Body_Variants))) {

  # Step One: Merge VCF Element with Correct GTF Columns
  # Isolate the element
  element <- Annotated_VCFs_Gene_Body_Variants[[element_name]] %>% as.data.frame()
  # Format the gtf - Calculate Sequence Length
  element_gtf <- Species_Canonical_Transcripts_gtf[[element_name]] 
  element_gtf$cDNA_length <- width(element_gtf)
  # Remove version number if present
  element_gtf$transcript_id <- sub("^(.*)[.].*", "\\1", element_gtf$transcript_id) 
  # Annotate with some information from gtf
  gtf_info <- data.frame(gene_id = element_gtf$gene_id,
                         transcript_id = element_gtf$transcript_id,
                         gene_name = element_gtf$gene_name,
                         gene_biotype = element_gtf$gene_biotype,
                         seqnames = seqnames(element_gtf),
                         start = start(element_gtf),
                         end = end(element_gtf),
                         cDNA_length = element_gtf$cDNA_length,
                         strand = strand(element_gtf))
  
  # Step Two: Merge with dNdS Estimates but including ALL genes
  # First need to isolate the species name, because this is how to pull out the dNdS element
  species <- str_split(element_name, "_", simplify = TRUE)[,4] %>% str_replace(., "^(.)(.*)$", "\\1_\\2")
  
  print(element_name)
  print(class(gtf_info))
  print(names(gtf_info))
  print(head(gtf_info))
  
  gtf_info <- merge(
                    gtf_info,
                    Species_dNdS_Ratios[[species]][,c("query_id","dNdS")],
                    by.x = "transcript_id",
                    by.y = "query_id",
                    all = TRUE)
  
  # Step Three: Merge GTF with Seq Length to Element
  element_info <- merge(
  element,
  gtf_info,
  by.x = c("GENE_ID","CHROM"),
  by.y = c("gene_id","seqnames"),
  all.x = FALSE,
  all.y = TRUE)

  element_info <- dplyr::select(
  element_info,
  -any_of(c("BIOTYPE","GENE_NAME","TRANSCRIPT_NAME","TRANSCRIPT_ID")))
  
  # Step Four: Replace elements in Annotated_VCFs_Raw with the updated one
  Annotated_VCFs_Gene_Body_Variants[[element_name]] <- subset(element_info, is.na(transcript_id) == FALSE)
  
  }



for (element_name in names(Annotated_VCFs_Gene_Body_Variants)) {

  # Step One: Isolate Components of Path to Load Reference Build(s)
  element <- Annotated_VCFs_Gene_Body_Variants[[element_name]] # Get the element data
  # Add Empty Columns to be filled with nucleotide compositions
  element$GC_Length <- ""
  element$AT_Length <- ""
  build <- str_split(element_name, "_", simplify = TRUE)[,5] # Isolate just the build (is a sub-directory in species folder)
  species <- str_split(element_name, "_", simplify = TRUE)[,4] %>% str_replace(., "^(.)(.*)$", "\\1_\\2")
  speciesANDbuild <- paste(species, build, sep = "_") 
  # Step Two: Calculate Nucleotide Composition per Transcript
  DNA_fasta <-
  Biostrings::readDNAStringSet(paste0("/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/", species, "/", build, "/", "sequences.fa"),
                             format = "fasta") 
  # Get Chromosome Names
  names(DNA_fasta) <- str_split(names(DNA_fasta), pattern = " ", simplify = TRUE)[,1]
  
      for (row in 1:nrow(element)) {
      #print(paste("Current Transcript is", i)) # This is i going forward
      # Retrieving Sequence Ranges
      iChrom <- element[row,"CHROM"]
      iStart <- element[row,"start"]
      iEnd <- element[row, "end"]
      # Next Thing
      element[row, "GC_Length"] <- 
      subseq(DNA_fasta[names(DNA_fasta) == iChrom, ], 
                start = as.integer(iStart), end = as.integer(iEnd)) %>% letterFrequency("GC", as.prob = FALSE) %>% as.numeric()
      element[row, "AT_Length"] <- 
      subseq(DNA_fasta[names(DNA_fasta) == iChrom, ], 
                start = as.integer(iStart), end = as.integer(iEnd)) %>% letterFrequency("AT", as.prob = FALSE) %>% as.numeric()
      }
  
  # Replace Element in Annotated_VCFs_Gene_Body_Variants with the updated one
  Annotated_VCFs_Gene_Body_Variants[[element_name]] <- element
  
}

Annotated_VCFs_Gene_Body_Variants <- lapply(
  Annotated_VCFs_Gene_Body_Variants,
  function(df) {
    df %>%
      group_by(GENE_ID) %>%
      mutate(UniqueVariants = n_distinct(Variant_ID[!is.na(Variant_ID)])) %>%
      mutate(UniqueVariantsPlusOne = UniqueVariants+0.5) %>%      
      mutate(UniqueVariantsPlusOnePerBP = UniqueVariantsPlusOne/cDNA_length) %>%
      ungroup()
  }
)


Annotated_VCFs_Gene_Body_Variants <- lapply(
  Annotated_VCFs_Gene_Body_Variants,
  function(df) {
df %>%
  group_by(GENE_ID) %>%
  mutate(
    TtoA_or_AtoT_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "T>A | A>T"]) %>% as.numeric,
        TtoA_or_AtoT_UniqueVariantsPlusOnePerBP = (TtoA_or_AtoT_UniqueVariants+1.0)/as.numeric(AT_Length),
    CtoG_or_GtoC_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "C>G | G>C"]) %>% as.numeric,
        CtoG_or_GtoC_UniqueVariantsPlusOnePerBP = (CtoG_or_GtoC_UniqueVariants+1.0)/as.numeric(GC_Length),
    CtoT_or_GtoA_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "C>T | G>A"]) %>% as.numeric,
        CtoT_or_GtoA_UniqueVariantsPlusOnePerBP = (CtoT_or_GtoA_UniqueVariants+1.0)/as.numeric(GC_Length),
    CtoA_or_GtoT_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "C>A | G>T"]) %>% as.numeric,
        CtoA_or_GtoT_UniqueVariantsPlusOnePerBP = (CtoA_or_GtoT_UniqueVariants+1.0)/as.numeric(GC_Length),
    TtoG_or_AtoC_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "T>G | A>C"]) %>% as.numeric,
        TtoG_or_AtoC_UniqueVariantsPlusOnePerBP = (TtoG_or_AtoC_UniqueVariants+1.0)/as.numeric(AT_Length),
    TtoC_or_AtoG_UniqueVariants = n_distinct(Variant_ID[MutationalCategory %in% "T>C | A>G"]) %>% as.numeric,
        TtoC_or_AtoG_UniqueVariantsPlusOnePerBP = (TtoC_or_AtoG_UniqueVariants+1.0)/as.numeric(AT_Length)) %>%
  ungroup()
  })


# 1 Mutation Count Distribution 
Genewise_Body_Mutation_Counts_Plots <- Map(
  function(df, name) {
    p1 <- ggplot(df, 
                 aes(x = reorder(GENE_ID, UniqueVariants), 
                     y = UniqueVariants, 
                     colour = gene_biotype)) + 
      geom_point(size = 2.8, alpha = 0.9) +
      labs(
        title = name, 
        x = paste(length(unique(df$GENE_ID)), "Unique Canonical Genes"), 
        y = "Gene Body Mutation Count",
        colour = NULL
      ) +
      theme(
        # Transparent backgrounds
        panel.background       = element_rect(fill = "transparent", colour = NA),
        plot.background        = element_rect(fill = "transparent", colour = NA),
        legend.background      = element_rect(fill = "transparent", colour = NA),
        legend.box.background  = element_rect(fill = "transparent", colour = NA),
        
        # Subtle grid
        panel.grid.major.y     = element_line(colour = "grey60", linewidth = 0.3),
        panel.grid.minor.y     = element_line(colour = "grey80", linewidth = 0.2),
        panel.grid.major.x     = element_blank(),
        panel.grid.minor.x     = element_blank(),
        
        # Axes and text styling
        axis.text.x   = element_blank(),
        axis.text.y   = element_text(colour = "white", size = 35),
        axis.title.x  = element_text(colour = "white", size = 14, face = "bold"),
        axis.title.y  = element_text(colour = "white", size = 14, face = "bold"),
        axis.ticks    = element_line(colour = "white"),
        axis.line     = element_line(colour = "white"),
        
        # Title
        plot.title    = element_text(hjust = 0.5, colour = "white", size = 16, face = "bold"),
        
        # Legend
        legend.position = "right",
        legend.text     = element_text(colour = "white", size = 20)
      ) +
      scale_x_discrete(expand = c(0.05, 0.01)) +
      scale_colour_manual(
        values = scales::hue_pal()(length(unique(df$gene_biotype)))
      )

    # --- Save individual plot ---
    n_biotypes <- length(unique(df$gene_biotype))
    total_width <- 30 + min(0.25 * n_biotypes, 20)

    ggsave(
      plot = p1,
      filename = paste0(
        "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
        name, "_Mutation_Count_Distribution_Plot.png"
      ),
      units = "cm",
      height = 10,
      width = total_width,
      dpi = 1500,
      bg = "transparent"
    )

    # --- Return plot for patchwork ---
    p1
  },
  Annotated_VCFs_Gene_Body_Variants,
  names(Annotated_VCFs_Gene_Body_Variants)
)

# --- Combine selected plots vertically ---
#Combined_Counts_Plots <- wrap_plots(Genewise_Body_Mutation_Counts_Plots[
#  c("Assaf_2017_MA_Dmelanogaster_BDGP6",
#    "SaxenaKonrad_2019_MA_Celegans_WS235",
#    "Behringer_2016_MA_Spombe_ASM294v2",
#    "Villalbadelapena_2023_MA_Ncrassa_NC12",
#    "Goldmann_2016_Trio_Hsapiens_hg38")
#], ncol = 1) &
#  theme(
#    plot.background       = element_rect(fill = "transparent", colour = NA),
#    panel.background      = element_rect(fill = "#eeeeee", colour = NA),
#    legend.background     = element_rect(fill = "transparent", colour = NA),
#    legend.box.background = element_rect(fill = "transparent", colour = NA)
#  )

# --- Save the combined figure ---
#ggsave(
#  filename = "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/Combined_Mutation_Count_Distribution_Plots.png",
#  plot = Combined_Counts_Plots,
#  width = 50,
#  height = 10 * length(Genewise_Body_Mutation_Counts_Plots),
#  units = "cm",
#  dpi = 800,
#  bg = "transparent"
#)


#stop("PAUSE HERE")



# 2 Mutation Rate Distribution 
Genewise_Body_Mutation_Rate_Plots <- Map(
  function(df, name) {
  p1 <-
  ggplot(df,
       aes(x = reorder(GENE_ID, UniqueVariantsPlusOnePerBP), y = UniqueVariantsPlusOnePerBP, colour = gene_biotype,
           label = GENE_ID)) + 
  geom_point(fill = "transparent") +
  labs(title = name,
       x = paste(length(unique(df$GENE_ID)),
                 "Unique Canonical Genes"), 
       y = "Gene Body Mutations (+1) per BP",
       colour = NULL) +
  theme(
        # transparent panel and plot backgrounds
        panel.background = element_rect(fill = NA, colour = NA),
        plot.background  = element_rect(fill = NA, colour = NA),
        # light-gray grid lines
        panel.grid.major = element_line(colour = "grey60", linewidth = 0.3),
        panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2),
        # white text, axes, and ticks
        axis.text.x   = element_blank(),
        axis.text.y   = element_text(colour = "white", size = 10),
        axis.title.x  = element_text(colour = "white", size = 14, face = "bold"),
        axis.title.y  = element_text(colour = "white", size = 14, face = "bold"),
        axis.ticks    = element_line(colour = "white"),
        axis.line     = element_line(colour = "white"),
        # title
        plot.title    = element_text(hjust = 0.5, colour = "white", size = 16, face = "bold"),
        # remove legend unless needed
        legend.position = "none") +
        scale_x_discrete(expand = c(0.05, 0.01)) 
    
    # Save the plot — change width/height as needed
    ggsave(plot = p1,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_Mutation_Rate_Distribution_Plot.png"),
      create.dir = TRUE,
      units = "cm",
             height = 10,
             width = 25,
      dpi = 300
    )
  },
  Annotated_VCFs_Gene_Body_Variants,
  names(Annotated_VCFs_Gene_Body_Variants)
)


# 3 Protein Coding Mutation Rate Distribution 
Genewise_Protein_Body_Mutation_Rate_Plots <- Map(
  function(df, name) {
    p1 <- ggplot(subset(df, gene_biotype == "protein_coding"), 
                 aes(
                   x = reorder(GENE_ID, UniqueVariantsPlusOnePerBP),
                   y = UniqueVariantsPlusOnePerBP,
                   colour = UniqueVariantsPlusOnePerBP,
                   label = GENE_ID
                 )) +
      geom_point(size = 4, fill = "transparent") +
      labs(
        title = name,
        x = paste(length(unique(df$GENE_ID)), "Unique Canonical Protein Coding Genes"), 
        y = "Gene Body Mutations (+1) per BP",
        colour = NULL
      ) +
      theme(
        plot.margin = margin(5, 5, 5, 30),
        panel.grid.major.y    = element_line(colour = "grey60", linewidth = 0.3),
        panel.grid.minor.y    = element_line(colour = "grey80", linewidth = 0.2),
        panel.grid.major.x    = element_blank(),
        panel.grid.minor.x    = element_blank(),
        axis.text.x   = element_blank(),
        axis.text.y   = element_text(colour = "white", size = 40),
        axis.title.x  = element_text(colour = "white", size = 14, face = "bold"),
        axis.title.y  = element_text(colour = "white", size = 14, face = "bold"),
        axis.ticks    = element_line(colour = "white"),
        axis.line     = element_line(colour = "white"),
        plot.title    = element_text(hjust = 0.5, colour = "white", size = 16, face = "bold"),
        legend.position = "none",
        legend.text     = element_text(colour = "white", size = 10)
      ) +
      scale_x_discrete(expand = c(0.05, 0.01)) +
      scale_y_continuous(position = "left") +
      scale_colour_gradient(low = "orange", high = "red") +
      theme(
  plot.background       = element_rect(fill = "transparent", colour = NA),
  panel.background      = element_rect(fill = "#eeeeee", colour = NA),
  legend.background     = element_rect(fill = "transparent", colour = NA),
  legend.box.background = element_rect(fill = "transparent", colour = NA))
    
    # --- Save individual plot ---
    ggsave(
      plot = p1,
      filename = paste0(
        "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
        name, "_Protein_Mutation_Count_Distribution_Plot.png"
      ),
      create.dir = TRUE,
      units = "cm",
      height = 10,
      width = 40,
      dpi = 1500,
      bg = "transparent"  # <--- preserve transparency on export
    )

    return(p1)
  },
  Annotated_VCFs_Gene_Body_Variants,
  names(Annotated_VCFs_Gene_Body_Variants)
)

# --- Combine all plots vertically ---
#Combined_Protein_Plots <- wrap_plots(Genewise_Protein_Body_Mutation_Rate_Plots[
#  c("Assaf_2017_MA_Dmelanogaster_BDGP6",
#    "SaxenaKonrad_2019_MA_Celegans_WS235",
#    "Behringer_2016_MA_Spombe_ASM294v2",
#    "Villalbadelapena_2023_MA_Ncrassa_NC12",
#    "Goldmann_2016_Trio_Hsapiens_hg38")], ncol = 1) &
#  theme(
#  plot.background       = element_rect(fill = "transparent", colour = NA),
#  panel.background      = element_rect(fill = "#eeeeee", colour = NA),
#  legend.background     = element_rect(fill = "transparent", colour = NA),
#  legend.box.background = element_rect(fill = "transparent", colour = NA))

# --- Save the combined figure ---
#ggsave(
#  filename = "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/Combined_Protein_Mutation_Count_Distribution_Plots.png",
#  plot = Combined_Protein_Plots,
#  width = 55,
#  height = 10 * length(Genewise_Protein_Body_Mutation_Rate_Plots),
#  units = "cm",
#  dpi = 900,
#  bg = "transparent")  # <--- critical for transparent composite




# 4 Protein dNdS Distribution
Genewise_Protein_dNdS_Plots <- Map(
  function(df, name) {
    p1 <-
    ggplot(
      subset(df, gene_biotype == "protein_coding" & !is.na(dNdS)),
      aes(x = reorder(GENE_ID, dNdS), y = dNdS, colour = dNdS)
    ) +
      geom_point(size = 5) +
      labs(
        title = name,
        x = paste(length(unique(df$GENE_ID)), "Unique Orthologs"),
        y = "dNdS Ratio"
      ) +
       theme(
        plot.margin = margin(5, 5, 5, 30),
        panel.grid.major.y    = element_line(colour = "grey60", linewidth = 0.3),
        panel.grid.minor.y    = element_line(colour = "grey80", linewidth = 0.2),
        panel.grid.major.x    = element_blank(),
        panel.grid.minor.x    = element_blank(),
        axis.text.x   = element_blank(),
        axis.text.y   = element_text(colour = "white", size = 15),
        axis.title.x  = element_text(colour = "white", size = 14, face = "bold"),
        axis.title.y  = element_text(colour = "white", size = 14, face = "bold"),
        axis.ticks    = element_line(colour = "white"),
        axis.line     = element_line(colour = "white"),
        plot.title    = element_text(hjust = 0.5, colour = "white", size = 16, face = "bold"),
        legend.position = "none",
        legend.text     = element_text(colour = "white", size = 10)
      ) +
      scale_x_discrete(expand = c(0.05, 0.01)) +
      scale_y_continuous(position = "left", n.breaks = 8) +
      scale_colour_gradient(low = "orange", high = "#00ff00") +
      theme(
  plot.background       = element_rect(fill = "transparent", colour = NA),
  panel.background      = element_rect(fill = "#eeeeee", colour = NA),
  legend.background     = element_rect(fill = "transparent", colour = NA),
  legend.box.background = element_rect(fill = "transparent", colour = NA))
    
 
# Save the plot with transparent background
ggsave(
  plot = p1,
  filename = paste0("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                    name, "_Protein_dNdS_Distribution_Plot.png"),
  bg = "transparent",   # ensures transparency in saved image
  units = "cm",
             height = 10,
             width = 16,
  dpi = 2500
)
  },
  Annotated_VCFs_Gene_Body_Variants,
  names(Annotated_VCFs_Gene_Body_Variants)
)


Annotated_VCFs_PC_Gene_Body_Variants <- lapply(
  Annotated_VCFs_Gene_Body_Variants,
  function(df) {
    df <- df %>%
      subset(gene_biotype == "protein_coding")
})


# 5 Mutation Rate vs dNdS Associations: Linear and Non-Linear
Species_dNdS_MutationRate_Spearman_Tests <- list()
Species_dNdS_MutationRate_Pearson_Tests <- list()

Genewise_MutationRate_by_dNdS_Plots <- Map(
  function(df, name) {
    
  # Run Spearman Correlations
  spearman_test <- 
    cor.test(log(df$dNdS),
             df$UniqueVariantsPlusOnePerBP,
             method = "spearman")
  # Run Pearson Correlations
  pearson_test <- 
    cor.test(log(df$dNdS),
             df$UniqueVariantsPlusOnePerBP,
             method = "pearson")
  
  # Save results into object
  Species_dNdS_MutationRate_Spearman_Tests[[name]] <<- spearman_test
  Species_dNdS_MutationRate_Pearson_Tests[[name]] <<- pearson_test
  
  # Create a Scatter Plot
  p1 <- 
  ggplot(
  data = subset(df, !is.na(dNdS)),
  aes(x = log(dNdS), y = UniqueVariantsPlusOnePerBP)
) +
  geom_point() +
  geom_smooth(
    method = "lm", 
    se = TRUE,
    colour = "#66ff00",                # regression line color
    fill = alpha("orange", 0.5)    # 50% transparent orange confidence band
  ) +
  geom_rug(alpha = 0.5, color = "orange") +          # axis rugs/carpets
  labs(
    title = paste(
      name, "\n",
      subset(df, !is.na(dNdS) & !is.na(UniqueVariantsPlusOnePerBP))$GENE_ID %>%
        unique() %>%
        length(),
      "Unique Protein Coding Genes with dNdS Estimate"
    ),
    subtitle = paste0(
      "Spearman Correlation: ",
      "p Value = ", signif(Species_dNdS_MutationRate_Spearman_Tests[[name]]$p.value, 5), ", ",
      "Coefficient = ", signif(Species_dNdS_MutationRate_Spearman_Tests[[name]]$estimate, 5),
      "\n",
      "Pearson Correlation: ",
      "p Value = ", signif(Species_dNdS_MutationRate_Pearson_Tests[[name]]$p.value, 5), ", ",
      "Coefficient = ", signif(Species_dNdS_MutationRate_Pearson_Tests[[name]]$estimate, 5)
    ),
    caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, colour = "white", hjust = .5),
    text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.subtitle = element_text(size = 8, color = "white", hjust = .5, face = "italic"),
    plot.caption = element_text(size = 6, color = "white", hjust = .5, face = "italic"),
    strip.text = element_text(color = "white"),
    panel.background = element_rect(fill = "#eeeeee", colour = "white"),  # transparent panel w/ black border
    plot.background = element_rect(fill = NA, colour = NA)         # transparent overall background
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10)

# Save the plot — transparent background, keep black border
ggsave(
  plot = p1,
  filename = paste0(
    "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
    name, "_dNdS_vs_Mutation_Rate_Plot.png"
  ),
  bg = "transparent",   # ensures background transparency in PNG
  width = 5,
  height = 4,
  dpi = 2500
)
    
    
    # Calculate some statistics beforehand
    # --- Fit GAM (nonparametric smooth) ---
    gam_fit <- gam(UniqueVariantsPlusOnePerBP ~ s(dNdS, k = 6), data = df)

    # --- Fit linear model first, then segmented ---
    lm_fit <- lm(UniqueVariantsPlusOnePerBP ~ dNdS, data = df)
    seg_fit <- segmented(lm_fit, seg.Z = ~ dNdS, npsi = 2)  # try 2 breakpoints

    # Create a prediction grid
    newdata <- data.frame(dNdS = seq(min(df$dNdS, na.rm = TRUE), max(df$dNdS, na.rm = TRUE), length.out = 200))

    # Predictions from both models
    newdata$GAM_fit <- predict(gam_fit, newdata, type = "response")
    newdata$SEG_fit <- predict(seg_fit, newdata, type = "response")
    breaks <- seg_fit$psi[, "Est."]

    p2 <-
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_line(data = newdata, aes(y = GAM_fit), color = "blue", size = 1.2) +
    geom_line(data = newdata, aes(y = SEG_fit), color = "orange", linetype = "dashed", size = 1) +
    geom_vline(xintercept = breaks, linetype = "dotted", color = "orange") +
    annotate("text", x = breaks, y = max(subset(df, is.na(dNdS) == FALSE)$UniqueVariantsPlusOnePerBP, na.rm = TRUE) * 0.5,
           label = paste0("Breakpoint\n", round(breaks, 3)),
           color = "orange", hjust = -0.1, size = 3.2) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = "Red = GAM smooth; Blue dashed = segmented regression",
      caption = "Blue line represents fit of generalised linear model: gam(UniqueVariantsPlusOnePerBP ~ s(dNdS, k = 6), data = df). \n
                Orange line represents fit of segmented linear regression: seg_fit <- segmented(lm_fit, seg.Z = ~ dNdS, npsi = 2)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "white", hjust = .5),
          text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.text = element_text(color = "white"),
          plot.subtitle = element_text(size = 8, color = "white", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 6, color = "white", hjust = .5, face = "italic"),
          strip.text = element_text(color = "white")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)

    # Save the plot — change width/height as needed
      ggsave(plot = p2,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_GAM_SegReg_dNdS_vs_Mutation_Rate_Plot.png"),
      create.dir = TRUE,
      width = 6,
      height = 6,
      dpi = 300)
    
  },
  Annotated_VCFs_PC_Gene_Body_Variants,
  names(Annotated_VCFs_PC_Gene_Body_Variants)
)


# 5 Mutation Rate vs dNdS Association of Genes with at least one MA Variant
Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly <- list()
Species_dNdS_MutationRate_Pearson_TestsMutatedOnly <- list()

MutatedOnly_Genewise_MutationRate_by_dNdS_Plots <- Map(
  function(df, name) {
  # Run Spearman Correlations
  spearman_test <- 
    cor.test(log(subset(df, UniqueVariants > 0)$dNdS),
             subset(df, UniqueVariants > 0)$UniqueVariantsPlusOnePerBP,
             method = "spearman")
  # Run Pearson Correlations
  pearson_test <- 
    cor.test(log(subset(df, UniqueVariants > 0)$dNdS),
             subset(df, UniqueVariants > 0)$UniqueVariantsPlusOnePerBP,
             method = "pearson")
  
  # Save results into object
  Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly[[name]] <<- spearman_test
  Species_dNdS_MutationRate_Pearson_TestsMutatedOnly[[name]] <<- pearson_test
  
  # Create a Scatter Plot
  p1 <- 
  ggplot(
  data = subset(df, !is.na(dNdS)),
  aes(x = log(dNdS), y = UniqueVariantsPlusOnePerBP)
) +
  geom_point() +
  geom_smooth(
    method = "lm", 
    se = TRUE,
    colour = "#66ff00",                # regression line color
    fill = alpha("#ff2600", 0.5)    # 50% transparent orange confidence band
  ) +
  geom_rug(alpha = 0.5, color = "#ff2600") +          # axis rugs/carpets
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(UniqueVariantsPlusOnePerBP) == FALSE &
                              UniqueVariants > 0)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate and >0 Variants"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_MutationRate_Pearson_TestsMutatedOnly[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_MutationRate_Pearson_TestsMutatedOnly[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(
    plot.title = element_text(size = 10, colour = "white", hjust = .5),
    text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.subtitle = element_text(size = 8, color = "white", hjust = .5, face = "italic"),
    plot.caption = element_text(size = 6, color = "white", hjust = .5, face = "italic"),
    strip.text = element_text(color = "white"),
    panel.background = element_rect(fill = "#eeeeee", colour = "white"),  # transparent panel w/ black border
    plot.background = element_rect(fill = NA, colour = NA)         # transparent overall background
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p1,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_Mutation_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 5,
      height = 4,
      dpi = 2500)

  },
  Annotated_VCFs_PC_Gene_Body_Variants,
  names(Annotated_VCFs_PC_Gene_Body_Variants)
)


# Block Loci into Equal Size Bins by dNdS
Annotated_VCFs_PC_Gene_Body_Variants <- lapply(
  Annotated_VCFs_PC_Gene_Body_Variants,
  function(df) {
    df %>%
      mutate(Overall_dNdS_Group = cut_number(dNdS + runif(length(dNdS), -1e-6, 1e-6), n = 5)) %>%
      mutate(Overall_dNdS_Group = factor(Overall_dNdS_Group, ordered = TRUE)) %>%
      mutate(Overall_dNdS_Group_More = cut_number(dNdS + runif(length(dNdS), -1e-6, 1e-6), n = 15)) %>%
      mutate(Overall_dNdS_Group_More = factor(Overall_dNdS_Group_More, ordered = TRUE))
  }
)



# 1 Mutation Rate vs dNdS Association 
Genewise_Binned_MutationRate_by_dNdS_Plots <- Map(
  function(df, name) {

# Subset to only loci with dNdS 
df <- subset(df, is.na(dNdS) == FALSE) 

# Plot Association of dNdS Bins with MutationsPlusOnePerBP 
p1 <-
  ggplot(
    data = subset(df, !is.na(Overall_dNdS_Group)),
    aes(
      x = Overall_dNdS_Group,
      y = UniqueVariantsPlusOnePerBP,
      fill = Overall_dNdS_Group
    )
  ) +
  theme_bw() +
  geom_boxplot(outliers = FALSE, colour = "black") +  # keep box borders black
  labs(
    title = paste(
      name, "\n",
      subset(df, !is.na(Overall_dNdS_Group) & !is.na(UniqueVariantsPlusOnePerBP))$GENE_ID %>%
        unique() %>%
        length(),
      "Unique Protein Coding Genes with dNdS Estimate"
    ),
    subtitle = paste(
      "Jonckheere =",
      jonckheere.test(
        df$UniqueVariantsPlusOne,
        df$Overall_dNdS_Group,
        alternative = "two.sided"
      )$statistic %>% round(3),
      "p =",
      jonckheere.test(
        df$UniqueVariantsPlusOne,
        df$Overall_dNdS_Group,
        alternative = "two.sided"
      )$p.value %>% round(3)
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 6, colour = "white", hjust = .5),
    text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text.x = element_text(size = 5, color = "white", vjust = .6, angle = 15),
    axis.text.y = element_text(size = 6, color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.subtitle = element_text(size = 7, color = "white", hjust = .5, face = "italic"),
    plot.caption = element_text(size = 6, color = "white", hjust = .5, face = "italic"),
    strip.text = element_text(color = "white"),
    panel.background = element_rect(fill = "#eeeeee", colour = "white"),  # transparent panel, black border
    plot.background = element_rect(fill = NA, colour = NA)         # transparent overall background
  ) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(values = c("#FFB703", "#FB8500", "#E3670E", "#D63B09", "#B21801")) +
  guides(fill = "none")

# Save the plot — transparent background, keep borders
ggsave(
  plot = p1,
  filename = paste0(
    "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
    name, "_Binned_dNdS_vs_Mutation_Rate_Plot.png"
  ),
  bg = "transparent",   # ensure transparent PNG background
  width = 2.75,
  height = 4,
  dpi = 2500
)


# Plot Association of dNdS Bins with MutationsPlusOnePerBP 
p2 <-
  ggplot(
    data = subset(df, !is.na(Overall_dNdS_Group)),
    aes(
      x = Overall_dNdS_Group_More,
      y = UniqueVariantsPlusOnePerBP,
      fill = Overall_dNdS_Group_More
    )
  ) +
  theme_bw() +
  geom_boxplot(outliers = FALSE, colour = "black") +  # keep box borders black
  labs(
    title = paste(
      name, "\n",
      subset(df, !is.na(Overall_dNdS_Group) & !is.na(UniqueVariantsPlusOnePerBP))$GENE_ID %>%
        unique() %>%
        length(),
      "Unique Protein Coding Genes with dNdS Estimate"
    ),
    subtitle = paste(
      "Jonckheere =",
      jonckheere.test(
        df$UniqueVariantsPlusOne,
        df$Overall_dNdS_Group_More,
        alternative = "two.sided"
      )$statistic %>% round(3),
      "p =",
      jonckheere.test(
        df$UniqueVariantsPlusOne,
        df$Overall_dNdS_Group_More,
        alternative = "two.sided"
      )$p.value %>% round(3)
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 6, colour = "white", hjust = .5),
    text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text.x = element_text(size = 4, color = "white", vjust = .6, angle = 15),
    axis.text.y = element_text(size = 6, color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.subtitle = element_text(size = 7, color = "white", hjust = .5, face = "italic"),
    plot.caption = element_text(size = 6, color = "white", hjust = .5, face = "italic"),
    strip.text = element_text(color = "white"),
    panel.background = element_rect(fill = "#eeeeee", colour = "white"),  # transparent panel, black border
    plot.background = element_rect(fill = NA, colour = NA)         # transparent overall background
  ) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(values = colorRampPalette(c("#FFB703", "#B21801"))(length(unique(df$Overall_dNdS_Group_More)))) +
  guides(fill = "none")

# Save the plot — transparent background, keep borders
ggsave(
  plot = p2,
  filename = paste0(
    "/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
    name, "_More_Binned_dNdS_vs_Mutation_Rate_Plot.png"
  ),
  bg = "transparent",   # ensure transparent PNG background
  width = 2.75,
  height = 4,
  dpi = 2500
)
      
  },
  Annotated_VCFs_PC_Gene_Body_Variants,
  names(Annotated_VCFs_PC_Gene_Body_Variants)
)

# 6 Mutation Rate vs dNdS Associations Facetted by SNP Category
# "T>A | A>T"
Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests <- list()
Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests <- list()
# "C>G | G>C"
Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests <- list()
Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests <- list()
# "C>T | G>A"
Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests <- list()
Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests <- list()
# "C>A | G>T"
Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests <- list()
Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests <- list()
# "T>G | A>C"
Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests <- list()
Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests <- list()
# "T>C | A>G"
Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests <- list()
Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests <- list()



Genewise_MutationRate_by_dNdS_by_Category_Plots <- Map(

  function(df, name) {
  # Run Spearman Correlations
  spearman_test_TtoA_or_AtoT <- 
    cor.test(df$dNdS,
             df$TtoA_or_AtoT_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  spearman_test_CtoG_or_GtoC <- 
    cor.test(df$dNdS,
             df$CtoG_or_GtoC_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  spearman_test_CtoT_or_GtoA <- 
    cor.test(df$dNdS,
             df$CtoT_or_GtoA_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  spearman_test_CtoA_or_GtoT <- 
    cor.test(df$dNdS,
             df$CtoA_or_GtoT_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  spearman_test_TtoG_or_AtoC <- 
    cor.test(df$dNdS,
             df$TtoG_or_AtoC_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  spearman_test_TtoC_or_AtoG <- 
    cor.test(df$dNdS,
             df$TtoC_or_AtoG_UniqueVariantsPlusOnePerBP,
             method = "spearman")
  
  # Run Pearson Correlations
  pearson_test_TtoA_or_AtoT <- 
    cor.test(df$dNdS,
             df$TtoA_or_AtoT_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  pearson_test_CtoG_or_GtoC <- 
    cor.test(df$dNdS,
             df$CtoG_or_GtoC_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  pearson_test_CtoT_or_GtoA <- 
    cor.test(df$dNdS,
             df$CtoT_or_GtoA_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  pearson_test_CtoA_or_GtoT <- 
    cor.test(df$dNdS,
             df$CtoA_or_GtoT_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  pearson_test_TtoG_or_AtoC <- 
    cor.test(df$dNdS,
             df$TtoG_or_AtoC_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  pearson_test_TtoC_or_AtoG <- 
    cor.test(df$dNdS,
             df$TtoC_or_AtoG_UniqueVariantsPlusOnePerBP,
             method = "pearson")
  
  # Save results into object
  Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests[[name]] <<- spearman_test_TtoA_or_AtoT
  Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests[[name]] <<- spearman_test_CtoG_or_GtoC
  Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests[[name]] <<- spearman_test_CtoT_or_GtoA
  Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests[[name]] <<- spearman_test_CtoA_or_GtoT
  Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests[[name]] <<- spearman_test_TtoG_or_AtoC
  Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests[[name]] <<- spearman_test_TtoC_or_AtoG
  #
  Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests[[name]] <<- pearson_test_TtoA_or_AtoT
  Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests[[name]] <<- pearson_test_CtoG_or_GtoC
  Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests[[name]] <<- pearson_test_CtoT_or_GtoA
  Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests[[name]] <<- pearson_test_CtoA_or_GtoT
  Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests[[name]] <<- pearson_test_TtoG_or_AtoC
  Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests[[name]] <<- pearson_test_TtoC_or_AtoG
  
  # Create a TtoA_or_AtoT Scatter Plot
  p_TtoA_or_AtoT <- 
  ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = TtoA_or_AtoT_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(TtoA_or_AtoT_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_TtoA_or_AtoT,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_TtoA_or_AtoT_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)
      
      
        
  # Create a CtoG_or_GtoC Scatter Plot
  p_CtoG_or_GtoC <- 
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = CtoG_or_GtoC_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(CtoG_or_GtoC_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_CtoG_or_GtoC,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_CtoG_or_GtoC_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)
      
      
      
      # Create a CtoT_or_GtoA Scatter Plot
    p_CtoT_or_GtoA <- 
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = CtoT_or_GtoA_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(CtoT_or_GtoA_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_CtoT_or_GtoA,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_CtoT_or_GtoA_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)
      
      
      
    # Create a CtoA_or_GtoT Scatter Plot
    p_CtoA_or_GtoT <- 
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = CtoA_or_GtoT_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(CtoA_or_GtoT_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_CtoA_or_GtoT,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_CtoA_or_GtoT_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)
      
      
      
    # Create a TtoG_or_AtoC Scatter Plot
    p_TtoG_or_AtoC <- 
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = TtoG_or_AtoC_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(TtoG_or_AtoC_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_TtoG_or_AtoC,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_TtoG_or_AtoC_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)
      
      
      
    # Create a TtoC_or_AtoG Scatter Plot
    p_TtoC_or_AtoG <- 
    ggplot2::ggplot(data = subset(df,
                 is.na(dNdS) == FALSE),
                    aes(x = dNdS, y = TtoC_or_AtoG_UniqueVariantsPlusOnePerBP)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste(name, "\n",
                       subset(df, is.na(dNdS) == FALSE & is.na(TtoC_or_AtoG_UniqueVariantsPlusOnePerBP) == FALSE)$GENE_ID %>% 
                       unique() %>% length(), "Unique Protein Coding Genes with dNdS Estimate"),
      subtitle = paste0("Spearman Correlation: ", 
                       "p Value = ", Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests[[name]]$estimate %>%
                         signif(., digits = 5),
                       "\n",
                       "Pearson Correlation: ", 
                       "p Value = ", Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests[[name]]$p.value %>% 
                         signif(., digits = 5),", ",
                       "Coefficient = ", Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests[[name]]$estimate %>% 
                         signif(., digits = 5)),
      caption = "Line of best fit represents linear regression line produced using: geom_smooth(method = 'lm', se = TRUE)") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, colour = "black", hjust = .5),
          text = element_text(color = "black"),
          axis.title.x = element_text(size = 10, color = "black", vjust = .6),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black"),
          plot.subtitle = element_text(size = 8, color = "black", hjust = .5, face = "italic"),
          plot.caption = element_text(size = 5, color = "black", hjust = .5, face = "italic"),
          strip.text = element_text(color = "black")) +
          scale_x_continuous(n.breaks = 10) +
          scale_y_continuous(n.breaks = 10)
  
    # Save the plot — change width/height as needed
      ggsave(plot = p_TtoC_or_AtoG,
      filename = paste("/home/ocdm0351/DPhil/R_Data/Plots/3_Mutation_and_dNdS_Association_Plots/",
                       name, "_dNdS_vs_TtoC_or_AtoG_Rate_MutatedOnly_Plot.png"),
      create.dir = TRUE,
      width = 4,
      height = 4,
      dpi = 300)

  },
  Annotated_VCFs_PC_Gene_Body_Variants,
  names(Annotated_VCFs_PC_Gene_Body_Variants)
)



# Saving Annotated Gene Body Variants
save(Annotated_VCFs_Gene_Body_Variants, 
     file = "/home/ocdm0351/DPhil/R_Data/Annotated_VCFs_Gene_Body_Variants")
save(Annotated_VCFs_PC_Gene_Body_Variants, 
     file = "/home/ocdm0351/DPhil/R_Data/Annotated_VCFs_PC_Gene_Body_Variants")

# Saving dNdS Mutation Rate Testing Statistics
save(Species_dNdS_MutationRate_Spearman_Tests, 
     file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_MutationRate_Spearman_Tests")
save(Species_dNdS_MutationRate_Pearson_Tests, 
     file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_MutationRate_Pearson_Tests")
save(Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly, 
     file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_MutationRate_Spearman_Tests_MutatedOnly")
save(Species_dNdS_MutationRate_Pearson_TestsMutatedOnly, 
     file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_MutationRate_Pearson_TestsMutatedOnly")

# Saving dNdS Mutation Category Association Testing Statistics
# "T>A | A>T"
save(Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoA_or_AtoT_Rate_Spearman_Tests") 
save(Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoA_or_AtoT_Rate_Pearson_Tests")
# "C>G | G>C"
save(Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoG_or_GtoC_Rate_Spearman_Tests")
save(Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoG_or_GtoC_Rate_Pearson_Tests")
# "C>T | G>A"
save(Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoT_or_GtoA_Rate_Spearman_Tests")
save(Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests , file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoT_or_GtoA_Rate_Pearson_Tests")
# "C>A | G>T"
save(Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoA_or_GtoT_Rate_Spearman_Tests") 
save(Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_CtoA_or_GtoT_Rate_Pearson_Tests") 
# "T>G | A>C"
save(Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoG_or_AtoC_Rate_Spearman_Tests") 
save(Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoG_or_AtoC_Rate_Pearson_Tests")
# "T>C | A>G"
save(Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoC_or_AtoG_Rate_Spearman_Tests") 
save(Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests, file = "/home/ocdm0351/DPhil/R_Data/Species_dNdS_TtoC_or_AtoG_Rate_Pearson_Tests") 
