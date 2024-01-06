### Main figures Bbsl manuscript
# December 4 2023
# lemieux@broadinstitute.org

# requires local path as borrelialr/notebooks

## to do
# alignemnt to ref, plasmid vs chromosome
# association with genus / species

# Read in libraries ####
library(tidyverse)
library(phytools)
library(ggtree)
library(ggpubr)
library(cowplot)
source("../scripts/Heatmap_PAM.R")
library(Biostrings)
library(treeio)
library(ggtext)

# Read in data ####
# read in trees and midpoint root


rundate <- "12_28_2023"

x <- read.newick(paste("../data/", rundate, "/core_alignment.nwk", sep=""))
x <- midpoint.root(x)
ggtree(x)

a <- read.newick(paste("../data/", rundate, "/accessory_binary_genes.fa.newick", sep=""))
a <- midpoint.root(a)
ggtree(a)

# read in metadata
metadata <- read_tsv("../data/borreliaceae_assembly_info.tsv") 
# need to correct some species-level identifications
metadata$`Organism Name`[metadata$`Assembly Name` == "ASM307659v1"] = "Borrelia salvatore" 
metadata$`Organism Name`[metadata$`Assembly Name` == "ASM2116577v2"] = "Borrelia mahuryensis" 
metadata$`Organism Name`[metadata$`Assembly Name` == "ASM1408402v1"] = "Borreliella maritima"
metadata$`Organism Name`[metadata$`Assembly Name` == "ASM3044027v1"] = "Borreliella maritima"
metadata$`Organism Name`[metadata$`Assembly Name` == "ASM3043611v1"] = "Borrelia rubricentralis"

metadata <- metadata %>% 
  mutate(assem_id = paste(`Assembly Accession`, "_",`Assembly Name`,"_", "genomic.fna", sep="")) %>% 
  mutate(species = str_split_fixed(`Organism Name`, " ", Inf)[,2]) %>% 
  mutate(genus = str_split_fixed(`Organism Name`, " ",Inf)[,1])
metadata <- cbind(gsub(" ", "_", metadata$assem_id), metadata)
names(metadata)[1] <- "label"

Bb_annotation <- read_csv("../data/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-02-14.csv", na = c("", "NA", "ND")) %>% 
  mutate(RST_Type = factor(RST_Type)) %>%
  mutate(Region = recode(Location, Nantucket = "US Northeast", NY = "US Northeast",
                         RI = "US Northeast", CT = "US Northeast", WI = "US Midwest",
                         Slovenia = "EU Slovenia"
  )) %>% 
  filter(!is.na(Original_Name)) %>%
  select(c(Rename_A, Original_Name, OspC_Type, RST_Type, Location))
names(Bb_annotation)[1] <- "Organism Infraspecific Names Strain"

metadata <- left_join(metadata, Bb_annotation)

# Create basic trees ####
p <- ggtree(x)

# read in core genome alignment
aln <- readDNAStringSet(paste("../data/", rundate, "/core_gene_alignment.aln", sep=""))
set.seed(11111)
# generate subsampled tree ####
metadata_tree <- metadata %>% 
  filter(label %in% names(aln)) %>%
  group_by(species) %>%
  filter(species != "burgdorferi" | species =="burgdorferi" & !is.na(OspC_Type)) %>%
  slice_sample(n = 25) %>% 
  distinct(genus, species, `Organism Infraspecific Names Strain`, .keep_all=TRUE)
aln_subset <- aln[metadata_tree$label]
writeXStringSet(aln_subset, paste("../data/", rundate, "/core_gene_alignment_subset.aln", sep=""))
system(paste("~/bin/fasttree/FastTree -nt -gtr ../data/", rundate, "/core_gene_alignment_subset.aln > ../data/", rundate, "/core_gene_alignment_subset_tree.nwk", sep=""))
y <- read.newick(paste("../data/", rundate, "/core_gene_alignment_subset_tree.nwk", sep=""))
y <- midpoint.root(y)

# generate B burgdorferi tree ####

metadata_bb <- metadata %>% 
  filter(label %in% names(aln)) %>%
  filter(species == "burgdorferi")

aln_subsetbb <- aln[metadata_bb$label]
writeXStringSet(aln_subsetbb, paste("../data/", rundate, "/core_gene_alignment_subsetbb.aln", sep=""))
system(paste("~/bin/fasttree/FastTree -nt -gtr ../data/", rundate, "/core_gene_alignment_subsetbb.aln > ../data/", rundate, "/core_gene_alignment_subsetbb_tree.nwk", sep=""))
system(paste("~/bin/iqtree-2.2.0-MacOSX/bin/iqtree2 -s ../data/", rundate, "/core_gene_alignment_subset.aln -B 1000 -T 10 --redo", sep=""))

z <- read.newick(paste("../data/", rundate, "/core_gene_alignment_subsetbb_tree.nwk", sep=""))
z <- midpoint.root(z)

iq_tree <- read.iqtree(paste("../data/", rundate, "/core_gene_alignment_subset.aln.contree", sep=""))
iq_tree@phylo <- midpoint.root(iq_tree@phylo)

# make B miyamotoi plot #### 

metadata_bm <- metadata %>% 
  filter(label %in% names(aln)) %>%
  filter(species == "miyamotoi")

aln_subsetbm <- aln[metadata_bm$label]
writeXStringSet(aln_subsetbm, paste("../data/", rundate, "/core_gene_alignment_subsetbm.aln", sep=""))
system(paste("~/bin/fasttree/FastTree -nt -gtr ../data/", rundate, "/core_gene_alignment_subsetbm.aln > ../data/", rundate, "/core_gene_alignment_subsetbm_tree.nwk", sep=""))
xx <- read.newick(paste("../data/", rundate, "/core_gene_alignment_subsetbm_tree.nwk", sep=""))
xx <- midpoint.root(xx)

metadata_treeplot <- metadata %>% 
  filter(label %in% y$tip.label)

# Make a variety of trees

# IQTree ####
ggtree(iq_tree, layout="fan", open.angle=5) %<+% 
  metadata_treeplot +
  #geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) + 
  geom_tiplab(aes(label = paste(species, `Organism Infraspecific Names Strain`,sep=" "), color = species), size = 1.6) + 
  geom_nodepoint(aes(label = round(UFboot,2), subset = UFboot > 90, vjust = 3), 
               ,alpha = 0.5, color = "blue", size = 1) + 
  geom_treescale(x = 0.2, y = 190, offset = 4, width = 0.1)+
  theme(legend.position = "none")

ggsave(paste("../data/", rundate, "/figures/Figure_1.jpg", sep=""), height =6, width = 6)


species_tree <- p %<+% metadata + 
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 3)
species_tree

metadata_filt <- metadata %>% filter(species %in% c("afzelii", "burgdorferi", "bavariensis", "garinii", "miyamotoi", 
                                                    "valaisiana", "crocidurae", "duttonii", "hermsii", "recurrentis", "turicatae")) 
species_tree_sub <- ggtree(y) %<+% metadata_filt + 
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 3) + 
  geom_treescale(x = 0.05, y = 190, offset = 4, width = 0.1)
species_tree_sub

species_tree_accessory <- ggtree(a) %<+% metadata + 
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 3)
species_tree_accessory


ggtree(x) %<+% metadata_filt +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4)


ggtree(y, layout = "circular") %<+% metadata_treeplot +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) + 
  geom_tiplab(aes(label = species), size = 1.6) + 
  theme(legend.position = "none")

ggsave(paste("../data/", rundate, "/figures/radial_tree_borrelia_annotated.jpg", sep=""), height =6, width = 6)

ggtree(z) %<+% metadata_bb +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) + 
  geom_tiplab(aes(label = species), size = 1.6) + 
  theme(legend.position = "none")

ggsave(paste("../data/", rundate, "/figures/bb_tree_borrelia.jpg", sep=""), height =6, width = 6)

bb_tree <- ggtree(z) %<+% metadata_bb +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) 
 
bm_tree <- ggtree(xx) %<+% metadata_bm +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) 


# generate figure 1:

genus_tree <- p %<+% metadata + 
  geom_tippoint(aes(color=`genus`, subset=!is.na(`genus`)), size = 3)
genus_tree


metadata_clean <- metadata %>% 
  filter(species %in% c("afzelii", "burgdorferi", "bavariensis", "garinii", "miyamotoi", 
                                "valaisiana", "crocidurae", "duttonii", "hermsii", "recurrentis", "turicatae") )
genus_tree_filt <- p %<+%  metadata_clean + 
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`)), size = 3)

# annotate tree with roary data and annotations
genotype <- as.data.frame(read_csv(paste("../data/", rundate, "/gene_presence_absence.csv", sep="")))
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-c(1:14)]
# convert to Rtab object
gene_mat <- mapply(function(x) ifelse(!is.na(x), 1, 0), genotype)
gene_mat <- data.frame(gene_mat, check.names = FALSE)
gene_mat$Gene <- rownames(genotype)
# read in annotation from roary
annotation <- read.csv(paste("../data/", rundate, "/gene_presence_absence.csv", sep=""), check.names=FALSE)
annotation <- annotation[,c(1,3)]
annotation$Gene_trunc <- gsub("group", "g", annotation$Gene)
annotation$Annotation <- gsub("hypothetical protein", "hyp", annotation$Annotation)


lipo_anno <- read_csv(paste("../data/", rundate, "/splip_out.csv", sep=""))
colnames(lipo_anno) <- c("Gene", "Lipo")
annotation <- left_join(annotation, lipo_anno)

# read in and annotate genome positions ####

pangenome_fa <- readDNAStringSet(paste("../data/", rundate, "/pan_genome_reference.fa", sep=""))
system(paste("minimap2 -x asm5 ../data/genomes/combined.fna ../data/", rundate, "/pan_genome_reference.fa > ../data/", rundate, "/pan_genome_alignment.paf", sep=""))
pangenome_mapping <- read_tsv(paste("../data/", rundate, "/pan_genome_alignment.paf", sep=""), col_names = FALSE) %>% 
  select(c(X1,X7)) 
names(pangenome_mapping)[1:2] <- c("Allele", "Mapping_size")
pangenome_mapping_unique <- pangenome_mapping %>% 
  distinct(Allele, .keep_all=TRUE)

pangenome_mapping_unique <- pangenome_mapping_unique %>%
  mutate(Chromosome = ifelse(Mapping_size > 800000, "Chromosome", "Plasmid"))

pangenome_db <- data.frame( Allele = str_split_fixed(names(pangenome_fa), " ", 2)[,1],
                            Gene = str_split_fixed(names(pangenome_fa), " ", 2)[,2])
# remove duplicate mapping
pangenome_db_unique <- pangenome_db %>% 
  distinct(Allele, .keep_all=TRUE)

pangenome_db_unique <- left_join(pangenome_db_unique, pangenome_mapping_unique)


annotation <- left_join(annotation, pangenome_db_unique)

gt.annotated <- merge(annotation, gene_mat) 
rownames(gt.annotated) <- paste(gt.annotated$Gene_trunc,":  ", gt.annotated$Annotation, gt.annotated$Allele, sep="")

gt.annotated <- gt.annotated[,-c(1:7)]

chrom.annotated <- as_tibble(annotation) %>% 
  filter(Chromosome == "Chromosome")

ch.annotated <- as.data.frame(left_join(chrom.annotated, gene_mat))
rownames(ch.annotated) <- paste(ch.annotated$Gene_trunc,":  ", ch.annotated$Annotation, sep="")
ch.annotated <- ch.annotated[,-c(1:7)]

plasmid.annotated <- as_tibble(annotation) %>% 
  filter(Chromosome == "Plasmid")

pl.annotated <- as.data.frame(left_join(plasmid.annotated, gene_mat))
rownames(pl.annotated) <- paste(pl.annotated$Gene_trunc,":  ", pl.annotated$Annotation, sep="")
pl.annotated <- pl.annotated[,-c(1:7)]


lipo.annotated <- as_tibble(annotation) %>% 
  filter(Lipo == "PROBABLE LIPOPROTEIN")

lp.annotated <- as.data.frame(left_join(lipo.annotated, gene_mat))
rownames(lp.annotated) <- paste(lp.annotated$Gene_trunc,":  ", lp.annotated$Annotation, sep="")
lp.annotated <- lp.annotated[,-c(1:7)]

# compute number of ORFS; might be nice to be able to plot these with a barplot next to the tree, but haven't figured out how to do this yet.
genome_length <- rowSums(t(gt.annotated))
genome_lengths <- tibble(`Number of ORFs` = genome_length, label = colnames(gt.annotated))
metadata <- left_join(metadata, genome_lengths, by = "label")
lipo_length <- rowSums(t(lp.annotated))
lipo_lengths <- tibble(`Number of Lipoproteins` = lipo_length, label = colnames(lp.annotated))
metadata <- left_join(metadata, lipo_lengths, by = "label")

#metadata[metadata$`Assembly Name` == "ASM2072076v1","Number of ORFs"] = NA

top_list <- names(table(metadata$species)[order(table(metadata$species), decreasing=TRUE)])

genus_plot <- metadata %>% #filter(`Assembly Level` %in%  c("Complete Genome", "Scaffold")) %>%
  filter(genus %in% c("Borrelia", "Borreliella")) %>%
  filter(`Number of ORFs` > 1000) %>% 
  #filter(species != "miyamotoi") %>%
ggboxplot(x = "genus", y = "Number of ORFs",
             add = "jitter", 
          add.params = list(alpha = 0.2), 
          outlier.shape = NA,
             alpha = 0.35) + 
  #stat_compare_means(comparisons = list(c("Borrelia", "Borreliella")))+ 
  rotate_x_text(angle = 90) + 
  theme(axis.text.x = element_text(face = "italic"))

species_plot <- metadata %>% 
  #filter(`Assembly Level` %in%  c("Complete Genome", "Scaffold")) %>%
  filter(species %in% top_list[1:25]) %>%
  filter(species != "sp.") %>%
  filter(`Number of ORFs` > 1000) %>%
  ggboxplot(x = "species", y = "Number of ORFs", color = "species",
               add = "jitter", 
            add.params = list(alpha = 0.25), 
            outlier.shape = NA) + 
  rotate_x_text(angle = 90) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(face = "italic"))

genus_plot_lipo <- metadata %>% #filter(`Assembly Level` %in%  c("Complete Genome", "Scaffold")) %>%
  filter(genus %in% c("Borrelia", "Borreliella")) %>%
  #filter(species != "miyamotoi") %>%
  #filter(`Number of ORFs` > 1000) %>%
  ggboxplot(x = "genus", y = "Number of Lipoproteins",
            add = "jitter", 
            add.params = list(alpha = 0.2), 
            outlier.shape = NA,
            alpha = 0.35) + 
  #stat_compare_means(comparisons = list(c("Borrelia", "Borreliella")))+ 
  rotate_x_text(angle = 90) + 
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

species_plot_lipo <- metadata %>% 
  #filter(`Assembly Level` %in%  c("Complete Genome", "Scaffold")) %>%
  filter(species %in% top_list[1:25]) %>%
  filter(species != "sp.") %>%
  filter(`Number of ORFs` > 1000) %>%
  ggboxplot(x = "species", y = "Number of Lipoproteins", color = "species",
            add = "jitter", 
            add.params = list(alpha = 0.25), 
            outlier.shape = NA) + 
  rotate_x_text(angle = 90) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(face = "italic"))


plot_grid(genus_plot, species_plot, genus_plot_lipo, species_plot_lipo, nrow = 2, rel_widths = c(1,2), labels = c("A", "B", "C", "D"))
ggsave(paste("../data/", rundate, "/figures/Figure_6.jpg", sep=""), height =8.4, width = 10)


# Generate MDS plots annotated with metadata ####
# Create core + accessory genome PAM ####

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 3 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
genus_tree_relab <- species_tree
p_hm <- heatMapPAM(genus_tree_relab,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_hm
ggsave(paste("../data/", rundate, "/figures/heatmap.jpg", sep=""), height =30, width = 50, limitsize = FALSE)

upper_threshold = 100 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 3 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
species_tree_accessory2 <- ggtree(a) %<+% metadata + 
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = RST_Type), size = 4)
p_hma <- heatMapPAM(species_tree_accessory2,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_hma
ggsave(paste("../data/", rundate, "/figures/heatmap_accessory.jpg", sep =""), height =30, width = 50, limitsize = FALSE)

# Generate MDS plots annotated with metadata ####
# Create core + accessory genome PAM ####

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 2 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
species_tree_relab <- species_tree_sub
p_hm2 <- heatMapPAM(species_tree_relab,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=7, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_hm2
ggsave(paste("../data/", rundate, "/figures/Figure_2.jpg", sep=""), height =5, width = 8, limitsize = FALSE)

# Create core + accessory genome PAM ####

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 2 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
p_hm3 <- heatMapPAM(bb_tree,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                    colnames_angle = -90,hjust =1, width=7, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_hm3
ggsave(paste("../data/", rundate, "/figures/heatmap_bb_subset.jpg", sep=""), height =5, width = 8, limitsize = FALSE)

bm_tree_annotated <- ggtree(xx) %<+% metadata_bm +
  geom_tippoint(aes(color=`species`, subset=!is.na(`species`), shape = genus), size = 2, alpha = 0.4) + 
  geom_tiplab(aes(label = `Organism Infraspecific Names Strain`), offset = 0.0005)

bm_tree_annotated

upper_threshold = 5 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 5 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
p_hm4 <- heatMapPAM(bm_tree_annotated,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                    colnames_angle = -90,width=6, font.size=4, cluster_cols=TRUE, 
                    offset = 0.005, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_hm4
ggsave(paste("../data/", rundate, "/figures/heatmap_bm_subset.jpg", sep=""), height =5, width = 20, limitsize = FALSE)


upper_threshold = 20 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 3 # remove from plotting anything in less than lower_threshold of isolates
lipo.trunc <- lp.annotated[rowSums(lp.annotated)> lower_threshold & (rowSums(lp.annotated) < (ncol(lp.annotated)-upper_threshold)),]
species_tree_relab2 <- species_tree_relab + theme(legend.text = element_text(face = "italic"))
p_lipo <- heatMapPAM(species_tree_relab,t(lipo.trunc), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_lipo
ggsave(paste("../data/", rundate, "/figures/Figure_3.jpg", sep=""), height =5, width = 7, limitsize = FALSE)

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 3 # remove from plotting anything in less than lower_threshold of isolates
chrom.trunc <- ch.annotated[rowSums(ch.annotated)> lower_threshold & (rowSums(ch.annotated) < (ncol(ch.annotated)-upper_threshold)),]
p_chrom <- heatMapPAM(species_tree_relab,t(chrom.trunc), colnames=FALSE, col_colours="blue", 
                     colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_chrom
ggsave(paste("../data/", rundate, "/figures/Figure_S2.jpg", sep=""), height =6, width = 8, limitsize = FALSE)

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 3 # remove from plotting anything in less than lower_threshold of isolates
plas.trunc <- pl.annotated[rowSums(pl.annotated)> lower_threshold & (rowSums(pl.annotated) < (ncol(pl.annotated)-upper_threshold)),]
p_plas <- heatMapPAM(species_tree_relab,t(plas.trunc), colnames=FALSE, col_colours="blue", 
                      colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_plas
ggsave(paste("../data/", rundate, "/figures/Figure_S3.jpg", sep=""), height =6, width = 8, limitsize = FALSE)

# correlate gene presence/absence with genotype ####

gt_mat <- t(gene_mat)
colnames(gt_mat) <- gt_mat[nrow(gt_mat),]
gt_mat <- as.data.frame(gt_mat[-c(nrow(gt_mat)),])
gt_mat <- gt_mat %>% mutate_if(is.character, as.numeric)
gt_mat <- gt_mat[,order(colSums(gt_mat))]

# construct outcome variables
gene_cors <- merge(metadata, gt_mat, by.x = "label", by.y = "row.names")
gene_cors$borrelia <- ifelse(gene_cors$genus != "Borreliella", 1, 0)
gene_cors$borreliella <- ifelse(gene_cors$genus == "Borreliella", 1, 0)
gene_cors$hermsii <- ifelse(gene_cors$species == "hermsii", 1, 0)
gene_cors$crocidurae <- ifelse(gene_cors$species == "crocidurae", 1, 0)
gene_cors$duttonii <- ifelse(gene_cors$species == "duttonii", 1, 0)


borreliaBeta <- rep(NA, ncol(gt_mat))
borreliellaBeta <- rep(NA, ncol(gt_mat))
hermsiiBeta <- rep(NA, ncol(gt_mat))
crociduraeBeta <- rep(NA, ncol(gt_mat))
duttoniiBeta <- rep(NA, ncol(gt_mat))
genenames <- rep(NA, ncol(gt_mat))

# need to get the hard coded number out of this...
for(i in 1:(ncol(gt_mat) - 24)){
  genenames[i] <- colnames(gene_cors)[(i+15)]
  borreliaBeta[i] <- glm(gene_cors$borrelia ~ gene_cors[,(i+15)] - 1, family = "binomial")$coefficients[1]
  borreliellaBeta[i] <- glm(gene_cors$borreliella ~ gene_cors[,(i+15)] -1, family = "binomial")$coefficients[1]
  hermsiiBeta[i] <- glm(gene_cors$hermsii ~ gene_cors[,(i+15)] - 1, family = "binomial")$coefficients[1]
  crociduraeBeta[i] <- glm(gene_cors$crocidurae ~ gene_cors[,(i+15)] -1, family = "binomial")$coefficients[1]
  duttoniiBeta[i] <- glm(gene_cors$duttonii ~ gene_cors[,(i+15)] -1, family = "binomial")$coefficients[1]
}

association_df <- data.frame(Gene = genenames, borreliaBeta = borreliaBeta, borreliellaBeta = borreliellaBeta,
                             hersmiiBeta = hermsiiBeta, crociduraeBeta = crociduraeBeta, duttoniiBeta = duttoniiBeta)

rosetters <- association_df %>% 
  filter(borreliaBeta > 12, borreliellaBeta < -12)
rosetters <- left_join(rosetters, lipo_anno) %>% 
  filter(Lipo == "PROBABLE LIPOPROTEIN")
rosetters_anno <- as.data.frame(left_join(rosetters, gene_mat))
rownames(rosetters_anno) <- paste(rosetters_anno$Gene,":  ", rosetters_anno$Annotation, sep="")
rosetters_anno <- rosetters_anno[,-c(1:6)]

p_ros <- heatMapPAM(species_tree_relab,t(rosetters_anno), colnames=FALSE, col_colours="blue", 
                     colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_ros

ggsave(paste("../data/", rundate, "/figures/Figure_S5.jpg", sep=""), height =6, width = 8, limitsize = FALSE)

rosetters <- left_join(rosetters, annotation)
barplot(table(rosetters$Chromosome))

pangenome_map <- pangenome_db %>% 
  mutate(genename = paste(Allele, Gene, sep=" "))

rosetters <- left_join(rosetters, pangenome_map)
pangenome_fa[rosetters$genename]

LDlipo <- association_df %>% 
  filter(borreliaBeta < -5, borreliellaBeta > 5)
LDlipo <- left_join(LDlipo, lipo_anno) %>% 
  filter(Lipo == "PROBABLE LIPOPROTEIN")
LDlipo_anno <- as.data.frame(left_join(LDlipo, gene_mat))
rownames(LDlipo_anno) <- paste(LDlipo_anno$Gene,":  ", LDlipo_anno$Annotation, sep="")
LDlipo_anno <- LDlipo_anno[,-c(1:6)]
metadata_filt <- metadata_filt %>% 
  mutate(OspC = ifelse(OspC_Type == "A", "A", "non-A"))

species_tree_relab2 <- ggtree(y) %<+% metadata_filt + 
  geom_tippoint(aes(color=species, subset=!is.na(species), shape = genus), size = 4)
p_LD <- heatMapPAM(species_tree_relab2,t(LDlipo_anno), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_LD

ggsave(paste("../data/", rundate, "/figures/Figure_S4A.jpg", sep=""), height =6, width = 8, limitsize = FALSE)


species_tree_relab3 <- ggtree(y) %<+% metadata_filt + 
  geom_tippoint(aes(color=RST_Type, subset=!is.na(RST_Type)), size = 2, alpha = 0.5)
p_LD_RST <- heatMapPAM(species_tree_relab3,t(LDlipo_anno), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_LD_RST

ggsave(paste("../data/", rundate, "/figures/Figure_S4B.jpg", sep=""), height =6, width = 8, limitsize = FALSE)

species_tree_relab4 <- ggtree(y) %<+% metadata_filt + 
  geom_tippoint(aes(color=OspC, subset=!is.na(OspC)), size = 2, alpha = 0.5)
p_LD_OspC <- heatMapPAM(species_tree_relab4,t(LDlipo_anno), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

p_LD_OspC
ß
ggsave(paste("../data/", rundate, "/figures/Figure_S4C.jpg", sep=""), height =6, width = 8, limitsßize = FALSE)
ßœ
LDlipo <- left_join(LDlipo, annotation)
#barplot(table(LDlipo$Chromosome))

# barplot of lipoprotein counts ####

annotation_accessory <- left_join(annotation, association_df) %>% 
  mutate(genus_weight = ifelse(borreliaBeta < -5 & borreliellaBeta > 5, "Borreliella", "Borrelia"))
annotation_accessory %>% 
  filter(!is.na(Chromosome), !is.na(genus_weight)) %>%
ggplot(aes(Chromosome, ..count..)) + geom_bar(aes(fill = Lipo), position = "dodge") + 
  labs(x = "SpLip Classification", y = "Count") + 
  scale_fill_discrete(name="") + 
  theme_bw() + 
  theme(legend.position = "top") + 
  facet_grid(~genus_weight)
  
ggsave(paste("../data/", rundate, "/figures/lipo_barplots.jpg", sep=""), height =4, width = 8)




