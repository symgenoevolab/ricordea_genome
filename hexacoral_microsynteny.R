################################################################################
############ Quantify microsynteny conservation across Hexacorallia ############
################################################################################

###### Code for manuscript Lewin, Sakagami, Piñon-Gonzalez et al. 2026, 
###### "Corallimorpharian genome supports the monophyly of Scleractinia and 
###### illuminates the evolution of coral calcification" 

###### Thomas D. Lewin
###### Sym Geno Evo Lab
###### Biodiversity Research Center, Academia Sinica
###### Created 14/02/2025
###### Last edited 22/09/2026

################################################################################
############################### PART 1: GENESPACE ##############################
################################################################################

###### This code follows the vignette from J T Lovell.
# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html

######################### STEP 0: Prepare input files ##########################

###### GENESPACE requires as input a directory with two subdirectories:
###### (1) "bed", containing bed files for each input genome
###### (2) "peptide", containing peptide files for each input genome

######################### STEP 1: Set up dependencies ##########################

# Set working directory 
setwd("/path/to/working/directory/")

# Load GENESPACE package
library(GENESPACE)

# Add OrthoFinder to path
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/path/to/Orthofinder/bin", sep = ":"
  )
)

# Add programmes from MCSCanX environment to path
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/path/to/MCScanX/bin", sep = ":"
  )
)

# Add MCScanX itself to path
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/path/to/MCScanX-master", sep = ":"
  )
)

# Add SciPy to path 
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/path/to/SciPy/bin", sep = ":"
  )
)

###################### Step 2: Set up the GENESPACE input ######################

# Set directory
genomeRepo <- "/path/to/input/directory/"

# Set directory
wd <- "/path/to/input/directory/"

# Set path to MCScanX directory
path2mcscanx <- "/path/to/MCScanX-master"

# Initiate the GENESPACE environment
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  blkSize = 5)
### For full list of parameters, see: https://rdrr.io/github/jtlovell/GENESPACE/man/init_genespace.html

############################# Run GENESPACE itself #############################

# Run GENESPACE
genespace_1 <- run_genespace(gpar, makePairwiseFiles = TRUE)

################################################################################
############################### PART 2: SYNTENET ###############################
################################################################################

###### This code follows the vignette from Bioconductor
# https://www.bioconductor.org/packages/release/bioc/vignettes/syntenet/inst/doc/syntenet.html#9_syntenet_as_a_synteny_detection_tool

######################### STEP 0: Prepare input files ##########################

###### SYNTENET requires as input two directories:
###### (1) "gtf_files", containing gtf files for each input genome
###### (2) "fasta_files", containing fasta files for each input genome

######################### STEP 1: Set up dependencies ##########################

# Load packages
library(syntenet)
library(GenomicRanges)
library(dplyr)
library(treeio)
library(ggtree)
library(ape)

######################### STEP 2: read in input files ##########################

###### fasta files

# Get fasta files into AAStringSet list object
fasta_dir <- file.path(getwd(), "fasta_files")
dir(fasta_dir)
aastringsetlist <- fasta2AAStringSetlist(fasta_dir)
proteomes <- aastringsetlist

# Modify gene names to remove transcript (e.g. g1000.t1 -> g1000)
# This only does so if it starts with 'g', such as from BRAKER
# also removes everything after the first space for all
proteomes <- lapply(proteomes, function(x) {
  names(x) <- sub("^(g[0-9]+)\\..*", "\\1", names(x))  # Remove suffix after dot if it starts with 'g'
  names(x) <- sub("\\s.*", "", names(x))  # Remove everything after the first space for all names
  return(x)
})

###### gtf files

#Get gtf files into GRangesList object
gtf_dir <- file.path(getwd(), "gtf_files")
dir(gtf_dir)

grangeslist <- gff2GRangesList(gtf_dir)
annotation <- grangeslist

##### Check the input data
check_input(proteomes, annotation)

############################ STEP 3: Run SYNTENET ##############################

# Process the data
pdata <- process_input(proteomes, annotation)

###### Run Diamond

data(blast_list)

if(diamond_is_installed()) {
  blast_list <- run_diamond(seq = pdata$seq)
}

# List names
names(blast_list)

###### Infer synteny network

net <- infer_syntenet(blast_list, pdata$annotation)

# Get a 2-column data frame of species IDs and names
id_table <- create_species_id_table(names(proteomes))

###### Phylogenomic profiling

clusters <- cluster_network(net)
write.csv(clusters, file = "clusters.csv")
profiles <- phylogenomic_profile(clusters)

# Set the order of species and their names
species_order <- setNames(
  # vector elements
  c(
    "Gve", "Dli", "Sca", "Nve", "Ppe", "Ros", "Rfl","Ryu", "Pve", "Mme", "Dcy", "Pcy", "Ssi" ,"Ami"
  ),
  # vector names
  c(
    "Gorgonia ventalina","Diadumene lineata", "Scolanthus callimorphus","Nematostella vectensis", "Plumapathes pennacea", "Rhodactis osculifera", "Ricordea florida", "Ricordea yuma", "Pocillopora verrucosa", "Meandrina meandrites","Dendrogyra cylindrus", "Porites cylindrica", "Siderastrea siderea", "Acropora millepora"
  )
)

# Annotate orders for each species
species_annotation <- data.frame(
  Species = species_order,
  Clade = c("Octocorallia", "Actiniaria",  "Actiniaria",  "Actiniaria",  "Antipatharia",  "Corallimorpharia",  "Corallimorpharia",  "Corallimorpharia" , "Scleractinia" , "Scleractinia" , "Scleractinia", "Scleractinia", "Scleractinia", "Scleractinia"
  )
)

# Plot heatmap
plot1 <- plot_profiles(
  profiles, 
  species_annotation,
  cluster_species = species_order
)

# Find clade-specific clusters
gs_clusters <- find_GS_clusters(profiles, species_annotation)
gs_clusters
write.csv(gs_clusters, file = "gs_clusters.csv")

# How many clade-specific clusters are there?
nrow(gs_clusters)

# Filter profiles matrix to only include group-specific clusters
idx <- rownames(profiles) %in% gs_clusters$Cluster
p_gs <- profiles[idx, ]

# Plot heatmap
plot2 <- plot_profiles(
  p_gs, species_annotation, 
  cluster_species = species_order, 
  cluster_columns = TRUE
)

# Visualize a network of first 5 GS-clusters
id <- gs_clusters$Cluster[1:5]
plot_network(net, clusters, cluster_id = id)

############## STEP 4: Microsynteny-based phylogeny construction ###############

bt_mat <- binarize_and_transpose(profiles)

included <- c("Gve", "Dli", "Nve", "Ppe", "Ros", "Rfl","Ryu", "Pve", "Mme", "Dcy", "Pcy", "Ssi" ,"Ami")
bt_mat <- bt_mat[rownames(bt_mat) %in% included, ]

# Remove non-variable sites
bt_mat <- bt_mat[, colSums(bt_mat) != length(included)]

write.csv(bt_mat, "bt_mat.csv")

# Add IQ-TREE into path 
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/path/to/iqtree/bin", sep = ":"
  )
)

# Check it's added correctly 
iqtree_is_installed()

# Run IQ-TREE
if(iqtree_is_installed()) {
  phylo <- infer_microsynteny_phylogeny(bt_mat, outgroup = "Gve", 
                                        threads = 30)
}

# Extract link to tree from phylo object
phylo

# Read tree in 
tree <- readLines("/path/to/treefile")

# Copy and paste tree into iToL for rendering
tree

################################################################################
############################### END OF SCRIPT ##################################
################################################################################
