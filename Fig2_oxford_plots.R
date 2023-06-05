library(ggpubr)
library(wesanderson)
library(tidyverse)
library(patchwork)

### Functions ###
read_busco <- function(buscoFile){
  busco <- read_tsv(buscoFile,
                    col_names = c("busco_id", "Status", "Sequence",
                                  "start", "end", "strand", "Score", "Length",
                                  "OrthoDB_url", "Description"),
                    col_types = c("ccciicdicc"),
                    comment = "#") %>%
    filter(Status == "Complete") %>%
    select(busco_id, Sequence, start, end)
  return(busco)
}

filter_buscos_x <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.x) %>%
    mutate(nGenes = n(),
           mxGpos = max(midpos_x)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}

filter_buscos_y <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.y) %>%
    mutate(nGenes = n(),
           mxGpos_y = max(midpos_y)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}

scale_lengths <- function(buscos){
  buscos$start <- buscos$start / 1000000
  buscos$end <- buscos$end / 1000000
  buscos$midpos <- (buscos$start + buscos$end)/2
  return(buscos)
}

make_pairwise_dataframe <- function(df1, df2){
  df1$midpos_x <- df1$midpos
  df2$midpos_y <- df2$midpos
  pairwise_df = merge(df1, df2, by="busco_id")
  return(pairwise_df)
}

### Main ###
setwd('/Users/cw22/Documents/R_work/Hyles_analysis/')


make_oxford_plot <- function(df){
  colour_palette <- c("black", "darkgrey", "#CC9088", "#6CE848", "#E6E6E0", "#D88E43", "#B5A4E2", "#6D7DE0",
                      "#DE4B44", "#E499D1", "#669F66", "#E054DC", "#93C656",
                      "#E1E9AA", "#C47AD7", "#A9E99E", "#C2EFDF", "#E5D265",
                      "#DB4DA2", "#D6ED49", "#DFBBD6", "#825ADF", "#78EAD3", "#63B9E6",
                      "#DA6F8F", "#608182", "#6A94D6", "#AEBDA0", "#E7BF9C", "#67E97D", "#73DAE0", "#AC33ED")
  
  
  the_plot <- ggplot() +  # scales=free, space=free is what you want to get boxes scaled by length :) 
    facet_grid(Sequence.x_f~Sequence.y_f, space = "free", scales="free")  + 
    geom_point(data = df, aes(x = midpos_y, y = midpos_x, color=Merian), size=0.1) + # was Sequence.x_f 
    theme_bw() + 
    theme(legend.position = "none", axis.ticks = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_rect(color = "#999999", size=0.3)) + 
    theme(panel.spacing = unit(0, "lines")) +
    theme(strip.background = element_blank()) +
    theme(strip.text.y = element_text(angle = 0,size=10)) +
    theme(strip.text.x = element_text(angle = 90,size=10)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    geom_point(data=df, aes(x=sequence_y_min, y=sequence_x_min), color="white") + # needed to scale each box
    geom_point(data=df, aes(x=sequence_y_max, y=sequence_x_max), color="white")  +   # needed to scale each box +
  scale_color_manual(values=colour_palette) +
    theme(strip.text = element_text(hjust = 0))
  return(the_plot)
}

Hyles_euphoribae <- read_busco('../Data/BUSCOs/All/Hyles_euphorbiae.tsv')
Hyles_vespertilio <- read_busco('../Data/BUSCOs/All/Hyles_vespertilio.tsv')
Deilephila_porcellus <- read_busco('../Data/BUSCOs/All/Deilephila_porcellus.tsv')
Manduca_sexta <- read_busco('../Data/BUSCOs/All/Manduca_sexta.tsv')
Laothoe_populi <- read_busco('../Data/BUSCOs/All/Laothoe_populi.tsv')

chr_conversions <- read.csv('../../Hyles_analysis/Table_S1_conversions.csv', header=FALSE, skip=4, sep=',')

# Reformat data
colnames(chr_conversions) <- c('Bmori_chr', 
                               'Hvespertilio_chr', 'Hvespertilio_accession', 'Hvespertilio_scaff','Hvespertilio_length', 
                               'Heuphorbia_chr', 'Heuphorbia_accession', 'Heuphorbia_scaff', 'Heuphorbia_length',
                               'Msexta_chr', 'Msexta_NCBI', 'Msexta_accession', 'Msexta_scaff', 'Msexta_length')

# Make conversion tables
Hyles_euphoribae_conversion <-chr_conversions[,c('Heuphorbia_chr','Heuphorbia_accession')]
Hyles_vespertilio_conversion <-chr_conversions[,c('Hvespertilio_chr','Hvespertilio_accession')]
Manduca_conversion <-chr_conversions[,c('Msexta_chr','Msexta_accession')]

colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr', 'Sequence')
colnames(Hyles_vespertilio_conversion) <- c('Relabelled_chr', 'Sequence')
colnames(Manduca_conversion) <- c('Relabelled_chr', 'Sequence')

Hyles_euphoribae <- scale_lengths(Hyles_euphoribae)
Hyles_vespertilio <- scale_lengths(Hyles_vespertilio)
Deilephila_porcellus <- scale_lengths(Deilephila_porcellus)
Manduca_sexta <- scale_lengths(Manduca_sexta)
Laothoe_populi <- scale_lengths(Laothoe_populi)

# Make pairwise dataframes
HylesE_vs_Deilephila <- make_pairwise_dataframe(Hyles_euphoribae,Deilephila_porcellus)
HylesE_vs_HylesV <- make_pairwise_dataframe(Hyles_euphoribae,Hyles_vespertilio)
HylesE_vs_Laothoe <- make_pairwise_dataframe(Hyles_euphoribae, Laothoe_populi)
HylesE_vs_Manduca <- make_pairwise_dataframe(Hyles_euphoribae,Manduca_sexta)


# species to be on x-axes needs to be specified second i.e. "x" actually is y
HylesE_vs_HylesV <- filter_buscos_x(HylesE_vs_HylesV,2) %>% filter_buscos_y(1)
HylesE_vs_HylesV <- merge(HylesE_vs_HylesV, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
HylesE_vs_HylesV$Merian <- factor(HylesE_vs_HylesV$Merian, levels = merian_order)

# Make columns to define limits of each facet
HylesE_vs_HylesV <- HylesE_vs_HylesV %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
HylesE_vs_HylesV <- HylesE_vs_HylesV %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
HylesE_vs_HylesV$sequence_y_min <- 0
HylesE_vs_HylesV$sequence_x_min <- 0

# Default order of chr is by size 
sequence_x_order <- HylesE_vs_HylesV %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- HylesE_vs_HylesV %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()
sequence_x_order <- rev(sequence_x_order)
sequence_y_order <- rev(sequence_y_order)

# sequence y order
sequence_y_order <- c("CM042833.1","CM042825.1","CM042808.1","CM042815.1",
                      "CM042818.1","CM042807.1","CM042805.1","CM042820.1",
                      "CM042813.1","CM042816.1","CM042812.1","CM042826.1",
                      "CM042814.1","CM042821.1","CM042809.1","CM042811.1",
                      "CM042806.1","CM042824.1","CM042822.1","CM042828.1",
                      "CM042817.1","CM042819.1","CM042810.1","CM042823.1",
                      "CM042831.1","CM042830.1","CM042827.1" ,"CM042832.1",
                      "CM042829.1","WUWR02000219.1")

HylesE_vs_HylesV$Sequence.x_f <- factor(HylesE_vs_HylesV$Sequence.x, levels = sequence_x_order)
HylesE_vs_HylesV$Sequence.y_f <- factor(HylesE_vs_HylesV$Sequence.y, levels = sequence_y_order)

# Specify sequences to flip to make orientations consistent
sequences_to_flip <- c("CM041050.1", "CM041032.1", "CM041041.1", 
                       "CM041037.1", "CM041046.1", "CM041031.1",
                       "CM041049.1", "CM041047.1", "CM041042.1",
                       "CM041048.1","CM041052.1", "CM041057.1",
                       "CM041054.1")

# Flip specified sequences
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_HylesV %>% filter(Sequence.x == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_HylesV <- HylesE_vs_HylesV %>%
    filter(Sequence.x != i) %>%
    bind_rows(subset_df)
}


# OPTIONAL: Just before plotting - replace Accession_IDs with relabelled chr name for H. euphorbiae
colnames(Hyles_vespertilio_conversion) <- c('Relabelled_chr_y', 'Sequence.y')
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr_x', 'Sequence.x')

HylesE_vs_HylesV <- merge(HylesE_vs_HylesV,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_HylesV <- merge(HylesE_vs_HylesV,Hyles_vespertilio_conversion, by="Sequence.y")
HylesE_vs_HylesV <- HylesE_vs_HylesV[,!names(HylesE_vs_HylesV) %in% c("Sequence.x", "Sequence.y")]
names(HylesE_vs_HylesV)[names(HylesE_vs_HylesV) == 'Relabelled_chr_x'] <- 'Sequence.x'
names(HylesE_vs_HylesV)[names(HylesE_vs_HylesV) == 'Relabelled_chr_y'] <- 'Sequence.y'

sequence_x_order <- HylesE_vs_HylesV %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull() %>% rev()
sequence_y_order <- HylesE_vs_HylesV %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull() %>% rev()


sequence_y_order <- c("HvChr1=Z","HvChr22#","HvChr5","HvChr12","HvChr15","HvChr4",
                      "HvChr2*","HvChr17","HvChr10","HvChr13","HvChr9","HvChr23",
                      "HvChr11","HvChr18","HvChr6","HvChr8","HvChr3","HvChr21",
                      "HvChr19","HvChr25","HvChr14","HvChr16","HvChr7","HvChr20",
                      "HvChr28#","HvChr27","HvChr24","HvChr29§","HvChr26**")



HylesE_vs_HylesV$Sequence.x_f <- factor(HylesE_vs_HylesV$Sequence.x, levels = sequence_x_order)
HylesE_vs_HylesV$Sequence.y_f <- factor(HylesE_vs_HylesV$Sequence.y, levels = sequence_y_order)


# NB: sequence X is actually on Y axis
HylesE_vs_HylesV_plot <- make_oxford_plot(HylesE_vs_HylesV) + 
  xlab("Hyles vespertilio") +   ylab("Hyles euphoribae") +
  theme(axis.title = element_text(face="italic"),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 
HylesE_vs_HylesV_plot

###

# Deilephila_porcellus vs Hyles_euphoribae
HylesE_vs_Deilephila <- filter_buscos_x(HylesE_vs_Deilephila,2) %>% filter_buscos_y(1)
HylesE_vs_Deilephila <- merge(HylesE_vs_Deilephila, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
HylesE_vs_Deilephila$Merian <- factor(HylesE_vs_Deilephila$Merian, levels = merian_order)

# Make columns to define limits of each facet
HylesE_vs_Deilephila <- HylesE_vs_Deilephila %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
HylesE_vs_Deilephila <- HylesE_vs_Deilephila %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
HylesE_vs_Deilephila$sequence_y_min <- 0
HylesE_vs_Deilephila$sequence_x_min <- 0


sequences_to_flip <- c("CM041030.1", "CM041037.1", "CM041051.1", "CM041046.1", "CM041031.1",
                       "CM041042.1", "CM041044.1", "CM041035.1", "CM041056.1",
                       "CM041052.1","CM041054.1" )
  
# Flip specified sequences
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Deilephila %>% filter(Sequence.x == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Deilephila <- HylesE_vs_Deilephila %>%
    filter(Sequence.x != i) %>%
    bind_rows(subset_df)
}

colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr', 'Sequence.x')
HylesE_vs_Deilephila <- merge(HylesE_vs_Deilephila,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Deilephila <- HylesE_vs_Deilephila[,!names(HylesE_vs_Deilephila) %in% c("Sequence.x")]
names(HylesE_vs_Deilephila)[names(HylesE_vs_Deilephila) == 'Relabelled_chr'] <- 'Sequence.x'
HylesE_vs_Deilephila$Sequence.x_f <- HylesE_vs_Deilephila$Sequence.x

# Default order of chr is by size 
sequence_x_order <- HylesE_vs_Deilephila %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull %>% rev()
sequence_y_order <- HylesE_vs_Deilephila %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull() %>% rev()


sequence_y_order <- c("LR999970.1","LR999971.1","LR999972.1","LR999975.1",
                     "LR999973.1","LR999974.1","LR999978.1","LR999980.1",
                    "LR999976.1","LR999977.1","LR999979.1","LR999981.1",
                   "LR999983.1","LR999984.1","LR999982.1","LR999985.1",
                  "LR999987.1","LR999986.1","LR999988.1","LR999989.1",
                 "LR999990.1","LR999991.1","LR999992.1","LR999993.1",
                "LR999995.1","LR999994.1","LR999996.1","LR999997.1",
                      "LR999998.1")
HylesE_vs_Deilephila$Sequence.x_f <- factor(HylesE_vs_Deilephila$Sequence.x, levels = sequence_x_order)
HylesE_vs_Deilephila$Sequence.y_f <- factor(HylesE_vs_Deilephila$Sequence.y, levels = sequence_y_order)

###
# Manduca_sexta vs Hyles_euphoribae
HylesE_vs_Manduca <- filter_buscos_x(HylesE_vs_Manduca,2) %>% filter_buscos_y(1)
HylesE_vs_Manduca <- merge(HylesE_vs_Manduca, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
HylesE_vs_Manduca$Merian <- factor(HylesE_vs_Manduca$Merian, levels = merian_order)

# Make columns to define limits of each facet
HylesE_vs_Manduca <- HylesE_vs_Manduca %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
HylesE_vs_Manduca <- HylesE_vs_Manduca %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
HylesE_vs_Manduca$sequence_y_min <- 0
HylesE_vs_Manduca$sequence_x_min <- 0


sequences_to_flip <- c("CM041058.1", "CM041050.1","CM041034.1","CM041030.1",
                       "CM041039.1", "CM041049.1", "CM041047.1",
                       "CM041042.1", "CM041035.1", "CM041055.1")

# Flip specified sequences
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Manduca %>% filter(Sequence.x == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Manduca <- HylesE_vs_Manduca %>%
    filter(Sequence.x != i) %>%
    bind_rows(subset_df)
}

colnames(Manduca_conversion) <- c('Relabelled_chr_y', 'Sequence.y')
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr_x', 'Sequence.x')
HylesE_vs_Manduca <- merge(HylesE_vs_Manduca,Manduca_conversion, by="Sequence.y")
HylesE_vs_Manduca <- merge(HylesE_vs_Manduca,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Manduca <- HylesE_vs_Manduca[,!names(HylesE_vs_Manduca) %in% c("Sequence.x", "Sequence.y")]
names(HylesE_vs_Manduca)[names(HylesE_vs_Manduca) == 'Relabelled_chr_y'] <- 'Sequence.y'
names(HylesE_vs_Manduca)[names(HylesE_vs_Manduca) == 'Relabelled_chr_x'] <- 'Sequence.x'
HylesE_vs_Manduca$Sequence.x_f <- HylesE_vs_Manduca$Sequence.x
HylesE_vs_Manduca$Sequence.y_f <- HylesE_vs_Manduca$Sequence.y

# Default order of chr is by size 
sequence_x_order <- HylesE_vs_Manduca %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull %>% rev()
sequence_y_order <- HylesE_vs_Manduca %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull() %>% rev()


sequence_y_order <- c("MsChr1=Z","MsChr22","MsChr15","MsChr12","MsChr5",
                      "MsChr4","MsChr2","MsChr17",
                      "MsChr10","MsChr13","MsChr9",
                      "MsChr23","MsChr11","MsChr18","MsChr6",
                      "MsChr8","MsChr21","MsChr19","MsChr3",
                      "MsChr25","MsChr16","MsChr14","MsChr7","MsChr20","MsChr28",
                      "MsChr27","MsChr24","MsChr26")

HylesE_vs_Manduca$Sequence.x_f <- factor(HylesE_vs_Manduca$Sequence.x, levels = sequence_x_order)
HylesE_vs_Manduca$Sequence.y_f <- factor(HylesE_vs_Manduca$Sequence.y, levels = sequence_y_order)

### Laothe_populi vs Hyles
HylesE_vs_Laothoe <- make_pairwise_dataframe(Hyles_euphoribae, Laothoe_populi)
HylesE_vs_Laothoe <- filter_buscos_x(HylesE_vs_Laothoe,2) %>% filter_buscos_y(1)
HylesE_vs_Laothoe <- merge(HylesE_vs_Laothoe, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
HylesE_vs_Laothoe$Merian <- factor(HylesE_vs_Laothoe$Merian, levels = merian_order)

# Make columns to define limits of each facet
HylesE_vs_Laothoe <- HylesE_vs_Laothoe %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
HylesE_vs_Laothoe <- HylesE_vs_Laothoe %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
HylesE_vs_Laothoe$sequence_y_min <- 0
HylesE_vs_Laothoe$sequence_x_min <- 0

sequences_to_flip <- c("CM041058.1", "CM041050.1","CM041033.1",
                       "CM041043.1","CM041032.1","CM041038.1",
                       "CM041037.1", "CM041045.1", "CM041046.1",
                       "CM041036.1", "CM041031.1", "CM041049.1",
                       "CM041047.1", "CM041053.1", "CM041048.1","CM041054.1")

# Flip specified sequences
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Laothoe %>% filter(Sequence.x == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Laothoe <- HylesE_vs_Laothoe %>%
    filter(Sequence.x != i) %>%
    bind_rows(subset_df)
}

# OPTIONAL: Just before plotting - replace Accession_IDs with relabelled chr name for H. euphorbiae
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr', 'Sequence.x')
HylesE_vs_Laothoe <- merge(HylesE_vs_Laothoe,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Laothoe <- HylesE_vs_Laothoe[,!names(HylesE_vs_Laothoe) %in% c("Sequence.x")]
names(HylesE_vs_Laothoe)[names(HylesE_vs_Laothoe) == 'Relabelled_chr'] <- 'Sequence.x'
HylesE_vs_Laothoe$Sequence.x_f <- HylesE_vs_Laothoe$Sequence.x

# Default order of chr is by size 
sequence_x_order <- HylesE_vs_Laothoe %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull %>% rev()
sequence_y_order <- HylesE_vs_Laothoe %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull() %>% rev()

sequence_y_order <- c("HG992149.1", "HG992146.1","HG992148.1","HG992155.1",
                      "HG992150.1","HG992151.1","HG992157.1","HG992156.1",
                      "HG992152.1","HG992153.1", "HG992154.1","HG992158.1",
                      "HG992147.1","HG992161.1","HG992159.1","HG992160.1",
                      "HG992164.1","HG992163.1","HG992165.1","HG992166.1",
                      "HG992168.1","HG992167.1","HG992169.1","HG992171.1",
                      "HG992170.1","HG992172.1","HG992173.1","HG992174.1")

HylesE_vs_Laothoe$Sequence.x_f <- factor(HylesE_vs_Laothoe$Sequence.x, levels = sequence_x_order)
HylesE_vs_Laothoe$Sequence.y_f <- factor(HylesE_vs_Laothoe$Sequence.y, levels = sequence_y_order)

HylesE_vs_HylesV_plot <- make_oxford_plot(HylesE_vs_HylesV) + 
  xlab("Hyles vespertilio") +   ylab("Hyles euphoribae") +
  theme(axis.title = element_text(face="italic", size=20),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 

HylesE_vs_Deilephila_plot <- make_oxford_plot(HylesE_vs_Deilephila) + 
  xlab("Deilephila porcellus") +   ylab("Hyles euphoribae") +
  theme(axis.title = element_text(face="italic", size=20),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 

HylesE_vs_Manduca_plot <- make_oxford_plot(HylesE_vs_Manduca) + 
  xlab("Manduca sexta") +   ylab("Hyles euphoribae") +
  theme(axis.title = element_text(face="italic", size=20),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 

HylesE_vs_Laothoe_plot <- make_oxford_plot(HylesE_vs_Laothoe) + 
  xlab("Laothoe populi") +   ylab("Hyles euphoribae") +
  theme(axis.title = element_text(face="italic", size=20),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 

HylesE_vs_HylesV_plot <- HylesE_vs_HylesV_plot + ggtitle('a)') + theme(title = element_text(face="bold"))
HylesE_vs_Deilephila_plot <- HylesE_vs_Deilephila_plot + ggtitle('b)') + theme(title = element_text(face="bold"))
HylesE_vs_Manduca_plot <- HylesE_vs_Manduca_plot + ggtitle('c)') + theme(title = element_text(face="bold"))
HylesE_vs_Laothoe_plot <- HylesE_vs_Laothoe_plot + ggtitle('d)') + theme(title = element_text(face="bold"))

ggsave('../../Hyles_analysis/revised_panelA.pdf', plot=HylesE_vs_HylesV_plot)
ggsave('../../Hyles_analysis/revised_panelB.pdf', plot=HylesE_vs_Deilephila_plot)
ggsave('../../Hyles_analysis/revised_panelC.pdf', plot=HylesE_vs_Manduca_plot)
ggsave('../../Hyles_analysis/revised_panelD.pdf', plot=HylesE_vs_Laothoe_plot)
