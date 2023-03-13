library(ggpubr)
library(wesanderson)
library(tidyverse)

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

make_oxford_plot <- function(df){
  the_plot <- ggplot(data = df, aes(x = midpos_x, y = midpos_y)) +  # scales=free, space=free is what you want to get boxes scaled by length :) 
    facet_grid(Sequence.x_f~Sequence.y_f, space = "free", scales="free")  + 
    geom_point(aes(color =color_f),size=0.1) + # was Sequence.x_f 
    scale_color_manual(values = wes_palette("Royal1")) +
    theme_bw() + 
    theme(legend.position = "none", axis.ticks = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_rect(color = "#999999", size=0.3)) + 
    theme(panel.spacing = unit(0, "lines")) +
    theme(strip.background = element_blank()) +
    theme(strip.text.y = element_text(angle = 0,size=6)) +
    theme(strip.text.x = element_text(angle = 90,size=6)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  return(the_plot)
}

### Main ###
setwd('/Users/cw22/Documents/R_work/Hyles_analysis/')

# Read in data
Hyles_euphoribae <- read_busco('Hyles_euphorbiae.tsv')
Hyles_vespertilio <- read_busco('Hyles_vespertilio.tsv')
Manduca_sexta <- read_busco('Manduca_sexta.tsv')
Deilephila_porcellus <- read_busco('Deilephila_porcellus.tsv')
Hemaris_fuciformis <- read_busco('Hemaris_fuciformis.tsv')
Mimas_tiliae <- read_busco('Mimas_tiliae.tsv')
Laothoe_populi <- read_busco('Laothoe_populi.tsv')
chr_conversions <- read.csv('Table_S1_conversions.csv', header=FALSE, skip=4, sep=',')

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

# Scale lengths to Mb
Hyles_euphoribae <- scale_lengths(Hyles_euphoribae)
Hyles_vespertilio <- scale_lengths(Hyles_vespertilio)
Manduca_sexta <- scale_lengths(Manduca_sexta)
Deilephila_porcellus <- scale_lengths(Deilephila_porcellus)
Hemaris_fuciformis <- scale_lengths(Hemaris_fuciformis)
Mimas_tiliae <- scale_lengths(Mimas_tiliae)
Laothoe_populi <- scale_lengths(Laothoe_populi)

# Make pairwise dataframes
HylesE_vs_HylesV <- make_pairwise_dataframe(Hyles_euphoribae, Hyles_vespertilio)
HylesE_vs_Deilephila <- make_pairwise_dataframe(Hyles_euphoribae, Deilephila_porcellus)
HylesE_vs_Manduca <- make_pairwise_dataframe(Hyles_euphoribae, Manduca_sexta)
HylesE_vs_Laothoe <- make_pairwise_dataframe(Hyles_euphoribae, Laothoe_populi)

# Filter BUSCOs (remove scaffs with < 3 BUSCOs)
HylesE_vs_HylesV <- filter_buscos_x(HylesE_vs_HylesV,2) %>% filter_buscos_y(2)
HylesE_vs_Deilephila <- filter_buscos_x(HylesE_vs_Deilephila,2) %>% filter_buscos_y(2)
HylesE_vs_Manduca <- filter_buscos_x(HylesE_vs_Manduca,2) %>% filter_buscos_y(2)
HylesE_vs_Laothoe <- filter_buscos_x(HylesE_vs_Laothoe,2) %>% filter_buscos_y(2)

### Hyles_euphoribae vs Laothoe_populi ###

# Default order of chr is by size 
sequence_x_order <- HylesE_vs_Laothoe %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- HylesE_vs_Laothoe %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()

# Manually specify order of chr to tidy plot
sequence_y_order <- c("HG992174.1","HG992173.1","HG992172.1","HG992170.1","HG992171.1","HG992169.1", "HG992167.1",
                      "HG992168.1","HG992166.1","HG992165.1","HG992163.1","HG992164.1","HG992160.1","HG992159.1",
                      "HG992161.1", "HG992147.1","HG992158.1","HG992154.1" ,"HG992153.1",
                      "HG992152.1","HG992156.1","HG992157.1","HG992151.1", "HG992150.1","HG992155.1",
                      "HG992148.1","HG992146.1","HG992149.1")

# Specify sequences to flip to make orientations consistent
sequences_to_flip <- c("HG992174.1","HG992171.1", "HG992166.1","HG992165.1","HG992163.1","HG992164.1","HG992160.1",
                       "HG992161.1","HG992154.1","HG992152.1","HG992156.1","HG992151.1", "HG992150.1",
                       "HG992147.1","HG992148.1","HG992146.1","HG992149.1","HG992147.1")

# Set manual order of chr as a factor
HylesE_vs_Laothoe$Sequence.x_f <- factor(HylesE_vs_Laothoe$Sequence.x, levels = sequence_x_order)
HylesE_vs_Laothoe$Sequence.y_f <- factor(HylesE_vs_Laothoe$Sequence.y, levels = sequence_y_order)

# Plot all non-rearranged chr in grey
HylesE_vs_Laothoe$color_f <- "grey"

# Plot rearranged chr in red
subset_df <- HylesE_vs_Laothoe %>% filter(Sequence.y == "HG992147.1")
subset_df$color_f <- "red"
HylesE_vs_Laothoe <- HylesE_vs_Laothoe %>%
  filter(Sequence.y != "HG992147.1") %>%
  bind_rows(subset_df)

# Flip specified sequences
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Laothoe %>% filter(Sequence.y == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Laothoe <- HylesE_vs_Laothoe %>%
    filter(Sequence.y != i) %>%
    bind_rows(subset_df)
}

# OPTIONAL: Just before plotting - replace Accession_IDs with relabelled chr name for H. euphorbiae
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr', 'Sequence.x')
HylesE_vs_Laothoe <- merge(HylesE_vs_Laothoe,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Laothoe <- HylesE_vs_Laothoe[,!names(HylesE_vs_Laothoe) %in% c("Sequence.x")]
names(HylesE_vs_Laothoe)[names(HylesE_vs_Laothoe) == 'Relabelled_chr.x'] <- 'Sequence.x'
HylesE_vs_Laothoe$Sequence.x_f <- HylesE_vs_Laothoe$Sequence.x

# Now need to reconfigure the order of x_sequences as defaults to alphabetical
sequence_x_order <- c("HeChr26","HeChr24","HeChr27","HeChr28","HeChr20","HeChr7","HeChr16","HeChr14",
                      "HeChr25","HeChr19","HeChr21","HeChr29",
                      "HeChr3","HeChr8","HeChr6","HeChr18",
                      "HeChr11","HeChr23","HeChr9","HeChr13","HeChr10","HeChr17","HeChr2","HeChr4",
                      "HeChr15","HeChr12","HeChr5","HeChr22","HeChr1=Z")
HylesE_vs_Laothoe$Sequence.x_f <- factor(HylesE_vs_Laothoe$Sequence.x, levels = sequence_x_order)


# Plot!
HylesE_vs_Laothoe_plot <- make_oxford_plot(HylesE_vs_Laothoe) + 
  xlab("Laothoe populi") + ylab("Hyles euphorbiae") + 
  theme(axis.title = element_text(face="italic"))
HylesE_vs_Laothoe_plot



### Hyles_euphoribae vs Deilephila_porcellus ###
sequence_x_order <- HylesE_vs_Deilephila %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- HylesE_vs_Deilephila %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()

sequence_y_order <- c("LR999998.1", "LR999997.1", "LR999996.1","LR999994.1","LR999995.1", "LR999993.1", "LR999992.1", "LR999991.1",
 "LR999990.1", "LR999989.1","LR999988.1", "LR999986.1", "LR999987.1","LR999985.1", "LR999982.1","LR999984.1",
 "LR999983.1", "LR999981.1", "LR999979.1","LR999977.1", "LR999976.1","LR999980.1","LR999978.1",
"LR999974.1","LR999973.1", "LR999975.1","LR999972.1", "LR999971.1", "LR999970.1")

sequences_to_flip <- c("LR999998.1","LR999996.1","LR999995.1", "LR999992.1", "LR999991.1", "LR999990.1",
                       "LR999987.1","LR999984.1", "LR999981.1", "LR999979.1","LR999978.1")
HylesE_vs_Deilephila$Sequence.x_f <- factor(HylesE_vs_Deilephila$Sequence.x, levels = sequence_x_order)
HylesE_vs_Deilephila$Sequence.y_f <- factor(HylesE_vs_Deilephila$Sequence.y, levels = sequence_y_order)

for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Deilephila %>% filter(Sequence.y == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Deilephila <- HylesE_vs_Deilephila %>%
    filter(Sequence.y != i) %>%
    bind_rows(subset_df)
}

HylesE_vs_Deilephila$color_f <- "grey"

# OPTIONAL: Just before plotting - replace Accession_IDs with relabelled chr name for H. euphorbiae
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr', 'Sequence.x')
HylesE_vs_Deilephila <- merge(HylesE_vs_Deilephila,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Deilephila <- HylesE_vs_Deilephila[,!names(HylesE_vs_Deilephila) %in% c("Sequence.x")]
names(HylesE_vs_Deilephila)[names(HylesE_vs_Deilephila) == 'Relabelled_chr.x'] <- 'Sequence.x'
HylesE_vs_Deilephila$Sequence.x_f <- HylesE_vs_Deilephila$Sequence.x

# Now need to reconfigure the order of x_sequences as defaults to alphabetical
sequence_x_order <- c("HeChr26","HeChr29","HeChr24","HeChr27","HeChr28","HeChr20","HeChr7","HeChr16","HeChr14",
                      "HeChr25","HeChr19","HeChr21",
                      "HeChr3","HeChr8","HeChr6","HeChr18",
                      "HeChr11","HeChr23","HeChr9","HeChr13","HeChr10","HeChr17","HeChr2","HeChr4",
                      "HeChr15","HeChr12","HeChr5","HeChr22","HeChr1=Z")
HylesE_vs_Deilephila$Sequence.x_f <- factor(HylesE_vs_Deilephila$Sequence.x, levels = sequence_x_order)


HylesE_vs_Deilephila_plot <- make_oxford_plot(HylesE_vs_Deilephila) + 
  xlab("Deilephila porcellus") + ylab("Hyles euphorbiae") + 
  theme(axis.title = element_text(face="italic"))
HylesE_vs_Deilephila_plot

### Hyles_euphoribae vs Manduca_sexta ###
sequence_x_order <- HylesE_vs_Manduca %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- HylesE_vs_Manduca %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()

# remove unscaffolded contigs
Manduca_sexta_unlocalised_contigs <- c("JACVES010001939.1", "JACVES010000056.1", "JACVES010000057.1")

HylesE_vs_Manduca <- HylesE_vs_Manduca[!HylesE_vs_Manduca$Sequence.y %in% Manduca_sexta_unlocalised_contigs, ]
sequence_y_order <- c("CM026232.1","CM026233.1","CM026234.1","CM026236.1","CM026235.1", "CM026237.1","CM026238.1","CM026240.1",
                     "CM026244.1", "CM026241.1","CM026247.1","CM026239.1","CM026245.1","CM026252.1","CM026242.1",
                     "CM026246.1","CM026250.1","CM026251.1","CM026253.1","CM026248.1","CM026243.1","CM026249.1",
                     "CM026254.1","CM026256.1","CM026255.1","CM026259.1","CM026257.1","CM026258.1")

sequences_to_flip <- c("CM026257.1","CM026258.1","CM026249.1","CM026246.1","CM026252.1","CM026239.1","CM026247.1", "CM026237.1","CM026238.1","CM026234.1")

HylesE_vs_Manduca$Sequence.x_f <- factor(HylesE_vs_Manduca$Sequence.x, levels = sequence_x_order)
HylesE_vs_Manduca$Sequence.y_f <- factor(HylesE_vs_Manduca$Sequence.y, levels = sequence_y_order)

for (i in sequences_to_flip){
  subset_df <- HylesE_vs_Manduca %>% filter(Sequence.y == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_Manduca <- HylesE_vs_Manduca %>%
    filter(Sequence.y != i) %>%
    bind_rows(subset_df)
}

HylesE_vs_Manduca$color_f <- "grey"

subset_df <- HylesE_vs_Manduca %>% filter(Sequence.y == "CM026259.1")
subset_df$color_f <- "red"
HylesE_vs_Manduca <- HylesE_vs_Manduca %>%
  filter(Sequence.y != "CM026259.1") %>%
  bind_rows(subset_df)

# OPTIONAL: Just before plotting - replace Accession_IDs with relabelled chr name for H. euphorbiae
colnames(Manduca_conversion) <- c('Relabelled_chr_y', 'Sequence.y')
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr_x', 'Sequence.x')
HylesE_vs_Manduca <- merge(HylesE_vs_Manduca,Manduca_conversion, by="Sequence.y")
HylesE_vs_Manduca <- merge(HylesE_vs_Manduca,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_Manduca <- HylesE_vs_Manduca[,!names(HylesE_vs_Manduca) %in% c("Sequence.x", "Sequence.y")]
names(HylesE_vs_Manduca)[names(HylesE_vs_Manduca) == 'Relabelled_chr_y'] <- 'Sequence.y'
names(HylesE_vs_Manduca)[names(HylesE_vs_Manduca) == 'Relabelled_chr_x'] <- 'Sequence.x'
HylesE_vs_Manduca$Sequence.x_f <- HylesE_vs_Manduca$Sequence.x
HylesE_vs_Manduca$Sequence.y_f <- HylesE_vs_Manduca$Sequence.y

sequence_x_order <- c("HeChr26","HeChr29","HeChr24","HeChr27","HeChr28","HeChr20","HeChr7","HeChr16","HeChr14",
                      "HeChr25","HeChr19","HeChr21",
                      "HeChr3","HeChr8","HeChr6","HeChr18",
                      "HeChr11","HeChr23","HeChr9","HeChr13","HeChr10","HeChr17","HeChr2","HeChr4",
                      "HeChr15","HeChr12","HeChr5","HeChr22","HeChr1=Z")

sequence_y_order <- c("MsChr26","MsChr24", "MsChr27", "MsChr28","MsChr20","MsChr7",
                      "MsChr16","MsChr14","MsChr25","MsChr19","MsChr21","MsChr3","MsChr8" ,"MsChr6",
                      "MsChr18","MsChr11","MsChr23","MsChr9","MsChr13","MsChr10",
                      "MsChr17","MsChr2","MsChr4","MsChr15","MsChr12","MsChr5","MsChr22","MsChr1=Z"
)


HylesE_vs_Manduca$Sequence.x_f <- factor(HylesE_vs_Manduca$Sequence.x, levels = sequence_x_order)
HylesE_vs_Manduca$Sequence.y_f <- factor(HylesE_vs_Manduca$Sequence.y, levels = sequence_y_order)

# Plot!
HylesE_vs_Manduca_plot <- make_oxford_plot(HylesE_vs_Manduca) + 
  xlab("Manduca sexta") + ylab("Hyles euphorbiae") + 
  theme(axis.title = element_text(face="italic"))

HylesE_vs_Manduca_plot

### Hyles_euphoribae vs Hyles_vespertilio ###
sequence_x_order <- HylesE_vs_HylesV %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- HylesE_vs_HylesV %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()

sequence_y_order <- c("CM042829.1","CM042832.1","CM042827.1","CM042830.1","CM042831.1" ,"CM042823.1","CM042810.1","CM042819.1" , 
                      "CM042817.1","CM042828.1","CM042822.1", "CM042824.1","CM042806.1"
                      ,"CM042811.1", "CM042809.1", "CM042821.1","CM042814.1","CM042826.1", "CM042812.1", "CM042816.1","CM042813.1","CM042820.1",
                      "CM042805.1", "CM042807.1","CM042818.1","CM042815.1" ,"CM042808.1", "CM042825.1","CM042833.1")

sequences_to_flip <- c("CM042829.1", "CM042832.1", "CM042827.1",  "CM042823.1", "CM042817.1", 
                       "CM042822.1", "CM042824.1","CM042806.1", "CM042821.1", "CM042812.1", "CM042816.1", "CM042807.1", "CM042825.1")
HylesE_vs_HylesV$Sequence.x_f <- factor(HylesE_vs_HylesV$Sequence.x, levels = sequence_x_order)
HylesE_vs_HylesV$Sequence.y_f <- factor(HylesE_vs_HylesV$Sequence.y, levels = sequence_y_order)
HylesE_vs_HylesV$color_f <- "grey"
for (i in sequences_to_flip){
  subset_df <- HylesE_vs_HylesV %>% filter(Sequence.y == i)
  max_pos <- max(subset_df$midpos_x)
  subset_df$midpos_x <- (subset_df$midpos_x -max_pos)*-1
  HylesE_vs_HylesV <- HylesE_vs_HylesV %>%
    filter(Sequence.y != i) %>%
    bind_rows(subset_df)
}

# Optional - just before plotting, replace NCBI accession numbers with relablled chr names for H.euphoribae vs H. vespertilio 
colnames(Hyles_vespertilio_conversion) <- c('Relabelled_chr_y', 'Sequence.y')
colnames(Hyles_euphoribae_conversion) <- c('Relabelled_chr_x', 'Sequence.x')

HylesE_vs_HylesV <- merge(HylesE_vs_HylesV,Hyles_euphoribae_conversion, by="Sequence.x")
HylesE_vs_HylesV <- merge(HylesE_vs_HylesV,Hyles_vespertilio_conversion, by="Sequence.y")
HylesE_vs_HylesV <- HylesE_vs_HylesV[,!names(HylesE_vs_HylesV) %in% c("Sequence.x", "Sequence.y")]
names(HylesE_vs_HylesV)[names(HylesE_vs_HylesV) == 'Relabelled_chr_x'] <- 'Sequence.x'
names(HylesE_vs_HylesV)[names(HylesE_vs_HylesV) == 'Relabelled_chr_y'] <- 'Sequence.y'
HylesE_vs_HylesV$Sequence.x_f <- HylesE_vs_HylesV$Sequence.x
HylesE_vs_HylesV$Sequence.y_f <- HylesE_vs_HylesV$Sequence.y

sequence_x_order <- c("HeChr26","HeChr29","HeChr24","HeChr27","HeChr28","HeChr20","HeChr7","HeChr16","HeChr14",
                      "HeChr25","HeChr19","HeChr21",
                      "HeChr3","HeChr8","HeChr6","HeChr18",
                      "HeChr11","HeChr23","HeChr9","HeChr13","HeChr10","HeChr17","HeChr2","HeChr4",
                      "HeChr15","HeChr12","HeChr5","HeChr22","HeChr1=Z")

sequence_y_order <-  c("HvChr26**","HvChr29§","HvChr24","HvChr27", "HvChr28#","HvChr20","HvChr7",
                       "HvChr16","HvChr14","HvChr25","HvChr19","HvChr21","HvChr3","HvChr8","HvChr6",
                       "HvChr18","HvChr11","HvChr23","HvChr9","HvChr13","HvChr10","HvChr17","HvChr2*","HvChr4","HvChr15","HvChr12","HvChr5","HvChr22#","HvChr1=Z")
HylesE_vs_HylesV$Sequence.x_f <- factor(HylesE_vs_HylesV$Sequence.x, levels = sequence_x_order)
HylesE_vs_HylesV$Sequence.y_f <- factor(HylesE_vs_HylesV$Sequence.y, levels = sequence_y_order)



HylesE_vs_HylesV_plot <- make_oxford_plot(HylesE_vs_HylesV) + xlab("Hyles vespertilio") + ylab("Hyles euphorbiae") + 
  theme(axis.title = element_text(face="italic"))
HylesE_vs_HylesV_plot

HylesE_vs_Laothoe_plot_noLabels <- HylesE_vs_Laothoe_plot + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
HylesE_vs_HylesV_plot_noLabels <- HylesE_vs_HylesV_plot + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
HylesE_vs_Manduca_plot_noLabels <- HylesE_vs_Manduca_plot + theme(strip.text.x = element_blank(), strip.text.y = element_blank())
HylesE_vs_Deilephila_plot_noLabels <- HylesE_vs_Deilephila_plot + theme(strip.text.x = element_blank(), strip.text.y = element_blank())

### Plot together ###
#HylesE_vs_HylesV_plot_legend <- HylesE_vs_HylesV_plot + theme(legend.position="left") +  
  #guides(color = guide_legend(override.aes = list(size = 5), 
  #                            title=expression(paste("Chromosome in ",italic(" H. euphorbiae")))))

# Extract the legend -  returns a gtable
#HylesE_legend <- get_legend(HylesE_vs_HylesV_plot_legend) %>% as_ggplot()  # Convert to a ggplot and print
#HylesE_legend

combined_plot <- ggarrange(HylesE_vs_HylesV_plot,HylesE_vs_Deilephila_plot, HylesE_vs_Manduca_plot,HylesE_vs_Laothoe_plot, labels=c("A","B","C","D"))
combined_plot_no_labels <- ggarrange(HylesE_vs_HylesV_plot_noLabels, HylesE_vs_Deilephila_plot_noLabels,HylesE_vs_Manduca_plot_noLabels,HylesE_vs_Laothoe_plot_noLabels)

ggsave(plot=combined_plot, 'Oxford_plots_redGrey_relablled_chr_130323.pdf', device='pdf', width = 9, height = 9, dpi = 300, units = "in")
ggsave(plot=combined_plot_no_labels, 'Oxford_plots_noFacetLabels_redGrey_relabelled_chr_130323.pdf', device='pdf', width = 9, height = 9, dpi = 300, units = "in")
