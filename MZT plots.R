# Juliet M Wong
# Wong et al. 2019 - Transcriptional profiles of early stage red sea urchins (Mesocentrotus
#franciscanus) reveal differential regulation of gene expression across
#development

install.packages("ggplot2")
library(ggplot2)
install.packages("viridis")
library(viridis)

# read in data for maternal degradation genes
smaug <- read.table(file = "smaug_data.txt", header = TRUE)
# create dummy matrix to hold info for reordering the categorical x axis, stage
foo = rep(0, nrow(smaug))
# define the order
foo[with(smaug, Stage == "EG")] = 1
foo[with(smaug, Stage == "CL")] = 2
foo[with(smaug, Stage == "MO")] = 3
foo[with(smaug, Stage == "BL")] = 4
foo[with(smaug, Stage == "GA")] = 5
foo[with(smaug, Stage == "PR")] = 6
foo[with(smaug, Stage == "PL")] = 7
# set the order
smaug$Stage = with(smaug, reorder(Stage, foo))
smaugtrace <- ggplot(data=smaug, aes(x=Stage, y=Count, group=Gene, colour=Gene))+geom_line()+
  labs(x="Stage", y="log2-counts per million (logCPM)") +
  scale_color_manual(values=cbPalette) +
  theme_classic()

# read in data for zygotic genes
zygotic <- read.table(file = "zygotic_data.txt", header = TRUE)
# create dummy matrix to hold info for reordering the categorical x axis, stage
zoo = rep(0, nrow(zygotic))
# define the order
zoo[with(zygotic, Stage == "EG")] = 1
zoo[with(zygotic, Stage == "CL")] = 2
zoo[with(zygotic, Stage == "MO")] = 3
zoo[with(zygotic, Stage == "BL")] = 4
zoo[with(zygotic, Stage == "GA")] = 5
zoo[with(zygotic, Stage == "PR")] = 6
zoo[with(zygotic, Stage == "PL")] = 7
# set the order
zygotic$Stage = with(zygotic, reorder(Stage, zoo))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
zygotictrace <- ggplot(data=zygotic, aes(x=Stage, y=Count, group=Gene, colour=Gene))+geom_line()+
  labs(x="Stage", y="log2-counts per million (logCPM)") +
  scale_color_manual(values=cbPalette) +
  theme_classic()

library("ggpubr")
ggarrange(smaugtrace, zygotictrace, labels = c("A", "B"), ncol = 1, nrow =2, align = "v")
    