# With this scripts we will make some plots describing publication bias in analyzed orders and barriers
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)

# Text size for both graphics
s <- 17

# Studied orders bias ####
# Reading tables of species used at each analyses
Hybrids <- read.delim("../../02_Genetic_Distances/figures/01_Hybrids_Cases_Table.tsv")
Barriers <- read.delim("../../05_Barriers_Calculation/figures/01_Merged_Table.tsv")

# Tidying tables
Hybrids <- gather(Hybrids, "Number", "Species", 3:4)
Hybrids <- Hybrids[,c(1,6)]

Barriers <- gather(Barriers, "Number","Species", 2:3)
Barriers <- Barriers[,c(1,9)]

# Removing duplicates
Hybrids <- unique(Hybrids)
Barriers <- unique(Barriers)

# Tabling
Hybrids <- as.data.frame(table(Hybrids$Order))
Barriers <- as.data.frame(table(Barriers$Order))

# Binding
Hybrids$Data <- "Hybrids"
Barriers$Data <- "Barriers"
Freqs <- rbind(Hybrids, Barriers)
rm(Barriers, Hybrids)
colnames(Freqs)[1] <- "Order"

# Now adding GBIF infotmation
Freqs <- rbind(Freqs, read.delim("../data/GBIF.tsv"))

# Factoring
Freqs$Data <- factor(Freqs$Data, levels = c("Species", "Hybrids", "Barriers"), labels = c("Total number\nof species", "Species reported\nto hybridize", "Species in reproductive\nbarriers studies"))

# Plotting
colors <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
p1 <- ggplot(Freqs) +
  geom_bar(aes(x=Data, y=Freq, fill=Order), position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) +
  labs(y="Relative Frequency") +
  theme_classic() +
  theme(text = element_text(family = "serif", size = s),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 20, lineheight = 0.75),
        legend.position =  c(1.01,.58),
        legend.spacing.y = unit(1.8, 'cm'),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.background = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))

# Reproductive barriers bias ####
Barriers <- read.delim("../../05_Barriers_Calculation/figures/01_Merged_Table.tsv")
Barriers <- Barriers[,c(1,4:8)]
#Barriers <- Barriers[Barriers$Order!="Odonata",]

# Splitting Barriers into orders
orders <- list()
for (i in unique(Barriers$Order)) {
  orders[[i]] <- Barriers[Barriers$Order==i,-1]
  orders[[i]][orders[[i]]=="?"] <- 0
  orders[[i]][orders[[i]]!="0"] <- 1
  orders[[i]] <- apply(orders[[i]], 2, as.numeric)
  orders[[i]] <- as.data.frame(colSums(orders[[i]]))
  orders[[i]]$Barrier <- row.names(orders[[i]])
}

# Binding rows
orders <- bind_rows(orders, .id="Order")
row.names(orders) <- 1:nrow(orders)
orders <- orders[,c(1,3,2)]
colnames(orders)[3] <- "Fx"

# Factoring
orders$Barrier <- factor(orders$Barrier, levels = rev(unique(orders$Barrier)), labels = rev(c("Assortative Mating +\nMechanical Isolation",
                                                                                     "Gametic or\nTactile Barriers",
                                                                                     "Hybrid\nInviability",
                                                                                     "Hybrid\nSterility",
                                                                                     "Partial Hybrid\nSterility")))

# Little changes:
orders <- orders[complete.cases(orders),]
orders[orders$Order=="Hymenoptera" & orders$Barrier=="Partial Hybrid\nSterility","Barrier"] <- "Hybrid\nSterility"

# Plotting
p2 <- ggplot(orders) +
  geom_bar(aes(y=Order, x=Fx, fill=Barrier), position = "fill", stat = "identity") +
  scale_fill_manual(values = c('#f6eff7','#bdc9e1','#67a9cf','#1c9099','#016c59'), guide = guide_legend(reverse = T)) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(labels = scales::percent, n.breaks = 6) +
  theme_classic() +
  labs(x="Relative Frequency") +
  theme(text = element_text(family = "serif", size = s),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid.major.x = element_line(color = "gray"),
        legend.position =  "bottom")

# Merging plots
png("../figures/01_Bias.png", width = 12.5984, height = 12.5984/2, units = "in", res = 300)
grid.arrange(p1, p2, nrow=1, widths = c(1, 2))
dev.off()

pdf("../figures/01_Bias.pdf", width = 12.5984, height = 12.5984/2)
grid.arrange(p1, p2, nrow=1, widths = c(1, 2))
dev.off()