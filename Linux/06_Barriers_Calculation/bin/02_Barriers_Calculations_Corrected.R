#!/usr/bin/Rscript
# 01_Barriers_Calculations.R
# Made by Luis Rodrigo Arce Valdés, to estimate RI using excel files and correlate it with genetic distance
rm(list = ls())

# Calling libraries
library(dplyr)
library(ggplot2)
library(tidyr)


# Four original orders ####

# Groups vector
groups <- c("Orthoptera", "Diptera", "Hymenoptera", "Lepidoptera")

# Barriers per group
barriers <- list()
for (i in groups) {
  barriers[[i]] <- read.table(paste0("../../05_Barriers_vs_Distance/data/01_raw/",i,"_Barriers.tsv"), sep = "\t", header = T)
}

# While estimating reproductive barriers we will do a different estimation for Hymenoptera (because they create hybrids of only one sex)
groups <- groups[groups!="Hymenoptera"]

# Estimating reproductive barriers
for (i in groups) {
  # Adding a column pasting the strings of each row
  barriers[[i]]$RI <- paste0(barriers[[i]]$Choice.Mating,
                             barriers[[i]]$Hybrid.Egg,
                             barriers[[i]]$Hybrid.Adult,
                             barriers[[i]]$One.Sex.Fertile.Hybrids,
                             barriers[[i]]$Fertile.Hybrids)
  # Looking for the last barrier each species can do gene flow
  RI <- vector()
  for (n in 1:nrow(barriers[[i]])) {
    pos <- unlist(gregexpr('N', barriers[[i]]$RI[n]))[1]
    RI <- append(RI, pos)
  }
  barriers[[i]]$pos <- RI
  # Removing -1s that have "?"
  barriers[[i]][barriers[[i]]$RI=="YYYYY","pos"] <- 6
  barriers[[i]] <- barriers[[i]][barriers[[i]]$pos!=-1,]
  # Estimating reproductive isolation barriers
  barriers[[i]]$RI <- 1-((barriers[[i]]$pos-1)/5)
  barriers[[i]] <- barriers[[i]][,-10]
}


# Now for hymenopterans with a small change (RI divided by 4 total barriers)
for (i in "Hymenoptera") {
  # Adding a column pasting the strings of each row
  barriers[[i]]$RI <- paste0(barriers[[i]]$Choice.Mating,
                             barriers[[i]]$Hybrid.Egg,
                             barriers[[i]]$Hybrid.Adult,
                             barriers[[i]]$One.Sex.Fertile.Hybrids)
  # Looking for the last barrier each species can do gene flow
  RI <- vector()
  for (n in 1:nrow(barriers[[i]])) {
    pos <- unlist(gregexpr('N', barriers[[i]]$RI[n]))[1]
    RI <- append(RI, pos)
  }
  barriers[[i]]$pos <- RI
  # Removing -1s that have "?"
  barriers[[i]][barriers[[i]]$RI=="YYYY","pos"] <- 5
  barriers[[i]] <- barriers[[i]][barriers[[i]]$pos!=-1,]
  # Estimating reproductive isolation barriers
  barriers[[i]]$RI <- 1-((barriers[[i]]$pos-1)/4)
  barriers[[i]] <- barriers[[i]][,-10]
}

# Creating groups vector
groups <- c(groups, "Hymenoptera")

# Adding column name
for (i in groups) {
  barriers[[i]]$Order <- i
}

# Merging dataframes
isolation <- bind_rows(barriers)

# Ordering columns
isolation <- isolation[,c(10,1,2,3,9)]

# Reading genetic distance data frame
COI <- read.delim("../../05_Barriers_vs_Distance/results/01_Genetic_Distances.tsv")

# Editing isolation species names
isolation$Cross <- paste0(isolation$Sp1,"_X_",isolation$Sp2)
isolation$Cross <- gsub(" ","_", isolation$Cross)
isolation <- isolation[,c(1,6,4,5)]

# Creating an additional dataframe for later
orders <- isolation

# Subseting odonates COIs
COI.odonata <- COI[COI$Order=="Odonata",]

# Merging isolation barriers and genetic distance
COI <- COI[match(intersect(COI$Cross, isolation$Cross),COI$Cross),]
isolation <- isolation[match(intersect(COI$Cross, isolation$Cross), isolation$Cross),]

# Checking if vectors are identical
identical(isolation$Cross, COI$Cross)

# Merging dataframes
isolation <- cbind(isolation, Gen=COI$Distance)
rm(COI)

# Phylogenetic correction
# First we will remove repeated comparisons
isolation <- separate(isolation, col=Cross, into=c("Sp1","Sp2"), sep = "_X_")

sp.order <- data.frame()
# Sorting species
for (i in 1:nrow(isolation)) {
  sp.order[i,c(1,2)] <- sort(as.character(isolation[i,c(2,3)]))
}

# Adding columns to isolation
isolation[,c(2,3)] <- sp.order[,c(1,2)]

# Removing duplicates
isolation <- unique(isolation)

# Ordering per species
isolation <- isolation[order(isolation$Sp1, isolation$Sp2),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()
order <- vector()

# For looping
for (i in unique(isolation$Sp1)) {
  ri <- mean(as.numeric(isolation[isolation$Sp1==i,"RI"]))
  gen <- mean(as.numeric(isolation[isolation$Sp1==i,"Gen"]))
  if (length(isolation[isolation$Sp1==i,4])==1) {
    s <- isolation[isolation$Sp1==i,3]
  }
  else {
    s <- paste0("AAA",i,"_average")
  }
  RI <- append(RI, ri)
  order <- append(order, isolation[isolation$Sp1==i,1][1])
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
isolation <- data.frame(Order=order,
                             Sp1=unique(isolation$Sp1),
                             Sp2=sp2,
                             RI=RI,
                             Gen=Gen)
rm(gen, Gen, i, n, order, pos, ri, RI, s, sp2)

# Sorting by second species
isolation <- isolation[order(isolation$Sp2, isolation$Sp1),c(1,3,2,4,5)]

# And repeating averages
sp.order <- data.frame()
for (i in 1:nrow(isolation)) {
  sp.order[i,c(1,2)] <- sort(as.character(isolation[i,c(2,3)]))
}

# Adding columns to isolation
isolation[,c(3,2)] <- sp.order[,c(1,2)]

# Ordering per species
isolation <- isolation[order(isolation$Sp2, isolation$Sp1),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()
order <- vector()

# For looping
for (i in unique(isolation$Sp2)) {
  ri <- mean(as.numeric(isolation[isolation$Sp2==i,"RI"]))
  gen <- mean(as.numeric(isolation[isolation$Sp2==i,"Gen"]))
  if (length(isolation[isolation$Sp2==i,4])==1) {
    s <- isolation[isolation$Sp2==i,3]
  }
  else {
    s <- "Average"
  }
  RI <- append(RI, ri)
  order <- append(order, isolation[isolation$Sp2==i,1][1])
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
isolation <- data.frame(Order=order,
                        Sp1=unique(isolation$Sp2),
                        Sp2=sp2,
                        RI=RI,
                        Gen=Gen)
rm(gen, Gen, i, order, ri, RI, s, sp2)

# Gsubing
isolation$Sp2 <- sub("^AAA.*","Average", isolation$Sp2)

# Function to get p value from linnear regressions
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# linear regressions
r2 <- vector()
p <- vector()
for (i in groups){
  x <- lm(formula = RI ~ Gen, data = isolation[isolation$Order==i,])
  p <- append(p, lmp(x))
  r2 <- append(r2, summary(x)$r.squared)
}

# Dataframing
regs <- data.frame(Order=groups, r2=r2, p=p)
rm(i,x,p,r2)

# Editing
regs$label <- paste0("r2 = ", round(regs$r2, 3), "; p = ", round(regs$p, 3))

# Plotting
dir.create("../figures_corrected", showWarnings = F)
png("../figures_corrected/01_Scatterplot.png", width = 18, height = 10, units = "cm", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, ncol = 2) +
  geom_smooth(aes(x=Gen, y=RI), formula = y~x, method = "lm", se = F) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  geom_label(data = regs, aes(x=0.12, y=1.15, label=label), size=2, family="serif") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "serif"),
        plot.margin = margin(5,10,5,5, "points"))
dev.off()

# Now transforming into categorical variables
# For non hymenopterans we will round to its closest 0.20
isolation[isolation$Order!="Hymenoptera","RI"] <- (2*round(isolation[isolation$Order!="Hymenoptera","RI"]*10/2))/10

# Now for hymenopterans we will round to the closest 0.25
isolation[isolation$Order=="Hymenoptera","RI"] <- (2.5*round(isolation[isolation$Order=="Hymenoptera","RI"]*10/2.5))/10

# First we will unify for hymenopterans
isolation[isolation==0.25] <- 0.40
isolation[isolation==0.50] <- 0.60
isolation[isolation==0.75] <- 0.80

# Now transforming ri into a factor
isolation <- isolation[order(isolation$RI, decreasing = T),]
isolation$RI <- factor(isolation$RI, levels = unique(isolation$RI),
                       labels = c("Choice_Mating", "Oviposition", "Hybrid_Survival",
                                  "Hybrid_Fertility","Haldanes_Rule","Complete_Gene_Flow"))
# Plotting
png("../figures_corrected/02_Barriers.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, ncol = 2) +
  geom_boxplot(aes(x=RI, y=Gen, color=Order)) +
  geom_point(aes(x=RI, y=Gen, color=Order)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 5))
dev.off()

# Statistical testing
sink("../figures_corrected/02_Barriers_testing.txt", append = F, split = T)
for(i in groups) {
  print("###########################")
  print(i)
  print(kruskal.test(Gen ~ RI, data = isolation[isolation$Order==i,]))
  for(c in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(isolation[isolation$Order==i,"Gen"], isolation[isolation$Order==i,"RI"], p.adjust.method = c))
  }
}
sink()

# Odonates ####

# Now we will include the odonate species
# Reading table of isolation in odonates
odonata <- read.delim("../../05_Barriers_vs_Distance/data/01_raw/Odonata_Barriers.tsv")

# Replacing "-1" by "0"
odonata[odonata==-1] <- 0

# Estimating RI
odonata$RI <- 1-(odonata$RI/6)

# Editing isolation species names
odonata$Cross <- paste0(odonata$Male,"_X_",odonata$Female)
odonata$Cross <- gsub(" ","_", odonata$Cross)
odonata.RI <- odonata[,c(10,9)]

# Creating additional df for later
odonata.extra <- odonata.RI

# Merging odonata barriers and genetic distance
COI.odonata <- COI.odonata[match(intersect(COI.odonata$Cross, odonata.RI$Cross),COI.odonata$Cross),]
odonata.RI <- odonata.RI[match(intersect(COI.odonata$Cross, odonata.RI$Cross), odonata.RI$Cross),]

# Checking if vectors are identical
identical(odonata.RI$Cross, COI.odonata$Cross)

# Merging dataframes
odonata.RI <- cbind(odonata.RI, Gen=COI.odonata$Distance)

# # Plotting
# png("../figures_corrected/03_Odonata_Scatterplot.png", width = 1800, height = 1000, units = "px", res = 300)
# ggplot(odonata.RI) +
#   geom_point(aes(x=Gen, y=RI)) +
#   theme_classic()
# dev.off()

# Now transforming ri into a factor
odonata.RI <- odonata.RI[order(odonata.RI$RI, decreasing = T),]
odonata.RI$RI <- factor(odonata.RI$RI, levels = unique(odonata.RI$RI),
                       labels = c("Sexual", "TandemAttempt", "Tandem",
                                      "Mating","Oviposition","Hybrid","CompleteGeneFlow"))

# # Plotting
# png("../figures_corrected/04_Odonata_Barriers.png", width = 2400, height = 1200, units = "px", res = 300)
# ggplot(odonata.RI) +
#   geom_boxplot(aes(x=RI, y=Gen)) +
#   geom_point(aes(x=RI, y=Gen)) +
#   geom_point(data = odonata.RI[odonata.RI$Cross=="Ischnura_elegans_X_Ischnura_graellsii" | odonata.RI$Cross=="Ischnura_graellsii_X_Ischnura_elegans",], aes(x=RI, y=Gen, color=Cross)) +
#   labs(caption = "Se retendrá la cruza recíproca que más anvance en la reproducción. Por ejemplo, Macho elegans con Hembra graellsii") +
#   theme_classic() +
#   theme(legend.position = "bottom")
# dev.off()


# Removing reciprocal crosses in Odonata ####

# To remove reciprocal crosses we will sort cross column alphabetically
odonata$code <- vapply(odonata$Cross, function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=''), '')

# Sorting first by code then by isolation
odonata <- odonata[order(odonata$code, odonata$RI),]

# Removing duplicated codes (species pairs)
odonata <- odonata[!duplicated(odonata$code),]

# Repeating analyses
odonata.RI <- odonata[,c(10,9)]

# Merging odonata barriers and genetic distance
COI.odonata <- COI.odonata[match(intersect(COI.odonata$Cross, odonata.RI$Cross),COI.odonata$Cross),]
odonata.RI <- odonata.RI[match(intersect(COI.odonata$Cross, odonata.RI$Cross), odonata.RI$Cross),]

# Checking if vectors are identical
identical(odonata.RI$Cross, COI.odonata$Cross)

# Merging dataframes
odonata.RI <- cbind(odonata.RI, Gen=COI.odonata$Distance)

# Phylogenetic correction in odonata ####
# First we will remove repeated comparisons
odonata.RI <- separate(odonata.RI, col=Cross, into=c("Sp1","Sp2"), sep = "_X_")

sp.order <- data.frame()
# Sorting species
for (i in 1:nrow(odonata.RI)) {
  sp.order[i,c(1,2)] <- sort(as.character(odonata.RI[i,c(1,2)]))
}

# Adding columns to odonata.RI
odonata.RI[,c(1,2)] <- sp.order[,c(1,2)]

# Removing duplicates
odonata.RI <- unique(odonata.RI)

# Ordering per species
odonata.RI <- odonata.RI[order(odonata.RI$Sp1, odonata.RI$Sp2),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()

# For looping
for (i in unique(odonata.RI$Sp1)) {
  ri <- mean(as.numeric(odonata.RI[odonata.RI$Sp1==i,"RI"]))
  gen <- mean(as.numeric(odonata.RI[odonata.RI$Sp1==i,"Gen"]))
  if (length(odonata.RI[odonata.RI$Sp1==i,"RI"])==1) {
    s <- odonata.RI[odonata.RI$Sp1==i,2]
  }
  else {
    s <- paste0("AAA",i,"_average")
  }
  RI <- append(RI, ri)
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
odonata.RI <- data.frame(Sp1=unique(odonata.RI$Sp1),
                        Sp2=sp2,
                        RI=RI,
                        Gen=Gen)
rm(gen, Gen, i, ri, RI, s, sp2)

# Sorting by second species
odonata.RI <- odonata.RI[order(odonata.RI$Sp2, odonata.RI$Sp1),c(2,1,3,4)]

# And repeating averages
sp.order <- data.frame()
for (i in 1:nrow(odonata.RI)) {
  sp.order[i,c(1,2)] <- sort(as.character(odonata.RI[i,c(1,2)]))
}

# Adding columns to odonata.RI
odonata.RI[,c(2,1)] <- sp.order[,c(1,2)]

# Ordering per species
odonata.RI <- odonata.RI[order(odonata.RI$Sp2, odonata.RI$Sp1),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()

# For looping
for (i in unique(odonata.RI$Sp2)) {
  ri <- mean(as.numeric(odonata.RI[odonata.RI$Sp2==i,"RI"]))
  gen <- mean(as.numeric(odonata.RI[odonata.RI$Sp2==i,"Gen"]))
  if (length(odonata.RI[odonata.RI$Sp2==i,4])==1) {
    s <- odonata.RI[odonata.RI$Sp2==i,2]
  }
  else {
    s <- "Average"
  }
  RI <- append(RI, ri)
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
odonata.RI <- data.frame(Sp1=unique(odonata.RI$Sp2),
                        Sp2=sp2,
                        RI=RI,
                        Gen=Gen)
rm(gen, Gen, i, ri, RI, s, sp2)

# Gsubing
odonata.RI$Sp2 <- sub("^AAA.*","Average", odonata.RI$Sp2)


# Creating regression results
x <- lm(RI~Gen, odonata.RI)
reg <- paste0("r2 = ", round(summary(x)$r.squared, 3), "; p = ", round(lmp(x), 3))

# Plotting
png("../figures_corrected/03_Odonata_Scatterplot_1D.png", width = 18, height = 10, units = "cm", res = 300)
ggplot(odonata.RI) +
  geom_smooth(aes(x=Gen, y=RI), formula = y~x, method = "lm", se = F) +
  geom_point(aes(x=Gen, y=RI)) +
  geom_label(data = regs, aes(x=0.25, y=0.90, label=reg), size=2, family="serif") +
  theme_classic()
dev.off()

# Now transforming ri into a factor
odonata.RI <- odonata.RI[order(odonata.RI$RI, decreasing = T),]

# Rounding RI to its closest 0.20
# For non hymenopterans we will round to its closest 0.20
odonata.RI$RI <- (2*round(odonata.RI$RI*10/2))/10

# Naming barriers
odonata.RI[odonata.RI$RI==0,"RI"] <- "CompleteGeneFlow"
odonata.RI[odonata.RI$RI=="0.2","RI"] <- "Hybrid"
odonata.RI[odonata.RI$RI=="0.4","RI"] <- "Oviposition"
odonata.RI[odonata.RI$RI=="0.6","RI"] <- "Mating"
odonata.RI[odonata.RI$RI=="0.8","RI"] <- "Tandem"
odonata.RI[odonata.RI$RI=="1","RI"] <- "TandemAttempt"

# factoring
odonata.RI$RI <- factor(odonata.RI$RI, levels = unique(odonata.RI$RI))

# Plotting
png("../figures_corrected/04_Odonata_Barriers_1D.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(odonata.RI) +
  geom_boxplot(aes(x=RI, y=Gen)) +
  geom_point(aes(x=RI, y=Gen)) +
  geom_point(data = odonata.RI[odonata.RI$Cross=="Ischnura_elegans_X_Ischnura_graellsii",], aes(x=RI, y=Gen), color="red") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

# Statistical testing
sink("../figures_corrected/06_Odonata_Barriers_1D_testing.txt", append = F, split = T)
print("Odonata: ")
print(kruskal.test(Gen ~ RI, data = odonata.RI))
  for(c in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(odonata.RI$Gen, odonata.RI$RI, p.adjust.method = c))
  }
sink()

# Merging Odonates and other insects ####
odonata.RI$Order <- "Odonata"
odonata.RI <- odonata.RI[,c(5,1:4)]

# Binding
isolation <- rbind(isolation, odonata.RI)

# Now unifing names of barriers
barriers <- data.frame(Original=unique(isolation$RI), New=c("Mating","Oviposition","Hybrid",
                                                "Complete_Gene_Flow","Complete_Gene_Flow","Complete_Gene_Flow",
                                                "Mating","Mating",
                                                "Hybrid","Complete_Gene_Flow"))

# Viewing barriers
barriers

# Factoring
isolation$RI <- factor(isolation$RI, levels = barriers$Original, labels = barriers$New)

# Plotting
# Plotting
png("../figures_corrected/07_Barriers_Merged.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, ncol = 2) +
  geom_boxplot(aes(x=RI, y=Gen, color=Order)) +
  geom_point(aes(x=RI, y=Gen, color=Order)) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()

# Statistical testing
sink("../figures_corrected/07_Barriers_Merged_testing.txt", append = F, split = T)
for(i in unique(isolation$Order)) {
  print("###########################")
  print(i)
  print(kruskal.test(Gen ~ RI, data = isolation[isolation$Order==i,]))
  for(c in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(isolation[isolation$Order==i,"Gen"], isolation[isolation$Order==i,"RI"], p.adjust.method = c))
  }
}
sink()

# Now scaterploting
isolation$RI <- factor(isolation$RI, levels = unique(isolation$RI), labels = c("1","0.666","0.333","0"))
isolation$RI <- as.numeric(as.vector(isolation$RI))

# linear regressions
r2 <- vector()
p <- vector()
for (i in unique(isolation$Order)){
  x <- lm(formula = RI ~ Gen, data = isolation[isolation$Order==i,])
  p <- append(p, lmp(x))
  r2 <- append(r2, summary(x)$r.squared)
}

# Dataframing
regs <- data.frame(Order=unique(isolation$Order), r2=r2, p=p)
rm(i,x,p,r2)

# Editing
regs$label <- paste0("r2 = ", round(regs$r2, 3), "; p = ", round(regs$p, 3))

# Plotting
png("../figures_corrected/08_Scatterplot_Merged.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, ncol = 2) +
  geom_smooth(aes(x=Gen, y=RI), formula = y~x, method = "lm", se = F) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  geom_label(data = regs, aes(x=0.225, y=0.10, label=label), size=2, family="serif") +
  theme_classic()
dev.off()

# Final plots (without COIs) ####
isolation <- orders
rm(orders)

# First we will unify for hymenopterans
isolation[isolation==0.25] <- 0.40
isolation[isolation==0.50] <- 0.60
isolation[isolation==0.75] <- 0.80

# Now transforming ri into a factor
isolation <- isolation[order(isolation$RI, decreasing = T),]
isolation$RI <- factor(isolation$RI, levels = unique(round(isolation$RI,1)),
                       labels = c("Choice_Mating", "Oviposition", "Hybrid_Survival",
                                  "Hybrid_Fertility","Haldanes_Rule","Complete_Gene_Flow"))

# Creating table
isolation.table <- as.data.frame(prop.table(table(isolation$RI, isolation$Order), margin = 2))
colnames(isolation.table)[1:2] <- c("Barrier","Order")

# Plotting
png("../figures_corrected/09_BarI.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  theme_classic() +
  theme(text = element_text(family = "serif"))
dev.off()

# Now without "complete gene flow"
isolation.table <- isolation[isolation$RI!="Complete_Gene_Flow",]  

# Creating table
isolation.table <- as.data.frame(prop.table(table(isolation.table$RI, isolation.table$Order), margin = 2))
colnames(isolation.table)[1:2] <- c("Barrier","Order")

# Plotting
png("../figures_corrected/10_BarII.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  theme_classic() +
  theme(text = element_text(family = "serif"))
dev.off()

# Now with Odonates
# Now transforming ri into a factor
odonata.extra <- odonata.extra[order(odonata.extra$RI, decreasing = T),]
odonata.extra$RI <- factor(odonata.extra$RI, levels = unique(odonata.extra$RI),
                        labels = c("Sexual", "TandemAttempt", "Tandem",
                                   "Mating","Oviposition","Hybrid","CompleteGeneFlow"))

# Creating table
isolation.table <- as.data.frame(prop.table(table(odonata.extra$RI)))
colnames(isolation.table)[1] <- c("Barrier")
isolation.table$Order <- "Odonata"

# Plotting
png("../figures_corrected/11_BarIII.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  theme_classic() +
  theme(text = element_text(family = "serif"))
dev.off()

# Binding
odonata.extra$Order <- "Odonata"
odonata.extra <- odonata.extra[,c(3,1,2)]
isolation <- rbind(isolation[,-3], odonata.extra)

# Viewing barriers
barriers

# Factoring
isolation$RI <- factor(isolation$RI, levels = barriers$Original, labels = barriers$New)

# Creating table
isolation.table <- as.data.frame(prop.table(table(isolation$RI, isolation$Order), margin = 2))
colnames(isolation.table)[1:2] <- c("Barrier","Order")

# Plotting
png("../figures_corrected/12_BarIV.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  theme_classic() +
  theme(text = element_text(family = "serif"))
dev.off()

# Now without "complete gene flow"
isolation.table <- isolation[isolation$RI!="Complete_Gene_Flow",]  

# Creating table
isolation.table <- as.data.frame(prop.table(table(isolation.table$RI, isolation.table$Order), margin = 2))
colnames(isolation.table)[1:2] <- c("Barrier","Order")

# Plotting
png("../figures_corrected/13_BarV.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  theme_classic() +
  theme(text = element_text(family = "serif"))
dev.off()
