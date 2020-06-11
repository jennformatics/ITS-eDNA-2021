#################
###  Startup  ###
#################

# Tidyverse, duh.
library("tidyverse")

# Graphics additions
library("RColorBrewer")
library("gplots"); library("cowplot"); library("egg")  # For "easy" plot stacking
library("grid"); library("gridExtra")                  # Also for plot combination
library("extrafont")                                   # (optional) For getting reasonable fonts on Windows

# Analysis additions
library("phyloseq")
library("vegan")
library("VennDiagram")
library("treemap")
library("Hmisc")          # Statistical additions

setwd("~/Box Sync/Projects/UNDERC-ITS/Analysis/04-new-plots-for-publication")
source("utility-functions.R")  # includes plot_with_hull; plot_alpha_diversity_by


#######################################
###  Load (highly customized) data  ###
#######################################

# Or skip to the "restart" section, below.

# plants.csv is a combined OTU table and taxon description table.
# It's separated in my local case at column 38. (Though I should probably use a dplyr select instead.)
rawdata <- read.table("plants.csv", sep=",", header=TRUE, row.names=1)
rawmetadata <- data.frame(rawdata[1:38])
rawmetadata$UnspacedSciName <- rownames(rawmetadata)
# This one will give errors about some entries being messy, with not exactly two components. That's fine.
inputtaxa <- data.frame(rawmetadata %>% separate(UnspacedSciName, into=c("Genus", "Species")) %>% replace_na(list(Species="sp.")))
# Actual taxon-vs-species read table. Includes controls, but they'll be taken out for lack of metadata shortly.
inputotus <- t(rawdata[39:dim(rawdata)[2]])

# Import sample metadata, and create "DistDepth" convenience column. Which I should probably rename "Compartment".
inputsamples <- read.table("sample-metadata.csv", sep=",", header=TRUE, row.names=1)
inputsamples$DistDepth <- paste(inputsamples$Distance, inputsamples$Depth, sep="-")

# Make PS
ot <- otu_table(inputotus, taxa_are_rows=FALSE)
tt <- tax_table(as.matrix(inputtaxa))
sd <- sample_data(inputsamples)
# sd$Run <- as.factor(sd$Run)    # May not be needed any more now that we have sets.

psraw <- phyloseq(ot, tt, sd)
# saveRDS(psraw, file="psraw.rds")   # Don't really need this one any more in later analysis.

# Remove all samples that have NAs in their metadata -- basically, the blanks.
blank_samples <- rownames(na.omit(sample_data(psraw)))
ps <- subset_samples(psraw, sample_names(psraw) %in% blank_samples)

ps_all <- dezero(ps)    # dezero is provided in utility-functions.R
ps_set_a <- dezero(subset_samples(ps, LakeSet=="Lake Set A"))
ps_set_b <- dezero(subset_samples(ps, LakeSet=="Lake Set B"))

save(file="ps-objects-set-a.RData", ps_set_a)
save(file="ps-objects-set-b.RData", ps_set_b)
save(file="ps-objects-all.RData", ps_all)

# From this point on, "ps" is the placeholder for whatever full or half dataset you're working with.
# Eventually we'll have: ps, psbig, psdd, pslake, psddten, pslaketen, pswet, pswetdd, pswetlake, pswetlakepct

# Make the new-primer six-lake set the primary data going forward.
ps <- ps_set_b


###############################
###  General restart point  ###
###############################

library("tidyverse"); library("phyloseq"); library("vegan"); library("grid"); library("gridExtra")
setwd("~/Box Sync/Projects/UNDERC-ITS/Analysis/04-new-plots-for-publication")
source("utility-functions.R")  # includes plot_with_hull; plot_alpha_diversity_by

# Pick one:
# currentfilename <- "ps-objects-set-a.RData"
# currentfilename <- "ps-objects-all.RData"
currentfilename <- "ps-objects-set-b.RData"

load(file=currentfilename)

# ps will now be a 6- or 12-lake set depending on your choice above.


##########################################
###  Rarefaction curve for inspection  ###
##########################################

pdf("plots/rarecurve.pdf", height=4, width=4)
rarecurve(t(otu_table(ps)), step=100, cex=0.5, label=FALSE)
dev.off()


############################
###  Raw "what's where"  ###
############################

# Generates basic bar graphs and heatmaps for internal consumption

# Basic sample aggregation levels
pslake <- merge_samples(ps, "Lake")     # 12 samples, 3
psdd <- merge_samples(ps, "DistDepth")

pstmp <- ps
# Remove all columns except the ones were about to glom on.
tax_table(pstmp) <- tax_table(pstmp)[,c("Vascularity", "WetlandSort", "WetlandStatus")]
pswet <- tax_glom(pstmp, "WetlandStatus")                 # 7 taxa
taxa_names(pswet) <- tax_table(pswet)[,"WetlandStatus"]
pswetlake <- merge_samples(pswet, "Lake")                 # 12 samples

# Get top n taxa for each lake, and overall.
n <- 20
get_tops <- function(lakename) {sort(t(otu_table(pslake)[lakename,]), decreasing=TRUE) %>% head(n) %>% rownames() }
bylaketops <- unlist(unlist(lapply(rownames(otu_table(pslake)), get_tops)))
overalltops <- sort(colSums(otu_table(pslake)), decreasing=TRUE) %>% head(n) %>% names()
overalltops %in% bylaketops   # All true at the moment
in_tops <- levels(factor(c(bylaketops, overalltops)))

pslaketops <- dezero(subset_taxa(pslake, rownames(tax_table(pslake)) %in% in_tops))
psddtops <- dezero(subset_taxa(psdd, rownames(tax_table(psdd)) %in% in_tops))

otutops <- data.frame(t(otu_table(pslaketops)), TaxName=row.names(t(otu_table(pslaketops))))
otutrans <- t(otu_table(pslake))
allsums <- data.frame(TaxTotal=rowSums(otutrans), TaxName=rownames(otutrans))
top_tax_table <- left_join(otutops, allsums)
rownames(top_tax_table) <- top_tax_table[,"TaxName"]
top_tax_table %>% arrange(TaxTotal)

# TaxName                 TopIn
# Nuphar.variegata	Bay, Hum, Ink
# Peridinium.wierzejskii	Mor
# Chrysophyceae.sp03	Long
# Synura.mammillosa	Ras
# Potamogeton	
# Gonyostomum.semen	
# Potamogeton.lucens	

# Plot absolute and relative read numbers per lake for vasculars and nonvasculars.
vascbar <- plot_bar_jd(pswetlake, fill="Vascularity")
pswetlakepct <- transform_sample_counts(pswetlake, function(x) x / sum(x))
vascbarpct <- plot_bar_jd(pswetlakepct, fill="Vascularity")
pslakebool <- decostand(otu_table(pslake), "pa") %>% otu_table() %>% phyloseq(tax_table(pslake), sample_data(pslake))
lakeboolbar <- plot_bar_jd(pslakebool, fill="Vascularity")

vascbars <- grid.arrange(vascbar, vascbarpct, lakeboolbar)
ggsave(vascbars, file="plots/vascular-bar.pdf", height=11, width=6)

pswetdd <- merge_samples(pswet, "DistDepth")
# pswetdd <- dezero(prune_taxa(taxa_names(pswetdd) != "UPL", pswetdd))
sample_data(pswetdd)[,"DistDepth"] <- sample_names(pswetdd)
sample_data(pswetdd) <- sample_data(pswetdd)[,"DistDepth"]

tax_descs <- c("Upland", "Fac. Upland", "Facultative", "Fac. Wetland", "Obligate Wetland", "Aquatic", "Algae")

tax_table(pswetdd) <-
    data.frame(tax_table(pswetdd), rowname=rownames(tax_table(pswetdd))) %>%
    arrange(desc(WetlandSort)) %>%
    cbind(tax_descs) %>%
    column_to_rownames() %>%
    as.matrix() %>%
    tax_table()

# taxa_names(pswetdd) <- tax_table(pswetdd)[,"tax_descs"]

sortorder <- tax_table(pswetdd) %>% as.data.frame() %>% arrange(WetlandSort) %>% select(tax_descs)
sortorder <- unname(sortorder[,1]) %>% as.character()

pswetddrel <- transform_sample_counts(pswetdd, function(x) { x/sum(x) })

# Apparently it doesn't like to do this with only six lakes; it puts NAs in places we can't use them.
lakewetheatmap <- plot_heatmap(pswetdd, taxa.order=sortorder, sample.order=c("Far-Deep", "Far-Shallow", "Near-Deep", "Near-Shallow"), low="blue", high="yellow")
ggsave(lakewetheatmap, file="plots/lake-by-wetland-heatmap.pdf", height=6, width=6)
lakewetheatmaprel <- plot_heatmap(pswetddrel, taxa.order=sortorder, sample.order=c("Far-Deep", "Far-Shallow", "Near-Deep", "Near-Shallow"), low="blue", high="yellow")
ggsave(lakewetheatmaprel, file="plots/lake-by-wetland-heatmap-rel.pdf", height=6, width=6)

# At this point we have: ps, psdd, pslake, psddten, pslaketen, pswet, pswetdd, pswetddrel, pswetlake, pswetlakepct
# ...but we're not going to save the ones we only used for the plots above.
save(file=currentfilename, ps, psdd, pslake, pswet)


#########################
###  Alpha Diversity  ###
#########################

# On not transforming before doing alpha diversity measures:
# https://github.com/joey711/phyloseq/issues/287
# https://bioconductor.statistik.tu-dortmund.de/packages/3.4/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

# On denoising and singletons/doubletons
# https://github.com/joey711/phyloseq/issues/445

# Waste not
# https://doi.org/10.1371/journal.pcbi.1003531

source("utility-functions.R")  # includes themes; plot_with_hull; plot_alpha_diversity_by
load(file=currentfilename)    # E.g., ps-objects-3738.RData
indices <- c("Shannon") # , "InvSimpson")  # "Simpson", "Chao1", 

theme_set(theme_minimal() + theme(plot.background = element_rect(fill="white", color="white")))

combined_alpha_plot <- function(inputps, namefrag) {
    p <- grid.arrange(plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices),
                     plot_alpha_diversity_by(psobj=inputps, var="Depth", measures=indices),
                     plot_alpha_diversity_by(psobj=inputps, var="Distance", measures=indices),
                     layout_matrix = rbind(c(1, 1), c(2, 3)))
    fname <- paste("alpha", namefrag, "pdf", sep=".")
    ggsave(p, filename=fname, width=8, height=8)
}

combined_alpha_plot(inputps=ps, namefrag="diversity")

# New specific plots for paper/presentation

inputps <- ps
lakemeta <- read.table("lake-metadata.txt", header=TRUE, sep="\t", colClasses=c("Lake"="character"))
lakemeta <- lakemeta[lakemeta["LakeSet"] == "B",]
# sample_data(ps) foo <- sample_data(as.matrix(left_join(sample_data(ps), lakemeta)))
# indices <- c("Shannon")

source("utility-functions.R")  # includes themes; plot_with_hull; plot_alpha_diversity_by

# Redefine some things in the plotting function.
plot_alpha_diversity_by <- function(var, psobj, measures=c("Shannon", "Simpson", "InvSimpson", "Chao1")) {
    p <- plot_richness(psobj, x=var, measures=measures, title="")
    p$layers <- p$layers[-1]  # Remove the default set of points since we can't manipulate them
    # Now, instead of individual points, do a boxplot.
    p <- p +   # theme(axis.text.x = element_text(face="italic")) +
        geom_point(position=position_jitter(w = 0.2, h = 0), aes(shape=DistDepth), size=2, alpha=0.5) + 
        scale_shape_manual(values=c("Far-Deep"=19, "Far-Shallow"=21, "Near-Deep"=15, "Near-Shallow"=22)) +
        labs(shape = "Compartment") + 
        stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
        labs(x=var, y="Shannon Index (H')", title=NULL, main=NULL)
    return(p)
}

library("cowplot")
library("egg")

p1 <- ggplot(lakemeta, mapping=aes(x=Lake, y=Area, group=1)) +
      geom_col(size=1, color="darkgreen", fill="darkgreen") +
      scale_fill_manual(values = "darkgreen") + labs(fill = "Lake Area") +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + ylab("Lake Area (ha)")

p2 <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices) +
      labs(shape = "Compartment") +
      theme(strip.text.x = element_blank(), axis.text.x = element_text(angle = -45, size=10)) 

ratio <- c(0.25, 0.75)
cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = ratio)
e <- egg::ggarrange(p1, p2, heights = ratio)
ggsave(e, filename="plots/complex-lake-alpha-div.pdf", h=4, w=8)

# l <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)
dep <- plot_alpha_diversity_by(psobj=inputps, var="Depth", measures=indices)
dis <- plot_alpha_diversity_by(psobj=inputps, var="Distance", measures=indices)

# ggsave(l, filename="plots/alpha.shannon.lake.pdf", height=4, width=6)
ggsave(dep, filename="plots/alpha.shannon.depth.pdf", height=4, width=4)
ggsave(dis, filename="plots/alpha.shannon.dist.pdf", height=4, width=4)



# Linear regression for alpha diversity vs. area. Will need to update for the overall Spearman thing.

alphadivmeans <- p2$data %>% group_by(Lake) %>% dplyr::summarize(mean=mean(value))
lakeareas <- lakemeta[c("Lake", "Area")]
mean_vs_area <- left_join(alphadivmeans, lakeareas)
mean_vs_area <- mutate(mean_vs_area, logarea=log(Area))
mean_vs_area <- mutate(mean_vs_area, squarea=sqrt(Area))

cor_test(mean_vs_area, "mean", "Area", method="spearman")    # rho=0.94, p = 0.005
cor_test(mean_vs_area, "mean", "Area", method="pearson")     # r=0.80, p=0.055
cor_test(mean_vs_area, "mean", "squarea", method="pearson")  # r=0.88, p=0.022

cor_test(mean_vs_area, "mean", "logarea", method="pearson")  # r=0.95, p=0.004

qplot(mean_vs_area$mean, mean_vs_area$logarea)   # R^2 = 0.8759, according to lm(formula = mean_vs_area$mean ~ mean_vs_area$logarea)
ggsave(filename="plots/area-diversity-dotplot.pdf", h=3, w=3)

# Generate Spearman correlations of Shannon diversity with lake variables

inputps <- ps
indices <- c("Shannon")

alphadiv <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)$data
lakeshannonmeans <- alphadiv %>% group_by(variable, Lake) %>% dplyr::summarize(mean=mean(value)) %>% spread(variable, mean)
lakemeta <- read.table("lake-metadata.txt", header=TRUE)
bylake <- left_join(lakemeta, lakeshannonmeans) %>%
          subset(lakemeta$Lake %in% sample_data(ps)$Lake) %>%
          select(-LakeAbbr, -LakeSet, -Date, -AlphaSort) %>%
          remove_rownames() %>% column_to_rownames("Lake")

          # mutate(ShannonRank=rank(Shannon), ShannonNorm=Shannon/sum(Shannon)) %>%

tart <- colorRampPalette(c("red3", "#FFFFBF", "#00441B"))
hm <- heatmap(cor(bylake, method="spearman"), symm=TRUE, Rowv=NA, Colv="Rowv", col=tart(100))

p <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures="Shannon") +
     scale_y_log10() + facet_grid(rows="DistDepth") + theme(axis.text.x=element_blank())

q <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures="Shannon") + scale_y_log10()
grid.arrange(p, q)

# ps_vasc <- dezero(subset_taxa(ps, Vascularity == "Vascular"))
# combined_alpha_plot(inputps=ps_vasc, namefrag="vasc")
# ps_nonvasc <- dezero(subset_taxa(ps, Vascularity == "Nonvascular"))
# combined_alpha_plot(inputps=ps_nonvasc, namefrag="nonvasc")

# Mann-Whitney, which is 2-sample Wilcoxson rank sum. Used in https://www.frontiersin.org/articles/10.3389/fpls.2016.02015/full

pval_mann_whitney <- function(indexname, inputps1, inputps2) {
    rich1 <- unlist(unname(estimate_richness(inputps1, measures=indexname)))
    rich2 <- unlist(unname(estimate_richness(inputps2, measures=indexname)))
    output <- wilcox.test(rich1, rich2, correct=FALSE)$p.value
    names(output) <- indexname
    return(output)
    }

inputps <- ps

print("Mann-Whitney test on alpha diversity by depth category:")
inputps1 <- subset_samples(inputps, Depth=="Shallow")
inputps2 <- subset_samples(inputps, Depth=="Deep")
rev(unlist(lapply(indices, function(x) { pval_mann_whitney(x, inputps1, inputps2) } ), recursive=TRUE))

print("Mann-Whitney test on alpha diversity by shore distance category:")
inputps1 <- subset_samples(inputps, Distance=="Far")
inputps2 <- subset_samples(inputps, Distance=="Near")
rev(unlist(lapply(indices, function(x) { pval_mann_whitney(x, inputps1, inputps2) } ), recursive=TRUE))

# Depth, for 3738 runs
    Chao1   Shannon 
0.8714475 0.7327865 
# Dist, for 3738 runs
     Chao1    Shannon 
0.07539991 0.00264300 

p3 <- plot_alpha_diversity_by(psobj=inputps, var="Depth", measures=indices) + labs(shape = "Compartment") + theme(strip.text.x = element_blank(), axis.text.x = element_text(angle = -45, size=10)) 

########################
###  Beta Diversity  ###
########################

# Lots of stuff about adonis is specific to the ps object, so wrap it in a function.
# The log+1 normalization with dEcoStand is baked in here.
    # "as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0"
    # Higher bases give less weight to quantities and more to presences, and logbase = Inf gives the presence/absence scaling.
    # Please note this is not log(x+1).

ps3233 <- subset_samples(ps, Run %in% c("32", "33"))
ps3738 <- subset_samples(ps, Run %in% c("37", "38"))
ps <- ps3738   # 3233

# Create "psbig", with only the taxa that appear in five or more samples, and only samples that have 100 or more reads.
taxa_under_five <- colnames(otu_table(ps))[unlist(lapply(colnames(otu_table(ps)), function(x) { sum(as.numeric(unname(otu_table(ps)[,x])>0)) })) < 5]
psbig <- subset_taxa(ps, !colnames(otu_table(ps)) %in% taxa_under_five)
psbig <- subset_samples(psbig, sample_sums(psbig) > 100)

pscomplement <- subset_taxa(ps, colnames(otu_table(ps)) %in% taxa_under_five)
table(tax_table(pscomplement)[,"Vascularity"])   # Breakdown of eliminated taxa.

save(file=currentfilename, ps, psbig, psdd, pslake, psddten, pslaketen, pswet, pswetdd, pswetlake, pswetlakepct)

bdiv <- function(psobj, frml, samblk="Lake", b=2) {
    sam <- data.frame(sample_data(dezero(psobj)))
    blkfactor <- as.factor(unname(as.matrix(sam[samblk])))
    onorm <- otu_table(decostand(data.frame(otu_table(dezero(psobj))), method="log", logbase=b, na.rm=TRUE), taxa_are_rows=FALSE)
    nperms <- 999; perm <- how(nperm = nperms); setBlocks(perm) <- with(sam, blkfactor)
    adonis(as.formula(paste("onorm ~ ", frml)), data=sam, permutations=perm)
}

bdivnoblk <- function(psobj, frml, b=2) {
    sam <- data.frame(sample_data(dezero(psobj)))
    onorm <- otu_table(decostand(data.frame(otu_table(dezero(psobj))), method="log", logbase=b, na.rm=TRUE), taxa_are_rows=FALSE)
    nperms <- 999; perm <- how(nperm = nperms)
    adonis(as.formula(paste("onorm ~ ", frml)), data=sam, permutations=perm)
}

==========

inputps <- ps
indices <- c("Shannon")

alphadiv <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)$data
lakeshannonmeans <- alphadiv %>% group_by(variable, Lake) %>% dplyr::summarize(mean=mean(value)) %>% spread(variable, mean)
lakemeta <- read.table("lake-metadata.txt", header=TRUE)
bylake <- left_join(lakemeta, lakeshannonmeans) %>%
          subset(lakemeta$Lake %in% sample_data(ps)$Lake) %>%
          select(-LakeAbbr, -LakeSet, -Date, -AlphaSort) %>%
          remove_rownames() %>% column_to_rownames("Lake")

          # mutate(ShannonRank=rank(Shannon), ShannonNorm=Shannon/sum(Shannon)) %>%

alpharanks <- bylake

# Do it manually instead. Sigh.
glarg <- data.frame(i=rep(0, 3))
lbls <- c("PDist", "PDep", "PDD")
for (k in rownames(alpharanks)) { tempps <- subset_samples(psbig, Lake==k);
    glarg <- data.frame(glarg, bdiv(tempps, "Distance+Depth+Distance:Depth")$aov.tab$"Pr(>F)"[1:3], row.names=lbls) }
colnames(glarg) <- c("i", as.character(rownames(alpharanks)))
glarg <- t(glarg[,2:dim(glarg)[2]])
pwithinlakes <- glarg %>% as.data.frame %>% rownames_to_column %>% rename(Lake=rowname)

glarg <- data.frame(i=rep(0, 3))
lbls <- c("R2Dist", "R2Dep", "R2DD")
for (k in rownames(alpharanks)) { tempps <- subset_samples(psbig, Lake==k);
    glarg <- data.frame(glarg, bdiv(tempps, "Distance+Depth+Distance:Depth")$aov.tab$R2[1:3], row.names=lbls) }
colnames(glarg) <- c("i", as.character(rownames(alpharanks)))
glarg <- t(glarg[,2:dim(glarg)[2]])
r2withinlakes <- glarg %>% as.data.frame %>% rownames_to_column %>% rename(Lake=rowname)

alpharanks <- alpharanks %>% rownames_to_column %>% rename(Lake=rowname)
moarcor <- alpharanks %>% left_join(pwithinlakes) %>% left_join(r2withinlakes) %>% as.data.frame
moarcor <- moarcor %>% select(-Lake)
# moarcor <- moarcor %>% select(-AlphaSort, -Date, -variable, -LakeAbbr, -Lake)
# tart <- colorRampPalette(c("red3", "#FFFFBF", "#00441B"))
# spec <- colorRampPalette(c("red3", "yellow", "#00441B", "blue"))
# spec2 <- colorRampPalette(brewer.pal(3, "Spectral"))
bluegold <- colorRampPalette(c("blue", "white", "gold"), breaks=c(-1, 0, 1))
heatmap(cor(moarcor, method="spearman"), symm=TRUE, Rowv=NA, Colv="Rowv", col=bluegold(100))
blackwhite <- colorRampPalette(c("black", "white", "black"))
heatmap(cor(moarcor, method="spearman"), symm=TRUE, Rowv=NA, Colv="Rowv", col=blackwhite(100))
# Nothing of interest.
# Or, wait; maybe effect of depth is related to area? and/or Chao richness?

lesscor <- as.tbl(moarcor) # %>% select(-starts_with("Norm"), -SimpsonValue, -ShannonRank, -SimpsonRank, -starts_with("Chao"), -starts_with("PD"), -R2DD)

library(gplots)
spears <- cor(lesscor, method="spearman")
colbreaks <- c(seq(min(spears), -0.01, length=50), 0, seq(0.01, max(spears), length=50))
heatmap.2(spears, symm=TRUE, Rowv=NA, Colv="Rowv", col=bluegold(100), breaks=colbreaks)

# Utilities from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software#correlation-matrix-with-significance-levels-p-value
library("Hmisc")

corrobj <- rcorr(as.matrix(lesscor))

flattenCorrMatrix <- function(corrobj) {
  cormat <- corrobj$r
  pmat <- corrobj$P
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

flatCorrMatrix <- function(dattab) {
  corrobj <- rcorr(as.matrix(lesscor))
  flattenCorrMatrix(corrobj)
}

ps3738 <- subset_samples(ps, Run %in% c("37", "38"))
inputps <- psbig
bigpermanova <- bdiv(inputps, "Lake+Distance+Depth+Distance:Depth")
runpermanova <- bdiv(inputps, "Run")

f <- file("permanova.txt", "w")
cat("By Lake, Distance, Depth, and Distance x Depth\n", file=f, sep="\n", append=TRUE)
write.table(as.data.frame(bigpermanova$aov.tab), file=f, sep="\t", na="", quote=FALSE, append=TRUE)
cat("\nBy sequencing run\n", file=f, sep="\n", append=TRUE)
write.table(as.data.frame(runpermanova$aov.tab), file=f, sep="\t", na="", quote=FALSE, append=TRUE)
close(f)

# bdiv(inputps, "Lake+Distance+Depth+DistDepth", b=Inf)
# bdiv(inputps, "Run", b=Inf)

# Data for taxon Venn
# install.packages("gplots")
# install.packages("VennDiagram")
library("gplots")
library("VennDiagram")

# install.packages('extrafont')
library(extrafont)
font_import(pattern="arial.*", prompt=FALSE)
font_import(pattern="times.*", prompt=FALSE)
loadfonts(device="win")

psbigdd <- merge_samples(psbig, "DistDepth")
psmed <- subset_taxa(ps, taxa_sums(ps) > 10)
psmeddd <- subset_taxa(psdd, taxa_sums(ps) > 10)
psvasc <- subset_taxa(ps, Vascularity=="Vascular")
psvascdd <- subset_taxa(psdd, Vascularity=="Vascular")
psmedvasc <- subset_taxa(psmed, Vascularity=="Vascular")
pstwovasc <- subset_taxa(psvasc, Org.Sample.Prs != "  1")
psbigvascdd <- subset_taxa(psbigdd, Vascularity=="Vascular")

cols=c("red", "blue", "orange", "purple")
cols=c("#999999", "#999999", "#999999", "#999999")
cols=c("#FF000066", "#0000FF66", "#FFFF0066", "#66009966")
cols=c("#FF000066", "#0000FF66", "#66000066", "#00003366")
cols=c("#FF660066", "#0000FF66", "#FFFF0066", "#00FF0066")

inputps <- psbig

psns <- dezero(subset_samples(inputps, DistDepth=="Near-Shallow"))
psnd <- dezero(subset_samples(inputps, DistDepth=="Near-Deep"))
psfs <- dezero(subset_samples(inputps, DistDepth=="Far-Shallow"))
psfd <- dezero(subset_samples(inputps, DistDepth=="Far-Deep"))

ns <- taxa_names(psns)
nd <- taxa_names(psnd)
fs <- taxa_names(psfs)
fd <- taxa_names(psfd)

a <- rownames(tax_table(inputps))
fs <- rownames(tax_table(psfs)); fd <- rownames(tax_table(psfd)); ns <- rownames(tax_table(psns)); nd <- rownames(tax_table(psnd))
faronly <- a[! a %in% c(ns, nd)]; nearonly <- a[! a %in% c(fs, fd)]; deeponly <- a[! a %in% c(ns, fs)]; shallowonly <- a[! a %in% c(nd, fd)]
fsonly <- a[! a %in% c(fd, ns, nd)]; fdonly <- a[! a %in% c(fs, nd, ns)]; nsonly <- a[! a %in% c(nd, fd, fs)]; ndonly <- a[! a %in% c(fs, fd, ns)]

venn.diagram(list("Far-Deep"=fd, "Far-Shallow"=fs, "Near-Deep"=nd, "Near-Shallow"=ns),
    filename="plots/venn-test.png", imagetype="png", fill=cols, col="white", cex=2, cat.cex=1.5, margin=.05, fontfamily="Arial")

faronly
nearonly
deeponly
shallowonly

# table(tax_table(dezero(subset_samples(psvasc, sample_names(psvasc)=="Near-Deep")))[,"WetlandStatus"])


##############
###  NMDS  ###
##############

# NMDS theme (defined in utility-functions.R)
theme_set(nmdstheme)

# ord.bray.k3 <- ordinate(ps, method="NMDS", distance="bray", k=3, trymax=100)
# ord.bray.k5 <- ordinate(ps, method="NMDS", distance="bray", k=5, trymax=100)
# ord.bray.k10 <- ordinate(ps, method="NMDS", distance="bray", k=10, trymax=100)
# ord.bray.k20 <- ordinate(ps, method="NMDS", distance="bray", k=20, trymax=100)
# saveRDS(ord.bray.k20, file="ord.bray.k20.rds")

ord.bray.k20 <- readRDS(file="ord.bray.k20.rds")

k=20
ord.nmds.bray <- eval(parse(text=paste("ord.bray.k", k, sep="")))

lakeplot <- plot_with_hull(ps, ord.nmds.bray, "Lake") + make_stress_label(ord.nmds.bray)
distplot <- plot_with_hull(ps, ord.nmds.bray, "Distance")
depthplot <- plot_with_hull(ps, ord.nmds.bray, "Depth")
runplot <- plot_with_hull(ps, ord.nmds.bray, "Run")
setplot <- plot_with_hull(ps, ord.nmds.bray, "LakeSet")
primerplot <- plot_with_hull(ps, ord.nmds.bray, "Primer") + scale_colour_manual(values = c("darkgreen", "purple"), aesthetics = c("colour", "fill"))

taxplot <- plot_with_hull(ps, ord.nmds.bray, "Vascularity", type="taxa")
taxplot <- taxplot + scale_colour_manual(values = c("green", "darkgreen"), aesthetics = c("colour", "fill"))

wetnessplot <- plot_with_hull_ordered(ps, ord.nmds.bray, "WetlandStatus", type="taxa")

ggsave(lakeplot, file="individual-nmdses/lakeplot-12lakes.pdf", width=6.25, height=5)
ggsave(distplot, file="individual-nmdses/distplot-12lakes.pdf", width=6.25, height=5)
ggsave(depthplot, file="individual-nmdses/depthplot-12lakes.pdf", width=6.25, height=5)
ggsave(runplot, file="individual-nmdses/runplot-12lakes.pdf", width=6.25, height=5)
ggsave(setplot, file="individual-nmdses/setplot-12lakes.pdf", width=6.25, height=5)
ggsave(primerplot, file="individual-nmdses/primerplot-12lakes.pdf", width=6.25, height=5)
ggsave(taxplot, file="individual-nmdses/taxplot-12lakes.pdf", width=6.25, height=5)
ggsave(wetnessplot, file="individual-nmdses/wetnessplot-12lakes.pdf", width=6.25, height=5)

# Example for removing k= and stress= text
primerplot$layers <- primerplot$layers[c(1,2)]
# ...or adding it
primerplot <- primerplot + make_stress_label(ord.nmds.bray)

# rowscols <- rbind(c(1,2),c(3,4))
# rowscols <- rbind(c(1,1,1,1,2,2),c(3,3,3,4,4,4))
# q1 <- arrangeGrob(lakeplot, runplot, depthplot, distplot, layout_matrix = rowscols)
# q1 <- arrangeGrob(q, top = make_nmds_title(paste("All Plant Categories, k=", k)))
q1 <- arrangeGrob(depthplot, distplot, layout_matrix = rbind(c(1,2)))

namefrag <- paste("k", k, sep="")
fname <- paste("plots/depthdist.", namefrag, ".pdf", sep="")
ggsave(q1, filename=fname, width=10, height=5)
ggsave(lakeplot, filename=paste("plots/lakeplot.", namefrag, ".pdf", sep=""), width=5, height=5)
ggsave(runplot, filename=paste("plots/runplot.", namefrag, ".pdf", sep=""), width=5, height=5)

q2 <- arrangeGrob(taxplot, wetnessplot, layout_matrix = rbind(c(1,2)))

namefrag <- paste("k", k, sep="")
fname <- paste("plots/taxa.", namefrag, ".pdf", sep="")
ggsave(q2, filename=fname, width=10, height=5)

# Just hulls, for overlays [not yet working]

p <- plot_ordination(ps, ord.nmds.bray, color="Primer")
DF <- p$data

x = colnames(DF)[1]
y = colnames(DF)[2]
ord_map = aes_string(x = x, y = y, color = "Primer", na.rm = TRUE)
q <- ggplot(DF, ord_map)
hulls <- plyr::ddply(na.omit(p$data), "Primer", find_hull)
q + geom_polygon(data = hulls, alpha = 0) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank())


, aes(fill=eval(as.name(paste(var))))) + 
             guides(fill=FALSE) + labs(color=var, fill=var) # + make_stress_label(ord)
    return(p)
    }

# Betadisper

d <- vegdist(otu_table(ps))
betadisper(d, "Depth")
# Error in x - c : non-conformable arrays


#################
###  Treemap  ###
#################

library("tidyverse")
library("phyloseq")
library("vegan")
library("grid"); library("gridExtra")
library("treemap")

setwd("~/Box Sync/Projects/UNDERC-ITS/Analysis/04-new-plots-for-publication")
source("utility-functions.R")  # includes plot_with_hull; plot_alpha_diversity_by

# Pick one:
currentfilename <- "ps-objects-set-a.RData"
currentfilename <- "ps-objects-set-b.RData"
currentfilename <- "ps-objects-all.RData"

load(file=currentfilename)


pstmp <- ps
# Simplify tax_table so it doesn't take forever to manipulate.
tax_table(pstmp) <- as.matrix(select(as.data.frame(tax_table(ps)[,c("Taxon.Category", "PlantCat", "Vascularity")]), Vascularity, PlantCat, Taxon.Category, everything()))

pscategories <- pstmp

twogreens=c("#61953D", "#A9D18E")

# General source data.frame for both read and count maps (counts will be added later)
tinput <- data.frame(Vasc=tax_table(pscategories)[,"Vascularity"], TaxCat=tax_table(pscategories)[,"Taxon.Category"], Reads=unname(rowSums(t(otu_table(pscategories)))))

# Thresholds, consolidation, and plotting for read map.
readthresh <- as.integer(sum(tinput[,"Reads"]) / 100)
readbig <- filter(tinput, Reads >= readthresh)
readother <- filter(tinput, Reads < readthresh) %>% group_by(Vascularity) %>% dplyr::summarize(Reads=sum(Reads)) %>% 
    cbind(Taxon.Category="other") %>% select(Vascularity, Taxon.Category, Reads)
fullreadinput <- rbind(readbig, readother) %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(Reads=sum(Reads))

jd_treemap(fullreadinput, "Reads", c(19, 13))
sum(fullreadinput$Reads)   # ~1.3M

# And that's the reads one.

# Thresholds, consolidation, and plotting for counts map.

t2input <- tinput %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(CatCount=n())
countthresh <- 10 # as.integer(n(t2input[,"CatCount"]) / 100)
countbig <- as.data.frame(filter(t2input, CatCount >= countthresh))
countother <- filter(t2input, CatCount < countthresh) %>% group_by(Vascularity) %>% dplyr::summarize(CatCount=sum(CatCount)) %>% 
    cbind(Taxon.Category="other") %>% select(Vascularity, Taxon.Category, CatCount)
fullcountinput <- rbind(countbig, countother) %>% group_by(Vascularity, Taxon.Category)

jd_treemap(fullcountinput, "CatCount", c(19, 13))
sum(fullcountinput$CatCount)   # 248

======================

pstmp <- ps

tax_table(pstmp) <- as.matrix(select(as.data.frame(tax_table(ps)[,c("Taxon.Acc", "Taxon.Category", "PlantCat", "Top.Taxon", "Vascularity")]), Vascularity, Top.Taxon, PlantCat, Taxon.Category, Taxon.Acc, everything()))

pscategories <- pstmp
taxa_names(pscategories) <- tax_table(pscategories)[,"Taxon.Acc"]

tinput <- data.frame(Vasc=tax_table(pscategories)[,"Vascularity"], TaxCat=tax_table(pscategories)[,"Taxon.Category"], Reads=unname(rowSums(t(otu_table(pscategories)))))

tinput %>% group_by(Vascularity, Taxon.Category) %>% summarize(CatCount=n())
