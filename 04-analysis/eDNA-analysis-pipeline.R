setwd("~/Box Sync/Projects/UNDERC-ITS/Analysis/ITS-eDNA-2020/04-analysis")
source("utility-functions.R")  # includes imports, functions, and other startup

currentfilename <- "data/ps.RData"
load(file=currentfilename)


#######################################
###  Load (highly customized) data  ###
#######################################

# Import sample metadata, and create "Compartment" column.
inputsamples <- read.table("data/sample-metadata.txt", sep="\t", header=TRUE, row.names=1) %>% rownames_to_column("ShortName")

# CHANGE HERE if going back to 12 lakes.
workingsamples <- inputsamples # %>% filter(!Run %in% c("32", "33"))

mungedsamples <- workingsamples %>%
  mutate(Depth = case_when(Depth == 'Shallow' ~ 'Surface', Depth == 'Deep' ~ 'Benthic', TRUE ~ 'NA')) %>%
  mutate(Distance = case_when(Distance == 'Near' ~ 'Nearshore', Distance == 'Far' ~ 'Offshore', TRUE ~ 'NA')) %>%
  mutate(Compartment = paste(Distance, Depth, sep="-")) %>%
  mutate(ShortCompartment = case_when(Compartment == 'Nearshore-Surface' ~ 'NearSurf',
                                      Compartment == 'Offshore-Surface' ~ 'OffSurf',
                                      Compartment == 'Nearshore-Benthic' ~ 'NearBenth',
                                      Compartment == 'Offshore-Benthic' ~ 'OffBenth',
                                      TRUE ~ 'NA')) %>%
  column_to_rownames("ShortName")

# plants.csv is a combined OTU table and taxon description table.
# It's separated in my local case between metadata field "Org.Lake.Prs" and sample name "MC-32".

rawdata <- read.table("data/plants-2020-05-27-ancestral-1st-manuscript.csv", sep=",", header=TRUE, row.names=1)
# rawdata <- read.table("data/plants.csv", sep=",", header=TRUE, row.names=1)

rawmetadata <- dplyr::select(rawdata, 1:Org.Lake.Prs)  # data.frame(rawdata[1:39])
rawmetadata$UnspacedSciName <- rownames(rawmetadata)
# This one will give errors about some entries being messy, with not exactly two components. That's fine.
inputtaxa <- data.frame(rawmetadata %>% separate(UnspacedSciName, into=c("Genus", "Species")) %>% replace_na(list(Species="sp.")))
# Actual taxon-vs-species read table. Includes controls, but they'll be taken out for lack of metadata shortly.

inputotucols <- select(rawdata, starts_with("MC"):last_col())
inputotus <- inputotucols %>% select(one_of(workingsamples$ShortName)) %>% t()

# Make PS
ot <- otu_table(inputotus, taxa_are_rows=FALSE)
tt <- tax_table(as.matrix(inputtaxa))
sd <- sample_data(mungedsamples)
sd$Run <- as.factor(sd$Run)    # To enable selection of samples by run

psraw <- phyloseq(ot, tt, sd)
# saveRDS(psraw, file="psraw.rds")   # Don't really need this one any more in later analysis.

# Remove all samples that have NAs in their metadata -- basically, the blanks.
nonblank_samples <- rownames(na.omit(sample_data(psraw)))
ps12old <- subset_samples(psraw, sample_names(psraw) %in% nonblank_samples)

currentfilename <- "data/ps.RData"
save(file=currentfilename, ps, ps6old, ps6new, ps12old, ps12new)  # plain ps is actually ps6new
load(file=currentfilename)

# From this point on, "ps" is the placeholder for whatever full or half dataset you're working with.
# Eventually we'll have: ps, psbig, psdd, pslake, psddten, pslaketen, pswet, pswetdd, pswetlake, pswetlakepct

# > sum(otu_table(ps12new))
# [1] 1819384
# > sum(otu_table(ps12old))
# [1] 1887671
# > sum(otu_table(ps6new))
# [1] 1325970
# > sum(otu_table(ps6old))
# [1] 1328139

##########################################
###  Rarefaction curve for inspection  ###
##########################################

pdf("plots/background/rarecurve.pdf", height=4, width=4)
rarecurve(t(otu_table(ps)), step=100, cex=0.5, label=FALSE)
dev.off()


############################
###  Raw "what's where"  ###
############################

# Generates basic bar graphs and heatmaps for internal consumption

pstmp <- ps6old

# Basic sample aggregation levels
pslake <- merge_samples(pstmp, "Lake")     # 12 samples, 3
psdd <- merge_samples(pstmp, "Compartment")

# Remove all columns except the ones we're about to glom on.
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
top_tax_table %>% arrange(TaxTotal) %>% tail(20)

# Plot absolute and relative read numbers per lake for vasculars and nonvasculars.
vascbar <- plot_bar_jd(pswetlake, fill="Vascularity")
pswetlakepct <- transform_sample_counts(pswetlake, function(x) x / sum(x))
vascbarpct <- plot_bar_jd(pswetlakepct, fill="Vascularity")

pstmp <- pslake
# Remove all columns except the ones we're about to glom on.
tax_table(pstmp) <- tax_table(pstmp)[,c("Vascularity", "WetlandSort", "WetlandStatus")]
pslakebool <- decostand(otu_table(pstmp), "pa") %>% otu_table() %>% phyloseq(tax_table(pstmp), sample_data(pstmp)) %>% tax_glom("Vascularity")
lakeboolbar <- plot_bar_jd(pslakebool, fill="Vascularity")

vascbars <- grid.arrange(vascbar, vascbarpct, lakeboolbar)
ggsave(vascbars, file="plots/background/vascular-bar.pdf", height=11, width=6)

# Eliminate all irrelevant sample metadata.
pswetdd <- dezero(merge_samples(pswet, "Compartment"))
# pswetdd <- dezero(prune_taxa(taxa_names(pswetdd) != "UPL", pswetdd))
sample_data(pswetdd)[,"Compartment"] <- sample_names(pswetdd)
sample_data(pswetdd) <- sample_data(pswetdd)[,"Compartment"]

# tax_descs <- c("Agricultural", "Upland", "Fac. Upland", "Facultative", "Fac. Wetland", "Obligate Wetland", "Aquatic", "Algae")
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

lakewetheatmap <- plot_heatmap(pswetdd, taxa.order=sortorder, sample.order=c("Far-Deep", "Far-Shallow", "Near-Deep", "Near-Shallow"), low="blue", high="yellow")
ggsave(lakewetheatmap, file="plots/background/lake-by-wetland-heatmap.pdf", height=6, width=6)
lakewetheatmaprel <- plot_heatmap(pswetddrel, taxa.order=sortorder, sample.order=c("Far-Deep", "Far-Shallow", "Near-Deep", "Near-Shallow"), low="blue", high="yellow")
ggsave(lakewetheatmaprel, file="plots/background/lake-by-wetland-heatmap-rel.pdf", height=6, width=6)

# At this point we have: ps, psdd, pslake, psddten, pslaketen, pswet, pswetdd, pswetddrel, pswetlake, pswetlakepct
# ...but we're not going to save the ones we only used for the plots above.
save(file=currentfilename, ps, psdd, pslake, pswet)


###############################
###  Alpha Diversity Plots  ###
###############################

# On not transforming before doing alpha diversity measures:
# https://github.com/joey711/phyloseq/issues/287
# https://bioconductor.statistik.tu-dortmund.de/packages/3.4/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

# On denoising and singletons/doubletons: https://github.com/joey711/phyloseq/issues/445
# Waste not: https://doi.org/10.1371/journal.pcbi.1003531

# Alpha diversity plots

theme_set(theme_minimal() + theme(plot.background = element_rect(fill="white", color="white")))

indices <- c("Shannon")  #, "InvSimpson", "Simpson", "Chao1")
# combined_alpha_plot(inputps=ps, namefrag="diversity")

# New specific plots for paper/presentation

inputps <- ps
lakemeta <- read.table("data/lake-metadata.txt", header=TRUE, sep="\t", colClasses=c("Lake"="character"))

p1 <- ggplot(lakemeta, mapping=aes(x=Lake, y=Area, group=1)) +
      geom_col(size=1, color="#0099FF", fill="#0099FF", width=0.5) +
      scale_y_continuous(breaks=c(0, 25, 50, 75), lim=c(0, 75)) + 
      theme(panel.grid.minor.y = element_blank(), plot.margin=unit(c(24, 0, 0, 0), "pt")) +
      scale_fill_manual(values = "#0099FF") + labs(fill = "Lake Area") +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + ylab("Lake Area (ha)")

p2 <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices) +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      theme(text = element_text(size=12))

p2 

ratio <- c(0.2, 0.8)
# cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = ratio)
e <- egg::ggarrange(p1, p2, heights = ratio)
ggsave(e, filename="plots/complex-lake-alpha-div.pdf", h=4, w=8)

scalefac <- 2
themetweak <- theme(legend.position="none", text = element_text(size=rel(2*scalefac)), axis.text.x = element_text(size=rel(2*scalefac), angle=0, hjust=0.5))
# l <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)
dep <- plot_alpha_diversity_by(psobj=inputps, var="Depth", measures=indices, scalefac=2) + themetweak
dis <- plot_alpha_diversity_by(psobj=inputps, var="Distance", measures=indices, scalefac=2) + themetweak

# ggsave(l, filename="plots/alpha.shannon.lake.pdf", height=4, width=6)
ggsave(dep, filename="plots/alpha.shannon.depth.pdf", height=4, width=4)
ggsave(dis, filename="plots/alpha.shannon.dist.pdf", height=4, width=4)


###############################
###  Alpha Diversity Stats  ###
###############################

# Linear regression for alpha diversity vs. area. Will need to update for the overall Spearman thing.

# Already done above:
lakemeta <- read.table("data/lake-metadata.txt", header=TRUE, sep="\t", colClasses=c("Lake"="character"))
indices <- c("Shannon")  # "InvSimpson", "Simpson", "Chao1"
inputps <- ps
pdata <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)$data  # Same as p2 above.

alphadivmeans <- pdata %>% group_by(Lake) %>% dplyr::summarize(mean=mean(value))
lakeareas <- lakemeta[c("Lake", "Area")]
mean_vs_area <- left_join(alphadivmeans, lakeareas)
mean_vs_area <- mutate(mean_vs_area, logarea=log(Area))
mean_vs_area <- mutate(mean_vs_area, squarea=sqrt(Area))

cor_test(mean_vs_area, "mean", "Area", method="spearman")    # rho=0.26, p = 0.658  --> rho=0.94? p=0.01
cor_test(mean_vs_area, "mean", "Area", method="pearson")     # r=0.14, p=0.796   0.8  0.056
cor_test(mean_vs_area, "mean", "squarea", method="pearson")  # r=0.16, p=0.76    0.88 0.022
cor_test(mean_vs_area, "mean", "logarea", method="pearson")  # r=0.2,  p=0.698   0.95 0.006

qplot(mean_vs_area$mean, mean_vs_area$logarea)   # R^2 = 0.8759, according to lm(formula = mean_vs_area$mean ~ mean_vs_area$logarea)
ggsave(filename="plots/background/area-diversity-dotplot.pdf", h=3, w=3)

# Generate Spearman correlations of Shannon diversity with lake variables

inputps <- ps
indices <- c("Shannon")

alphadiv <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)$data
lakeshannonmeans <- alphadiv %>% group_by(variable, Lake) %>% dplyr::summarize(mean=mean(value)) %>% spread(variable, mean)
bylake <- left_join(lakemeta, lakeshannonmeans) %>%
          subset(lakemeta$Lake %in% sample_data(ps)$Lake) %>%
          select(-LakeAbbr, -Primer, -Date, -AlphaSort, -NearBelowMixed, -FarBelowMixed, -SetSort) %>%
          remove_rownames() %>% column_to_rownames("Lake")

          # mutate(ShannonRank=rank(Shannon), ShannonNorm=Shannon/sum(Shannon)) %>%

tart <- colorRampPalette(c("red3", "#FFFFBF", "#00441B"))
hm <- heatmap(cor(bylake, method="spearman"), symm=TRUE, Rowv=NA, Colv="Rowv", col=tart(100))

p <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures="Shannon") +
     scale_y_log10() + facet_grid(rows="Compartment") + theme(axis.text.x=element_blank())

q <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures="Shannon") + scale_y_log10()
grid.arrange(p, q)
ggsave("plots/background/alternate-alpha-plots.pdf", h=8, w=8)

# ps_vasc <- dezero(subset_taxa(ps, Vascularity == "Vascular"))
# combined_alpha_plot(inputps=ps_vasc, namefrag="vasc")
# ps_nonvasc <- dezero(subset_taxa(ps, Vascularity == "Nonvascular"))
# combined_alpha_plot(inputps=ps_nonvasc, namefrag="nonvasc")

# Mann-Whitney, which is 2-sample Wilcoxson rank sum. Used in https://www.frontiersin.org/articles/10.3389/fpls.2016.02015/full

inputps <- ps

print("Mann-Whitney test on alpha diversity by depth category:")
inputps1 <- subset_samples(inputps, Depth=="Surface")
inputps2 <- subset_samples(inputps, Depth=="Benthic")
rev(unlist(lapply(indices, function(x) { pval_mann_whitney(x, inputps1, inputps2) } ), recursive=TRUE))

print("Mann-Whitney test on alpha diversity by shore distance category:")
inputps1 <- subset_samples(inputps, Distance=="Offshore")
inputps2 <- subset_samples(inputps, Distance=="Nearshore")
rev(unlist(lapply(indices, function(x) { pval_mann_whitney(x, inputps1, inputps2) } ), recursive=TRUE))

# But now...0.74 and 0.002. Distance is significant again.


########################
###  Beta Diversity  ###
########################

# Lots of stuff about how adonis is specific to the ps object, so wrap it in a function.
# The log+1 normalization with dEcoStand is baked in here.
    # "as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0"
    # Higher bases give less weight to quantities and more to presences, and logbase = Inf gives the presence/absence scaling.
    # Please note this is not log(x+1).

inputps <- ps

# Create "psbig", with only the taxa that appear in five or more samples, and only samples that have 100 or more reads.
taxa_under_five <- colnames(otu_table(inputps))[unlist(lapply(colnames(otu_table(inputps)), function(x) { sum(as.numeric(unname(otu_table(inputps)[,x])>0)) })) < 5]
psbig <- subset_taxa(inputps, !colnames(otu_table(inputps)) %in% taxa_under_five)
psbig <- subset_samples(psbig, sample_sums(psbig) > 100)

pscomplement <- subset_taxa(ps, colnames(otu_table(ps)) %in% taxa_under_five)
table(tax_table(pscomplement)[,"Vascularity"])   # Breakdown of eliminated taxa.

save(file=currentfilename, ps, psbig, psdd, pslake, psddtops, pslaketops, pswet, pswetdd, pswetlake, pswetlakepct)

inputps <- psbig    # or psbig
indices <- c("Shannon")

alphadiv <- plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices)$data
lakeshannonmeans <- alphadiv %>% group_by(variable, Lake) %>% dplyr::summarize(mean=mean(value)) %>% spread(variable, mean)
lakemeta <- read.table("data/lake-metadata.txt", header=TRUE)
bylake <- left_join(lakemeta, lakeshannonmeans) %>%
          subset(lakemeta$Lake %in% sample_data(ps)$Lake) %>%
          select(-LakeAbbr, -Date, -AlphaSort) %>%
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
moarcor <- moarcor %>% select(-Lake, -Primer, -FarBelowMixed, -NearBelowMixed, -SetSort)
bluegold <- colorRampPalette(c("blue", "white", "gold"))
# heatmap(cor(moarcor, method="spearman"), symm=TRUE, Rowv=NA, Colv="Rowv", col=bluegold(100))

lesscor <- as.tbl(moarcor) %>% select(-starts_with("Norm"), -starts_with("Chao"), -starts_with("PD"), -R2DD)

spears <- cor(lesscor, method="spearman")
colbreaks <- c(seq(min(spears), -0.01, length=50), 0, seq(0.01, max(spears), length=50))
pdf("plots/sophisticated-spearman-heatmap.pdf")
heatmap.2(spears, symm=TRUE, Rowv=NULL, Colv="Rowv", col=bluegold(100), breaks=colbreaks,
          trace="none", density.info="none", dendrogram="none",
          colsep=5,
          rowsep=5,
          sepcolor="black",
          sepwidth=c(0.1,0.1) )
dev.off()

# Utilities from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software#correlation-matrix-with-significance-levels-p-value

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

fcm <- flatCorrMatrix()

arrange(fcm, abs(cor)) %>% filter(p<0.1)
# row   column       cor          p
# 1     Area  Shannon 0.7995966 0.05621806
# 2 WestRank MaxDepth 0.8231026 0.04417122
# 3 MaxDepth  Clarity 0.8353977 0.03841103


###################
###  PERMANOVA  ###
###################

inputps <- psbig
bigpermanova <- bdiv(inputps, "Lake+Distance+Depth+Distance:Depth")
runpermanova <- bdiv(inputps, "Run")

f <- file("plots/permanova.txt", "w")
cat("By Lake, Distance, Depth, and Distance x Depth\n", file=f, sep="\n", append=TRUE)
write.table(as.data.frame(bigpermanova$aov.tab), file=f, sep="\t", na="", quote=FALSE, append=TRUE)
cat("\nBy sequencing run\n", file=f, sep="\n", append=TRUE)
write.table(as.data.frame(runpermanova$aov.tab), file=f, sep="\t", na="", quote=FALSE, append=TRUE)
close(f)

# bdiv(inputps, "Lake+Distance+Depth+Compartment", b=Inf)
# bdiv(inputps, "Run", b=Inf)


##############
###  NMDS  ###
##############

# NMDS theme (defined in utility-functions.R)
theme_set(nmdstheme)

# ord.bray.k2 <- ordinate(ps, method="NMDS", distance="bray", k=2, trymax=100)
# ord.bray.k3 <- ordinate(ps, method="NMDS", distance="bray", k=3, trymax=100)
# ord.bray.k4 <- ordinate(ps, method="NMDS", distance="bray", k=4, trymax=100)
# ord.bray.k5 <- ordinate(ps, method="NMDS", distance="bray", k=5, trymax=100)
# ord.bray.k10 <- ordinate(ps, method="NMDS", distance="bray", k=10, trymax=100)
# ord.bray.k20 <- ordinate(ps, method="NMDS", distance="bray", k=20, trymax=100)
# saveRDS(ord.bray.k20, file="data/ord.bray.k20.rds")

inputps <- dezero(ps)
ordfile <- "data/ord.bray.k20.rds"
if(file.exists(ordfile)) {
    ord.bray.k20 <- readRDS(file=ordfile)
} else {
    k=20
    ord.bray.k20 <- ordinate(inputps, method="NMDS", distance="bray", k=20, trymax=100)
    saveRDS(ord.bray.k20, file=ordfile)
}
ord.nmds.bray = ord.bray.k20

lakeplot <- plot_with_hull(inputps, ord.nmds.bray, "Lake") + make_stress_label(ord.nmds.bray)
distplot <- plot_with_hull(inputps, ord.nmds.bray, "Distance")
depthplot <- plot_with_hull(inputps, ord.nmds.bray, "Depth")
compplot <- plot_with_hull(inputps, ord.nmds.bray, "Compartment")
runplot <- plot_with_hull(inputps, ord.nmds.bray, "Run")
primerplot <- plot_with_hull(inputps, ord.nmds.bray, "Primer") + scale_color_manual(values = c("darkgreen", "purple"), aesthetics = c("color", "fill"))

vascplot <- plot_with_hull(inputps, ord.nmds.bray, "Vascularity", type="taxa")
vascplot <- vascplot + scale_color_manual(values = c("green", "darkgreen"), aesthetics = c("color", "fill"))

wetnessplot <- plot_with_hull_ordered(inputps, ord.nmds.bray, "WetlandStatus", type="taxa")

speciesplot <- taxplot_with_names(inputps, ord.nmds.bray, "WetlandStatus", "NCBISciName") +
    scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow"))

ggsave(lakeplot, file="plots/nmds-lake.pdf", width=6.25, height=5)
ggsave(distplot, file="plots/nmds-dist.pdf", width=6.25, height=5)
ggsave(depthplot, file="plots/nmds-depth.pdf", width=6.25, height=5)
ggsave(runplot, file="plots/background/nmds-run.pdf", width=6.25, height=5)
ggsave(primerplot, file="plots/background/nmds-primer.pdf", width=6.25, height=5)
ggsave(vascplot, file="plots/background/nmds-vascularity.pdf", width=6.25, height=5)
ggsave(wetnessplot, file="plots/background/nmds-wetness.pdf", width=6.25, height=5)

ggsave(speciesplot, file="plots/background/nmds-species.pdf", width=6.25, height=5)
ggsave(speciesplot, file="plots/background/nmds-species-large.pdf", width=18, height=18)

# Example for removing k= and stress= text
# primerplot$layers <- primerplot$layers[c(1,2)]
# ...or adding it
# primerplot <- primerplot + make_stress_label(ord.nmds.bray)


# ============== Experimental =============== #

pssp2 <- subset_samples(ps, Lake=="Hummingbird")
ord.bray.sp2 <- ordinate(pssp2, method="NMDS", distance="bray", k=4, trymax=100)
sp2b <- plot_with_hull(pssp2, ord.bray.sp2, "Depth") + make_stress_label(ord.bray.sp2)
sp2b

pssp6 <- subset_samples(ps, Lake=="Raspberry")
ord.bray.sp6 <- ordinate(pssp6, method="NMDS", distance="bray", k=4, trymax=100)
sp6b <- plot_with_hull(pssp6, ord.bray.sp6, "Distance") + make_stress_label(ord.bray.sp6)
sp6b

# or...
depthlakes <- depthplot
depthlakes$layers <- depthlakes$layers[1]
depthlakes <- depthlakes + facet_wrap("Lake", scales="fixed") +
  geom_hline(yintercept=0, color="grey", size=0.5) +
  geom_vline(xintercept=0, color="grey", size=0.5) +
  theme(panel.spacing = unit(1.5, "lines"))

distlakes <- distplot
distlakes$layers <- distlakes$layers[1]
distlakes <- distlakes + facet_wrap("Lake", scales="fixed") +
  geom_hline(yintercept=0, color="grey", size=0.5) +
  geom_vline(xintercept=0, color="grey", size=0.5) +
  theme(panel.spacing = unit(1.5, "lines"))

complakes <- compplot
complakes$layers <- complakes$layers[1]
complakes <- complakes + facet_wrap("Lake", scales="fixed") +
  geom_hline(yintercept=0, color="grey", size=0.5) +
  geom_vline(xintercept=0, color="grey", size=0.5) +
  theme(panel.spacing = unit(1.5, "lines"))

pssp1 <- subset_samples(ps, Lake=="Bay")
ord.bray.sp1 <- ordinate(pssp1, method="NMDS", distance="bray", k=4, trymax=100)
sp1 <- taxplot_with_names(pssp1, ord.bray.sp1, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
  make_stress_label(ord.bray.sp1)
sp1

sp1b <- plot_with_hull(inputps, ord.bray.sp1, "Compartment") + make_stress_label(ord.nmds.bray)
sp1b

pssp2 <- subset_samples(ps, Lake=="Hummingbird")
ord.bray.sp2 <- ordinate(pssp2, method="NMDS", distance="bray", k=4, trymax=100)
sp2 <- taxplot_with_names(pssp2, ord.bray.sp2, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
  make_stress_label(ord.bray.sp2)
sp2

sp2b <- plot_with_hull(inputps, ord.bray.sp2, "Compartment") + make_stress_label(ord.nmds.bray)
sp2b

pssp3 <- subset_samples(ps, Lake=="Inkpot")
ord.bray.sp3 <- ordinate(pssp3, method="NMDS", distance="bray", k=4, trymax=100)
sp3 <- taxplot_with_names(pssp3, ord.bray.sp3, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
  make_stress_label(ord.bray.sp3)
sp3

pssp4 <- subset_samples(ps, Lake=="Long")
ord.bray.sp4 <- ordinate(pssp4, method="NMDS", distance="bray", k=4, trymax=100)
sp4 <- taxplot_with_names(pssp4, ord.bray.sp4, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) + 
  make_stress_label(ord.bray.sp4)
sp4

pssp5 <- subset_samples(ps, Lake=="Morris")
ord.bray.sp5 <- ordinate(pssp5, method="NMDS", distance="bray", k=4, trymax=100)
sp5 <- taxplot_with_names(pssp5, ord.bray.sp5, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
  make_stress_label(ord.bray.sp5)
sp5

pssp6 <- subset_samples(ps, Lake=="Raspberry")
ord.bray.sp6 <- ordinate(pssp6, method="NMDS", distance="bray", k=4, trymax=100)
sp6 <- taxplot_with_names(pssp6, ord.bray.sp6, "WetlandStatus", "SimpleSciName") +
  scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
  make_stress_label(ord.bray.sp6)
sp6


# And why doesn't this work?
bylake_text_nmds <- function(inlake, k) {
  lakeps <- subset_samples(ps, Lake==inlake)
  ord.bray <- ordinate(lakeps, method="NMDS", distance="bray", k=k, trymax=100)
  p <- taxplot_with_names(lakeps, ord.bray, "WetlandStatus", "SimpleSciName") +
    scale_color_manual(values=c("#00AA0033", "turquoise", "blue", "purple", "red", "orange", "yellow")) +
    make_stress_label(ord.bray)
  return(p)
}

bylake_text_nmds(inlake="Raspberry", k=4)
# Error in eval(e, x, parent.frame()) : object 'inlake' not found

# Betadisper

d <- vegdist(otu_table(ps))
betadisper(d, sample_data(ps)[,"Depth"], type="centroid")
# Error in names(y) <- nm[!nas] : 
#   'names' attribute [2] must be the same length as the vector [1]
betadisper(d, sample_data(ps)[,"Depth"], type="median")
# Error in x - c : non-conformable arrays

levels(as.factor(is.na(sample_data(ps)[,"Depth"])))
# [1] "FALSE"

# ============== End experimental =============== #


#######################
###  Venn Diagrams  ###
#######################

font_import(pattern="arial.*", prompt=FALSE)
font_import(pattern="times.*", prompt=FALSE)
loadfonts(device="win")

# Data for taxon Venn
psbigdd <- merge_samples(psbig, "Compartment")
psmed <- subset_taxa(ps, taxa_sums(ps) > 10)
psmeddd <- subset_taxa(psdd, taxa_sums(ps) > 10)
psvasc <- subset_taxa(ps, Vascularity=="Vascular")
psvascdd <- subset_taxa(psdd, Vascularity=="Vascular")
psmedvasc <- subset_taxa(psmed, Vascularity=="Vascular")
pstwovasc <- subset_taxa(psvasc, Org.Sample.Prs != "  1")
psbigvascdd <- subset_taxa(psbigdd, Vascularity=="Vascular")

cols=c("#FF660066", "#0000FF66", "#FFFF0066", "#00FF0066")

inputps <- psbig

psns <- dezero(subset_samples(inputps, Compartment=="Nearshore-Surface"))
psnd <- dezero(subset_samples(inputps, Compartment=="Nearshore-Benthic"))
psfs <- dezero(subset_samples(inputps, Compartment=="Offshore-Surface"))
psfd <- dezero(subset_samples(inputps, Compartment=="Offshore-Benthic"))

ns <- taxa_names(psns)
nd <- taxa_names(psnd)
fs <- taxa_names(psfs)
fd <- taxa_names(psfd)

a <- rownames(tax_table(inputps))
fs <- rownames(tax_table(psfs)); fd <- rownames(tax_table(psfd)); ns <- rownames(tax_table(psns)); nd <- rownames(tax_table(psnd))
faronly <- a[! a %in% c(ns, nd)]; nearonly <- a[! a %in% c(fs, fd)]; deeponly <- a[! a %in% c(ns, fs)]; shallowonly <- a[! a %in% c(nd, fd)]
fsonly <- a[! a %in% c(fd, ns, nd)]; fdonly <- a[! a %in% c(fs, nd, ns)]; nsonly <- a[! a %in% c(nd, fd, fs)]; ndonly <- a[! a %in% c(fs, fd, ns)]

venn.diagram(list("Far-Deep"=fd, "Far-Shallow"=fs, "Near-Deep"=nd, "Near-Shallow"=ns),
    filename="plots/compartment-venn.png", imagetype="png", fill=cols, col="white", cex=2, cat.cex=1.5, margin=.05, fontfamily="Arial")

faronly
nearonly
deeponly
shallowonly

# table(tax_table(dezero(subset_samples(psvasc, sample_names(psvasc)=="Near-Deep")))[,"WetlandStatus"])


#################
###  Treemap  ###
#################

pstmp <- ps
# Simplify tax_table so it doesn't take forever to manipulate.
tax_table(pstmp) <- as.matrix(select(as.data.frame(tax_table(ps)[,c("Taxon.Category", "PlantCat", "Vascularity")]), Vascularity, PlantCat, Taxon.Category, everything()))

pscategories <- pstmp

# General source data.frame for both read and count maps (counts will be added later)
tinput <- data.frame(Vasc=tax_table(pscategories)[,"Vascularity"], TaxCat=tax_table(pscategories)[,"Taxon.Category"], Reads=unname(rowSums(t(otu_table(pscategories)))))

# Thresholds, consolidation, and plotting for read map.
readthresh <- as.integer(sum(tinput[,"Reads"]) / 25)
readbig <- filter(tinput, Reads >= readthresh)
readother <- filter(tinput, Reads < readthresh) %>% group_by(Vascularity) %>% dplyr::summarize(Reads=sum(Reads)) %>% 
    cbind(Taxon.Category="other") %>% select(Vascularity, Taxon.Category, Reads)
fullreadinput <- rbind(readbig, readother) %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(Reads=sum(Reads))

pdf(file="plots/treemap-reads.pdf", h=2.5, w=7.5)
# dev.new(h=2.5, w=7.5)
treemap_jd(fullreadinput, "Reads", c(19, 13))
dev.off()
sum(fullreadinput$Reads)   # ~1.3M

svg(file="plots/treemap-reads.svg", h=2.5, w=7.5)
# dev.new(h=2.5, w=7.5)
treemap_jd(fullreadinput, "Reads", c(19, 13))
dev.off()
sum(fullreadinput$Reads)   # ~1.3M

# And that's the reads one.

# Thresholds, consolidation, and plotting for counts map.

t2input <- tinput %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(CatCount=n())
countthresh <- 10 # as.integer(n(t2input[,"CatCount"]) / 100)
countbig <- as.data.frame(filter(t2input, CatCount >= countthresh))
countother <- filter(t2input, CatCount < countthresh) %>% group_by(Vascularity) %>% dplyr::summarize(CatCount=sum(CatCount)) %>% 
    cbind(Taxon.Category="other") %>% select(Vascularity, Taxon.Category, CatCount)
fullcountinput <- rbind(countbig, countother) %>% group_by(Vascularity, Taxon.Category)

pdf(file="plots/treemap-counts.pdf", h=2.5, w=7.5)
treemap_jd(fullcountinput, "CatCount", c(19, 13))
dev.off()
sum(fullcountinput$CatCount)   # 248

# Tabular version

pstmp <- ps

tax_table(pstmp) <- as.matrix(select(as.data.frame(tax_table(ps)[,c("Taxon.Acc", "Taxon.Category", "PlantCat", "Top.Taxon", "Vascularity")]), Vascularity, Top.Taxon, PlantCat, Taxon.Category, Taxon.Acc, everything()))

pscategories <- pstmp
taxa_names(pscategories) <- tax_table(pscategories)[,"Taxon.Acc"]

tinput <- data.frame(Vasc=tax_table(pscategories)[,"Vascularity"], TaxCat=tax_table(pscategories)[,"Taxon.Category"], Reads=unname(rowSums(t(otu_table(pscategories)))))
counttab <- tinput %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(CatCount=n()) %>% as.data.frame()
readtab <- tinput %>% group_by(Vascularity, Taxon.Category) %>% dplyr::summarize(SumReads=sum(Reads)) %>% as.data.frame()
outtab <- left_join(counttab, readtab)

f <- file("plots/background/taxon-category-counts.txt", "w")
write.table(as.data.frame(outtab), file=f, sep="\t", na="", quote=FALSE, append=TRUE)
close(f)

