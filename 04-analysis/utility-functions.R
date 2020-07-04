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
library("Hmisc"); library("rstatix")                   # Statistical additions

# General functions and settings

dezero <- function(psobj) {
    psnew <- prune_taxa(taxa_sums(psobj) > 0, psobj)
    psnew <- prune_samples(sample_sums(psnew) > 0, psnew)
    return(psnew)
    }


combined_alpha_plot <- function(inputps, namefrag) {
    p <- grid.arrange(plot_alpha_diversity_by(psobj=inputps, var="Lake", measures=indices),
                      plot_alpha_diversity_by(psobj=inputps, var="Depth", measures=indices),
                      plot_alpha_diversity_by(psobj=inputps, var="Distance", measures=indices),
                      layout_matrix = rbind(c(1, 1), c(2, 3)))
    fname <- paste("plots/background/alpha", namefrag, "pdf", sep=".")
    ggsave(p, filename=fname, width=8, height=8)
}

four_compartments_ordered <- c("Offshore-Benthic", "Offshore-Surface", "Nearshore-Benthic", "Nearshore-Surface")
four_shapes_ordered <- c(19, 21, 15, 22)
light="deepskyblue"; dark="blue"
four_colors_ordered <- c(dark, light, dark, light)
shapedf <- reshape2::melt(data.frame(k=four_compartments_ordered, v=four_shapes_ordered), variable.name="k", value.name="v")
shapes <- setNames(as.character(shapedf$v), shapedf$k)
colordf <- reshape2::melt(data.frame(k=four_compartments_ordered, v=four_colors_ordered), variable.name="k", value.name="v")
colors <- setNames(as.character(colordf$v), colordf$k)

plot_alpha_diversity_by <- function(var, psobj, measures=c("Shannon", "Simpson", "InvSimpson", "Chao1"), scalefac=1) {
    p <- plot_richness(psobj, x=var, measures=measures, title="")
    p$layers <- p$layers[-1]  # Remove the default set of points since we can't manipulate them
    # Now, instead of individual points, do a boxplot.
    p <- p +   # theme(axis.text.x = element_text(face="italic")) +
        geom_point(position=position_jitter(w = 0.2, h = 0), aes(shape=Compartment, color=Compartment), size=scalefac, stroke = 1.2*scalefac) + 
        stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5, show.legend=FALSE) +
        labs(x=var, y="Shannon Index (H')", title=NULL, main=NULL) +
        scale_discrete_manual("Compartment", values=as.integer(shapes), labels=four_compartments_ordered, aesthetics="shape") +
        scale_color_manual("Compartment", values=colors, labels=four_compartments_ordered, aesthetics=c("color", "fill")) +
        labs(shape = "Compartment") + 
        theme(text = element_text(size=12), strip.text.x = element_blank())
    return(p)
}

plot_with_hull <- function(psobj, ord, var, type="samples") {
    if (var=="Lake") {
        p <- plot_ordination(psobj, ord, color=var, type=type)  # + guides(col = guide_legend(ncol=1))
    } else {
        p <- plot_ordination(psobj, ord, color=var, shape=var, type=type)
    }
    p$data <- na.omit(p$data)
    hulls <- plyr::ddply(na.omit(p$data), var, find_hull)
    p <- p + geom_polygon(data = hulls, alpha = 0.1, aes(fill=eval(as.name(paste(var))))) + 
             guides(fill=FALSE) + labs(color=var, fill=var) # + make_stress_label(ord)
    return(p)
    }

plot_with_hull_ordered <- function(psobj, ord, var, type="samples") {
    psobj <- subset_taxa(psobj, WetlandStatus != "UPL")
    p <- plot_ordination(psobj, ord, color=var, shape=var, type=type)
    p$data <- na.omit(p$data)
    hulls <- plyr::ddply(na.omit(p$data), var, find_hull)
    p$data$WetlandStatus <- factor(p$data$WetlandStatus, levels=c("ALG", "AQU", "OBL", "FACW", "FAC", "FACU"))
    p <- p + geom_polygon(data = hulls, alpha = 0.1, aes(fill=eval(as.name(paste(var))))) +
             guides(fill=FALSE) + labs(color=var, fill=var) # + make_stress_label(ord)
    return(p)
    }

taxplot_with_names <- function(psobj, ord, var, textfield, type="taxa") {
    p <- plot_ordination(psobj, ord, color=var, type=type, label=textfield) + guides(col = guide_legend(ncol=1))
    return(p)
}

lakeplot_with_hull_ordered <- function(psobj, ord, var, type="samples") {
    p <- plot_ordination(psobj, ord, color=var)
    p$data <- na.omit(p$data)
    hulls <- plyr::ddply(na.omit(p$data), var, find_hull)
    sixDs <- c("#DDDDDD", "#DDDDDD", "#DDDDDD", "#DDDDDD", "#DDDDDD", "#DDDDDD")
    bolder <- c("#F60000", "#FF6D01", "#FFC404", "#46AF3B", "#4198FF", "#4815AA")
    lakecols <- c(sixDs, bolder)    # Reverse these two to make the version with the other six grey
    p$data$Lake <- factor(p$data$Lake, levels=lakesordered)  # lakesordered is defined outside the function
    p <- p + geom_polygon(data = hulls, alpha = 0.1, aes(fill=Lake)) +
        labs(color=var, fill=var) + # make_stress_label(ord) +
        scale_colour_manual(values = lakecols, aesthetics = c("colour", "fill"))
    return(p)
    }

    # This is what goes in geom_polygon to make the lake hulls come out in pretty colors. fill=var doesn't work; why?
    # aes(fill=Lake) or aes(fill=as.name(var))

plot_bar_jd <- function (physeq, x = "Sample", y = "Abundance", ylab = "Abundance", fill = NULL, title = NULL, facet_grid = NULL) 
{
    theme_set(nmdstheme)

    mdf = psmelt(physeq)
    p = ggplot(mdf, aes_string(x = reorder(mdf$Sample, mdf$Primer), y = y, fill = fill))
    p = p + geom_bar(stat = "identity", position = "stack" ) + scale_colour_manual(values = twogreens, aesthetics = c("colour", "fill"))
    p = p + theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=0))
    p = p + ylab(ylab) + xlab(NULL)
    if (!is.null(facet_grid)) {
        p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

treemap_jd <- function(fulldata, fieldname, fontsizes=c(20, 16), colors=twogreens) {
    fulldata <- fulldata %>% mutate("SortCol"=rev(rank(Vascularity)))
    tm2 <- treemap(fulldata, index=c("Vascularity", "Taxon.Category"), vSize=fieldname, type="categorical",
                   title="", position.legend="none",
                   fontsize.labels=fontsizes, fontcolor.labels=c("white", "black"),
                   fontface.labels=c(4, 1), bg.labels=c("transparent"),
                   align.labels=list(c("left", "top"), c("right", "bottom")), overlap.labels=0.5, inflate.labels=F,
                   border.lwds=c(5, 2), border.col=c("#385723", "white"),
                   vColor="Vascularity", palette=colors,
                   aspRatio=3, algorithm="pivotSize", sortID="SortCol", draw=T)
    return(tm2)
}

pval_mann_whitney <- function(indexname, inputps1, inputps2) {
    rich1 <- unlist(unname(estimate_richness(inputps1, measures=indexname)))
    rich2 <- unlist(unname(estimate_richness(inputps2, measures=indexname)))
    output <- wilcox.test(rich1, rich2, correct=FALSE)$p.value
    names(output) <- indexname
    return(output)
}

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

# Specific utility functions and definitions used by the above

nmdstheme <- theme_set(theme_minimal(base_size = 12)) +
             theme_update(plot.background = element_rect(color="black", fill="white", size=0.25),
                          plot.margin = unit(c(.5, .5, .5, .5), "cm"))

make_stress_label <- function(ord) {
    thing <- annotate("text", label=paste("k=", ord$ndim, "; Stress=", as.character(format(ord$stress, digits=5)), sep=""),
                      x=-Inf, y = Inf, vjust=1.5, hjust=0)
    return(thing) }

make_nmds_title <- function(datasetname) {
    return(textGrob(label = paste("NMDS on ", datasetname, sep=""),
    x = unit(0, "lines"), y = unit(0, "lines"), hjust = 0, vjust = 0, gp = gpar(fontsize = 14))) }

make_alpha_title <- function(datasetname) {
    return(textGrob(label = paste("Alpha diversity for ", datasetname, sep=""),
    x = unit(0, "lines"), y = unit(0, "lines"), hjust = 0, vjust = 0, gp = gpar(fontsize = 14))) }

# Function to generate "hulls" around different sets of points in NMDS data.
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

a_lakes <- c("Bergner", "Bolger", "Brown", "Crampton", "Cranberry", "Reddington")
b_lakes <- c("Bay", "Hummingbird", "Inkpot", "Long", "Morris", "Raspberry")
lakesordered <- c(a_lakes, b_lakes)

twogreens=c("#A9D18E", "#61953D")
