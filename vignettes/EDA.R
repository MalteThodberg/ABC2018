## ------------------------------------------------------------------------
library(ABC2018)
library(ggplot2)
theme_set(theme_bw())
library(RColorBrewer)
library(pheatmap)
library(edgeR)

data(zebrafish)

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(zebrafish$Expression > 1)
trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(y=dge, log=TRUE)

## ------------------------------------------------------------------------
# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)

## ------------------------------------------------------------------------
ggplot(data.frame(pca$x, zebrafish$Design), aes(x=PC1, y=PC2, 
                            color=gallein, label=rownames(zebrafish$Design))) +
    geom_text()

## ------------------------------------------------------------------------
ggplot(data.frame(pca$x, zebrafish$Design), aes(x=PC1, y=PC3, 
                            color=gallein, label=rownames(zebrafish$Design))) +
    geom_text()

## ---- eval=FALSE---------------------------------------------------------
#  pheatmap(mat=EM, annotation_col=zebrafish$Design, scale="row")
#  # Running this takes waaaaay to long!

## ---- tidy=TRUE----------------------------------------------------------
pheatmap(mat=EM, kmeans_k=10, annotation_col=zebrafish$Design, scale="row")

## ---- tidy=TRUE----------------------------------------------------------
pheatmap(mat=EM, kmeans_k=100, annotation_col=zebrafish$Design, scale="row")

## ---- tidy=TRUE----------------------------------------------------------
pheatmap(mat=EM, kmeans_k=1000, annotation_col=zebrafish$Design, scale="row")

## ------------------------------------------------------------------------
# Transpose for samples-wise distances
distmat <- dist(t(EM))
distmat

## ---- tidy=TRUE----------------------------------------------------------
pheatmap(mat=as.matrix(distmat), color=brewer.pal(name="RdPu", n=9), clustering_distance_rows=distmat, clustering_distance_cols=distmat, annotation_col=zebrafish$Design, annotation_row=zebrafish$Design)

## ---- tidy=TRUE----------------------------------------------------------
class(UScitiesD)
UScitiesD

## ---- tidy=TRUE----------------------------------------------------------
# Call cmdscale
mds <- cmdscale(UScitiesD, k=2)
mds

## ---- tidy=TRUE----------------------------------------------------------
ggplot(as.data.frame(mds), aes(x=V1, y=V2, label=rownames(mds))) + geom_text()

## ---- tidy=TRUE----------------------------------------------------------
mds <- cmdscale(distmat, k=2)

## ---- tidy=TRUE----------------------------------------------------------
ggplot(data.frame(mds, zebrafish$Design), aes(x=X1, y=X2, color=gallein,
                                              label=rownames(mds))) + geom_text()

## ---- tidy=TRUE----------------------------------------------------------
plotMDS(EM, top=1000)

## ---- tidy=TRUE----------------------------------------------------------
plotMDS(EM, top=10)

## ---- tidy=TRUE----------------------------------------------------------
plotRLDF(EM, nprobes=100, labels.y=zebrafish$Design$gallein, trend=TRUE, robust=TRUE)

## ---- tidy=TRUE----------------------------------------------------------
# Perfrom PLSDA
plsda <- DiscriMiner::plsDA(variables=scale(t(EM)), group=zebrafish$Design$gallein, autosel=FALSE, comps=2)

# Lots of output, we are interested in the components
summary(plsda)

## ------------------------------------------------------------------------
ggplot(data.frame(plsda$components, zebrafish$Design), aes(x=t1, y=t2, 
        color=gallein, label=rownames(zebrafish$Design))) + geom_text()

## ------------------------------------------------------------------------
data(pasilla)

## ------------------------------------------------------------------------
data(tissues)

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(pasilla$Expression > 1)
trimmed_em <- subset(pasilla$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(y=dge, log=TRUE)

## ------------------------------------------------------------------------
plotDensities(EM)

## ------------------------------------------------------------------------
# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Save data for plotting
P <- data.frame(pasilla$Design, pca$x)

# Inspect components
summary(pca)

## ------------------------------------------------------------------------
qplot(data=P, x=PC1, y=PC2, geom=c("text"), 
			color=type, label=rownames(P))

## ------------------------------------------------------------------------
qplot(data=P, x=PC1, y=PC2, geom=c("text"), 
			color=condition, label=rownames(P))

## ------------------------------------------------------------------------
# Perform MDS, but only save output
mds <-plotMDS(EM, top=100, gene.selection="common", plot=FALSE)

# Save as a data.frame
P <- data.frame(pasilla$Design, mds$cmdscale.out)

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, geom=c("text"), 
			color=type, label=rownames(P))

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, geom=c("text"), 
			color=condition, label=rownames(P))

## ---- tidy=TRUE----------------------------------------------------------
# LDF with only 100 genes
ldf <- plotRLDF(EM, nprobes=100, labels.y=pasilla$Design$condition, trend=TRUE, robust=TRUE, plot = FALSE)

# Save as a data.frame
P <- data.frame(pasilla$Design, ldf$training)

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, geom=c("text"), 
			color=type, label=rownames(P))

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, geom=c("text"), 
			color=condition, label=rownames(P))

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(tissues$Expression > 2)
trimmed_em <- subset(tissues$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(y=dge, log=TRUE)

## ------------------------------------------------------------------------
plotDensities(EM, legend = FALSE)

## ------------------------------------------------------------------------
# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Save data for plotting
P <- data.frame(tissues$Design, pca$x)

# Inspect components
summary(pca)

## ------------------------------------------------------------------------
qplot(data=P, x=PC1, y=PC2, geom="text", 
			color=gender, label=tissue.type)

