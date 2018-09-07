# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
# 	* Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

library(forcats)

#### Rat neuron ####

# Get data
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData")
load(file=con)
close(con)
bot <- hammer.eset

rat <- list(Design=pData(bot),
									 Annotation=fData(bot),
									 Expression=as.matrix(exprs(bot)))

names(rat$Expression) <- NULL

# fix typo
rat$Design$Time <- gsub(x=rat$Design$Time, pattern="2months", replacement="2 months")

devtools::use_data(rat, overwrite = TRUE)

#### Asmann Tissues ####

# Get data
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bot = bodymap.eset

tissues <- list(Design=pData(bot),
											Annotation=fData(bot),
											Expression=as.matrix(exprs(bot)))

names(tissues$Expression) <- NULL

devtools::use_data(tissues, overwrite = TRUE)

#### NEW ??? ####

load("http://duffel.rail.bio/recount/v2/SRP051825/rse_gene.Rdata")

# Trim
SE <- rse_gene
SE <- SE[rowSums(rpkm(assay(SE), gene.length = rowData(SE)$bp_length) >= 1) >= 6,]

SE$Group <- c(rep("NEC", 9), rep("Ctrl", 5)) %>% factor
SE$Gender <- c("M", "M", "F", "M", "F", "F", "M", "F", "F",
               "F", "F", "M", "M", "F")
SE$GA <- do.call(rbind, SE$characteristics) %>% as.data.frame %>% use_series(V2) %>% str_remove("gestational age:")

# Normalize
dge <- assay(SE) %>% DGEList(genes=rowData(SE)) %>% calcNormFactors()
EM <- cpm(dge, log=TRUE)
plotDensities(EM, group=SE$Group)
plotMDS(EM, top=1000, col=as.integer(SE$Group), pch=19)

mod <- model.matrix(~Gender+Group, data=colData(SE))
SVs <- wsva(EM, design=mod, n.sv = 5, weight.by.sd = TRUE)
mod <- cbind(mod, SVs)

# edgeR
disp <- estimateDisp(dge, design=mod, robust=TRUE)
fit <- glmQLFit(y=disp, design=mod, robust=TRUE)
res <- glmQLFTest(fit, coef=3)
summary(decideTests(res))

v <- voomWithQualityWeights(dge, design=mod, plot=TRUE)
#v <- voom(dge, design=mod)
#SE$Weight <- v$targets$sample.weights
fit <- lmFit(v, design = mod)
eb <- eBayes(fit, trend=FALSE, robust=TRUE)
dt <- decideTests(eb)
summary(dt)
topTable(eb, coef=3)

pca <- prcomp(t(EM), scale=TRUE)
data.frame(colData(SE), pca$x) %>%
    ggplot(aes(x=PC1, y=PC2, color=Gender)) +
    geom_point()




load(url("http://duffel.rail.bio/recount/v2/SRP055438/rse_gene.Rdata"))

# Trim
SE <- rse_gene
SE <- SE[rowSums(rpkm(assay(SE), gene.length = rowData(SE)$bp_length) >= 1) >= 3,]

# Groups
raw <- do.call(rbind, SE$characteristics) %>% as.data.frame

SE$Group <- ifelse(raw$V3 == "disease: Crohn's Disease", "Crohns", "Ctrl") %>%
    factor %>% relevel("Ctrl")
SE$Gender <- str_split(raw$V2, pattern = " ", simplify = TRUE)[,2] %>% factor
SE$Age  <- str_split(raw$V1, pattern = " ", simplify = TRUE)[,2] %>% as.integer
SE$Montreal  <- str_split(raw$V6, pattern = " ", simplify = TRUE)[,3] %>% factor(exclude="NA") %>% fct_infreq()
SE$Stage <- str_split(raw$V4, pattern = " ", simplify = TRUE)[,2] %>% factor %>% relevel("NA")

# Normalize
dge <- assay(SE) %>% DGEList(genes=rowData(SE)) %>% calcNormFactors()
EM <- cpm(dge, log=TRUE, prior.count = 5)
plotDensities(EM, group=SE$Group)
MDS <- plotMDS(EM, top=5000, plot=FALSE)
data.frame(colData(SE), MDS$cmdscale.out) %>%
    ggplot(aes(x=X1, y=X2, color=Stage)) +
    geom_point(aes(size=Weight))

mod <- model.matrix(~Gender+Stage, data=colData(SE))

v <- voomWithQualityWeights(dge, design=mod, plot=TRUE)
#v <- voom(dge, design=mod)
SE$Weight <- v$targets$sample.weights
fit <- lmFit(v, design = mod)
eb <- eBayes(fit, trend=FALSE, robust=TRUE)
dt <- decideTests(eb)
summary(dt)
topTable(eb, coef=3)



mod <- model.matrix(~Gender+Group, data=colData(SE))
SVs <- wsva(EM, design=mod, n.sv = 5, weight.by.sd = TRUE)
mod <- cbind(mod, SVs)


# edgeR
disp <- estimateDisp(dge, design=mod, robust=TRUE)
fit <- glmQLFit(y=disp, design=mod, robust=TRUE)
res <- glmQLFTest(fit, coef=3)
summary(decideTests(res))

pca <- prcomp(t(EM), scale=TRUE)
data.frame(colData(SE), pca$x) %>%
    ggplot(aes(x=PC1, y=PC2, color=Stage, size=Weight)) +
    geom_point()

