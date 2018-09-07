library(ThodbergMisc)
collections("Genomics")
collections("DE")
collections("Hadley")
library("tidyverse")

#### IBD RNA-Seq ####

load(url("http://duffel.rail.bio/recount/v2/SRP042228/rse_gene.Rdata"))

# Trim
SE <- rse_gene

# Groups
raw <- do.call(rbind, SE$characteristics) %>% as.data.frame

SE$Group <- str_split(raw$V5, pattern = " ", simplify = TRUE)[,2] %>%
    factor %>%
    fct_collapse(Ctrl=c("Not", "not")) %>%
    fct_relevel("Ctrl")

SE$Type <- str_split(raw$V6, pattern = " ", simplify = TRUE)[,3] %>%
    factor %>%
    fct_collapse(Ctrl=c("Not", "not")) %>%
    fct_relevel("Ctrl")

SE$Gender  <- str_split(raw$V2, pattern = " ", simplify = TRUE)[,2] %>% factor
SE$Age  <- str_split(raw$V3, pattern = " ", simplify = TRUE)[,4] %>% as.numeric()

#SE <- subset(SE, select=Type != "iCD" & Age <= 9)
SE <- SE[rowSums(rpkm(assay(SE), gene.length = rowData(SE)$bp_length) >= 5) >= 8,]

#SE$Montreal  <- str_split(raw$V6, pattern = " ", simplify = TRUE)[,3] %>% factor(exclude="NA") %>% fct_infreq()

# Normalize
dge <- assay(SE) %>% DGEList(genes=rowData(SE)) %>% calcNormFactors()
EM <- cpm(dge, log=TRUE, prior.count = 5)
plotDensities(EM, group=SE$Group)
# MDS <- plotMDS(EM, top=100, plot=FALSE)
# data.frame(colData(SE), MDS$cmdscale.out) %>%
#     ggplot(aes(x=X1, y=X2, color=Group)) +
#     geom_point(aes(size=Weight))

mod <- model.matrix(~Gender+Type, data=colData(SE))

#v <- voomWithQualityWeights(dge, design=mod, plot=TRUE)
#v <- voom(dge, design=mod)
v <- EM
#SE$Weight <- v$targets$sample.weights
fit0 <- lmFit(v, design = mod)
aw <- SE$Weight <- arrayWeightsQuick(v, fit)
fit <- lmFit(v, design=mod, weights = aw)
eb <- eBayes(fit, trend=TRUE, robust=TRUE)
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
    ggplot(aes(x=PC3, y=PC4, color=Type, size=Weight)) +
    geom_point()

data.frame(colData(SE), pca$x) %>%
    as_tibble() %>%
    gather(key="PC", value="Score", starts_with("PC")) %>%
    filter(PC %in% paste0("PC", 1:9)) %>%
    ggplot(aes(x=Type, y=Score, weight=Weight)) +
    geom_boxplot() +
    facet_wrap(~PC)


