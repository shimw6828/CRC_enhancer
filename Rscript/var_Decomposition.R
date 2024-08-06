#0. Set working directory, change as needed
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1008-y
setwd("/home/shimw/project/enhancer_map/conserved/")
library(lme4)
library(optimx)

datasets = as.data.frame(scan("Stanford_datasets.txt",list(setname="",seqBatch="",species="",tissue=""),sep="\t"))
rawCounts <- as.matrix(read.table('Stanford_datasets_fpkmMat.txt',header=FALSE,sep='\t'))
colnames(rawCounts) <- datasets$setname
geneDetails <- as.data.frame(scan("ortholog_GC_table.txt",skip=1,list(mouse_name="",mouse_GC = 0.0,human_name = "",human_GC=0.0)))
rownames(rawCounts) <- geneDetails$human_name 
rownames(datasets) <- datasets$setname
rownames(geneDetails) <- geneDetails$human_name 


# mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# annotLookup <- getBM(mart = mart,
#                      attributes = c("external_gene_name","ensembl_gene_id", "gene_biotype"),
#                      filters = "external_gene_name", values = geneDetails$human_name,
#                      uniqueRows=TRUE)
# 
# names(annotLookup)[2] = "Gene.stable.ID"
# H0M_logFC1 = dplyr::left_join(H0M_logFC1, annotLookup)
mat <- rawCounts[apply(rawCounts,1,max) >= .1,]
mat.vd <- log10(mat + 0.01)
genes.filtered <- geneDetails$human_name [apply(rawCounts,1,max) >= .1]


tissue <- unname(sapply(colnames(rawCounts), function(x){
  stringr::str_split(x, pattern = " \\(")[[1]][1]
}))
species <- unname(sapply(colnames(rawCounts), function(x){
  if (stringr::str_split(x, pattern = " \\(")[[1]][2]=="h)") {
    return("human")
  }else{
    return("mouse")
  }
}))

var.explain <- matrix(, nrow = length(genes.filtered), ncol = 2)
for(i in 1:length(genes.filtered)){
  model.dat <- data.frame(geneExp = unname(unlist(mat.vd[i,])), species = species, tissue = tissue)
  mixed.lmer <- lmer(geneExp ~ (1|species) + (1|tissue), data = model.dat, REML = F,
                     control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  var.explain[i,] <- (c(as.data.frame(VarCorr(mixed.lmer))$vcov[1], as.data.frame(VarCorr(mixed.lmer))$vcov[2]) / sum(as.data.frame(VarCorr(mixed.lmer))$vcov))
}

percent <- function(x){
  if(as.numeric(x["V1"]) > 0.60){
    return("High across species")
  }else if(as.numeric(x["V2"]) > 0.60){
    return("High across tissues")
  }else{
    return("None")
  }
}
var.explain = as.data.frame(var.explain)
var.explain$state = apply(var.explain, 1, percent)
var.explain$gene_id <- genes.filtered
colnames(var.explain) <- c("tissues", "species", "state", "gene_id")

ggplot(var.explain, aes(x = tissues, y = species, col = state)) + geom_point()

#var.explain[var.explain$state=="High across species",]
readr::write_csv(var.explain[var.explain$state=="High across species",],
                 "/home/shimw/project/enhancer_map/conserved/species-specific.csv")
readr::write_csv(var.explain,
                 "/home/shimw/project/enhancer_map/conserved/var.csv")


mat[,1:13]
mat[,14:26]
library(TissueEnrich)
library(SummarizedExperiment)
se<-SummarizedExperiment(assays = SimpleList(as.matrix(mat[,1:13])),rowData = row.names(mat[,1:13]),colData = colnames(mat[,1:13]))
output<-teGeneRetrieval(se)
head(assay(output))
aaa = as.data.frame(assay(output))
aaa$Tissue

test <- load("/home/shimw/project/combine-expression.rda")
test
tissue_specific = dataset[["ENCODE Dataset"]][["tissueSpecificGenes"]]

tissue_specific$Tissue=="Intestine"
tissue_specific$Group=="Tissue-Enriched"
colon_tissue_specific <- tissue_specific[tissue_specific$Tissue=="Intestine"&tissue_specific$Group=="Tissue-Enriched",]
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]
colon_tissue_specific[colon_tissue_specific$Gene%in%HOM_one2one$Mouse.gene.stable.ID,]
readr::write_csv(colon_tissue_specific[colon_tissue_specific$Gene%in%HOM_one2one$Mouse.gene.stable.ID,],
                 "/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
##human
tissue_specific = dataset[["Protein-Atlas"]][["tissueSpecificGenes"]]
unique(tissue_specific$Tissue)
tissue_specific$Tissue=="Colon"
tissue_specific$Group=="Tissue-Enriched"
colon_tissue_specific <- tissue_specific[tissue_specific$Tissue=="Colon"&tissue_specific$Group=="Tissue-Enriched",]
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]
colon_tissue_specific[colon_tissue_specific$Gene%in%HOM_one2one$Gene.stable.ID,]
readr::write_csv(colon_tissue_specific[colon_tissue_specific$Gene%in%HOM_one2one$Mouse.gene.stable.ID,],
                 "/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
