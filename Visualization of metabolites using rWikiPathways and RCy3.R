##Visualization of metabolite accumulation using pathway taken from Wikipathways
##Ken Kamiya, wrote at 2019/06/25

#prepare these packages
install.packages(C("rWikipathways", "RCy3"))
library(rWikiPathways)
library(RCy3)

#get graphIDs from certain pathway from WikiPathway
xrefs <- getXrefList(pathway="WP3622", systemCode="Ck")
xrefs.graphIds = NULL
for(i in 1:length(xrefs)){
  pathways <- findPathwaysByXref(xrefs[[i]], systemCode="Ck") 
  targeted.pathway <- pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")][[1]]
  xrefs.graphIds[[i]] <- as.character(targeted.pathway$fields$graphId$values)
}


#prepare MetaboAnalystR
install.packages("pacman")
library(pacman)
pacman::p_load(Rserve, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, xcms, CAMERA, fgsea, MSnbase, BiocParallel, metap, reshape2, scales)

install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
library(MetaboAnalystR)


#prepare metabolite data
library(DiffCorr)
data("AraMetLeaves") #可視化用データ用意

##get metabolite ids lists using MetaboAnalystR
mSet<-InitDataObjects("NA", "utils", FALSE)
cmpd.vec<- rownames(AraMetLeaves)
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name", T, T, T, T, T);
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "Inositol-1-phosphate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Inositol-1-phosphate", "Inositol phosphate");
mSet<-PerformDetailMatch(mSet, "D-Glucose-6-phosphate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "D-Glucose-6-phosphate", "Beta-D-Glucose 6-phosphate");

##get graphIDs from prepared metabolite data
mSet.map.table = mSet$dataSet$map.table
mSet.graphIds = NULL
for(i in 1:length(mSet.map.table[,"ChEBI"])){
  pathways <- findPathwaysByXref(mSet.map.table[,"ChEBI"][[i]], systemCode="Ce")
  targeted.pathway = pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")]
  if(length(targeted.pathway)==0){
    mSet.graphIds[[i]] = "NA"
  }
  if(!length(targeted.pathway)==0){
    mSet.graphIds[[i]] <- targeted.pathway[[1]]$fields$graphId$values
  }
}

##Calculate log2FC from Col0.1 and mto1.1
log2FC = log2(AraMetLeaves[,"mto1.1"] / AraMetLeaves[,"Col0.1"])
min.log2FC = min(log2FC,na.rm=TRUE)
max.log2FC = max(log2FC,na.rm=TRUE)
abs.log2FC = max(abs(min.log2FC),max.log2FC)
data.values = c(-abs.log2FC,0,abs.log2FC)

##make matrix
mSet.map.table.graphId.log2FC = data.frame(mSet.map.table, mSet.graphIds, log2FC)


#connect to Cytoscape and visualize content of metabolites
cytoscapePing()
RCy3::commandsRun('wikipathways import-as-pathway id=WP3622') 
toggleGraphicsDetails()
loadTableData(mSet.map.table.graphId.log2FC, data.key.column = "mSet.graphIds", table.key.column = "GraphID")
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2FC", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
