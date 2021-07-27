####
# miRSCAPE - Inferring miRNA expression in single-cell clusters
####
# Prepare scmRNA data
#
# To load the Seurat Object
# require(Seurat)
# example = readRDS('example.rds')
#
# clustId = c('Ductal cell type 1', 'Ductal cell type 2')
# denem = modifySeuratObject(pbmc = example, clusterId = clustId)
modifySeuratObject <-
  function(pbmc,
           clusterId = 'all',
           scaled = FALSE) {
    if (scaled) {
      dat = pbmc[["RNA"]]@scale.data
    }
    else{
      dat = data.frame(example[["RNA"]]@data)
    }
    metaDat <- pbmc@meta.data
    
    exp <- list()
    if (clusterId == 'all')
    {
      clusterId = unique(metaDat$cluster)
    }
    for (i in 1:length(clusterId)) {
      aa = which(metaDat$cluster %in% clusterId[i])
      smp = row.names(metaDat[aa, ])
      aa = which(colnames(dat) %in% smp)
      exp[[i]] = rowSums(dat[, aa])
    }
    exp = do.call(cbind, exp)
    colnames(exp) = clusterId
    exp = scale(exp)
    return(exp)
  }

# Prepare bulk data
#
# bulkk = bulkTransform(mrna)
bulkTransform <- function(bulkmRNA) {
  bulkmRNA = bulkmRNA + 1
  bulkmRNA = log(bulkmRNA)
  bulkmRNA = scale(bulkmRNA)
  return(bulkmRNA)
}

# Predict miRNA expression
#
# pred = miRSCAPE(bulkmRNA = bulkk, bulkmiRNA = bulk_miRNA, scmRNA = denem, nrnds = 20)
miRSCAPE <-
  function(bulkmRNA,
           bulkmiRNA,
           scmRNA,
           bstr = 'gbtree',
           objt = "reg:linear",
           mdpth = 4,
           ett = 0.3,
           nrnds = 20,
           echoIn = 10,
           esr = 3) {
    commGenes = intersect(rownames(bulkmRNA), rownames(scmRNA))
    bulk_in = match(commGenes, rownames(bulkmRNA))
    sc_in = match(commGenes, rownames(scmRNA))
    
    bulkmRNA = bulkmRNA[bulk_in, ]
    scmRNA = scmRNA[sc_in, ]
    
    require(xgboost)
    resultt = data.frame()
    for (i in 1:dim(bulkmiRNA)[1]) {
      dtrain <- xgb.DMatrix(data = t(bulkmRNA), label = t(bulkmiRNA[i, ]))
      dtest = scmRNA
      bst <- xgboost(
        data = dtrain,
        booster = bstr,
        objective = objt,
        max.depth = mdpth,
        eta = ett,
        nrounds = nrnds,
        print_every_n = echoIn,
        early_stop_round = esr
      )
      pr <- predict(bst, t(dtest))
      resultt  = rbind(resultt, t(pr))
    }
    row.names(resultt) = rownames(bulkmiRNA)
    return(resultt)
  }