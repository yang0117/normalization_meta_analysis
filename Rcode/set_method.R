setMethod("rma", "SnpCnvFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            sql <- paste("SELECT fsetid, fid FROM pmfeatureCNV")
            probeInfo <- dbGetQuery(db(object), sql)
            probeInfo <- probeInfo[order(probeInfo[["fsetid"]]),]
            pms <- exprs(object)[probeInfo[["fid"]],,drop=FALSE]
            exprs <- basicRMA(pms, probeInfo[["fsetid"]], normalize, background)
            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
          })