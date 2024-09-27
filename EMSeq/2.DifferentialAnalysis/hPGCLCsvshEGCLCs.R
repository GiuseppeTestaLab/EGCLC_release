library(DSS)

#to modify
Type_1 <- "hEGCLCs"
Type_2 <- "hPGCLCs"

#to not modify

bsseq_obj <- readRDS("~/DataDir/2.DifferentialAnalysis/Input/bsseq_obj_sharedby75ofall.rds") #getting bsseq object where sample selection and CpG filtering were already performed 

group_1 <- pData(bsseq_obj)[pData(bsseq_obj)$Type %in% Type_1, "SampleID"]
group_2 <- pData(bsseq_obj)[pData(bsseq_obj)$Type %in% Type_2, "SampleID"]

dmlTest <- DSS::DMLtest(bsseq_obj, group1=group_1,
                        group2=group_2, smoothing = TRUE)  #performing differential methylation test

output_file_name <- paste0("~/DataDir/2.DifferentialAnalysis/Output/", Type_1, "vs", Type_2, ".rds")

saveRDS(dmlTest, file = output_file_name) #saving output of the test

