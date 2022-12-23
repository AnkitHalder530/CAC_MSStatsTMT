library(MSstatsTMT)
raw.pd <- read.delim("PD_CAC.txt")
head(raw.pd)
annotation.pd <- read.csv(file="PD_Annotation_file.csv", header=TRUE)
head(annotation.pd)
input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd,rmPSM_withMissing_withinRun = TRUE,useNumProteinsColumn = TRUE,useUniquePeptide = TRUE)
head(input.pd)
write.csv(input.pd,"Input.pd.csv")
quant.msstats <- proteinSummarization(input.pd,method = "msstats",global_norm = TRUE,reference_norm = TRUE)
head(quant.msstats)
x<-quant.msstats
write.csv(x,"Quant_msstats.csv")
dataProcessPlotsTMT(data=quant.msstats,
                    type='ProfilePlot',
                    width = 35,
                    height = 10)
dataProcessPlotsTMT(data=quant.msstats,
                    type='QCPlot',
                    width = 35,
                    height = 10)
x<-quant.msstats$ProteinLevelData
write.csv(x,"Protein_datawoMBimpute.csv")
head(x)
data<-quant.msstats
test.pairwise_1 = groupComparisonTMT(quant.msstats,
                                     contrast.matrix = "pairwise",
                                     moderated = TRUE,
                                     adj.method = "BH",
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE,
                                     save_fitted_models = TRUE,
                                     use_log_file = TRUE,
                                     append = FALSE,
                                     verbose = TRUE,
                                     
)
head(test.pairwise_1$ComparisonResult)
Y<-test.pairwise_1$ComparisonResult
write.csv(Y,"Comparison_result_1.csv")