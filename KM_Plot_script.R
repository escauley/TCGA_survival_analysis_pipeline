library(survival)
library (survminer)

study = 'BRCA'

gene = 'PYGM'
  
data_file <- read.csv(file = paste('TCGA-',study,'/',gene,'.csv',sep = ''), header = TRUE, stringsAsFactors = FALSE)

sample_size = nrow(data_file)/2

surv_object <- Surv(time = data_file$survival, event = data_file$outcome)

fit1 <- survfit(surv_object ~ group, data = data_file)

ggsurvplot(fit1, data = data_file,
           pval = TRUE,
           title = gene,
           font.title = c(12, "bold", "black"),
           subtitle = paste('TCGA-',study,sep = ''),
           font.subtitle = c(10, "plain", "black"),
           legend.title =c(""),
           legend.labs = c(paste('High expression N = ',sample_size,sep = ''),paste('Low expression N = ',sample_size,sep = '')),
           xlab = "Time (Days)"
           )
