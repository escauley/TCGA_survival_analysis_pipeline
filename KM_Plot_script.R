library(survival)
library(survminer)

# Define the TCGA study name and gene name to generate a KM Plot.
study = 'UCEC'
gene = 'GALNT15'

# Define the Veridis color codes to use. 
veridis_high = "#440154FF"
veridis_med = "#238A8DFF"
veridis_low = "#FDE725FF"

# Define the size of text in the KM Plots
kmplot_theme = theme_bw(base_size = 18)

# Load the data (must be in same working directory as R project) 
data_file <- read.csv(file = paste('TCGA-',study,'/',gene,'.csv',sep = ''), header = TRUE, stringsAsFactors = FALSE)

# Automatically gather the sample size from the data
sample_size = nrow(data_file)/2

# Perform the survival analysis step 1
surv_object <- Surv(time = data_file$survival, event = data_file$outcome)

# Perform the survival analysis step 2
fit1 <- survfit(surv_object ~ group, data = data_file)

# Make the KM Plot
ggsurvplot(fit1, data = data_file,
  pval = TRUE,
  title = gene,
  font.title = c(20, "bold", "black"),
  subtitle = paste('TCGA-',study,sep = ''),
  font.subtitle = c(18, "plain", "black"),
  legend.title =c(""),
  legend.labs = c(paste('High expression N = ',sample_size,sep = ''),paste('Low expression N = ',sample_size,sep = '')),
  xlab = "Time (Days)",
  palette = c(veridis_high,veridis_low),
  ggtheme = kmplot_theme
  )

