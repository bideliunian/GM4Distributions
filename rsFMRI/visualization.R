#########################################
# visualize the functional signals
#############################################

library(ggplot2)
library(tidyr)

file_path <- "D:/Research/GM4Distributions/processed_data/ADHD200_CC200_TCs_filtfix/NYU"
reading_path <- "D:/Research/GM4Distributions/processed_data/NYU_processed"

id <- "0010001"
node <- '100'
load(paste(reading_path, "/", id, ".RData", sep = ""))

signal_df <- as.data.frame(t(processed_list[[node]]))
signal_mean <- read.delim(paste(file_path, "/", id, "/sfnwmrda", 
                                id,"_session_1_rest_1_cc200_TCs.1D", sep=""))
signal_mean_test <- signal_mean[paste('Mean_', node, sep = "")]

names(signal_df) <- gsub('V', 'Signal', names(signal_df))
signal_df <- data.frame(SignalMean = rowMeans(signal_df), Time = seq(0, length = nrow(signal_df)) * 2, signal_df)


# check processed mean signal is equal to the one provided by Athena pipline
# all.equal(signal_df$SignalMean, signal_mean_test$Mean_40)

signal_df_long <- signal_df %>%
  pivot_longer(starts_with("Signal"), names_to = "SignalName", values_to = "Value")


p <- ggplot(signal_df_long, aes(x = Time, y = Value, color = SignalName)) +
  geom_line() +
  xlab("Time") +
  ylab("Signal Value") +
  ggtitle("BOLD Signal Visualization") +
  theme_minimal() + 
  theme(legend.position = "none")

p + geom_line(data = signal_df_long[signal_df_long$SignalName == "SignalMean",], size = 1.5, color = 'black') 


