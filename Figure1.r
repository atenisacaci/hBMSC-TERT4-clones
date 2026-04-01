##Figure 1C
# Barplot of fraction of bone formed for each cell line 
Implants <- read.xlsx("Implants_all_clones.xlsx", 
                      sheetIndex = 1, header=TRUE)
names(Implants)[2] <- "Data"

library(dplyr)
df.summary <-Implants %>%
  group_by(Clones) %>%
  summarise(
    sd = sd(Data, na.rm = TRUE),
    Data = mean(Data)
  )
df.summary

library(ggplot2)
# Default bar plot
conditions <- c("AD10", "DD8", "CB", "CD")

ggplot(Implants, aes(Clones, Data)) + scale_x_discrete(limits = conditions)+ 
  geom_bar(stat = "identity", data = df.summary,
           fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.2),
               color = "black") + 
  geom_errorbar(
    aes(ymin = Data-sd, ymax = Data+sd),
    data = df.summary, width = 0.2) 
	
# Final aesthetics, such as colors and line drawing, were done in Illustrator
