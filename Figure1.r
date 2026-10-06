##Figure 1C
#Load the following objects from https://osf.io/wxpgn/
# Bar plot of fraction of bone formed for each cell line 
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
##Figure 1E
#Cell morphology
HBF <- read_delim("HBF_nuclei.csv",delim = ";",locale = locale(encoding = "UTF-8"),show_col_types = FALSE)

LBF <-read_delim("LBF_nuclei.csv",delim = ";",locale = locale(encoding = "UTF-8"),show_col_types = FALSE)

HBF$group <- "HBF"
LBF$group   <- "LBF"

df <- bind_rows(HBF, LBF)
df$length <- df$Major
df$width  <- df$Minor
df$size   <- df$Area
df$elongation <- df$length / df$width
df$Roundness <- df$width/df$length     
df$Circularity <- df$Circ


df <- df %>%filter(size > 20,width > 2,length > 2)

#Cell size comparison
ggplot(df, aes(x = group, y = size, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  labs(title = "Cell Size Comparison",
       y = "Area") +
  theme_minimal()

# Cell width comparison
ggplot(df, aes(x = group, y = width, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  labs(title = "Cell Width Comparison",
       y = "Width") +
  theme_minimal()

#Cell length comparison
ggplot(df, aes(x = group, y = length, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  labs(title = "Cell Length Comparison",
       y = "Length") +
  theme_minimal()

#Roundness comaparision 
ggplot(df, aes(x = group, y = Roundness, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  labs(title = "Cell Shape (Roundness) Comparison",
       y = "Width / Length") +
  theme_minimal()

##Statistical test
df %>%
  wilcox_effsize(size ~ group)

wilcox.test(size ~ group, data=df)
wilcox.test(length ~ group, data=df)
wilcox.test(width ~ group, data=df)
wilcox.test(Roundness ~ group, data=df)	
# Final aesthetics, such as colors and line drawing, were done in Illustrator
