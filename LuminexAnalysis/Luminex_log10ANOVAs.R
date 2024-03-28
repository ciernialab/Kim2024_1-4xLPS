#make heatmaps for luminex data
#1/30/2020 Annie Ciernia

library(tidyverse)
library(xlsx)
library(ComplexHeatmap)
library(dendsort)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubr)
library(rstatix)
library(car)
library(broom)

################################################################################################
#read in luminex results
################################################################################################
#CB
CB <- read.xlsx("cytokine_chemokineData.xlsx", sheetIndex = 1)
CB$region <- c("cerebellum")
#HC
HC <- read.xlsx("cytokine_chemokineData.xlsx", sheetIndex = 2)
HC$region <- c("hippocampus")
#serum
serum <- read.xlsx("cytokine_chemokineData.xlsx", sheetIndex = 3)
serum$region <- c("serum")
#FC
FC <- read.xlsx("cytokine_chemokineData.xlsx", sheetIndex = 4)
FC$region <- c("frontal cortex")
#STR
STR <- read.xlsx("cytokine_chemokineData.xlsx", sheetIndex = 5)
STR$region <- c("striatum")

#################################################################################"###############
# Compare different regions for each gene
################################################################################################

#get common gene sets
common <- rbind(FC,HC,STR,CB,serum)
common$sample <- rep(c(1,2,3,4),(nrow(common)/4))
#reformat
common2 <- common %>% 
  group_by(region, cytokine.chemokine, sample) %>%
  gather(treatment,pg_per_ml,2:6)

#fix group names for treatment
common2$treatment <- gsub("X","",common2$treatment)
#fix order for treatment
common2$treatment <- factor(common2$treatment,levels = c("PBS"  , "1xLPS" ,"2xLPS" ,"3xLPS","4xLPS"))


#fix order for regions
common2$region <- factor(common2$region,levels = c("serum", "frontal cortex", "hippocampus"  ,  "striatum"   ,    "cerebellum"))

################################################################################
#statistics
################################################################################

#log-transformed the raw data 
library(nlme)
library(lsmeans)

logpgml <- as.data.frame(common2)
logpgml$logpgml <- log(logpgml$pg_per_ml)
logpgml <- logpgml[complete.cases(logpgml),]

df <- logpgml %>% select(-pg_per_ml) %>%
  group_by(region, sample, cytokine.chemokine) %>%
  spread(treatment, logpgml)

write.csv(df,"log.pg.per.ml.luminexdata.csv")


#test for outliers
outlier <- logpgml %>%
  group_by(cytokine.chemokine,region,treatment) %>%
  identify_outliers(logpgml) %>%
  filter(is.extreme == TRUE)

#The normality assumption can be checked by computing Shapiro-Wilk test for each outcome variable at each level of the grouping variable. If the data is normally distributed, the p-value should be greater than 0.05.
logpgml <- as.data.frame(logpgml)

logpgml$logpgml <- as.numeric(logpgml$logpgml)

mydata <-logpgml
mydata$inter<-with(mydata,interaction(cytokine.chemokine,region,treatment))
mydata1<-mydata[mydata$inter %in% names(which(table(mydata$inter) > 3)), ]  
library(plyr) 
ddply(mydata1, .(inter), summarize, n=length(logpgml))


#can only test for conditions with sufficient varience
out <- NULL
for (i in unique(mydata1$inter)) {
  tmp <- mydata1 %>% filter(inter ==i)
  
  if(length(unique(tmp$logpgml))>1) {
    
    SW <- tmp %>%
      group_by(inter) %>%
      shapiro_test(logpgml) 
    
    out <- rbind(out, SW)
  }
}

out$pass <- ifelse(out$p < 0.05, "significant", "ns")

table(out$pass)

#48 conditions v violated the SW test > KW test


################################################################################
#statistics: Kruskal-Wallis rank sum stat (no repeated measure)
#for use with non-parameteric data
################################################################################

#log-transformed the raw data 
logpgml <- as.data.frame(common2)
logpgml$logpgml <- log(logpgml$pg_per_ml)
logpgml <- logpgml[complete.cases(logpgml),]


#mixed model for brain region and treatment for each cytokine/chemokine
#loop through each gene and calculate kruskal wallis test
library(FSA)
kw_fullmodel2 <- NULL
dunn_fullmodel2 <- NULL

genes <- unique(as.character(logpgml$cytokine.chemokine))
brainregion <- unique(as.character(logpgml$region))

for(i in genes){
  
  print(i)
  #subset data for matching each sample/condition combo
  tmp <- logpgml %>%
    filter(cytokine.chemokine== i) 
  
  kw_fullmodel <- NULL
  dunn_fullmodel <- NULL
  for(n in brainregion){
    tmp2 <- tmp %>% filter(region == n)
    
    if(nrow(tmp2)>2) { #only run if there is data for this gene/region combo
      #full model kw:
      model <- kruskal.test(logpgml ~ treatment, data = tmp2) 
      
      kwresults <- data.frame(
        method = model$method,
        pvalue =  model$p.value,
        statistic =  model$statistic,
        df =  model$parameter,
        cytokine.chemokine = paste(i),
        region = paste(n)
      )

      #posthocs: dunn's test
      
      dunnTest <- dunnTest(logpgml ~ treatment,data=tmp2,method="none")
      dunnTest <- as.data.frame(print(dunnTest))
      
      #select comparisons of interest: PBS vs LPS
      outsum <- dunnTest[c(7:10),c(1:3)]
      
      #correct for multiple comparisons with dunn's test
      outsum$padjust <- p.adjust(outsum$P.unadj, method = "BH")
      
      outsum$gene <- paste(i)
      outsum$region <- paste(n)
      
      kw_fullmodel <- rbind(kw_fullmodel,kwresults)
      dunn_fullmodel <- rbind(dunn_fullmodel,outsum)
    }
  }
  
  kw_fullmodel2 <- rbind(kw_fullmodel2,kw_fullmodel)
  dunn_fullmodel2 <- rbind(dunn_fullmodel2,dunn_fullmodel)
}


kw_fullmodel2$Significant <- ifelse(kw_fullmodel2$pvalue < 0.05, "*", "ns")
dunn_fullmodel2$Significant <- ifelse(dunn_fullmodel2$padjust < 0.05, "*", "ns")


library(openxlsx)
write.xlsx(kw_fullmodel2,file="KWper_gene_fullmodel.xlsx")
write.xlsx(dunn_fullmodel2,file="Dunn_BHposthocs_per_gene_fullmodel.xlsx")

#################################################################################"###############
#boxplots graphs
#################################################################################"###############

pdf("cytokine_boxplot.pdf", height = 50, width =8.5)       

cbPalette <- c("lightblue", "yellow", "red","blue","purple")
ggplot(common2, aes(x=treatment, y=log(pg_per_ml)),group=treatment) + 
  facet_wrap(~cytokine.chemokine*region, scales = "free", ncol=5)+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=treatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=treatment),size=1) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Treatment") + #,limits = order2)+
  scale_y_continuous(name = "log(pg/ml)")

dev.off();

#filter for PBS vs LPS comparissons only
stat.test <- logpgml %>%
  group_by(region,cytokine.chemokine) %>%
  dunn_test(logpgml ~ treatment, p.adjust.method="none") %>%
  filter(group1=="PBS") %>%
  ungroup() %>%
  group_by(region,cytokine.chemokine) %>%
  adjust_pvalue(method = "BH") %>% #corrected only within each measure
  add_significance("p.adj") 


loc <- logpgml %>%
  group_by(region,cytokine.chemokine) %>%
  dunn_test(logpgml ~ treatment, p.adjust.method="none") %>%
  filter(group1=="PBS") %>%
  adjust_pvalue(method = "BH") %>% #corrected only within each measure
  add_significance("p.adj") %>%
  add_xy_position(x = "treatment")

stat.test$y.position <- loc$y.position


test <- dunn_fullmodel2 %>% 
  separate(Comparison, into=c("group1","group2"), sep=" - ") %>%
  mutate(.y. = c("logpgml")) %>%
  mutate(n1=4, n2=4) %>%
  select(gene,region,.y., group1, group2,n1,n2,Z,P.unadj,padjust,Significant)

colnames(test) <- c("cytokine.chemokine", "region", ".y.", "group1", "group2",    "n1",    "n2", "statistic","p", "p.adj", "p.adj.signif")

test$region <- factor(test$region,levels = c("serum", "frontal cortex", "hippocampus"  ,  "striatum"   ,    "cerebellum"))


test <- as.tibble(test)

test <- test %>% add_xy_position(x = "treatment")
  
#plot
cbPalette <- c("lightblue", "yellow", "red","blue","purple")

pdf("cytokine_boxplot_dunns.pdf", height = 50, width =8.5)       

bxp <- ggboxplot(logpgml, x = "treatment", y = "logpgml", fill = "#00AFBB",
                 facet.by = c("cytokine.chemokine", "region")
) 

bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns = TRUE, tip.length = 0)+
  scale_y_continuous(name = "log(pg/ml)", expand = expansion(mult = c(0.05, 0.10))) +
  geom_point(position = position_dodge(width = 0.90),aes(group=treatment),size=1) + 
  scale_fill_manual(values = cbPalette) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj.signif") +
 # theme_cowplot(font_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Treatment") 
  
dev.off()



pdf("cytokine_boxplot_dunns.pdf", height = 50, width =8.5)       

cbPalette <- c("lightblue", "yellow", "red","blue","purple")
ggplot(logpgml, aes(x=treatment, y=logpgml,group=treatment)) + 
  facet_wrap(~cytokine.chemokine*region, scales = "free", ncol=5)+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=treatment))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  geom_point(position = position_dodge(width = 0.90),aes(group=treatment),size=1) + 
    stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj.signif") +
  scale_fill_manual(values = cbPalette) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj.signif") +
  theme_cowplot(font_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Treatment") +
  scale_y_continuous(name = "log(pg/ml)", expand = expansion(mult = c(0.05, 0.10)))

dev.off();

bxp <- ggboxplot(df, x = "supp", y = "len", fill = "#00AFBB")
bxp + 
  stat_pvalue_manual(stat.test, label = "p") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


################################################################################
#complex heatmap data wrangling
################################################################################

#calculate the mean pg/ml across replicates for each region/treatment/cytokine:
Means <- common2 %>% group_by(region,cytokine.chemokine, treatment) %>%
  summarise(mean_pg_per_ml =  mean(pg_per_ml)) %>%
  group_by(region,cytokine.chemokine, treatment) %>%
  spread(treatment, mean_pg_per_ml) %>%
  arrange(region,cytokine.chemokine)

Means <- as.data.frame(Means)

#fold change = experiment/control
Means$FC1xLPSvsPBS <- Means$`1xLPS`/Means$PBS
Means$FC2xLPSvsPBS <- Means$`2xLPS`/Means$PBS
Means$FC3xLPSvsPBS <- Means$`3xLPS`/Means$PBS
Means$FC4xLPSvsPBS <- Means$`4xLPS`/Means$PBS

#put in log scale
Means[8:11] <- log(Means[8:11])
colnames(Means)[8:11] <- paste("log(",colnames(Means[8:11]),")",sep="")

#put adjusted pvalues in same format as input data for heatmaps
colnames(dunn_fullmodel2)[5] <- c("cytokine.chemokine")
pval <- dunn_fullmodel2 %>% dplyr::select(region, cytokine.chemokine, Comparison,padjust) 
pval <- as.data.frame(pval)
pval$Comparison <- as.character(pval$Comparison)
pval2 <- pval %>%
  spread(Comparison,padjust,3:4)
colnames(pval2)[3:6] <- paste("BHpvalue",colnames(pval2)[3:6],sep=".")


#merge in data:
MeansPval <- merge(Means, pval2 , by=c("region","cytokine.chemokine"), all=T)

write.xlsx(MeansPval,file="heatmapdata.xlsx")

#log FC data:
mat <- MeansPval[,grepl("FC",colnames(MeansPval))]
mat$region <- MeansPval$region
mat$cytokine.chemokine <- MeansPval$cytokine.chemokine

#fix col names 
colnames(mat) <- gsub("log\\(FC","",colnames(mat))
colnames(mat) <- gsub("\\)","",colnames(mat))


#pvalue data:
adjpval <- MeansPval[,grepl("BHpvalue",colnames(MeansPval))]
adjpval$region <- MeansPval$region
adjpval$cytokine.chemokine <- MeansPval$cytokine.chemokine


#make NA =1
adjpval[is.na(adjpval)] <- 1

#split by region for FC:
PFCdf <- mat %>% filter(region == "frontal cortex")
rownames(PFCdf) <- PFCdf$cytokine.chemokine
PFCdf$region <- NULL
PFCdf$cytokine.chemokine <- NULL
PFCdf <- as.matrix(PFCdf)

serumdf <- mat %>% filter(region == "serum")
rownames(serumdf) <- serumdf$cytokine.chemokine
serumdf$region <- NULL
serumdf$cytokine.chemokine <- NULL
serumdf <- as.matrix(serumdf)

HCdf <- mat %>% filter(region == "hippocampus")
rownames(HCdf) <- HCdf$cytokine.chemokine
HCdf$region <- NULL
HCdf$cytokine.chemokine <- NULL
HCdf <- as.matrix(HCdf)

STRdf <- mat %>% filter(region == "striatum")
rownames(STRdf ) <- STRdf$cytokine.chemokine
STRdf$region <- NULL
STRdf$cytokine.chemokine <- NULL
STRdf <- as.matrix(STRdf)

CBdf <- mat %>% filter(region == "cerebellum")
rownames(CBdf) <- CBdf$cytokine.chemokine
CBdf$region <- NULL
CBdf$cytokine.chemokine <- NULL
CBdf <- as.matrix(CBdf)


#split by region for adjpval:
PFCadjpval <- adjpval %>% filter(region == "frontal cortex")
rownames(PFCadjpval) <- PFCadjpval$cytokine.chemokine
PFCadjpval$region <- NULL
PFCadjpval$cytokine.chemokine <- NULL

serumadjpval <- adjpval %>% filter(region == "serum")
rownames(serumadjpval) <- serumadjpval$cytokine.chemokine
serumadjpval$region <- NULL
serumadjpval$cytokine.chemokine <- NULL


HCadjpval <- adjpval %>% filter(region == "hippocampus")
rownames(HCadjpval) <- HCadjpval$cytokine.chemokine
HCadjpval$region <- NULL
HCadjpval$cytokine.chemokine <- NULL


STRadjpval <- adjpval %>% filter(region == "striatum")
rownames(STRadjpval ) <- STRadjpval$cytokine.chemokine
STRadjpval$region <- NULL
STRadjpval$cytokine.chemokine <- NULL

CBadjpval <- adjpval %>% filter(region == "cerebellum")
rownames(CBadjpval) <- CBadjpval$cytokine.chemokine
CBadjpval$region <- NULL
CBadjpval$cytokine.chemokine <- NULL


################################################################################
#complex heatmap making
################################################################################
#set colors for heatmap
library(circlize)

max <- max(mat[,c(1:4)], na.rm=T)
min <- min(mat[,c(1:4)], na.rm=T)

col_fun = colorRamp2(c(-5,0, 10), c("steelblue","white",  "red"))


#scale by row: 
#dataRowNorm <- t(apply(mat, 1, function(x)  (x - mean(x)) / sd(x)))

#PFCdfNS <- apply(PFCdfNS, 2, as.numeric )
#distfunc <- function(x) daisy(x,metric="gower")
#d <- distfunc(PFCdfNS)
#dend = dendsort(hclust(dist(PFCdfNS)))

#PFC

ht1 = Heatmap(PFCdf, name = "log(Fold Change pg/ml)", column_title = "Frontal Cortex", 
              cluster_rows = FALSE,
              cluster_row_slices =FALSE,
              column_order = colnames(PFCdf),  border = TRUE,
              col = col_fun, 
              # row_split = siggois$PFCpattern,
              #row_split = siggois$PFCpattern,
              #row_split = 3,#number of clusters
              #top_annotation = ha,
              #add pvalues
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(PFCadjpval[i, j] <0.05)
                  #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
                  grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              }
              )

draw(ht1)

#HC
ht2 = Heatmap(HCdf, name = "log(Fold Change pg/ml)", column_title = "Hippocampus", 
              cluster_rows = FALSE,
              cluster_row_slices =FALSE,
              column_order = colnames(HCdf),  border = TRUE,col = col_fun, 
              #add pvalues
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(HCadjpval[i, j] <0.05)
                  #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
                  grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              })

draw(ht2)

#STR
ht3 = Heatmap(STRdf, name = "log(Fold Change pg/ml)", column_title = "Striatum", 
              cluster_rows = FALSE,
              cluster_row_slices =FALSE,
              column_order = colnames(STRdf),  border = TRUE,col = col_fun,
              #add pvalues
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(STRadjpval[i, j] <0.05)
                  #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
                  grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              })

draw(ht3)

#CB
ht4 = Heatmap(STRdf, name = "log(Fold Change pg/ml)", column_title = "Cerebellum", 
              cluster_rows = FALSE,
              cluster_row_slices =FALSE,
              column_order = colnames(CBdf),  border = TRUE,col = col_fun,
              #add pvalues
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(CBadjpval[i, j] <0.05)
                  #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
                  grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              })

draw(ht4)

#Serum
ht5 = Heatmap(serumdf, name = "log(Fold Change pg/ml)", column_title = "Serum", 
              cluster_rows = FALSE,
              cluster_row_slices =FALSE,
              column_order = colnames(serumdf),  border = TRUE,col = col_fun,
              #add pvalues
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(serumadjpval[i, j] <0.05)
                  #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
                  grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              })

draw(ht5)

ht_list = ht5+ ht1 + ht2 + ht3 +ht4

pdf("log10FCheat.luminex.pdf", width=8, height=6)

draw(ht_list, row_title = "Luminex Cytokines & Chemokines", row_title_gp = gpar(col = "black")
)

dev.off()







