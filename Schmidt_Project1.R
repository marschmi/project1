library(ggplot2)


phyla <- read.csv(file = "Phyla_Table.csv")
phyla <- subset(phyla,  select=c("sample", "trophicstate",  "Phylum", "abundPerPhylum", "filter"))

class <- read.csv(file = "Class_Table.csv")
class <- subset(class,  select=c("sample", "trophicstate",  "Phylum", "Class", "AbundPerClass", "filter"))

#Setting up our data frame
#To include the sub-phyla of the Proteobacteria merge the alpha, beta, delta, epsilon, gamma, 
#1.  Subset Proteo
proteo <- subset(class, Phylum == "Proteobacteria")
#2.  Subtract Proteos from phyla table
phyla <- subset(phyla, Phylum != "Proteobacteria")
#3.  Make class and phyla table the same -> a.  Remove Phylum column in class table
proteo$Phylum = NULL
#b.  Renmae Class column to be named Phylum to match phyla table
colnames(proteo) <- c("sample", "trophicstate", "Phylum", "SeqAbundance", "filter")
colnames(phyla) <- c("sample", "trophicstate", "Phylum", "SeqAbundance", "filter")
#4.  Rbind class and phyla tables.
data <- rbind(proteo, phyla)
#I normalized/rarefied at 8500 seqs, so to get relative abundance I am dividing by 8500
data$RelAbundance <- data$SeqAbundance/8500


#Summary stats for ratio analysis
##  Function below from this website:  http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
##  Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


summarystats <- summarySE(data, measurevar = "SeqAbundance", groupvars = c("filter", "Phylum"))
names(summarystats)[4]<-"Mean_SeqAbund"
free <- subset(summarystats, filter == "Free")
particle <- subset(summarystats, filter == "Particle")

#Calculate the log 2 ratio and create a new data frame
oddsratio <- log2(particle$Mean_SeqAbund/free$Mean_SeqAbund) #creates a vector
odds_ratio <- data.frame(free$Phylum, oddsratio) #combine the Phylum names with it
odds_ratio <- rename(odds_ratio, c("free.Phylum" = "Phylum"))
#  http://learnr.wordpress.com/2009/06/01/ggplot2-positioning-of-barplot-category-labels/
odds_ratio$colour <- ifelse(odds_ratio$oddsratio < 0, "firebrick1","steelblue")
odds_ratio$hjust <- ifelse(odds_ratio$oddsratio > 0, 1.3, -0.3)
#odds_ratio$Phylum <- factor(odds_ratio$Phylum, levels = odds_ratio$Phylum[order(odds_ratio$oddsratio)])

#Let's delete Phyla from our data frame that are Inf, -Inf, or NaN
odds_ratio <- odds_ratio[-31, ] #Remove Dictoglomi
odds_ratio <- odds_ratio[-26, ]# Candidate_division_TM7
odds_ratio <- odds_ratio[-24, ]# Candidate_division_OP8
odds_ratio <- odds_ratio[-22, ]# Candidate_division_OP11
odds_ratio <- odds_ratio[-15, ]# Thermotogae

#Order by oddsratio!
#odds_ratio <- odds_ratio[with(odds_ratio, order(-oddsratio)),] #doesn't work right!!!!
odds_ratio$Phylum <- factor(odds_ratio$Phylum, levels = odds_ratio$Phylum[order(odds_ratio$oddsratio)])


#  http://learnr.wordpress.com/2009/06/01/ggplot2-positioning-of-barplot-category-labels/
ratio_plot <- ggplot(odds_ratio, aes(y=oddsratio, x = Phylum, label = Phylum, hjust=hjust, fill = colour)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  #geom_bar(stat = "identity", aes(fill = colour)) +
  geom_text(aes(y=0, fill = colour)) + theme_bw() + 
  coord_flip() + ggtitle("Particle-Association vs. Free-Living") +
  labs(y = "Log2 Abundance", x = "") + scale_x_discrete(breaks = NULL) +
  scale_fill_manual(name  ="", breaks=c("firebrick1", "steelblue"), 
                    labels=c("More Abundant in Free-Water", "More Abundant in Particles"),
                    values = c("royalblue","violetred")) +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=0, colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        legend.justification=c(1,0), legend.position=c(1,0.6)); ratio_plot  




#Use a non-parametric test pairwise test = Pairwise Wilcoxon Test.
#Null Hypothesis:  The sequence abundance of particle and free living are the same (= to 0).
#Alternative Hypothesis:  The seq abundance of particle and free are different.
# If P is less than 0.05 then we reject the null hypothesis.
#wilcox.test(particle, free)

#To run a wilcox.test we can't just compare means -> we have to compare the distributions of our data
#To get the distributions:
dist <- subset(data, select = c("Phylum", "SeqAbundance", "filter", "sample")) #Get distributions from original file

#Let's make sure that all of the sample names are the same for particles as they are for frees
#This way we can check for what samples are missing, I first have to change the names of the samples to be the same
dist$sample2 <- gsub("3um", "", dist$sample) #create new column with substitute nothing for "3um"
dist$sample = NULL   #delete old column
names(dist)[4]<-"sample" #give original name again
as.factor(dist$sample)
#Have the same sample name now

#Get the distributions of the data for each factor.
free_dist <- subset(dist, filter == "Free") #Make a free-living dist
part_dist <- subset(dist, filter == "Particle") #Make a particle-associated dist
dim(free_dist); dim(part_dist) #They have different dims - time to fix this!

#What samples are missing from each data set
sux <- setdiff(free_dist$sample, part_dist$sample); sux
this_sux <- setdiff(part_dist$sample, free_dist$sample); this_sux

#Remove the above non-duplicate and duplicate samples that we found with sux and this_sux
dim(free_dist); dim(part_dist) 
free_dist <- subset(free_dist, sample != "BASE2" & sample != "LONE1" & sample != "LONE2" & sample != "LONH2")
part_dist<- subset(part_dist, sample != "BSTE2" & sample !="GULH1" & sample !="LEEE1" & sample !="SIXE1" & sample !="SIXH2")
dim(free_dist); dim(part_dist) # I'm so happy!!!

#We should also remove the phyla that gave us Inf, -Inf, or NaN from our log2 fold ratio
free_dist <- subset(free_dist, Phylum != "Thermotogae" & Phylum != "Candidate_division_OP11" & Phylum != "Candidate_division_TM7" & Phylum != "Candidate_division_OP8" & Phylum != "Dictyoglomi")
part_dist <- subset(part_dist, Phylum != "Thermotogae" & Phylum != "Candidate_division_OP11" & Phylum != "Candidate_division_TM7" & Phylum != "Candidate_division_OP8" & Phylum != "Dictyoglomi")
dim(free_dist); dim(part_dist) # Good

##Make a vector of all the names of the Phyla in our samples
phylumlist_free<-as.character(unique(free_dist$Phylum))
phylumlist_part<-unique(part_dist$Phylum)
setdiff(phylumlist_free, phylumlist_part)
phylumlist <- phylumlist_free

#Make data frame filled with NAs to hold results inside of
wresults<-data.frame(matrix(NA,length(phylumlist),2))
names(wresults) <- c("Phylum", "P_Value")

for(i in 1:length(phylumlist)){
  free_sub<- subset(free_dist, Phylum == phylumlist[i])
  part_sub<- subset(part_dist, Phylum == phylumlist[i])
  wilcox <- wilcox.test(part_sub$SeqAbundance, free_sub$SeqAbundance, paired = TRUE)
  wresults[i,1]<-phylumlist[i]
  wresults[i,2]<-wilcox$p.value
}
# warnings()

#We have to bonferroni correct because we ran 37 tests --> some are significant when they shouldn't be!
num_tests <- length(unique(wresults$Phylum))
wresults$Bonferroni_P<-wresults$P_Value/num_tests
P_cutoff <- 0.05/37
length(wresults$Bonferroni_P[wresults$Bonferroni_P<P_cutoff])
sigresults<-subset(wresults,Bonferroni_P<P_cutoff)



































#############################  CYANOBACTERIA  ######################
#let's try cyanobacteria
free_cyanos <- subset(free_dist, filter == "Free" & Phylum == "Cyanobacteria")
part_cyanos <- subset(part_dist, filter == "Particle" & Phylum == "Cyanobacteria")
dim(free_cyanos); dim(part_cyanos) 
length(unique(free_cyanos$sample)); length(unique(part_cyanos$sample))

#Now we can do the test: thank god!
wilcox_cyano <- wilcox.test(part_cyano$SeqAbundance, free_cyano$SeqAbundance, paired = TRUE)
#And we reject our Null Hypothesis!  PARTY!


#############################  VERRUCOMICROBIA  ######################
########
free_verrucos <- subset(free_dist, filter == "Free" & Phylum == "Verrucomicrobia")
part_verrucos  <- subset(part_dist, filter == "Particle" & Phylum == "Verrucomicrobia")
<- subset(verrucos, filter == "Particle")
dim(free_verrucos); dim(part_verrucos) 
length(unique(free_verrucos$sample)); length(unique(part_verrucos$sample))

#Now we can do the test: thank god!
wilcox_verrucos <- wilcox.test(part_verrucos$SeqAbundance, free_verrucos$SeqAbundance, paired = TRUE); wilcox_verrucos
#Barely non-significant, so close that I need a beer


#############################  PLANCTOMYCETES  ######################
########
free_planctos <- subset(free_dist, filter == "Free" & Phylum == "Planctomycetes")
part_planctos  <- subset(part_dist, filter == "Particle" & Phylum == "Planctomycetes")
dim(free_planctos); dim(part_planctos) 
length(unique(free_planctos$sample)); length(unique(part_planctos$sample))

wilcox_planctos <- wilcox.test(part_planctos$SeqAbundance, free_planctos$SeqAbundance, paired = TRUE); wilcox_planctos
wilcox_planctos2 <- wilcox.test(part_planctos$SeqAbundance, free_planctos$SeqAbundance); wilcox_planctos2


#Reject the null - woooo


#############################  BACTERIODETES  ######################
########
free_bacteriods <- subset(free_dist, filter == "Free" & Phylum == "Bacteroidetes")
part_bacteriods  <- subset(part_dist, filter == "Particle" & Phylum == "Bacteroidetes")
dim(free_bacteriods); dim(part_bacteriods) 
length(unique(free_bacteriods$sample)); length(unique(part_bacteriods$sample))

wilcox_bacteriods <- wilcox.test(part_bacteriods$SeqAbundance, free_bacteriods$SeqAbundance, paired = TRUE); wilcox_bacteriods
#Accept null boooooooooo

#############################  BETAPROTEOBACTERIA  ######################
########
free_betaproteos <- subset(free_dist, filter == "Free" & Phylum == "Betaproteobacteria")
part_betaproteos  <- subset(part_dist, filter == "Particle" & Phylum == "Betaproteobacteria")
dim(free_betaproteos); dim(part_betaproteos) 
length(unique(free_betaproteos$sample)); length(unique(part_betaproteos$sample))

wilcox_betaproteos <- wilcox.test(part_betaproteos$SeqAbundance, free_betaproteos$SeqAbundance, paired = TRUE); wilcox_betaproteos
wilcox_betaproteos2 <- wilcox.test(part_betaproteos$SeqAbundance, free_betaproteos$SeqAbundance); wilcox_betaproteos2
#Accept null boooooooooo








#Remove the Phyla that cause infs, -Inf, and Nans
#means <- means[-31, ] #Remove Dictoglomi
#means <- means[-26, ]# Candidate_division_TM7
#means <- means[-24, ]# Candidate_division_OP8
#means <- means[-22, ]# Candidate_division_OP11
#means <- means[-15, ]# Thermotogae



ddply(means, function(x) with(x,wilcox.test(Particle_MeanAbund,Free_MeanAbund))$p.value)






group.pairs <- combn(unique(df$group),2,simplify=FALSE)
# this loops over the 2nd margin - the columns - of df and makes each column
# available as x
apply(df[-1], 2, function(x)
  # this loops over the list of group pairs and makes each such pair
  # available as an integer vector y
  lapply(group.pairs, function(y)
    wilcox.test(x[df$group %in% y],df$group[df$group %in% y])))
















### This stuff is junk that I'm afraid to let go

apply(df[1],2,function(x) wilcox.test(x,df$group))

group.pairs <- combn(unique(means$Particle_MeanAbund),2,simplify=FALSE)

apply(means[-1], 2, function(x)
  # this loops over the list of group pairs and makes each such pair
  # available as an integer vector y
  lapply(group.pairs, function(y)
    wilcox.test(x[df$?wilcoxgroup %in% y],df$group[df$group %in% y])))

apply

wilcox.test(41.575, 21.38462)

wilcox.test(means$Particle_MeanAbund, means$Free_MeanAbund)



apply(df[-1],2,function(x) kruskal.test(x,df$group))

#### Time to test for significance using pairwise wilcoxon tests between the means of the particle and free living for each phyla
free_test <- subset(free, select = c("filter", "Phylum", "Mean_SeqAbund")); 
particle_test <- subset(particle, select = c("filter", "Phylum", "Mean_SeqAbund")); 
#Change the name of the Mean_SeqAbund column to include free or particle
names(free_test)[3]<-"Free_MeanAbund"
names(particle_test)[3]<-"Particle_MeanAbund"
#Check to make sure the data frames are in the same order --> yes!
dim(particle_test); dim(free_test); particle_test$Phylum; free_test$Phylum
#Make a new data frame with 3 columns: 1. Phyla name, 2. Particle Abund, 3. Free Abundance
means <- particle_test
means$filter = NULL
means$Free_MeanAbund <- free_test$Free_MeanAbund

#Remove the same infs, -Inf, and Nans
means <- means[-31, ] #Remove Dictoglomi
means <- means[-26, ]# Candidate_division_TM7
means <- means[-24, ]# Candidate_division_OP8
means <- means[-22, ]# Candidate_division_OP11
means <- means[-15, ]# Thermotogae




dist <- subset(data, Phylum == "Acidobacteria" & Phylum == "Actinobacteria" &
                 Phylum == "Alphaproteobacteria" & Phylum == "Betaproteobacteria" &
                 Phylum == "Chlamydiae" & Phylum == "Cyanobacteria" &
                 Phylum == "Deferribacteres" & Phylum == "Deltaproteobacteria" &
                 Phylum == "Elusimicrobia" & Phylum == "Epsilonproteobacteria" &
                 Phylum == "Gammaproteobacteria" & Phylum == "Gemmatimonadetes" &
                 Phylum == "SPOTSOCT00m83" & Phylum == "TA18" &
                 Phylum == "unclassified" & Phylum == "Armatimonadetes" &
                 Phylum == "Bacteroidetes" & Phylum == "BD1" &
                 Phylum == "Candidate_division_BRC1" & Phylum == "Candidate_division_OD1" &
                 Phylum == "Candidate_division_OP3" & Phylum == "Candidate_division_SR1" &
                 Phylum == "Candidate_division_WS3" & Phylum == "Chlorobi" &
                 Phylum == "Chloroflexi" & Phylum == "Deinococcus-Thermus" &
                 Phylum == "Fibrobacteres" & Phylum == "Firmicutes" &  
                 Phylum == "Fusobacteria" & Phylum == "Lentisphaerae" &
                 Phylum == "NPL-UPA2" & Phylum == "Planctomycetes" &
                 Phylum == "Spirochaetae" & Phylum == "Tenericutes" &
                 Phylum == "TM6" & Phylum == "Verrucomicrobia" & Phylum == "WCHB1-60")
