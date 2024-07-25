############################### READ ME #########################################

# MS Proteomics Pipeline - DIA Unlabeled - Written for DIA-NN report.tsv files ONLY
  #DIA-NN:  Demichev, V., et al. (2020). https://doi.org/10.1038/s41592-019-0638-x

#Written by Sarah Garcia, Lea Barny, Jake Hermanson - 2024

#Set Working Directory - unique to each system
#setwd("path/to/your/report.tsv/file")

setwd("C:/Users/nyanc/OneDrive/Documents/All Files - personal/Work_Vandy/R_misc/Generic_MS_pipeline/Version1.0_Generic")

#Goal: Analyze Mass Spectrometry Data (Un-labeled DIA)
#       and perform statistical testing

#LOAD ALL FUNCTIONS AND LIBRARIES THEN START AT STEP ONE

############################ Load Tools ########################################

#BEFORE BEGINNING: Load ALL listed libraries and functions below:
#                    (Collapse functions for easier viewing)

#Load libraries
library(tidyverse)        #used in almost every function
library(data.table)       #used in almost every function
library(diann)            #DIA only (peptide to protein roll-up)
library(broom)            #Used for t.tests
library(ggfortify)        #Used for the PCA plot
library(writexl)
library(readxl)

### (OPTIONAL LIBRARIES) ### - Load for each optional step if desired
library(cleaver)          #used for sequence coverage
library(EnhancedVolcano)  #used for EnhancedVolcano plot
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')


library(VennDiagram)      #used for making VennDiagrams
library(NbClust)          ####
library(clusterCrit)      # ^
library(factoextra)       #Used for k-means clustering and heatmap
library(cluster)          # V
library(fpc)              # V
library(pheatmap)         ####

#Load functions
load_file <- function(file){
  as.data.frame(fread(file, stringsAsFactors = FALSE))
}                         #works with .csv and .tsv; not .xlsx
diann.dia.qc <- function(raw.data, cont.name){
  
  cat("Cleaning up raw data...\n")                  #print message to console
  
  before <- length(unique(raw.data$Protein.Group))  #store size for comparison
  
  #Data Clean-up - Parameters set at: Q.value <= 0.01; 1+ Peptide per Protein
  raw.data.filtered <- raw.data %>%                                       #Takes our imported raw data
    group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
    filter(!str_detect(Protein.Group, "^Cont_")) %>%                      #filter out groups that match a given contaminant
    filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
    group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
    filter(length(unique(Stripped.Sequence))>1)                           #Filters for lengths > 1 of pasted sequence 

  
  after <- length(unique(raw.data.filtered$Protein.Group))                #store size for comparison
  
  removed <- (1-(after/before))*100                                       #get size comparison
  
  cat("\nRemoved", removed, "% of the data due to q.value > 0.01 \n and/or only 1 peptide per protein observed \n and/or determined to be a contaminant.\n")
  
  export.file <- raw.data.filtered %>%                                    #take our qc filtered data                   
    dplyr::select(Run, Protein.Group, Genes, Stripped.Sequence,           #select relevant columns
                  Precursor.Id, Precursor.Normalised) %>%   
    mutate(Precursor.Normalised = na_if(Precursor.Normalised, 0)) %>%     #replace any zeros with NA
    tidyr::unite(Protein.Info, Genes, Protein.Group, sep = "*")           #combine two columns, separated by '*'
  
}       #import diann file and contaminants file; apply quality control filters
assign.group <- function(df, groups = treatments){
  
  sample.names <- unique(df$Run)                                     #store unique sample names
  
  temp <- sample.names %>% data.frame(Run = .) %>%                   #make names a data.frame
    arrange(Run) %>%                                                 #arrange names alphabetically
    mutate(group = NA, sample = NA) %>% data.frame()                 #make two empty columns
  
  cat("Groups:", paste0(groups[1:length(groups)]))                   #print message to console
  
  for(i in 1:nrow(temp)){                                            #for each row in 'temp'
    user_input <- readline(prompt = paste0(                          #prompt the user to give
      "Assign ", temp[i,"Run"], " to group [1 - ",                   ## each sample a group number
      length(groups), "] : "))
    temp[i,"group"] <- groups[as.numeric(user_input)]                #store the sample's group assignment
  }
  
  name.fix <- function(df){                                          #function definition
    df <- df %>% data.frame()                                        #make df a data.frame
    max_parts <- max(str_count(df[ ,1], "_"))                        #count max pieces after separation
    separated_df <- df %>%                                           #take data.frame
      separate_wider_delim(col = 1,                                  #separate wider
                           names = paste0("X_", 1:(max_parts + 1)),  ### make names dynamic for varying sizes
                           delim = "_", too_few = "align_start") %>% 
      dplyr::select(where(~length(unique(.)) > 1)) %>%               #select columns with at least 2 different entries
      unite(Run, starts_with("X_"), sep = "_", na.rm = TRUE) %>%     #combine Run column with columns starting with 'X_'
      arrange(Run)                                                   #arrange by the 'Run' column alphabetically
  } 
  
  temp[ ,"sample"] <- name.fix(temp$Run)                             #run previously define function on sample names
  
  sample.info <- temp %>% mutate(                                    #take our almost finished data.frame
    sample = paste0(temp$sample, '_', group)                         #combine sample name with associated treatment group
  )
  
}   #Assign treatment groups
medNorm <- function(df, samples = Run, 
                    values = Precursor.Normalised) { 
  df %>% data.frame() %>%                                #take your specified data.frame
    mutate(global_median = median({{ values }},          #calculate the global median
                                  na.rm = TRUE)) %>% 
    group_by({{ samples }}) %>%                          #group by sample name
    mutate(x_median = median({{ values }},               #calculate each sample's median
                             na.rm=TRUE)) %>% 
    mutate(norm_factor = global_median / x_median) %>%   #calculate the norm_factor
    mutate(Norm_abundance = norm_factor * {{ values }}   #multiply values by the norm_factor
    ) %>%
    select(-global_median, -x_median, -norm_factor)      #remove irrelevant columns
} #median normalize data
max_LFQ <- function(df, values = "Norm_abundance") {
  diann_maxlfq(df, sample.header = "Run", group.header = "Protein.Info",
               id.header = "Precursor.Id", quantity.header = values,
  )
} #simpler wrapper for diann_maxLFQ
filter.NA <- function(df, points = 3){
  
  test <- list()                                         #initialize empty list
  
  df <- df %>% rownames_to_column(var = "Protein.Info")
  
  for(i in treatments){                                  #for each treatment specified
    sample.number <- df %>%                              #take our wide data.frame
      dplyr::select(Protein.Info, ends_with(i)) %>%      #select treatment columns
      ncol(.) - 1                                       #store number of treatment columns
    NA.amount <- sample.number - points                  #number of treatment columns minus minimum number of points
    test[[i]] <- df %>%                                  #take our wide data.frame
      dplyr::select(Protein.Info, ends_with(i)) %>%      #select treatment columns  
      filter(rowSums(is.na(.)) <= NA.amount)             #filter for columns that meet the criteria
  }
  
  combined_df <- Reduce(function(x, y) full_join(x, y, by = c("Protein.Info")), test)
  sample.amount <- ncol(combined_df) - 1
  combined_df_filter <- combined_df %>% distinct() %>% filter(rowSums(is.na(.)) < sample.amount)
}               #filter for minimum # of points per observation
info.group <- function(df){
  group.size <- data.frame(df) %>%                       #get group size (rows)
    group_by(Run) %>% group_size()
  group.name <- data.frame(df) %>%                       #get group labels
    group_by(Run) %>% group_keys()
  group.info <- bind_cols(group.name, group.size)        #create group info matrix
  names(group.info) <- c("Run", "Count")                 #rename column names
  group.info <- group.info                               #export file to global environment 
}                          #function to obtain sample sizes
average.FC <- function(data){
  test <- list()
  for(i in treatments){
    test[[i]] <- data %>% mutate(AverageFC = rowMeans(dplyr::select(., ends_with(i)), na.rm = TRUE)) %>% 
      dplyr::select(AverageFC) %>% rownames_to_column() %>% setNames(c("Protein.Info", paste0("Avg.", i)))
  }
  Avg.FC.log2 <- test %>% purrr::reduce(full_join, by = "Protein.Info")
}                        #function to calculate treatment avg.Log2FC

### (OPTIONAL FUNCTIONS) ### - Load for each optional step if desired
digest.one.species <- function(species.list, dir) {
  
  chomp.protein<-function(species.list,dir){  
    
    proteome <- read_delim(dir, delim = "\n", col_names = FALSE, col_types = "c") #Will get one file 
    proteome.list <- proteome %>% na.omit() %>% data.frame()
    
    i=0                                             #Initialize empty objects for faster processing time
    flag=0
    y<-data.frame()
    protein.size<-list()
    sequence.all<-list()
    alist<-list()
    species.name<-species.list                      #this will be user supplied, must match how it looks in FASTA
    b<-''
    c<-''
    
    cat(paste("Now importing:",dir, "\n"))        #basic progress checker
    for(i in proteome.list[ ,1]){                 #for each line of inputted proteomeFASTA
      if(str_detect(i,"^>")){                     #if row starts with >
        if(flag==1){                              #flag check to concatenate sequence lines
          sequence<-str_c(alist,collapse = "")    #collapse
          aasize<-nchar(sequence)                 #get amino acid size
          protein.size<-rbind.data.frame(protein.size,aasize)  #store size
          sequence.all<-rbind.data.frame(sequence.all,sequence)  #store sequence
          alist<-list()                           #empty temp variable
          flag=0                                  #reset flag
        }
        flag=0                                    #reset
        pattern.species<-paste0("OS=",species.list)#make species name pattern
        pattern.gn<-paste0("GN=\\S+")             #make gene name pattery
        pattern.type<-paste0("^>\\S+")            #make >tr or >sp pattern 
        a<-str_extract(i,pattern = pattern.species) #find and store species name (user supplies species name)
        b<-str_extract(i,pattern = pattern.gn)    #find and store gene name
        c<-str_extract(i,pattern = pattern.type)  #find and store pattern type (>tr or >sp, should be only two options)
        abc<-c(a,b,c)                             #make above identifiers into easy to bind object
        y<-rbind(y,abc)                           #rbind to fill rows of new df
      }
      else{
        flag=1                                    #Raise the flag!
        b<-str_trim(i,side = c("right"))          #Cut off new line character
        alist[[i]]<-b                             #Storing all sequence lines of current protein
      }
    }
    if(flag==1){                                  #this catches the last protein's sequence VVV
      sequence<-str_c(alist,collapse = "")        #
      aasize<-nchar(sequence)                     #
      protein.size<-rbind.data.frame(protein.size,aasize) #
      sequence.all<-rbind.data.frame(sequence.all,sequence) #
      alist<-list()                               #
      flag=0                                      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    }
    y[,4]<-protein.size                             #add sizes to df
    y[,5]<-sequence.all                             #add sequences to df
    colnames(y)<-c("Species","Gene.Name","UniProt Type","AA Length","Sequence") #add colnames to object 'y'
    y$Species <- str_replace_all(y$Species, "OS=", "")  #remove OS= from species column
    y$`Gene.Name` <- str_replace_all(y$`Gene.Name`, "GN=", "")  #remove GN= from gene name column
    proteomes<<-y %>% separate(`UniProt Type`, into = c("discard", "UniProt.AC","discard2"), sep = "\\|") %>%
      dplyr::select(-discard,-discard2)              #Drop the 'discard' columns
  }   #Load 'chomp' function to import proteome files
  chomp.protein(species.list,dir)                  #Generate 'Proteomes" object (both species combined)
  
  species.A <<- species.list[[1]]                  #Store species names
  
  proteomes.A <- proteomes %>%                     #Take our Proteomes object
    filter(Species == species.A)
  
  #If missing a gene name, fill in name with UniProt Accession number
  proteomes.A$Gene.Name <- ifelse(is.na(proteomes.A$Gene.Name), 
                                  proteomes.A$UniProt.AC, proteomes.A$Gene.Name)
  
  cat("Now cleaving peptides...\n")
  
  #Create filtered peptide list
  tryp.digest.A <- proteomes.A$Sequence %>%       #Looks at protein sequences one at a time
    cleave(enzym = "trypsin", missedCleavages = 0:1, unique = TRUE) %>% #Specify trypsin enzyme, 0-1 missed cleavages, and only unique peptides (per protein)
    lapply(data.frame) %>%                        #Makes all list entries into data frames
    setNames(str_c(proteomes.A$Gene.Name, proteomes.A$UniProt.AC, proteomes.A$Species, sep = '*')) %>% #Rename data frames with the associated protein name
    lapply(bind_rows) %>%                         #Collapses the peptides into single protein entries (all peptides are grouped under their associated protein)
    lapply(function(df) {                         #Make a small function to apply across the list
      filtered_df <- filter(df, str_length(X..i..) > 6 & str_length(X..i..) < 31) #filter for peptides with charLength 6 < x < 31
      colnames(filtered_df) <- "Peptide"          #Rename the column names of each protein data frame
      return(filtered_df)                         #Gives the edited filtered_df back to the parent pipe sequence
    }) %>% discard(~ nrow(.) == 0)                #discard proteins with zero peptides (after filtering for length)
  
  cat("Performing species identification...\n")
  
  #Assemble Species.A Peptides
  species.A.peptides <- tryp.digest.A %>%         #Takes our tryptic digestion
    bind_rows(.id = "source") %>%                 #Collapses protein data frames into a single data frame
    dplyr::select(source, Peptide) %>%            #Selects for specified columns in specified order
    separate_wider_delim(cols = source, delim = '*', names = c("Gene.Name", "UniProt.Id", "Species")) #Separates 'source' column into three columns (descriptors)
  
  #Group by similar peptides, filter for groups with only one entry (proteotypic)
  proteotypic.A <- species.A.peptides %>% group_by(Peptide) %>% filter(n() == 1) %>% mutate(P=1)
  
  #Group by similar peptides, filter for groups with one+ entry (razor)
  multi.A <- species.A.peptides %>% group_by(Peptide) %>% filter(n() > 1) %>% mutate(P=0)
  
  #We now have a master list of all peptides with Proteotypic Status
  all.species.all.peptides<<-bind_rows(proteotypic.A, multi.A)
  
  cat("Done!\n")
  
}  #function to generate library for sequence coverage
protein.cov <- function(df){
  
  cat("Gathering peptide information...\n")
  
  #Subset peptides
  peptides.stripped <- df %>% data.frame() %>%
    dplyr::select(Run, Protein.Info, Stripped.Sequence)
  
  #get protein sequence information
  proteomes.x <- proteomes %>% unite(Protein.Info, Gene.Name, UniProt.AC, sep = "*") %>%
    dplyr::select(Protein.Info, `AA Length`, Sequence)
  
  #get peptide positions
  peptide.positions <- peptides.stripped %>% data.frame() %>% 
    inner_join(proteomes.x, by = "Protein.Info") %>%
    distinct() %>% group_by(Run) %>%
    mutate(Pep.aa = nchar(Stripped.Sequence),
           Start = map2(Stripped.Sequence, Sequence, ~ gregexpr(.x, .y)[[1]])) %>%
    unnest(Start) %>% mutate(End = (Start + Pep.aa) - 1) %>%
    filter(Start > 0 & End > 0)
  
  cat("Mapping peptides...\n")
  
  coverage_data <- peptide.positions %>%
    group_by(Run, Protein.Info, Sequence, `AA Length`) %>%
    summarise(
      coverage_map = list(
        Reduce(`|`, map2(Start, End, ~ {
          coverage <- rep(FALSE, unique(`AA Length`))
          coverage[.x:.y] <- TRUE
          coverage
        }))
      ),
      .groups = 'drop'
    ) %>% mutate(protSum = map_int(coverage_map, sum), 
                 perCov = (protSum / `AA Length`) * 100) %>%
    select(Run, Protein.Info, protSum, perCov)
  
  protein.coverage <<- peptide.positions %>%
    left_join(coverage_data, by = c("Run", "Protein.Info")) %>%
    dplyr::select(Run, Protein.Info, perCov) %>% distinct()
  
  cat("Done! 'protein.coverage' created.\n")
  
}                         #functino to calculate protein sequence coverage
volcano.curve <- function(result.volcano,                #function to plot curvature volcano
                          c.95 = 0.4,             #Define curvature
                          c.99 = 0.8,             #Define high-confidence curvature
                          x0 = 0.5,               #Define minimum fold change
                          adjpCutoff = 0.01){
  
  #adding curve cutoff
  curve_cutoff <- function(lfc, c, x0, sigma) {
    return(c * sigma / (lfc - x0))
  }
  
  mirrored_function <- function(x, c, x0, sigma) {
    ifelse(abs(x) > x0, curve_cutoff(abs(x), c, x0, sigma), NA)
  }
  
  # Calculate standard deviation of log2 fold changes
  sigma <- sd(result.volcano$estimate, na.rm = TRUE)
  
  # Calculate curve cutoff values
  result.volcano$curve_cutoff_95 <- curve_cutoff(abs(result.volcano$estimate), c.95, x0, sigma)
  result.volcano$curve_cutoff_99 <- curve_cutoff(abs(result.volcano$estimate), c.99, x0, sigma)
  result.volcano$significant_med <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_95 &
    abs(result.volcano$estimate) >= x0
  result.volcano$significant_high <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_99 &
    abs(result.volcano$estimate) >= x0
  
  #Calculate graph sizes
  min.x <- as.integer(min(result.volcano$estimate, na.rm = TRUE)) - 1
  max.x <- as.integer(max(result.volcano$estimate, na.rm = TRUE)) + 1
  
  min.y <- 0
  max.y <- as.integer(max(-log10(result.volcano$p.adj), na.rm = TRUE)) + 1
  
  # ggplot volcano plot with specific labels
  volcano_plot <<- ggplot(result.volcano, aes(x = estimate, y = -log10(p.adj)#, label = rownames(result.volcano)
  )) +
    geom_point(alpha = 0.2, colour = 'grey30') +
    stat_function(fun = function(x) mirrored_function(x, c = c.95, x0 = x0, sigma = sigma), color = "black") +
    stat_function(fun = function(x) mirrored_function(x, c = c.99, x0 = x0, sigma = sigma), color = "blue") +
    theme_bw() +
    labs(title = "Volcano Plot with Curve Cutoff Based on Standard Deviation",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted p-value") +
    scale_x_continuous(breaks = seq(min.x, max.x, by = 1), limits = c(min.x, max.x))+
    scale_y_continuous(breaks = seq(min.y, max.y, by = 1), limits = c(min.y, max.y)) +
    geom_point(data = subset(result.volcano, significant_med), aes(color = "med-conf"), size = 2) +
    geom_point(data = subset(result.volcano, significant_high), aes(color = "high-conf"), size = 2)
  
  return(volcano_plot)
} #Define adj.p-value cutoff

############### Step One: Input and Reformat Data ##############################

#Load data files - .csv or .tsv (doesn't like .xlsx)
#raw.data <- load_file("your.diann.report.tsv")
raw.data <- load_file("report.tsv")

#Clean up data files
#raw.data.qc <- diann.dia.qc(raw.data, cont.name = "your.contaminants.file.fasta")
raw.data.qc <- diann.dia.qc(raw.data, cont.name = "contaminants.fasta") #Will work regardless of contaminant presence or absence
#This works regardless of if a contaminant file was used to process the data 

################################EASY METHOD For Data Import without Searching without Contaminent File: 
import.diann.org <- function(file = "report.tsv"){cat("Importing report.tsv file...\n")
  
  #Import .tsv file
  raw.data <<- data.frame(diann::diann_load(file)) #Import the output .tsv file from DIA-NN:
  
  #Add list of contaminants to be removed from data before protein quantification: 
  contaminants <- c("P00761","Q32MB2","P19013","Q7RTT2","P15636","P09870","Q9R4J5",
                    "P0C1U8","P00766","P13717","Q9U6Y5","P21578","O76009","O76011",
                    "O76013","O76014","O76015","P08779","Q14525","Q14532","Q15323",
                    "Q92764","Q14533","Q9NSB4","P78385","Q9NSB2","P78386","O43790",
                    "Q6IFU5","Q9UE12","Q8IUT8","Q6NT21","Q6ISB0","Q6NTB9","Q6IFU6",
                    "P04264","P13647","P35908","P13645","P35527","A3EZ79","P02533",
                    "P02538","P48668","P04259","A3EZ82","Q2KIG3","Q0VCM5","Q3SZ57",
                    "Q9N2I2","Q3SZH5","P28800","Q1A7A4","P41361","Q2YDI2","Q3Y5Z3",
                    "P81644","Q2KJ83","Q2KIT0","A2I7N3","Q3SZV7","Q2KJC7","Q3SZR3",
                    "Q28107","P02672","Q1RMN8","Q58D62","P06868","Q2KJF1","P02584",
                    "P02777","Q3SX14","P17697","Q6T181","P34955","P21752","Q32PJ2",
                    "Q28194","P00978","Q5XQN5","Q32PI4","Q9TTE1","Q2KIU3","P67983",
                    "Q28065","Q862S4","Q2KIF2","Q3SX28","Q0V8M9","Q148H6","Q29RQ1",
                    "Q95M17","P07224","Q2HJF0","Q2KIH2","Q04695","A2I7N0","P12763",
                    "P17690","P02769","P02676","P50448","P01030","P01966","P00735",
                    "Q03247","Q3ZBS7","Q2UVX4","Q9TT36","Q28085","Q3SX09","Q3ZBD7",
                    "Q3MHN2","Q9TRI1","P15497","Q95121","Q05443","P02070","Q2KIS7",
                    "Q3MHH8","Q3T052","Q3KUS7","Q1RMK2","Q2TBQ1","Q05B55","A2I7N1",
                    "P04258","Q2KJ62","Q0IIK2","Q3MHN5","P02662","P02663","P02666",
                    "P02668","P31096","P02754","P00711","P62894","Q29443","P19001",
                    "A2AB72","Q8VED5","Q61726","Q3ZAW8","P50446","Q497I4","Q9D312",
                    "Q922U2","Q8BGZ7","A2A4G1","Q9QWL7","Q6IME9","Q6NXH9","A2VCT4",
                    "P07744","Q6IFZ6","Q6IFX2","Q9R0H5","Q3TTY5","Q0VBK2","Q61782",
                    "A2A5Y0","Q99PS0","Q9D646","P05784","Q9DCV7","Q9Z2K1","P07477",
                    "P05787","Q7Z794","Q9BYR9","Q9BYQ5","Q9BYR8","Q9BYQ7","Q3LI72",
                    "Q9BYR4","Q9BYQ8","P60413","P19012","Q2M2I5","O95678","Q01546",
                    "Q99456","Q9H552","P35900","Q3SY84","Q8N1A0","Q5XKE5","P12035",
                    "Q9C075","P08729","Q7Z3Y8","Q7RTS7","Q7Z3Y9","Q7Z3Z0","Q7Z3Y7",
                    "P08727","Q3KNV1","Q86YZ3","P20930","Q5D862") #MaxQuant contaminant list
  
  cat("Cleaning up raw data...\n")
  
  before <- length(unique(raw.data$Protein.Group))
  
  #Data Clean-up - Parameters set at: Q.value <= 0.01; 1+ Peptide per Protein
  raw.data.filtered <- raw.data %>%                                      #Takes our imported raw data
    group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
    filter(sum(grepl(Protein.Group[1], contaminants)) == 0) %>%           #filter out groups that match a given contaminant
    filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
    group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
    filter(length(unique(paste0(Modified.Sequence, Precursor.Charge)))>1) #Filters for lengths > 1 of pasted sequence 
  
  #1 peptide per protein: 
  #raw.data.filtered <- raw.data %>%                                      #Takes our imported raw data
  #group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
  #filter(sum(grepl(Protein.Group[1], contaminants)) == 0) %>%           #filter out groups that match a given contaminant
  #filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
  #group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
  #filter(paste0(Precursor.Charge)>1)
  
  after <- length(unique(raw.data.filtered$Protein.Group))
  
  removed <- (1-(after/before))*100
  
  cat("\nRemoved", removed, "% of the data due to q.value > 0.01 and/or \n only 1 peptide per protein observed.")
  
  raw.data.filtered <- raw.data.filtered %>%                              #take our qc filtered data                   
    dplyr::select(Run, Protein.Group, Genes, Stripped.Sequence,           #select relevant columns
                  Precursor.Id, Precursor.Normalised) %>% 
    mutate(Precursor.Normalised = na_if(Precursor.Normalised, 0)) %>%
    tidyr::unite(Protein.Info, Genes, Protein.Group, sep = "*")  
  
  raw.data <<- raw.data.filtered
  
}  
raw.data.qc <- import.diann.org(file= "report.tsv")
########################################################################################################

#Define treatment groups
#treatments <- c("treatmentA", "treatmentB", "treatmentC", "etc.")
treatments <- c("DMSO", "6H", "18H")
treatments <- c("GFP", "WT", "34", "47", "169", "360", "400", "453")  #Add different constructs for ADA2
treatments <- c("2711_Rep4", "2711_Rep9", "3363_Rep2", "3363_R148W", "1_C2_C1")

#Generate treatment groups - requires user input
sample.groups <- assign.group(raw.data.qc, groups = treatments)

#Replace sample names with reformatted ones
raw.data.qc <- raw.data.qc %>% inner_join(sample.groups, by = "Run") %>%
  mutate(Run = sample) %>% select(-sample)

    ### (OPTIONAL) Export peptide groups ###
    peptides.raw.qc.f <- raw.data.qc.f %>% dplyr::select(Run, Protein.Info, Precursor.Id, Precursor.Normalised) %>%
      pivot_wider(names_from = Run, values_from = Precursor.Normalised)
    
    write_xlsx(peptides.raw.qc, "peptides.raw.qc.xlsx")

### Visualization ###

    #Boxplot 1 - A. Raw abundances across samples grouped by treatment
    box.precursor.int <- ggplot(raw.data.qc.f, aes(x= Run, y=Precursor.Normalised, fill = group)) +
      geom_hline(yintercept = median(raw.data.qc.f$Precursor.Normalised, na.rm = TRUE),  #this adds a line at the global median
                 color = "red4", linetype = "solid", size = 1) +
      geom_boxplot()+
      labs(title= "A. Distribution of Precursor Intensities Grouped by Treatment", 
           x= "Sample",
           y= "Unnormalized Abundance") +
      scale_y_log10() +
      theme(axis.text.x= element_text(angle = 45, hjust =1)) +
      facet_grid(cols = vars(group), scales = "free_x")
    
    box.precursor.int
    
    ggsave("PrecursorIntentsity_BoxPlot.pdf", plot = box.precursor.int, width = 12, height = 6, dpi = 300, device = cairo_pdf)

############## Step Two: Normalization #########################################

#Global Median Normalization (DIA: Performed on the peptide level) 
# Specify: samples = column_samplenames    values = column_values
peptides.medNorm <- medNorm(raw.data.qc.f, samples = Run, values = Precursor.Normalised)

#Visualize of precursor intensity after normalization
box.precursor.int.norm <- ggplot(peptides.medNorm, aes(x= Run, y=Norm_abundance, fill = group)) +
  geom_hline(yintercept = median(raw.data.qc.f$Precursor.Normalised, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "B. Distribution of Precursor Intensities Grouped by Treatment", 
       x= "Sample",
       y= "Normalized Abundance") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

box.precursor.int.norm

ggsave("PrecursorIntentsity_BoxPlot_Norm.pdf", plot = box.precursor.int.norm, width = 12, height = 6, dpi = 300, device = cairo_pdf)


#Protein Roll-up using maxLFQ algorithm [~1min.]
proteins.medNorm <- max_LFQ(peptides.medNorm, values = "Norm_abundance")

proteins.medNorm.export <- proteins.medNorm %>% data.frame() %>% rownames_to_column(var= "Protein.Info")

write_xlsx(proteins.medNorm.export, "ProteinAbund_Nolog2.xlsx")

    ### (OPTIONAL) Calculate Sequence Coverage for each Protein [~5min.]  ###

    #Define species list and proteome directory
    species.list<-c("Homo sapiens")                          #define species
    #Folder of the location of the single species proteome, do not add the name of the .FASTA
    dir<- "C:/Users/barny/OneDrive/Desktop/ADA2_Minsoo/DIANN_Interactomics_DIA_ADA2/IncludingTechRep/Proteome/UniProt_Human_reviewed_cont_03-25-2014_ADA2_myc_his.fasta"
    #Input UniProt fasta file name, full location and the name of the fasta file be sure that "C:/" 
    file.exists(dir) #To check that the needed file is accessed in the assigned directory 
    
    #Run functions
    digest.one.species(species.list, dir)
    protein.cov(df = peptides.medNorm)
    
    protein.cov.wide <- protein.coverage %>%
      pivot_wider(names_from = Run, values_from = perCov)
    
    write_xlsx(protein.cov.wide, "ProteinCoverage.xlsx")
    
    print(mean(protein.coverage$perCov))
    print(median(protein.coverage$perCov))
    
    #### Testing ####
    
    print(length(unique(peptides.medNorm$Protein.Info)))
    #Proteins: 7,121
    
    protein.cov.wide <- protein.coverage %>% pivot_wider(names_from = "Run", values_from = "perCov")
    print(length(unique(protein.cov.wide$Protein.Info)))
    #Proteins: 7,048
    
    proteins.medNorm.edit <- proteins.medNorm %>% data.frame() %>% rownames_to_column(var = "Protein.Info")
    print(length(unique(proteins.medNorm.edit$Protein.Info)))
    #Proteins: 7,121
    
    #See which proteins don't have a coverage calculation
    temp <- data.frame(setdiff(proteins.medNorm.edit$Protein.Info, protein.cov.wide$Protein.Info))
    
#Log2 Transformation
proteins.medNorm.log2 <- proteins.medNorm %>% data.frame() %>% log2()
    
proteins.medNorm.export <- proteins.medNorm.log2 %>% data.frame() %>% rownames_to_column(var= "Protein.Info")
    
write_xlsx(proteins.medNorm.export, "ProteinAbund_Nolog2.xlsx")

##Prepare the data for visualization 
sample.groups.v2 <- sample.groups %>% 
  mutate(Run = ifelse(grepl("^[0-9]", sample), paste0("X", sample), sample)) %>% 
  select(-sample)

### Visualization of Protein Level Data ### - (1) Distribution of protein abundances by injection and (2) number of proteins identified per sample

boxplot.proteins.log2 <- proteins.medNorm.log2 %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
  inner_join(sample.groups.v2, by = "Run") %>% na.omit()

#Boxplot 1 - Normalized Log2 Protein Abundances grouped by Treatment
boxplot.proteins.log2 <- ggplot(boxplot.proteins.log2, aes(x= Run, y=Log2_medNorm, fill = group)) +
  geom_hline(yintercept = median(boxplot.proteins.log2$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "A. Distribution of Protein Abundances Grouped by Treatment", 
       x= "Sample",
       y= "Log2 Median Normalized Protein Abundance") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

boxplot.proteins.log2

ggsave("boxplot.proteins.log2.pdf", plot = boxplot.proteins.log2, width = 12, height = 6, dpi = 300, device = cairo_pdf)

#Barplot 2 - Unique Protein IDs Identified per Sample
protein.count <- proteins.medNorm.log2 %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>% 
  na.omit() %>% info.group() %>% inner_join(sample.groups.v2, by = "Run") %>% na.omit()

protein.count.plot <- ggplot(protein.count, aes(x= Run, y=Count, fill = group)) +
  geom_col()+
  labs(title= "A. Unique Protein IDs Identified per Sample", 
       x= "Sample",
       y= "Count") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(group), scales = "free_x")

protein.count.plot

ggsave("Proteins.Ided.PerInj.pdf", plot = protein.count.plot, width = 12, height = 6, dpi = 300, device = cairo_pdf)

#PCA Plot: Data validation at the protein level 
#Copy our data
data <- proteins.medNorm.log2

#Replace all missing values with 1
data[is.na(data)] <- 1

#Transpose the data
data.t <- data %>% t() %>% data.frame()

#Store our sample names
sample_names <- data.t %>% rownames()

#Run the PCA analysis
data.prcomp <- data.t %>% select(where(~ var(.) != 0)) %>%  #select columns without zero variance
  prcomp(scale = TRUE)

#Label our samples with their treatment groups
data.groupInfo <- data.t %>%
  rownames_to_column(var = "Run") %>%
  inner_join(sample.groups.v2, by = "Run")

#Plot the PCA analysis
PCA <-   autoplot(data.prcomp, data=data.groupInfo, label= FALSE, size=4, colour = 'group') +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  labs(title = "A. Principal Component Analysis Plot")

PCA

ggsave("PCA.pdf", plot = PCA, width = 14, height = 6, dpi = 300, device = cairo_pdf)



 ### (OPTIONAL) Bait Normalization ### - For Interactomics Applications:
      
      #Define bait using UniProt Accession ID (Ex. Human GAPDH)
      bait <- "your.bait.UniProt.accession.ID"
      
      bait <- "Q9NZK5" #ADA2
      
      #Retrieve bait protein information
      bait.quant <- proteins.medNorm.log2 %>% rownames_to_column(var = "Protein.Info") %>%
        separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
        filter(Accession %in% bait) %>% 
        pivot_longer(cols = 3:ncol(.), names_to = "Run", values_to = "bait") %>% 
        dplyr::select(Run, bait)
      
      #Barplot 1 - Visualize Bait Levels Across Samples
      bait.quant.exp <- bait.quant %>% inner_join(sample.groups.v2, by = "Run") %>%
        ggplot(aes(x= Run, y=bait, fill = group)) +
        geom_col()+
        labs(title= "A. Bait Protein Abundances", 
             x= "Sample",
             y= "Log2 Normalized Abundance") +
        theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
        facet_grid(cols = vars(group), scales = "free_x")
      
      bait.quant.exp
      
      ggsave("BaitAbundances_Overall.pdf", plot = bait.quant.exp, width = 12, height = 6, dpi = 300, device = cairo_pdf)
      
      #Normalize to bait protein levels
      proteins.log2.baitNorm <- proteins.medNorm.log2 %>% rownames_to_column(var = "Protein.Info") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2" ) %>%
        inner_join(bait.quant, by = "Run") %>% mutate(Log2norm = Log2 - bait) %>%
        dplyr::select(Protein.Info, Run, Log2norm) %>%
        pivot_wider(names_from = "Run", values_from = "Log2norm") %>%
        column_to_rownames(var = "Protein.Info")
      
      boxplot.proteins.log2.baitnorm <- proteins.log2.baitNorm %>%
        rownames_to_column(var = "Protein.Info") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
        inner_join(sample.groups.v2, by = "Run") %>% na.omit()

      #Boxplot 2 - Bait Normalized Log2 Abundances grouped by Treatment
      global.bait.norm <- ggplot(boxplot.proteins.log2, aes(x= Run, y=Log2_medNorm, fill = group)) +
        geom_hline(yintercept = median(boxplot.proteins.log2$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
                   color = "red4", linetype = "solid", size = 1) +
        geom_boxplot()+
        labs(title= "A. Distribution of Values Grouped by Treatment Normalized to the Bait Protein", 
             x= "Sample",
             y= "Log2 Abundance Normalzied to Bait") +
        theme(axis.text.x= element_text(angle = 45, hjust =1)) +
        facet_grid(cols = vars(group), scales = "free_x")
      
      global.bait.norm
      
      ggsave("BaitAbundances_Overall_Norm.pdf", plot = global.bait.norm, width = 12, height = 6, dpi = 300, device = cairo_pdf)

      
##################If performing statistics in Prism- code to make a heatmap: 
#Import a list of differentially expressed proteins for all variants being studied: 
diff.prot  <- read_xlsx("ProteinAbund_log2.xlsx",  sheet= "Sig_5") #Call the column Gene
sig_genes <- diff.prot %>% data.frame() %>% unique()

diff.prot.baitnorm <- boxplot.proteins.log2.baitnorm %>%
  separate(Protein.Info, into = c("Gene", "Accession"), sep = "\\*") %>%
  filter(Gene %in% sig_genes$Significant) %>%
  unite(Protein.Info, Gene, Accession, sep = "*") %>%
  select(-group) %>%
  pivot_wider(., names_from = Run, values_from = Log2_medNorm ) %>%
  column_to_rownames(., var = "Protein.Info")


#Define the control group
control <- treatments[2]  #WT in the case of ADA2

#Get Control information
data.control <- diff.prot.baitnorm %>% dplyr::select(ends_with(control)) %>% 
  mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE))

#Calculate Fold Changes for Original Data Frame
FC.log2 <- diff.prot.baitnorm %>% mutate(Mean = data.control$mean, StDev = data.control$sd) %>%
  mutate_at(vars(-Mean, -StDev), ~ (.-Mean)) %>% dplyr::select(-Mean, -StDev) %>% filter(rowSums(!is.na(.)) > 0)

#This block calculates the average fold change of each protein per treatment group
Avg.FC.log2 <- FC.log2 %>% average.FC() %>% filter(rowSums(!is.na(.)) > 2) %>%
  separate_wider_delim(cols = Protein.Info, delim = '*', names = c("Protein", "UniProt.Id"))


write_xlsx(Avg.FC.log2, "HeatMap_SigEnriched_BaitNormalized.xlsx")

#Make a heatmap in prism or continue below to clustering

############## Step Three: Statistics ##########################################

#Protein Filtering - If no NA filtering is needed, set points = 1
#Filters for proteins with at least x observations per treatment group
proteins.filtered <- filter.NA(proteins.medNorm.log2, points = 2) #wide format

### Options: ###
# Multiple t.tests w/ p.adj(method = "fdr")  -  For comparison of 2 groups
# ANOVA / Mixed Model  -  For comparison of >= 3 groups

#Perform Multiple unpaired t.tests with equal variance - Visualize via Volcano Plot [Ex. DIA]
comparison <- c("treatment", "control") #define in order of [treatment - control]

comparison <- c("18H", "DMSO")

comparison <- c("360", "GFP")

#Note on factors: factors are used for categorical data to designate the levels, however the dataframe is not changed when factors are added

proteins.log2.long.stat <- proteins.filtered %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
  na.omit() %>% inner_join(sample.groups.v2, by = "Run") %>%
  filter(group %in% comparison) %>%
  mutate(group = factor(group, levels = comparison)) %>%
  group_by(Protein.Info, group) %>% filter(n() > 2) %>% ungroup() %>%
  group_by(Protein.Info) %>% filter(n_distinct(group) == 2) %>%
  do(tidy(t.test(Log2 ~ group, data = ., var.equal = TRUE)))

#Adjust p-values
result.volcano <- proteins.log2.long.stat %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr", n = length(p.value))) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  column_to_rownames(var = 'Protein.Info') %>% data.frame()

  ###   Simple Volcano Plot - EnhancedVolcano package   ###
volcano <- EnhancedVolcano(result.volcano, lab = rownames(result.volcano), 
                    x = "estimate", y = "p.adj", pointSize = 3.0,
                    FCcutoff = 0.5,  #set fold change cutoff
                    pCutoff = 0.01   #set p.adj cutoff
                    )
    
    #Filter for significant proteins - Enhanced Volcano
    sig.prt.EV <- result.volcano %>% filter(abs(estimate) > 0.5 & p.adj < 0.01) %>%
      rownames_to_column(var = "Protein.Info")
  
  ###                                                  ###
  
  ###   Volcano plot using curvature function          ###
    volcano.curve(result.volcano)
  
  #Filter for significant proteins - Curvature Volcano
    volcano_plot.data <- volcano_plot$data %>% 
      filter(significant_med == TRUE | significant_high == TRUE) %>%
      rownames_to_column(var = "Protein.Info")
    
  ###                                                 ###
  
  ###   Venn Diagram Comparing Both Volcano Plots ###
    venn.plot <- venn.diagram(
      x = list(Set1 = sig.prt.EV$Protein.Info, Set2 = volcano_plot.data$Protein.Info),
      category.names = c("EV", "Curve"),
      filename = NULL,
      output = FALSE,
      fill = c("blue", "yellow"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.col = c("black", "black")
    )
    
    grid.newpage()
    grid.draw(venn.plot)
  
  ###                                            ###

#Perform One-way ANOVA
  temp <- proteins.filtered %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
    inner_join(sample.groups.v2, by = "Run") %>%
    filter(is.finite(Log2))
  
  anova_result <- aov(Log2 ~ group, data = temp)  
  
  summary(anova_result)
  
  tukey_result <- TukeyHSD(anova_result)
  
  tukey_df <- as.data.frame(tukey_result$group)
  
  tukey_df_sig <- filter(tukey_df, `p adj` < 0.01 )

############### Additional Analysis ############################################

#Step Seven: Calculate Fold Changes normalized to user-defined control #The following is carried out in wide format

#Define the control group
control <- treatments[1]  

#Get Control information
data.control <- proteins.medNorm.log2 %>% dplyr::select(ends_with(control)) %>% 
  mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE))

#Calculate Fold Changes for Original Data Frame
FC.log2 <- proteins.medNorm.log2 %>% mutate(Mean = data.control$mean, StDev = data.control$sd) %>%
  mutate_at(vars(-Mean, -StDev), ~ (.-Mean)) %>% dplyr::select(-Mean, -StDev) %>% filter(rowSums(!is.na(.)) > 0)

#Plot Fold Changes per sample
FC.log2 %>% rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2FC") %>%
  inner_join(sample.groups.v2, by = "Run") %>%
  ggplot(aes(x = Run, y = Log2FC, fill = group)) + 
  geom_boxplot(outliers = TRUE) + theme_bw() +
  labs(title = "A. Boxplot of Average Log2 Fold Changes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

#This block calculates the average fold change of each protein per treatment group
Avg.FC.log2 <- FC.log2 %>% average.FC() %>% filter(rowSums(!is.na(.)) > 2) %>%
  separate_wider_delim(cols = Protein.Info, delim = '*', names = c("Protein", "UniProt.Id"))

#This block reformats the data into long format, then plots a boxplot of the average fold changes per treatment
Avg.FC.log2 %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Treatment", values_to = "AvgFC") %>%
  ggplot(aes(x = Treatment, y = AvgFC, fill = Treatment)) + 
  geom_boxplot(outliers = TRUE) + theme_bw() +
  labs(title = "A. Boxplot of Average Log2 Fold Changes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value")

#Step Eight: Clustering (Using K-Means Clustering)

clusterDF <- Avg.FC.log2
clusterDF[is.na(clusterDF)] <- 0

clusterDF.mean <- clusterDF %>% #filter(Avg.50 != 0 | Avg.60 != 0) %>%
  unite(Protein.Info, Protein, UniProt.Id, sep = "*") %>% column_to_rownames(var = "Protein.Info")

######## Run this when deciding the number of kmeans clusters ###############
#Change this to change the limit to amount of clusters you calculate for
clusLength <- 20

# #remove the normalized condtion column that is all 0
# clusterDF.mean <- clusterDF.mean.all[ ,-1]

result <- fviz_nbclust(clusterDF.mean, kmeans, method = "wss", k.max = clusLength)
print(result)
elbow <- result$data

#start for the silhouette score calculation
silPoints <- 1
df.silhouette <- data.frame(k = 2:clusLength)

for(i in 1:20){  for(k in 2:clusLength){
  # Perform k-means clustering
  kmeans_result <- kmeans(clusterDF.mean, centers = k)
  # Compute silhouette information using silhouette function
  sil_info <- silhouette(kmeans_result$cluster, dist(as.matrix(clusterDF.mean)))
  silPoints[k-1] <- summary(sil_info)$avg.width
}
  silhouetteDf <- data.frame(k = 2:clusLength, Score = silPoints)
  df.silhouette <- left_join(df.silhouette, silhouetteDf, by = 'k')
}

row_sds <- apply(df.silhouette, 1, sd)
row_mean <- apply(df.silhouette, 1, median)

df.silMean <- data.frame(k = df.silhouette$k, means = row_mean, stdev = row_sds)

silhouttePlot <- ggplot(silhouetteDf, aes(x = k, y = Score)) +
  geom_line() +
  theme_classic()
silhouttePlot

silhouttePlot <- ggplot(df.silMean, aes(x = k, y = means)) +
  geom_line() +
  theme_classic() + labs(title = "Average of 100 Silhouette Scores")
silhouttePlot

clusAmount = 4 #Change clusAmount to be the amount of clusters you decide on, based on the above plots

d <- clusterDF.mean

# Perform hierarchical clustering
dist_matrix <- dist(d)  # Compute distance matrix
hc <- hclust(dist_matrix, method = "complete")  # Perform hierarchical clustering

# Perform hierarchical clustering
clustering_rows <- hclust(dist(d)) 

# Create the heatmap
set.seed(52)   
pheatmap.kmean<-pheatmap(d, kmeans_k = clusAmount, cluster_rows = TRUE, cluster_cols = FALSE,
                         clustering_method = "complete",  # You can change the method as needed
                         cutree_rows = 1,  # Number of clusters for rows, estimated from dendrogram above
                         main = "Heatmap of Proteins with Hierarchical Clustering using K-means", 
                         labels_row = d$Protein.Group,
                         color = hcl.colors(50, "Temps"), fontsize_row = 10, display_numbers = TRUE)

cluster.values <- pheatmap.kmean$kmeans$cluster %>% data.frame() %>% rownames_to_column() %>%
  setNames(c("Protein.Info","Cluster"))

clusterDF.mean.fixed <- clusterDF.mean %>% 
  rownames_to_column(var = "Protein.Info") %>% full_join(cluster.values, by = "Protein.Info") %>%
  separate_wider_delim(col = Protein.Info, names = c("Protein", "UniProt.AC"), delim = "*")

#Export proteins with cluster assignment
#write.csv(clusterDF.mean.fixed, "Avg.FC.log2.A.Kmean.csv")

