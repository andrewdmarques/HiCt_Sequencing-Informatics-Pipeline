library(gsheet)
library(kableExtra)
library(staplr)

################################################################################
# Define Parameters
################################################################################

minAmpliconLength <- 280    # Minimum amplicon length 
maxAmpliconLength <- 300    # Maximum amplicon length
minCoverage <- 50           # Minimum coverage for a variant to be called
minAF <- 2/3                # Minimum allele frequency for a variant to be called
minMQ <- 30                 # Minimum mapping quality for a variant to be called 
minQUAL <- 0                # Minimum base quality for a variant to be called
QCThreshold <- 95           # Minimum percent 50x coverage for a specimen to have a consensus sequence made.

################################################################################
# Define Functions
################################################################################

run_conda_command <- function(env, command){
  # Write a bash script that will run the blast.
  filepath01 <- '.'
  p <- c('#!/Bin/bash', 'source ~/anaconda3/etc/profile.d/conda.sh', paste0('conda activate ', env), command, 'conda deactivate')
  writeLines(p, file.path(filepath01, 'temp.script'))
  system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
  comm <- paste0(file.path(filepath01, 'temp.script'))
  print(command)
  
  # Run the bash script.
  system(comm)
}

make_consensus <- function(i, minAF, minCoverage, ref){
  vsp <- ref$sample_id[i]

  # Determine if there are enough reads to make a consensus.
  if(file.size(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv')) == 0){
    ref$error[i] <- paste0(ref$error[i],'|no consensus made - too few reads|')
    print('|no consensus made - too few reads|')
  }
  if(file.size(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv')) > 0){
  
  # Open the samtools pileup
  pil1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv'), sep='\t',header=FALSE)
  colnames(pil1) <- c('CHROM','POS','REF','COVERAGE','READ','QAUL')
  
  # Determine the occurs of base per position.
  pil2 <- pil1
  pil2$Consensus <- NA
  pil2$A_count <- 0
  pil2$C_count <- 0
  pil2$G_count <- 0
  pil2$T_count <- 0
  # Reorder columns
  pil2 <- pil2[ , c(1,2,3,4,7,8,9,10,11,5,6)]  
  for(j in 1:length(pil2$POS)){
    read <- pil2$READ[j]
    # Handle forward and reverse reads in the same way 
    read <- gsub('.',',',read, fixed = T)
    read <- gsub('a','A', read)
    read <- gsub('c','C', read)
    read <- gsub('g','G', read)
    read <- gsub('t','T', read)
    # Remove the markers for the start and ends of reads
    read <- gsub('^','',read)
    read <- gsub('$','',read)
    # Remove the deletions/insertions 
    
    # Count the number of bases present for each nucleotide.
    pil2$A_count[j] <- lengths(regmatches(read, gregexpr("A", read)))
    pil2$C_count[j] <- lengths(regmatches(read, gregexpr("C", read)))
    pil2$G_count[j] <- lengths(regmatches(read, gregexpr("G", read)))
    pil2$T_count[j] <- lengths(regmatches(read, gregexpr("T", read)))
    if(pil2$REF[j] == 'A'){pil2$A_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
    if(pil2$REF[j] == 'C'){pil2$C_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
    if(pil2$REF[j] == 'G'){pil2$G_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
    if(pil2$REF[j] == 'T'){pil2$T_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
    
    # Determine which nucleotide to use as the consensus.
    a <- pil2$A_count[j]
    c <- pil2$C_count[j]
    g <- pil2$G_count[j]
    t <- pil2$T_count[j]
    
    # Determine which nucleotide has the maximum
    max <- max(a,c,g,t)
    
    # Check that the max is above the threshold.
    if(pil2$COVERAGE[j] == 0){max <- -1}else{
    if(max/pil2$COVERAGE[j] < minAF){max <- -1}}
    
    # Assign the consensus. 
    if(pil2$A_count[j] == max){pil2$Consensus[j] <- 'A'}
    if(pil2$C_count[j] == max){pil2$Consensus[j] <- 'C'}
    if(pil2$G_count[j] == max){pil2$Consensus[j] <- 'G'}
    if(pil2$T_count[j] == max){pil2$Consensus[j] <- 'T'}
    if(max == -1){pil2$Consensus[j] <- 'N'} # If there is no clear maximum, then place an N for consensus.
    if(pil2$COVERAGE[j] < minCoverage){pil2$Consensus[j] <- 'N'} # If the coverage is too low, then place an N for consensus.
  }
  # Make the consensus sequence into a data frame.
  col <- c('fasta')
  consensus <- data.frame(matrix(NA, nrow = 2, ncol = length(col)))
  colnames(consensus) <- col
  consensus$fasta[1] <- paste0('>', vsp,'_2')
  consensus$fasta[2] <- paste(t(pil2$Consensus), collapse = '')
  
  # Save the consensus sequence.
  file <- paste0('./03_Consensus/',vsp,'_consensus2.fasta')
  write.table(consensus, file=file, quote=FALSE, sep='\t', col.names = F, row.names = F)
  
  # Save the data frame with nucleotide positions.
  write.csv(pil2, paste0('./02_Variant-Calling/',vsp,'_base-calls.csv'))
  }
  return(ref)
}
  
major_variants <- function(vsp, minCoverage, minAF, minMQA, minQUAL){
  # vsp <- ref$sample_id[i]
  vcf1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'), comment.char = '#',sep='\t',header=FALSE)
  colnames(vcf1) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENOME')
  # Separate the columns with semicolons to be individual columns 
  # Determine the colnames of the vcf file.
  temp <- strsplit(as.character(data.frame(vcf1$INFO[1])), ";")
  temp <- t(sapply(temp[as.logical(lengths(temp))], function(a) c(a, rep("",max(lengths(temp))-length(a)))))
  col <- as.character(gsub("=.*","",temp))
  # Add the column names to the vcf data frame.
  for(j in 1:length(col)){
    vcf1[col[j]] <- NA
  }
  # Populate the new column from the INFO
  for(j in 1:length(vcf1$POS)){
    temp <- strsplit(as.character(data.frame(vcf1$INFO[j])), ";")
    temp <- data.frame(gsub(".*=","",t(sapply(temp[as.logical(lengths(temp))], function(a) c(a, rep("",max(lengths(temp))-length(a)))))))
    colnames(temp) <- col
    vcf1[j,11:(11+length(col))] <- temp[1,]
  }
  # Save the variant calls as tsv file.
  write.table(vcf1, paste0('./02_Variant-Calling/',vsp,'_allVariants.tsv'), sep='\t',row.names = F, col.names = T)
  vcf1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.tsv'), sep='\t',header=T)
  vcf2 <- vcf1
  
  # Filter the variant calls by minimum coverage, allele frequence, mapping quality, and base quality 
  vcf2$MQA<- as.numeric(vcf2$MQS)/as.numeric(vcf2$DP)/as.numeric(vcf2$AF) # Get the average mapping quality score.
  vcf2 <- subset(vcf2, as.numeric(vcf2$DP) >= minCoverage)
  vcf2 <- subset(vcf2, as.numeric(vcf2$AF) >= minAF)
  vcf2 <- subset(vcf2, as.numeric(vcf2$MQA) >=minMQ)
  # vcf2 <- subset(vcf2, vcf2$QUAL >= minQUAL) # Commented out because the qual score is not correct.
  write.table(vcf2, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.tsv'), sep='\t',row.names = F, col.names = T)
  vcf2 <- vcf2[,1:10]
  
  # Prepare the headers to save the vcf file
  # vcf3 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'),sep=';',header=FALSE)
  vcf3 <- read.csv2(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'), header = F,sep = "!", quote="")
  vcf3 <- subset(vcf3,grepl('#',vcf3$V1))
  temp <- data.frame(rbind(rep(NA, ncol(vcf2)),rep(NA, ncol(vcf2))))
  n <- length(vcf3$V1)
  temp <- as.data.frame(lapply(temp, rep, n))[1:n,]
  colnames(temp) <- colnames(vcf2)
  vcf4 <- rbind(temp,vcf2)
  vcf4$CHROM[1:n] <- vcf3$V1
  vcf4[is.na(vcf4)] <-  ''
  
  # Save the filtered vcf file.
  write.table(vcf4, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), sep = '\t', quote = F, row.names = F, col.names = F)
  # vcf5 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), header = F, sep = '!')
  vcf5 <- read.csv2(paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), header = F,sep = "!", quote="")
  vcf5 <- data.frame(gsub('\t\t\t\t\t\t\t\t\t','',vcf5$V1))
  write.table(vcf5, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), sep = '', quote = F, row.names = F, col.names = F)
}

rename_consensus <- function(vsp){
  # vsp <- ref$sample_id[i]
  fasta <- read.csv2(paste0('./03_Consensus/',vsp,'_consensus.fasta'), header = F,sep = "!", quote="")
  fasta$V1[1] <- paste0('>',vsp)
  write.table(fasta, paste0('./03_Consensus/',vsp,'_consensus.fasta'), sep = '', quote = F, row.names = F, col.names = F)
}
  
extract_stat <- function(i, ref){
  # Read in the mapping stat file.
  stat <- read.csv('temp.txt', header = F, sep = '!')
  
  # Get the percent reads mapped
  temp <- subset(stat,grepl('Percent mapped',stat$V1))
  temp <- sub('.*:', '', temp)
  temp <- gsub(' ', '', temp, fixed = TRUE)
  temp <- gsub('\t', '', temp, fixed = TRUE)
  ref$percent_aligned[i] <- as.numeric(temp)
  
  # Get the average coverage
  temp <- subset(stat,grepl('Average coverage:',stat$V1))
  temp <- sub('.*:', '', temp)
  temp <- gsub(' ', '', temp, fixed = TRUE)
  temp <- gsub('\t', '', temp, fixed = TRUE)
  ref$average_coverage[i] <- as.numeric(temp)
  
  # Get the read 50x coverage for the 550bp target region.
  temp <- read.csv(paste0('./02_Variant-Calling/',ref$sample_id[i],'_basecov.txt'),sep='\t')
  temp <- subset(temp, temp$Coverage >= 50)
  ref$percent50x_coverage[i] <- round(length(temp$Coverage)/550*100,2)
  
  return(ref)
}

make_report <- function(i, ref){
  ## Make the cover sheet.
  sample <- data.frame(t(ref[i,]))
  colnames(sample) <- c(' ')

  table <- sample %>%
    kbl(caption = ' ', 
        font_size = 30)  %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    kable_styling(latex_options="scale_down")
  table <- add_header_above(table, c('HiCt Sample Report' = 2), font_size = 30)
  
  save_kable(table, paste0('temp_page1report.pdf'))

  ## Make the coverage sheet.
  vsp <- ref$sample_id[i]
  cov1 <- read.csv(paste0('./02_Variant-Calling/',ref$sample_id[i],'_basecov.txt'),sep='\t')
  cov1$Pos <- cov1$Pos + 1
  # Make sure that data is present for all 550 positions
  cov2 <- data.frame(matrix(NA, nrow = 550, ncol = 1))
  colnames(cov2) <- c('Pos')
  cov2$Pos <- as.numeric(seq(1:550))
  cov2 <- merge(cov1, cov2, by = 'Pos', all.y = T)
  cov2[is.na(cov2)] <-  0
  
  # Generate the report pdf pages.
  # pdf('temp_page2report.pdf', paper = 'a4', height = 0, width = 0) 
  pdf('temp_page2report.pdf', height = 17.83, width = 13.78)
  par(mfrow=c(4,1),xpd = T, mar = c(8,8,8,8))
  # Plot the whole coverage.
  barplot(cov2$Coverage, border = c('grey40'),col = c('grey40'), xaxs = 'i', yaxs = 'i', cex.sub = 2, cex.axis = 1.8, cex.lab = 2, cex.names=1.62, 
          xlab ='Genome position', ylab = 'Coverage', names.arg = seq(1:550)+22772) # 22772 is the start position for this 550bp fragment.
  
  # Plot the 50x coverage.
  temp <- cov2
  temp$Coverage <- ifelse(temp$Coverage>50,50,temp$Coverage) # Limit coverage to 50
  barplot(temp$Coverage, border = c('grey40'),col = c('grey40'), xaxs = 'i', yaxs = 'i', ylim = c(0,50), cex.axis = 1.8, cex.lab = 2, cex.names=1.62,
          xlab = '', ylab = 'Coverage', names.arg = seq(1:550)+22772) # 22773 is the start position for this 550bp fragment.
  
  # If a consensus was made, plot the coverage map.
  if(!grepl('no consensus', ref$error[i])){
  # Load in the base call data.
  alle1 <- read.csv(paste0('./02_Variant-Calling/',ref$sample_id[i],'_base-calls.csv'),sep=',')
  alle2 <- alle1[,1:10]
  
  # Count the consensus nucleotides
  alle2$Consensus_count <- 0
  alle2$Total <- 0
  for(j in 1:length(alle2$POS)){
    alle2$Total[j] <- sum(alle2$A_count[j]+alle2$C_count[j]+alle2$G_count[j]+alle2$T_count[j])
    alle2$A_count[j] <- alle2$A_count[j]/alle2$Total[j]
    alle2$C_count[j] <- alle2$C_count[j]/alle2$Total[j]
    alle2$G_count[j] <- alle2$G_count[j]/alle2$Total[j]
    alle2$T_count[j] <- alle2$T_count[j]/alle2$Total[j]
    if(alle2$REF[j] == 'A'){
      alle2$Consensus_count[j] <- alle2$A_count[j]
      alle2$A_count[j] <- 0
    }
    if(alle2$REF[j] == 'C'){
      alle2$Consensus_count[j] <- alle2$C_count[j]
      alle2$C_count[j] <- 0
    }
    if(alle2$REF[j] == 'G'){
      alle2$Consensus_count[j] <- alle2$G_count[j]
      alle2$G_count[j] <- 0
    }
    if(alle2$REF[j] == 'T'){
      alle2$Consensus_count[j] <- alle2$T_count[j]
      alle2$T_count[j] <- 0
    }
  }
  alle3 <- alle2[,7:11]
  alle4 <- data.frame(t(alle3))
  colnames(alle4) <- as.character(alle2[,3]+22772)
  
  # Pick nucleotide colors
  labels <- c('A', 'C', 'G', 'T', 'Ref')
  colors <- c('#f54242','#a0e655','#55b3e6', '#e155e6', 'grey40')
  # Plot the variance for each position.
  barplot(as.matrix(alle4),
          ylab = 'Allele Frequency', border = NA, 
          ylim = c(0,1), xaxs = 'i', yaxs = 'i',  col = colors, cex.axis = 1.2, cex.lab = 2, cex.names=1.62)
  legend(1,1.3,legend = labels, colors, horiz = T, cex = 1.8)
  }
  dev.off()
  
  ## Prepare the third page of the report.
  # Only prepare the third page if a consesus sequence was made.
  if(!grepl('no consensus', ref$error[i])){
  # Determine the lineage most likely associated with the mutation.
  mut1 <- read.table(paste0('./02_Variant-Calling/',vsp,'_filteredVariants.tsv'), sep='\t',header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  mut1$POS <- as.numeric(mut1$POS)
  mut1$POS <- mut1$POS + 22772
  mut_ref <- read.csv('Bin/unique-mutations-ref.csv')
  colnames(mut_ref)[1] <- 'num'
  lin_ref <- read.csv('Bin/prop-of-lineage-prevelancy.csv')
  colnames(lin_ref)[1] <- 'mutation'
  col <- c('nt_mutation', 'gene', 'aa_mutation','top_lineages')
  mut2 <- data.frame(matrix('Unknown', nrow = length(mut1$POS), ncol = length(col)))
  colnames(mut2) <- col

  for(j in 1:length(mut1$POS)){
    mut2$nt_mutation[j] <- paste0(mut1$REF[j], mut1$POS[j], mut1$ALT[j])
    # Determine the gene and aa mutation that are affected.
    temp <- subset(mut_ref, mut_ref$nt_mutation == mut2$nt_mutation[j])
    # Determine the top lineages associated with a given mutation.
    if(length(temp$num) > 0){
      mut2$gene[j] <- temp$gene[1]
      mut2$aa_mutation[j] <- temp$aa_mutation[1]
    }
    temp2 <- subset(lin_ref, lin_ref$mutation == mut2$aa_mutation[j])
    if(length(temp2$mutation) > 0){
      top1 <- data.frame(t(temp2))
      colnames(top1)[1] <- 'frequency'
      top1$lineage <- rownames(top1)
      top1 <- top1[2:length(top1$frequency),]
      top1 <- top1[order(top1$frequency, decreasing = TRUE),]
      top1$frequency <- as.character(round(as.numeric(top1$frequency),2))
      mut2$top_lineages[j] <- paste0(top1$lineage[1], ' (',top1$frequency[1],'); ',top1$lineage[2], ' (',top1$frequency[2],'); ',top1$lineage[3], ' (',top1$frequency[3],')')
    }
  }
  
  # Save this mutation report.
  write.csv(mut2, paste0('./05_Summaries/temp_',vsp,'mutations.csv'), row.names = F)
  
  # Rename Columns
  colnames(mut2) <- c('nt_mutation', 'gene', 'aa_mutation','top_lineages (propotion)')
  
  # Save the output table.
  table <- mut2 %>%
    kbl(caption = 'Mutations') %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    kable_styling(latex_options="scale_down")
  
  save_kable(table, paste0('temp_page3report.pdf'))
  }
  
  ## Combine the PDFs together.
  staple_pdf(input_files = list.files(full.names=TRUE,pattern="report.pdf"), output_filepath = paste0('./04_Reports/',vsp,'_report.pdf'), overwrite = TRUE)
  system('rm -r temp_page*')
}

# Run this command to generate a reference file of nucleotide mutations for coding mutations.
unique_mutations <- function(){
  mut1 <- read.csv('Bin/mutations.csv', header = T)
  # Determine the nucleotide mutation
  mut1$nt_mutation <- paste0(mut1$REF,mut1$POS,mut1$ALT)
  # Determine the aa mutation
  mut1$mutation <- paste0(mut1$gene,'_', mut1$type)
  # If it is an intergenic mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(mut1$mutation == 'intergenic_',paste0(mut1$gene,'_', mut1$REF,mut1$POS,mut1$ALT),mut1$mutation)
  # If it is a silent mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(mut1$type == 'silent',paste0(mut1$gene,'_',mut1$type,'_', mut1$REF,mut1$POS,mut1$ALT),mut1$mutation)
  # If it is a deletion mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(grepl('del',mut1$type),paste0(mut1$gene,'_',mut1$POS,'_', mut1$ALT),mut1$mutation)
  # Determine the unique mutations.
  mut2 <- cbind(mut1$nt_mutation,mut1$genes,mut1$mutation)
  colnames(mut2) <- c('nt_mutation', 'gene', 'aa_mutation')
  mut2 <- mut2[!duplicated(mut2),]
  write.csv(mut2, 'Bin/unique-mutations-ref.csv')
}

# This command can be run to generate a new file that determines the prevelancy of a lineage for each particular mutation.
lineage_prevelancy <- function(){
  mut1 <- read.csv('Bin/mutations.csv', header = T)
  met1 <- read.csv('Bin/genomeMetaData.csv')
  # Determine the aa mutation
  mut1$mutation <- paste0(mut1$gene,'_', mut1$type)
  # If it is an intergenic mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(mut1$mutation == 'intergenic_',paste0(mut1$gene,'_', mut1$REF,mut1$POS,mut1$ALT),mut1$mutation)
  # If it is a silent mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(mut1$type == 'silent',paste0(mut1$gene,'_',mut1$type,'_', mut1$REF,mut1$POS,mut1$ALT),mut1$mutation)
  # If it is a deletion mutation, then include the nucleotide mutation.
  mut1$mutation <- ifelse(grepl('del',mut1$type),paste0(mut1$gene,'_',mut1$POS,'_', mut1$ALT),mut1$mutation)
  # Determine the lineages
  mut1$lineage <- ''
  for(j in 1:length(mut1$VSP)){
    temp <- subset(met1,mut1$VSP[j] == met1$lab_id)
    mut1$lineage[j] <- temp$lineage[1]
  }
  # Determine the frequency of a lineage for a given mutation.
  mutation <- unique(mut1$mutation)
  lineage <- unique(mut1$lineage)
  # Make data frame that will have the mutations as column names and each column will be the prevelance of that mutation for a given row (rows are lineages).
  mut2 <- data.frame(matrix(0, nrow = length(mutation), ncol = length(lineage)))
  colnames(mut2) <- lineage
  rownames(mut2) <- mutation 
  # Iterate through mutations and record the relevancy.
  for(j in 1:length(mutation)){
    moi <- mutation[j] # Determine the mutation of interest
    temp <- subset(mut1, mut1$mutation == moi)
    temp <- data.frame(table(temp$lineage))
    temp$Freq <- temp$Freq/sum(temp$Freq)
    for(k in 1:length(lineage)){
      temp2 <- subset(temp, temp$Var1 == lineage[k])
      if(length(temp2$Var1) > 0){mut2[j,k] <- temp2$Freq[1]}
    }
  }
  # Save the data frame.
  write.csv(mut2,'Bin/prop-of-lineage-prevelancy.csv')
}

################################################################################
# Initialization
################################################################################

# Make the necessary file architecture.
directories <- c('00_Raw-Data','temp_Sam-Files','01_Bam-Files','02_Variant-Calling','03_Consensus','04_Reports','05_Summaries')
for(i in 1:length(directories)){
  if(!file.exists(directories[i])){dir.create(directories[i])}
}

# Gather all of the files in the directory.
files <- list.files("./00_Raw-Data/")

# Filter files to only the fastq files. 
files1 <- data.frame(files)
files1$filter <- "n"
for(i in 1:length(files1$files)){
  if(grepl("fastq",files1$files[i]) && grepl("R1",files1$files[i])){
    files1$filter[i] <- "y"
  }
}
ref <- subset(files1, files1$filter == "y")
rownames(ref) <- NULL
ref$files <- gsub('R1','RX',ref$files)

# Save the VSP of the sample as the sample name.
ref$sample_id <- gsub("_.*",'',ref$files)

# Pull the metadata.
metadata <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1QORRoqIGk_2UgtxG1SNj0Ex09KgqURenQvr3SfcWsmk/edit?usp=sharing')
# Determine the sample dates.
ref$sample_date <- NA
for(i in 1:length(ref$sample_date)){
  temp <- subset(metadata,metadata$sample_id == ref$sample_id[i])
  if(length(temp$VSP)>0){ref$sample_date[i] <- temp$sampleCollection_date[1]}
}
# Determine the VSP.
ref$vsp <- gsub("-.*",'',ref$sample_id )
# Determine the report generation date.
ref$report_date <- format(Sys.Date(), "%-m/%-d/%Y")
ref <- ref[,c(5,3,4,6,1,2)] # Reorder reference sheet to be vsp, sample_id, date, filter, files

# Initialize the reference file
ref$reads_aligned <- 0
ref$percent_aligned <- 0
ref$average_coverage <- 0
ref$percent50x_coverage <- 0
ref$error <- ''
ref$command_sam <- ''
ref$command_bam <- ''
ref$command_filt1 <- ''
ref$command_qual <- ''
ref$command_sort <- ''
ref$command_index <- ''
ref$command_pileup1 <- ''
ref$command_pileup2 <- ''
ref$command_var1 <- ''
ref$command_var2 <- ''
ref$command_consensus <- ''
ref$version_bwa <- ''
ref$version_samtools <- ''
ref$version_bcftools <- ''
ref$version_bbmap <- ''

# Make an index file for the reference.
system('cp ./Bin/nc_045512.2_rbd.fasta ./')
system('bwa index nc_045512.2_rbd.fasta')

## Recording Version History.
# BWA
system('bwa 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bwa <- gsub('Version: ','', version$V1[2])
# Samtools
system('samtools 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_samtools <- gsub('Version: ','', version$V1[2])
# Bcftools
system('bcftools 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bcftools <- gsub('Version: ','', version$V1[2])
# BBmap
system('./Bin/bbmap/bbmap.sh -version 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bbmap <- gsub('BBMap version ','', version$V1[2])

################################################################################
# Execute commands for each sample
################################################################################

# Run the commands for each file.
for(i in 1:length(ref$files)){
# for(i in 1:10){
  # Align the reads to the reference.
  ref$command_sam[i] <- paste0('bwa mem -M nc_045512.2_rbd.fasta 00_Raw-Data/',gsub("RX", "R1", ref$files[i]),' 00_Raw-Data/', gsub("RX", "R2", ref$files[i]), ' >  ./temp_Sam-Files/', ref$sample_id[i],'_genome.sam')
  system(ref$command_sam[i])

  # Convert sam to bam files.
  ref$command_bam[i] <- paste0('samtools view -S -b -h ./temp_Sam-Files/', ref$sample_id[i],'_genome.sam',' > ./01_Bam-Files/', ref$sample_id[i],'_genome.bam')
  system(ref$command_bam[i])

  # Filter reads based on size
  ref$command_filt1[i] <- paste0('samtools view -F 8 -h ./01_Bam-Files/', ref$sample_id[i],'_genome.bam', " | awk 'length($10) > ", as.character(minAmpliconLength), " && length($10) < ", as.character(maxAmpliconLength), " || $1 ~ /^@/' | ", 'samtools view -b -h > ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.bam')
  system(ref$command_filt1[i])

  # Filter reads by mapping quality
  ref$command_qual[i] <- paste0('samtools view -h -q ',minMQ,' -b ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.bam > ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.bam')
  system(ref$command_qual[i])

  # Determine the number of reads that aligned.
  lengthAlignedReadsIDs <- length(system(paste0('samtools view ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.bam | cut -f 1 | uniq'), intern = TRUE))

  # Report the number of reads that aligned.
  ref$reads_aligned[i] <- lengthAlignedReadsIDs

  # Sort the sample reads
  ref$command_sort[i] <- paste0('samtools sort -o ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.sorted.bam ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.bam')
  system(ref$command_sort[i])

  # Index the sample
  ref$command_index[i] <- paste0('samtools index ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.sorted.bam')
  system(ref$command_index[i])

  # Make the pileup through bbmap
  ref$command_pileup1[i] <- paste0('./Bin/bbmap/pileup.sh in=./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.sorted.bam', ' ref=nc_045512.2_rbd.fasta out=./02_Variant-Calling/',ref$sample_id[i],'_covstats.txt basecov=./02_Variant-Calling/',ref$sample_id[i],'_basecov.txt countgc=f overwrite=t 2> temp.txt')
  system(ref$command_pileup1[i])
  
  # Extract the stats from bbmap and record them to the stat file.
  ref <- extract_stat(i, ref)

  # Make the pileup through samtools
  ref$command_pileup2[i] <- paste0('samtools mpileup ./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.sorted.bam -f ./nc_045512.2_rbd.fasta > ./02_Variant-Calling/',ref$sample_id[i],'_pileup.tsv')
  system(ref$command_pileup2[i])
  
  # Make a consensus file if there are enough reads.
  if(ref$percent50x_coverage[i] < QCThreshold){ref$error[i] <- paste0(ref$error[i],'|no consensus made - coverage too low|')}
  if(!grepl('too low coverage', ref$error[i])){ref <- make_consensus(i, minAF, minCoverage, ref)}

  # Call variants
  ref$command_var1[i] <- paste0('./Bin/bbmap/callvariants.sh in=./01_Bam-Files/',ref$sample_id[i],'_genome.filt.qual.sorted.bam ref=nc_045512.2_rbd.fasta out=./02_Variant-Calling/',ref$sample_id[i],'_allVariants.vcf shist=./02_Variant-Calling/',ref$sample_id[i],'_variantQualityHisto.txt rarity=0 overwrite=t clearfilters')
  system(ref$command_var1[i])
  
  # If there are enough reads, then attempt to make a consensus.
  if(!grepl('no consensus made',ref$error[i])){ # Check to make sure a concensus should be made.
    # Filter variant calls
    major_variants(ref$sample_id[i], minCoverage, minAF, minMQA, minQUAL)
    
    # Compress the variant calls and index them
    ref$command_var2[i] <- paste0('bcftools view ./02_Variant-Calling/',ref$sample_id[i],'_filteredVariants.vcf -Oz -o ./02_Variant-Calling/',ref$sample_id[i],'_filteredVariants.vcf.gz')
    system(ref$command_var2[i])
    ref$command_var2[i] <- paste0('bcftools index ./02_Variant-Calling/',ref$sample_id[i],'_filteredVariants.vcf.gz')
    system(ref$command_var2[i])
    
    # # Get the consensus sequence.
    ref$command_consensus[i] <- paste0('cat nc_045512.2_rbd.fasta | bcftools consensus ./02_Variant-Calling/',ref$sample_id[i],'_filteredVariants.vcf.gz > ./03_Consensus/',ref$sample_id[i],'_consensus.fasta')
    system(ref$command_consensus[i])
    
    # Rename the consensus sequence. 
    rename_consensus(ref$sample_id[i])
  }
  
  # Generate the sample report.
  make_report(i, ref)
}

# Record the reference file.
write.csv(ref, './05_Summaries/reference.csv', row.names = F)

################################################################################
# Cleanup
################################################################################
system('rm -r ./nc_045512.2*')
system('rm -r temp*')
