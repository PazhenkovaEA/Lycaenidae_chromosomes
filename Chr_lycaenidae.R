library(pafr)
library(dplyr)
require(RIdeogram)
library(syntenyPlotteR)

ali <- read_paf("Lbel_Lcor.paf")

# filtration

long_ali <- subset(ali, alen >= 2e3 & mapq >= 60)

# dotplot
a <- dotplot(long_ali, label_seqs=TRUE, xlab = "L. bellargus", ylab = "L. coridon", dashes = T, order_by = "provided", line_size	= 2,
             ordering = list(c(1:44, "Z"), c(1:89, "Z"))) +
  theme_bw()+
  theme(axis.title =element_text(size =50, face = "italic"),axis.text = element_blank())

a$layers[[6]]$aes_params$size <- 8 # to change chromosome names size
a$layers[[5]]$aes_params$size <- 8
a$layers[[4]]$aes_params$colour = 'blue'
a$layers[[3]]$aes_params$colour = 'blue'


# Create long syntenic blocks
long_ali = long_ali[grep("chromosome", long_ali$qname),]
long_ali = long_ali[grep("chromosome", long_ali$tname),]
# filter out missmappings
long_ali <- long_ali[order(long_ali[,"tname"], long_ali[,"tstart"],long_ali["qname"] ),]
long_ali <- long_ali %>% mutate(keep = ifelse((qname == lag(qname)) | (qname == lead(qname)), "keep", "remove"))
long_ali <- long_ali[long_ali$"keep" == "keep",]

filtered_ali <- long_ali %>% group_by(tname, qname, strand) %>% mutate(count = n()) %>% filter(count >= 100)
filtered_ali <- filtered_ali  %>% arrange(tname, qname, strand)
a <- filtered_ali %>% group_by(tname, qname, strand)%>% 
  mutate(start = dplyr::first(tstart), end = dplyr::last(tend), tarSt =dplyr::first(qstart), tarEnd = dplyr::last(qend))

# Format for syntenyPlotteR
a$chr <- str_split_fixed(a$tname, "_", 3)[,3]
a$tarChr <- str_split_fixed(a$qname, "_", 3)[,3]
a$tarSpecies <- "Lysandra coridon"
a <- ungroup(a)
a <- a %>% select(chr, start, end, tarChr,tarSt, tarEnd, strand, tarSpecies)
a <- unique(a)
write.table(a, "synt_cor_bell1.txt", row.names = F, col.names = F, quote =F, sep = "\t")
# Draw macrosynteny plot
draw.pairwise1("synt_cor_bell1.txt", "test5", "bellargus_chr.txt", "coridon_chr.txt", "Lysandra bellargus","Lysandra coridon")




