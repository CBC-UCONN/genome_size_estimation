library(tidyverse) 

# The purpose of the script is to demonstrate in principle how
# distributions of k-mer frequencies can be used to 
# infer genome size in a really simple way. Actual 
# empirical data are much messier, so methods implemented
# in software such as genomescope are preferable. 

# The script works by simulating a random diploid genome 
# and sequencing reads from that genome. 


# first specify genome length (don't pick a realistic eukaryotic number)
G <- 100000

# In diploid organisms, there are usually differences
# between homozygous chromosomes. Simulate set the frequency here. 

# mutations:
	# set an arbitrary translation vector to "mutate" the homologous chromosome
	mut <- c(A="T",C="G",G="C",T="A")
	# heterozygosity (fraction of heterozygous sites)
	hz <- 0.00
	# expected number of variable sites
	nvar <- round(hz * G)


# create random genome sequence
seq1 <- sample(c("A","C","G","T"),size=G,replace=TRUE) 
# create homologous sequence (with heterozygous sites)
seq2 <- seq1
	if(nvar > 0){
		s <- sample(1:G,nvar)
		seq2[s] <- mut[seq2[s]]
	}

# concatenate bases into sequences
seq1 <- paste(seq1, collapse="")
seq2 <- paste(seq2, collapse="")

# put both sequences into a list object
seqlist <- list(seq1,seq2)

# now we have a diploid genome with heterozygous sites
# simulate sequencing data with error + chopping up k-mers

# expected per-base sequencing coverage
ec <- 60

# read length
rl <- 100

# number of reads needed to achieve coverage
rn <- round(G/rl * ec)

# random read start positions:
rstart <- sample(1:(G-rl+1),rn,replace=TRUE)

# which chromosome should the read come from?
rchrom <- sample(1:2,rn,replace=TRUE)

# how many errors should the read have?
	# per read error rate
	e <- 0.000 * rl
	# add error to this many reads
	rerror <- round(e * rn)
	# place error in these positions:
	rerrpos <- sample(1:rl,rn,replace=TRUE)

# initialize read vector of length 'rn'
reseq <- character(rn)

# extract first 'rerror' sequences, add errors
if(rerror > 0){
	for(i in 1:rerror){
		# extract read
		reseq[[i]] <- substr(
			seqlist[[rchrom[i]]],
			rstart[i],
			(rstart[i] + rl - 1))
		# add error
		substr(reseq[[i]],rerrpos[i],rerrpos[i]) <- mut[substr(reseq[[i]],rerrpos[i],rerrpos[i])]
		if((i %% 1000) == 0){print(i)}
	}
}

# extract remaining sequences without errors
for(i in (rerror+1):rn){
	# extract read
	reseq[[i]] <- substr(
		seqlist[[rchrom[i]]],
		rstart[i],
		(rstart[i] + rl - 1))
	if((i %% 1000) == 0){print(i)}
}


# decompose reads into k-mers

# k-mer length
kl <- 20

# number of k-mers in a read
kn <- rl - kl + 1

# initialize kmer vector of length 'kn * rn'
kmers <- character(kn * rn)

# extract 'kn' k-mers from each of 'rn' reads
for(i in 1:rn){
	for(j in 1:kn){
		k <- (i - 1) * kn + j
		kmers[k] <- substr(reseq[[i]],j,(j+kl-1))
	}
	if((i %% 1000) == 0){print(i)}
}

# collect frequency of each k-mer (note, "non-canonical" k-mer tally)
cc <- table(kmers)

# collect frequencies of k-mer frequencies (k-mer spectrum)
ff <- table(cc) %>% data.frame()
	ff[,1] <- ff[,1] %>% as.character() %>% as.numeric()

# plot k-mer spectrum
plot(ff,pch=20,type="b",xlab="k-mer frequency",ylab="frequency of k-mer frequencies")


# How can we infer the size of the genome from the k-mer spectrum?

# with no error, no heterozygosity, and an unbiased sequencing 
# process, we could calculate it as the total number of k-mers
# divided by the average k-mer frequency:

# using average k-mer frequency:
Kmean <- sum(ff[,1]*ff[,2])/sum(ff[,2])
Ktotal <- sum(ff[,1]*ff[,2])
Ktotal / Kmean + kl - 1
# this won't be correct if e > 0

# sequencing error complicates things, so we want to try to exclude 
# those very low frequency k-mers that are probably errors
# exclude k-mers with frequency < 5
ex <- 1:4

Kmean <- sum(ff[-ex,1]*ff[-ex,2])/sum(ff[-ex,2])
Ktotal <- sum(ff[-ex,1]*ff[-ex,2])
Ktotal / Kmean + kl - 1
# this won't be correct if hz > 0

# heterozygosity complicates things further. the average
# k-mer frequency is no longer useful. now our plot has 2
# peaks. the left-most peak is composed of heterozygous k-mers. 
# we want to use the mode of the right-most peak. 
ex <- 1:4
Kmode <- ff[-ex,][which.max(ff[-ex,2]),1]
Ktotal <- sum(ff[-ex,1]*ff[-ex,2])

Ktotal/Kmode + kl - 1


