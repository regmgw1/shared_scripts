# hacked by Gareth Wilson - original function genomeVector forms part of MEDIPS package by Lukas Chavez

genomeVectorStrand <- function (data = NULL, extend = 400, bin_size = 50,strandSpec = NULL) 
{
    if (class(data) != "MEDIPSset") 
        stop("Must specify a MEDIPSset object.")
    chr = regions_chr(data)
    start = regions_start(data)
    stop = regions_stop(data)
    strand = regions_strand(data)
    chromosomes = chr_names(data)
    chr_lengths = chr_lengths(data)
    no_chr_windows = ceiling(chr_lengths/bin_size)
    supersize_chr = cumsum(no_chr_windows)
    cat("Create the genome vector...\n")
    genomeVec_chr = vector(length = supersize_chr[length(chromosomes)], 
        mode = "character")
    genomeVec_pos = vector(length = supersize_chr[length(chromosomes)], 
        mode = "numeric")
    total = length(chromosomes)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for (i in 1:length(chromosomes)) {
        setTxtProgressBar(pb, i)
        if (i == 1) {
            genomeVec_chr[1:no_chr_windows[i]] = chromosomes[i]
            genomeVec_pos[1:no_chr_windows[i]] = seq(1, chr_lengths[i], 
                bin_size)
        }
        if (i > 1) {
            genomeVec_chr[(supersize_chr[i - 1] + 1):(supersize_chr[i - 
                1] + no_chr_windows[i])] = chromosomes[i]
            genomeVec_pos[(supersize_chr[i - 1] + 1):(supersize_chr[i - 
                1] + no_chr_windows[i])] = seq(1, chr_lengths[i], 
                bin_size)
        }
    }
    genomeVec_signal = vector(length = supersize_chr[length(chromosomes)], 
        mode = "numeric")
    cat("\nDistribute reads over genome...\n")
    for (i in 1:length(chromosomes)) {
        setTxtProgressBar(pb, i)
                genomeVec_signal[genomeVec_chr == chromosomes[i]] = MEDIPS.distributeReads(start[chr == 
            chromosomes[i] & strand == strandSpec], stop[chr == chromosomes[i] & strand == strandSpec], strand[chr == 
            chromosomes[i] & strand == strandSpec], genomeVec_pos[genomeVec_chr == chromosomes[i]], 
            extend)
      
    }
    cat("\n")
    MEDIPSsetObj = new("MEDIPSset", genome_chr = genomeVec_chr, 
        genome_pos = genomeVec_pos, genome_raw = genomeVec_signal, 
        extend = extend, bin_size = bin_size, sample_name = sample_name(data), 
        genome_name = genome_name(data), regions_chr = regions_chr(data), 
        regions_start = regions_start(data), regions_stop = regions_stop(data), 
        regions_strand = regions_strand(data), number_regions = number_regions(data), 
        chr_names = chr_names(data), chr_lengths = chr_lengths(data))
    return(MEDIPSsetObj)
}


