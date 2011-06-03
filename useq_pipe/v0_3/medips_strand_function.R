genomeVectorStrand <- function (data = NULL, extend = 400, bin_size = 50,strandSpec = NULL) 
{
    if (class(data) != "MEDIPSset") 
        stop("Must specify a MEDIPSset object.")
    #if (strandSpec == NULL) 
    #   stop("Must specify a strand (+ or -).")
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

exportWigStrand <- function (data = NULL, file = NULL, raw = FALSE, descr = "", pattern.density = FALSE, sample_name = "", color = "255,0,0") 
{
    if (class(data) != "MEDIPSset") 
        stop("Must specify a MEDIPSset object.")
    if (is.null(file)) {
        stop("Must specify an output file.")
    }
    genome_chr = genome_chr(data)
    bin_size = bin_size(data)
    chr_names = chr_names(data)
    if (pattern.density) {
        output_data = genome_CF(data)
        sample_name = paste("Pattern density ", distFunction(data), 
            ", ", fragmentLength(data), sep = "")
        descr = sample_name
    }
    else {
        if (!raw) {
            output_data = genome_norm(data)
            sample_name = paste(sample_name, "_normalized", sep = "")
        }
        else {
            output_data = genome_raw(data)
            output_data = output_data/(number_regions(data)/1e+06)
            sample_name = paste(sample_name, "_raw", sep = "")
        }
    }
    maxLimit = max(output_data)
    maxLimit = maxLimit * 0.75
    header_out = paste("track type=wiggle_0 name=\"", sample_name, 
        "\" description=\"", descr, "\" visibility=full viewLimits=0:", 
        maxLimit, " autoScale=on color=",color," maxHeightPixels=100:50:20 graphType=bar priority=20", 
        sep = "")
    write.table(header_out, file = file, sep = "", quote = F, 
        row.names = F, col.names = F)
    for (i in 1:length(chr_names)) {
        cat(paste("Writing data for ", chr_names[i], "...\n", 
            sep = ""))
        chr_header = paste("fixedStep chrom=", chr_names[i], 
            " start=1 step=", bin_size, " span=", bin_size, sep = "")
        write.table(chr_header, file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
        temp_data = output_data[genome_chr == chr_names[i]]
        write.table(temp_data, file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
    }
}
