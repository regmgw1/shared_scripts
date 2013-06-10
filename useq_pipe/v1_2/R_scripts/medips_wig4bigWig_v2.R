# hacked by Gareth Wilson - original function MEDIPS.exportWIG forms part of MEDIPS package by Lukas Chavez. This is for MEDIPS v1.10

exportWigTrim <- function (Set = NULL, file = NULL, format = "rpkm", descr = "") 
{
    if (is.null(file)) {
        stop("Must specify an output file!")
    }
    if (format != "count" & format != "rpkm" & format != "pdensity") {
        stop("Format pareeter must be either count, rpkm or pdensity!")
    }
    if (format == "pdensity") {
        if (class(Set) != "COUPLINGset") {
            stop("For exporting pattern densities, must specify a COUPLINGset object!")
        }
        window_size = window_size(Set)
        chr_names = chr_names(Set)
        chr_lengths = chr_lengths(Set)
        chromosomes = chr_names(Set)
        output_data = genome_CF(Set)
        header_out = paste("track type=wiggle_0 name=\"", descr, 
            "\" description=\"", descr, "\" visibility=full autoScale=on color=0,200,100 maxHeightPixels=100:50:20 graphType=bar priority=20", 
            sep = "")
    }
    else if (format == "count") {
        if (class(Set) != "MEDIPSset") {
            stop("For exporting count intensities, must specify a MEDIPSset object!")
        }
        window_size = window_size(Set)
        chr_names = chr_names(Set)
        chr_lengths = chr_lengths(Set)
        chromosomes = chr_names(Set)
        output_data = genome_count(Set)
        header_out = paste("track type=wiggle_0 name=\"", descr, 
            "\" description=\"", descr, "\" visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", 
            sep = "")
    }
    else if (format == "rpkm") {
        if (class(Set) != "MEDIPSset") {
            stop("For exporting rpm intensities, must specify a MEDIPSset object!")
        }
        window_size = window_size(Set)
        chr_names = chr_names(Set)
        chr_lengths = chr_lengths(Set)
        chromosomes = chr_names(Set)
        output_data = (genome_count(Set) * 10^9)/(window_size * 
            number_regions(Set))
        header_out = paste("track type=wiggle_0 name=\"", descr, 
            "\" description=\"", descr, "\" visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", 
            sep = "")
    }
    no_chr_windows = ceiling(chr_lengths/window_size)-1
    supersize_chr = cumsum(no_chr_windows)
    genome_chr = as.vector(seqnames(MEDIPS.GenomicCoordinates(supersize_chr, 
        no_chr_windows, chromosomes, chr_lengths, window_size)))
    write.table(header_out, file = file, sep = "", quote = F, 
        row.names = F, col.names = F)
    for (i in 1:length(chr_names)) {
        cat(paste("Writing data for ", chr_names[i], "...\n", 
            sep = ""))
        chr_header = paste("fixedStep chrom=", chr_names[i], 
            " start=1 step=", window_size, " span=", window_size, 
            sep = "")
        write.table(chr_header, file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
        temp_data = output_data[genome_chr == chr_names[i]]
        write.table(temp_data, file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
    }
}

