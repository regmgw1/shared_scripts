# hacked by Gareth Wilson - original function MEDIPS.exportWIG forms part of MEDIPS package by Lukas Chavez

exportWigTrim <- function (data = NULL, file = NULL, raw = FALSE, descr = "", pattern.density = FALSE) 
{
    if (class(data) != "MEDIPSset") 
        stop("Must specify a MEDIPSset object.")
    if (is.null(file)) {
        stop("Must specify an output file.")
    }
    sample_name = sample_name(data)
    genome_chr = genome_chr(data)
    bin_size = bin_size(data)
    chr_names = chr_names(data)
    if (pattern.density) {
        output_data = genome_CF(data)
        sample_name = paste("Pattern density ", distFunction(data), ", ", fragmentLength(data), sep = "")
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
        maxLimit, " autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", 
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
        end<-length(temp_data) - 1
        write.table(temp_data[1:end], file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
       }
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
	end<-length(temp_data) - 1
        write.table(temp_data[1:end], file = file, sep = "", quote = F, 
            row.names = F, col.names = F, append = T)
    }
}
