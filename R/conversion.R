#' @name biocrates
#' 
#' @title Convert Biocrates xlsx output to \code{SummarizedExperiment} object
#' 
#' @description 
#' The function \code{biocrates} will create a \code{SummarizedExperiment} 
#' from a Biocrates xlsx file. 
#' The function \code{biocrates} takes as input the path to a .xlsx file 
#' (Biocrates output) and additional parameters given to the \code{read.xlsx} 
#' function from the \code{openxlsx} package (e.g. specifying the sheet name or 
#' index by \code{sheet}).
#' 
#' @details 
#' The column "Sample Identification" has to contain unique identifiers
#' (no duplications).
#' 
#' @param file \code{character}
#' @param sheet \code{character} or \code{numeric}, the name or index of the 
#' sheet to read data from
#' @param ... additional parameters given to \code{read.xlsx}
#' 
#' @examples
#' file <- "data/biocrates_test_file.xlsx"
#' biocrates(file = file, sheet = 1)
#' 
#' @usage biocrates(file, sheet, ...)
#'
#' @return 
#' \code{SummarizedExperiment} object
#' 
#' @export
#' 
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr select
#' @importFrom SummarizedExperiment SummarizedExperiment
biocrates <- function(file, sheet, ...) {
    
    xls <- openxlsx::read.xlsx(file, sheet = sheet, ...)
    
    ## colnames is in the first row, assign to colnames attribute 
    ## and remove the first row
    colnames(xls) <- xls[1, ]
    xls <- xls[-1, ]
    
    ## remove columns that are NA
    xls <- xls[, !is.na(colnames(xls))]
    
    ## check if Choline is the last column (this is typically the case)
    if (ncol(xls) != which(colnames(xls) == "Choline")) {
        print("Please check 'file': Column 'Choline' is not the last column in 'file' as expected.")
        print("I do not expect that there are feature columns containing intensities beyond the column 'Choline'.")
        print("I will select all columns until the column 'Choline' and continue.")
    }
    ## unncessesary step but keep it for clarity
    xls <- xls[, seq(1, which(colnames(xls) == "Choline"))]
    
    ## find the columns that contain the metabolites, row 1 contains class,
    ## row 2 contains LOD (row 1 and 2 is NA for columns not containing the 
    ## metabolites)
    inds_met <- !is.na(xls[1, ])
    ## set the first TRUE value to FALSE since it contains the label of the row
    inds_met[which(inds_met)[1]] <- FALSE
    ## set the column C0 if it exists to TRUE (for some older versions of files
    ## this is necessary)
    if ("C0" %in% colnames(inds_met))
        inds_met[1, "C0"] <- TRUE
    
    ## find the rows that contain the samples
    ## find from the back the first FALSE entry, set all following TRUEs to FALSE
    inds_name <- !is.na(xls[, 1])
    inds_name <- inds_name[length(inds_name):1]
    first_false <- which(!inds_name)[1]
    if (length(first_false) == 1)
        inds_name[which(!inds_name)[1]:length(inds_name)] <- FALSE
    inds_name <- inds_name[length(inds_name):1]
    
    ## combine all HMDB ids into a character string for each feature
    inds_name_hmdb <- !inds_name
    inds_name_hmdb[1] <- FALSE
    hmdb <- apply(xls[inds_name_hmdb, inds_met], 2, function(cols_i) {
        cols_i_hmdb <- grep(cols_i, pattern = "HMDB", value = TRUE)
        paste(cols_i_hmdb, collapse = ", ")
    })
    
    ## create rowData 
    rD <- data.frame(feature = make.names(colnames(xls)[inds_met]), 
        feature_original = colnames(xls)[inds_met],
        class = as.character(xls[1, inds_met]),
        HMDB_ids = hmdb)
    rownames(rD) <- rD[["feature"]]
    
    ## create colData
    ## rename column "Sample Identification" to "name" and move to the beginning
    ## of cD
    colnames(xls) <- colnames(xls) |>
        make.names()
    cD <- xls[inds_name, seq_len(min(which(inds_met)) - 1)]
    cD_tmp <- cD
    colnames(cD_tmp) <- tolower(colnames(cD_tmp))
    cD <- data.frame(name = make.names(cD_tmp[, "sample.identification"]), 
        name_original = cD_tmp[, "sample.identification"], cD)
    rownames(cD) <- cD[["name"]]
    
    ## create assay, set values of 0 to NA
    a <- xls[inds_name, inds_met]
    a <- as.matrix(a)
    mode(a) <- "numeric"
    a[a == 0] <- NA
    rownames(a) <- cD[["name"]]

    ## create SummarizedExperiment and return
    SummarizedExperiment::SummarizedExperiment(assays = t(a), 
        rowData = rD, colData = cD)
}


#' @name metaboscape
#' 
#' @title Convert MetaboScape xlsx output to \code{SummarizedExperiment} object
#' 
#' @description 
#' The function \code{metaboscape} will create a \code{SummarizedExperiment} 
#' from a MetaboScape xlsx file. 
#' The function \code{metaboscape} takes as input the path to a .xlsx file 
#' (MetaboScape output) and additional parameters given to the \code{read.xlsx} 
#' function from the \code{openxlsx} package (e.g. specifying the sheet name or 
#' index by \code{sheet}).
#' 
#' @details 
#' Heuristics are run that select the sample columns based on expected 
#' metadata columns. Data upload may lead to insufficient results.
#' 
#' @param file \code{character}
#' @param sheet \code{character} or \code{numeric}, the name or index of the 
#' sheet to read data from
#' @param ... additional parameters given to \code{read.xlsx}
#' 
#' @examples
#' file <- "data/metaboscape_test_file.xlsx"
#' metaboscape(file = file, sheet = 1)
#' 
#' @usage metaboscape(file, sheet, ...)
#'
#' @return 
#' \code{SummarizedExperiment} object
#' 
#' @export
#' 
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr select
#' @importFrom SummarizedExperiment SummarizedExperiment
metaboscape <- function(file, sheet, ...) {
    
    xls <- openxlsx::read.xlsx(file, sheet = sheet, colNames = FALSE, ...)
    
    ## find the first row of xls that contains complete observations, assume
    ## that this row contains the header/colnames, assign colnames and 
    ## remove preceding rows
    nas <- apply(xls, 1, function(rows_i) sum(!is.na(rows_i)))
    ind_colnames <- which(nas == ncol(xls))[1]
    colnames(xls) <- xls[ind_colnames, ]
    xls <- xls[-seq_len(ind_colnames), ]
    
    ## find the columns that contain the features
    cols <- colnames(xls)
    
    ## make all characters lower-case for colnames(.f) to grep small differences
    ## in orthography of colnames
    .cols <- tolower(cols) |>
        make.names()
    .xls <- xls
    colnames(.xls) <- .cols
    
    ## UNICODE \u00E5 refers to "small letter a with a ring above", 
    ## UNICODE \u03B4 refers to delta
    cols_rD <- c("rt..min.", "ccs..\u00E5..", "\u03B4ccs....", "m.z.meas.", 
        "m.meas.", "ions", "ms.ms", "qc.rsd....")
    inds_samp <- which(!.cols %in% cols_rD)
    cols_samp <- cols[inds_samp]
    
    
    ## create rowData
    rD <- data.frame(feature = paste("feature", rownames(.xls), sep = "_"))
    if ("rt..min." %in% .cols) rD$rt_min <- .xls[, "rt..min."] |>
        as.numeric()
    ## UNICODE refers to "small letter a with a ring above"
    if ("ccs..\u00E5.." %in% .cols) rD$ccs_a2 <- .xls[, "ccs..\u00E5.."] |>
        as.numeric() 
    ## UNICODE refers to delta
    if ("\u03B4ccs...."  %in% .cols) rD$deltaccs_percent <- .xls[, "\u03B4ccs...."] |> 
        as.numeric()
    if ("m.z.meas."  %in% .cols) rD$mz <- .xls[, "m.z.meas."] |>
        as.numeric()
    if ("m.meas."  %in% .cols) rD$molecular_mass <- .xls[, "m.meas."] |>
        as.numeric()
    if ("ions"  %in% .cols) rD$ions <- .xls[, "ions"]
    if ("ms.ms"  %in% .cols) rD$msms <- .xls[, "ms.ms"] |>
        as.logical()
    if ("qc.rsd...."  %in% .cols) rD$qc_rsd_percent <- .xls[, "qc.rsd...."] |>
        as.numeric()
    rownames(rD) <- rD[["feature"]]
    
    ## create colData
    cD <- data.frame(name = make.names(cols_samp), name_original = cols_samp)
    rownames(cD) <- cD[["name"]]
    
    ## create assay, set values of 0 to NA
    a <- xls[, cols_samp]
    a <- as.matrix(a)
    mode(a) <- "numeric"
    a[a == 0] <- NA
    colnames(a) <- cD[["name"]]
    rownames(a) <- rownames(rD)
    
    ## create SummarizedExperiment and return
    SummarizedExperiment::SummarizedExperiment(assays = a, 
        rowData = rD, colData = cD)
}


#' @name maxquant
#'  
#' @title Convert MaxQuant txt, tsv, or xlsx output to 
#' \code{SummarizedExperiment} object 
#' 
#' @description 
#' The function \code{maxquant} will create a \code{SummarizedExperiment} from a
#' MaxQuant tsv, txt, or xlsx file.
#' The function \code{maxquant} takes as input the path to a .tsv, .txt, or 
#' .xlsx file (MaxQuant output). Additional parameters can be given to the 
#' \code{read.xlsx} function from the \code{openxlsx} package (e.g. 
#' specifying the sheet name or index by \code{sheet}) or the \code{read.table} 
#' function from the \code{utils} (depending on the \code{type} argument) .
#' 
#' @details 
#' The argument \code{intensity} will specify if the \code{iBAQ} or 
#' \code{LFQ} values are taken. If \code{intensity} is set to \code{"none"}, 
#' heuristics are run that select the sample columns based on expected 
#' metadata columns. In the latter case, data upload may lead to insufficient
#' results.
#' 
#' The argument \code{type} will specify if the data is loaded from 
#' \code{tsv}, \code{txt}, or \code{xlsx} files.
#'  
#' @param file \code{character}
#' @param intensity \code{character}, either \code{"iBAQ"}, \code{"LFQ"}, or 
#' \code{"none"}
#' @param sheet \code{character} or \code{numeric}, the name or index of the 
#' sheet to read data from
#' @param type \code{character}, either \code{"tsv"},  \code{"txt"}, or \code{"xlsx"}
#' @param ... additional parameters given to \code{read.xlsx} (for 
#' \code{type = "xlsx"}) or \code{read.table}
#' (for \code{type = "tsv"}/\code{type = "txt"}) 
#'
#' @examples
#' file <- "data/maxquant_test_file.xlsx"
#' maxquant(file = file, intensity = "LFQ", type = "xlsx", sheet = 1)
#'
#' @return 
#' \code{SummarizedExperiment} object
#'
#' @export
#' 
#' @importFrom openxlsx read.xlsx
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment SummarizedExperiment
maxquant <- function(file, intensity = c("iBAQ", "LFQ", "none"), sheet, 
    type = c("tsv", "txt", "xlsx"), ...) {
    
    intensity <- match.arg(intensity)
    type <- match.arg(type)
    
    if (type == "xlsx")
        f <- openxlsx::read.xlsx(file, sheet = sheet, ...)
    if (type %in% c("txt", "tsv"))
        f <- utils::read.table(file, sep = "\t", dec = ".", header = TRUE, ...)
    
    ## names of proteins is in the first col, assign and remove the first col
    rownames(f) <- f[, 1]
    f <- f[, -1]
    
    ## find the columns that contain the features
    cols <- colnames(f)
    
    ## make all characters lower-case for colnames(.f) to grep small differences
    ## in orthography of colnames
    .cols <- tolower(cols)
    .f <- f
    colnames(.f) <- .cols
    
    if (intensity %in% c("iBAQ", "LFQ")) {
        
        inds_samp <- grep(pattern = intensity, cols)
        cols_samp <- cols[inds_samp]
        
        ## remove the column that only contains "intensity"
        inds_samp <- inds_samp[cols_samp != intensity]
        cols_samp <- cols_samp[cols_samp != intensity]    
    } 
    if (intensity == "none") {
        cols_rD <- c("best.ms.ms", "charges", "count", "evidence.ids",
            "fasta.headers", "first.protein.description", "genes", "gene.names", 
            "id", "length", "majority.protein.ids", "mass", "missed.cleavages",
            "mod..peptide.ids", "mol..weight..kda.", "ms.ms.ids",
            "number.of.proteins", "only.identified.by.site", 
            "oxidation..m..site.ids", "oxidation..m..site.positions",
            "peptide.counts..all.", "peptide.counts..razor.unique.", 
            "peptide.counts..unique.", "peptides", "peptide.ids",
            "peptide.is.razor", "potential.contaminant", "protein.group",
            "protein.ids", "protein.names", "proteins", "q.value", 
            "razor...unique.peptides", "reverse", "sequence", 
            "sequence.coverage....", "sequence.length", "sequence.lengths", 
            "unique.peptides", "unique..proteins.", 
            "unique...razor.sequence.coverage....", 
            "unique.sequence.coverage....")
        inds_samp <- which(!.cols %in% cols_rD)
        cols_samp <- cols[inds_samp]
    }
    
    ## create rowData
    rD <- data.frame(feature = rownames(.f))
    if ("best.ms.ms" %in% .cols) rD$best_MS_MS <- .f[, "best.ms.ms"]
    if ("charges" %in% .cols) rD$charges <- .f[, "charges"]
    if ("count" %in% .cols) rD$count <- .f[, "count"]
    if ("evidence.ids" %in% .cols)
        rD$evidence_ids <- .f[, "evidence.ids"]
    if ("fasta.headers" %in% .cols) rD$fasta_header <- .f[, "fasta.headers"]
    if ("first.protein.description" %in% .cols) 
        rD$first_protein_description <- .f[, "first.protein.description"]
    if ("genes" %in% .cols) rD$genes <- .f[, "genes"]
    if ("gene.names" %in% .cols) rD$gene_name <- .f[, "gene.names"]
    if ("id" %in% .cols) rD$id <- .f[, "id"]
    if ("length" %in% .cols) rD$length <- .f[, "length"]
    if ("majority.protein.ids" %in% .cols)
        rD$majority_protein_ids <- .f[, "majority.protein.ids"]
    if ("mass" %in% .cols) rD$mass <- .f[, "mass"]
    if ("missed.cleavages" %in% .cols) 
        rD$missed_cleavages <- .f[, "missed.cleavages"]
    if ("mod..peptide.ids" %in% .cols)
        rD$mod_peptide_ids <- .f[, "mod..peptide.ids"]
    if ("mol..weight..kda." %in% .cols)
        rD$mol_weight_kDa <- .f[, "mol..weight..kda."]
    if ("ms.ms.ids" %in% .cols) rD$ms_ms_ids <- .f[, "ms.ms.ids"]
    if ("number.of.proteins" %in% .cols) 
        rD$number_of_proteins <- .f[, "number.of.proteins"]
    if ("only.identified.by.site" %in% .cols)
        rD$only_identified_by_site <- .f[, "only.identified.by.site"]
    if ("oxidation..m..site.ids" %in% .cols)
        rD$oxidation_m_site_ids <- .f[, "oxidation..m..site.ids"]
    if ("oxidation..m..site.positions" %in% .cols)
        rD$oxidation_m_site_positions <- .f[, "oxidation..m..site.positions"]
    if ("peptide.counts..all." %in% .cols)
        rD$peptide_counts_all <- .f[, "peptide.counts..all."]
    if ("peptide.counts..razor.unique." %in% .cols)
        rD$peptide_counts_razor_unique <- .f[, "peptide.counts..razor.unique."]
    if ("peptide.counts..unique." %in% .cols)
        rD$peptide_counts_unique <- .f[, "peptide.counts..unique."]
    if ("peptides"  %in% .cols) rD$peptides <- .f[, "peptides"]
    if ("peptide.ids" %in% .cols) rD$peptide_ids <- .f[, "peptide.ids"]
    if ("peptide.is.razor" %in% .cols)
        rD$peptide_is_razor <- .f[, "peptide.is.razor"]
    if ("potential.contaminant" %in% .cols)
        rD$potential_contaminant <- .f[, "potential.contaminant"]
    if ("protein.group" %in% .cols) rD$protein_group <- .f[, "protein.group"]
    if ("protein.ids" %in% .cols) rD$protein_ids <- .f[, "protein.ids"]
    if ("protein.names" %in% .cols) rD$protein_names <- .f[, "protein.names"]
    if ("proteins" %in% .cols) rD$proteins <- .f[, "proteins"]
    if ("q.value" %in% .cols) rD$q_value <- .f[, "q.value"]
    if ("razor...unique.peptides" %in% .cols)
        rD$razor_unique_peptides <- .f[, "razor...unique.peptides"]
    if ("reverse" %in% .cols) rD$reverse <- .f[, "reverse"]
    if ("sequence" %in% .cols) rD$sequence <- .f[, "sequence"]
    if ("sequence.coverage...." %in% .cols)
        rD$sequence_coverage <- .f[, "sequence.coverage...."]
    if ("sequence.length" %in% .cols)
        rD$sequence_length <- .f[, "sequence.length"]
    if ("sequence.lengths" %in% .cols)
        rD$sequence_lengths <- .f[, "sequence.lengths"]
    if ("unique.peptides" %in% .cols) 
        rD$unique_peptides <- .f[, "unique.peptides"]
    if ("unique..proteins." %in% .cols) 
        rD$unique_proteins <- .f[, "unique..proteins."]
    if ("unique...razor.sequence.coverage...." %in% .cols)
        rD$unique_razor_sequence_coverage <- .f[, "unique...razor.sequence.coverage...."]
    if ("unique.sequence.coverage...." %in% .cols)
        rD$unique_sequence_coverage <- .f[, "unique.sequence.coverage...."]
    rownames(rD) <- rD[["feature"]]

    ## create colData
    cD <- data.frame(name = make.names(cols_samp), name_original = cols_samp)
    name_cut <- gsub(paste0("^", intensity), "", cD$name)
    name_cut <- gsub("^[. _]", "", name_cut)
    cD$name_cut <- name_cut
    rownames(cD) <- cD[["name"]]
    
    ## create assay, set values of 0 to NA
    a <- f[, cols_samp]
    a <- as.matrix(a)
    mode(a) <- "numeric"
    a[a == 0] <- NA
    colnames(a) <- cD[["name"]]
    
    ## create SummarizedExperiment and return
    SummarizedExperiment::SummarizedExperiment(assays = a, 
        rowData = rD, colData = cD)
}

#' @name diann
#' 
#' @title Convert DIA-NN tsv output to \code{SummarizedExperiment} object
#' 
#' @description 
#' The function \code{diann} will create a \code{SummarizedExperiment} from a
#' DIA-NN tsv file. 
#' The function \code{diann} takes as input the path to a .tsv 
#' file (DIA-NN output). Additional parameters can be given to the 
#' \code{read.table} function from the \code{utils} package.
#' 
#' @details 
#' \code{intensity} is set by default to \code{"none"}, 
#' Heuristics are run that select the sample columns based on expected 
#' metadata columns. Data upload may lead to insufficient
#' results.
#' 
#' \code{type} is set by default to \code{"tsv"}.
#'  
#' The function \code{diann} is a wrapper function with pre-set arguments 
#' \code{intensity} and \code{type} of the function \code{maxquant}. For further
#' information see also the help page of \code{maxquant}.
#'  
#' @param file \code{character}
#' @param ... additional parameters given to \code{read.table}
#'
#' @examples
#' file <- "data/diann_test_file.tsv"
#' diann(file = file)
#'
#' @return 
#' \code{SummarizedExperiment} object
#'
#' @export
diann <- function(file, ...) {
    .type <- "tsv"
    .intensity <- "none"
    maxquant(file, intensity = .intensity, type = .type, quote = "", ...)
}

#' @name spectronaut
#' 
#' @title Convert Spectronaut xlsx output to \code{SummarizedExperiment} object
#' 
#' @description 
#' The function \code{spectronaut} will create a \code{SummarizedExperiment} 
#' from a Spectronaut xlsx file. 
#' The function \code{spectronaut} takes as input the path to a .xlsx file 
#' (Spectronaut output).
#' 
#' @details 
#' The function requires that the intensity values are stored in the sheet
#' \code{sheetIntensities} and the sample annotations in the sheet
#' \code{sheetAnnotation}.
#' 
#' The sample names are taken from the column \code{"SAMPLE_IDs"} from the 
#' sheet \code{sheetAnnotation}.
#' 
#' @param file \code{character}
#' @param sheetIntensities \code{character} or \code{numeric}, name or index of the 
#' sheet where the intensities are stored
#' @param sheetAnnotation \code{character} or \code{numeric}, name or index of the 
#' sheet where the annotations are stored 
#' @param ... additional parameters given to \code{read.xslx}
#'
#' @examples
#' file <- "data/spectronaut_test_file.xlsx"
#' spectronaut(file = file, sheetIntensities = 1, 
#'     sheetAnnotation = 2)
#' 
#' @usage spectronaut(file, sheetIntensities, sheetAnnotation, ...)
#'
#' @return 
#' \code{SummarizedExperiment} object
#'
#' @export
#' 
#' @importFrom openxlsx read.xlsx
#' @importFrom SummarizedExperiment SummarizedExperiment
spectronaut <- function(file, sheetIntensities = 1, sheetAnnotation = 2, ...) {
    
    xls <- openxlsx::read.xlsx(file, sheet = sheetIntensities, ...)
    cD <- openxlsx::read.xlsx(file, sheet = sheetAnnotation, ...)
    
    ## get the name of the samples
    samps <- cD[, "Sample_IDs"]
    samps <- make.names(samps)
    
    ## names of proteins is in the first col, assign and remove the first col
    if ("PG.ProteinGroups" %in% colnames(xls))
        rownames(xls) <- xls[, "PG.ProteinGroups"]
    else 
        stop("column 'PG ProteinGroups' not present in file")

    ## find the columns that contain the metabolites
    colnames(xls) <- make.names(colnames(xls))
    a <- dplyr::select(xls, samps)

    ## create rowData
    rD <- data.frame(feature = rownames(xls))
    if ("PG.Molecularweight" %in% colnames(xls))
        rD$PG_Molecularweight <- xls[, "PG.Molecularweight"]
    if ("PG.Genes" %in% colnames(xls))
        rD$PG_Genes <- xls[, "PG.Genes"]
    if ("PG.CellularComponent"  %in% colnames(xls))
        rD$PG_CellularComponent <- xls[, "PG.CellularComponent"]
    if ("PG.BiologicalProcess" %in% colnames(xls)) 
        rD$PG_BiologicalProcess <- xls[, "PG.BiologicalProcess"]
    if ("PG.MolecularFunction" %in% colnames(xls)) 
        rD$PG_MolecularFunction <- xls[, "PG.MolecularFunction"]
    rownames(rD) <- rD[["feature"]]

    ## create colData
    colnames(cD)[colnames(cD) == "Sample_IDs"] <- "name"
    cD$name <- make.names(samps)
    cD$name_original <- samps
    rownames(cD) <- cD[["name"]]

    ## create assay, set values of 0 to NA
    a <- as.matrix(a)
    mode(a) <- "numeric"
    a[a == 0] <- NA
    colnames(a) <- cD[["name"]]

    ## create SummarizedExperiment
    SummarizedExperiment::SummarizedExperiment(assays = a, 
        rowData = rD, colData = cD)
}
