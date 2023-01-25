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
#' file <- "path/to/biocrates/object"
#' \donttest{biocrates(file = file, sheet = 1)}
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
    
    ## colnames is in the first row, assign and remove the first row
    colnames(xls) <- make.names(xls[1, ])
    xls <- xls[-1, ]
    
    ## take the columns until Choline
    xls <- xls[, seq(1, which(colnames(xls) == "Choline"))]
    
    ## find the columns that contain the metabolites, row 1 contains class,
    ## row 2 contains LOD (row 1 and 2 is NA for columns not containing the 
    ## metabolites)
    inds_met <- !is.na(xls[1, ])
    ## set the first TRUE value to FALSE since it contains the label of the row
    inds_met[which(inds_met)[1]] <- FALSE
    
    ## create rowData 
    rD <- data.frame(feature = colnames(xls)[inds_met], 
                class = as.character(xls[1, inds_met]))
    rownames(rD) <- rD[["feature"]]
    
    ## find the rows that contain the samples
    inds_name <- !is.na(xls[, 1])
    
    ## create colData
    ## rename column "Sample Identification" to "name" and move to the beginning
    ## of cD
    cD <- xls[inds_name, seq_len(min(which(inds_met)) - 1)]
    cD <- data.frame(name = cD[, "Sample.Identification"], cD)
    cD <- dplyr::select(cD, -c("Sample.Identification"))
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


#' @name maxQuant
#' 
#' @title Convert MaxQuant xlsx or txt output to \code{SummarizedExperiment} object
#' 
#' @description 
#' The function \code{maxQuant} will create a \code{SummarizedExperiment} from a
#' MaxQuant xlsx or txt file. 
#' The function \code{maxQuant} takes as input the path to a .xlsx or .txt file 
#' (MaxQuant output) and additional parameters given to the \code{read.xlsx} 
#' function from the \code{openxlsx} package (e.g. specifying the sheet name 
#' or index by \code{sheet}).
#' 
#' @details 
#' The argument \code{intensity} will specify if the \code{iBAQ} or 
#' \code{LFQ} values are taken. 
#' 
#' The argument \code{type} will specify if the data is loaded from \code{txt} 
#' or \code{xlsx} files.
#'  
#' @param file \code{character}
#' @param intensity \code{character}, either \code{"iBAQ"} or \code{"LFQ"}
#' @param sheet \code{character} or \code{numeric}, the name or index of the 
#' sheet to read data from
#' @param type \code{character}, either \code{"txt"} or \code{"xlsx"}
#' @param ... additional parameters given to \code{read.xlsx} (for 
#' \code{type = "xlsx"})
#'
#' @examples
#' file <- "path/to/maxQuant/object.txt"
#' \donttest{maxQuant(file = file, intensity = "iBAQ", type = "txt")}
#'
#' @return 
#' \code{SummarizedExperiment} object
#'
#' @export
#' 
#' @importFrom openxlsx read.xlsx
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment SummarizedExperiment
maxQuant <- function(file, intensity = c("iBAQ", "LFQ"), sheet, 
    type = c("txt", "xlsx"), ...) {
    
    intensity <- match.arg(intensity)
    type <- match.arg(type)
    
    if (type == "xlsx")
        f <- openxlsx::read.xlsx(file, sheet = sheet, ...)
    if (type == "txt")
        f <- utils::read.table(file, sep = "\t", dec = ".", header = TRUE)
    
    ## names of proteins is in the first col, assign and remove the first col
    rownames(f) <- f[, 1]
    f <- f[, -1]
    
    ## find the columns that contain the features
    cols <- colnames(f)
    inds_samp <- grep(pattern = intensity, cols)
    cols_samp <- cols[inds_samp]

    ## remove the column that only contains "intensity"
    inds_samp <- inds_samp[cols_samp != intensity]
    cols_samp <- cols_samp[cols_samp != intensity]
    
    ## make all characters lower-case for colnames(.f) to grep small differences
    ## in orthography of colnames
    .cols <- tolower(cols)
    .f <- f
    colnames(.f) <- .cols
    
    ## create rowData
    rD <- data.frame(feature = rownames(.f))
    if ("best.ms.ms" %in% .cols) rD$best_MS_MS <- .f[, "best.ms.ms"]
    if ("charges" %in% .cols) rD$charges <- .f[, "charges"]
    if ("count" %in% .cols) rD$count <- .f[, "count"]
    if ("evidence.ids" %in% .cols)
        rD$evidence_ids <- .f[, "evidence.ids"]
    if ("fasta.headers" %in% .cols) rD$fasta_header <- .f[, "fasta.headers"]
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
    cD <- data.frame(name = cols_samp)
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
#' file <- "path/to/spectronaut/object"
#' \donttest{spectronaut(file = file, sheetIntensitities = 1, 
#'     sheetAnnotation = 2, ...)}
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
    cD$name <- samps
    rownames(cD) <- cD[["name"]]

    ## create assay, set values of 0 to NA
    a <- as.matrix(a)
    mode(a) <- "numeric"
    a[a == 0] <- NA

    ## create SummarizedExperiment
    SummarizedExperiment::SummarizedExperiment(assays = a, 
        rowData = rD, colData = cD)
}
