## function biocrates
test_that("biocrates", {
    expect_is(biocrates, "function")
    expect_error(biocrates(file = ""), "File does not exist")
    
    file <- system.file("extdata", "biocrates_test_file.xlsx", 
        package = "MatrixQCvisUtils")
    se <- biocrates(file = file, sheet = 1)
    cols_samples <- c("sample_1", "sample_2", "sample_3", "sample_4", 
        "sample_5", "sample_6", "sample_7", "sample_8", "sample_9", "sample_10")
    expect_equal(colnames(se), cols_samples)
    expect_equal(rownames(se)[1:5], 
        c("C0", "C2", "C3", "C3.DC..C4.OH.", "C3.OH"))
    expect_equal(dim(se), c(630, 10))
    
    ## check colData
    expect_equal(rownames(colData(se)), cols_samples)
    expect_equal(se$name, cols_samples)
    expect_equal(se$name_original, cols_samples)
    expect_equal(se$Plate.Bar.Code, rep("plate_bar_code_1", 10))
    expect_equal(se$Sample.Bar.Code, c("sample_bar_code_1", "sample_bar_code_2",
        "sample_bar_code_3", "sample_bar_code_4", "sample_bar_code_5", 
        "sample_bar_code_6", "sample_bar_code_7", "sample_bar_code_8", 
        "sample_bar_code_9", "sample_bar_code_10"))
    expect_equal(se$Sample.Identification, cols_samples)
    expect_equal(se$Submission.Name, cols_samples)
    expect_equal(se$Species, rep("human", 10))
    expect_equal(se$Material, rep("cells", 10))
    expect_equal(se$Cell.Number, c("4000000", "4000000", "2000000", "2000000",
        "2000000", "5000000", "5000000", "5000000", "5000000", "2000000"))
    expect_equal(se$Cell.Extraction.Volume..µl., rep("60", 10))
    expect_equal(se$Org..Info, rep("human", 10))
    expect_equal(se$Measurement.Time, rep("2024.03.08", 10))
    
    ## check rowData
    expect_equal(rownames(rowData(se)), rownames(se))
    expect_equal(rowData(se)$feature[1:5], 
        c("C0", "C2", "C3", "C3.DC..C4.OH.", "C3.OH"))
    expect_equal(rowData(se)$feature_original[1:5], 
                 c("C0", "C2", "C3", "C3-DC (C4-OH)", "C3-OH"))
    expect_equal(rowData(se)$class[1:5], rep("Acylcarnitines", 5))
    expect_equal(unique(rowData(se)$class), c("Acylcarnitines", "Alkaloids",
        "Amine Oxides", "Aminoacids", "Aminoacids Related", "Bile Acids",
        "Biogenic Amines", "Carboxylic Acids", "Ceramides",
        "Cholesterol Esters", "Cresols", "Diacylglycerols", "Dihydroceramides",
        "Fatty Acids", "Glycosylceramides", "Hormones", "Indoles Derivatives",
        "Nucleobases Related", "Phosphatidylcholines", "Sphingomyelins",     
        "Sugars", "Triacylglycerols", "Vitamins & Cofactors"))
    expect_equal(rowData(se)$HMDB_ids[1:5], 
        c( "HMDB0000062", "HMDB0000201", "HMDB0000824", "HMDB0002095", 
            "HMDB0013125"))

    ## check assay
    expect_equal(colnames(assay(se)), cols_samples)
    expect_equal(rownames(assay(se))[1:5], 
        c("C0", "C2", "C3", "C3.DC..C4.OH.", "C3.OH"))
    expect_equal(sum(is.na(assay(se))), 0)
    expect_equal(sum(assay(se)), 360713.422, tolerance = 1e-03)
})

## function metaboscape
test_that("metaboscape", {
    expect_is(metaboscape, "function")
    expect_error(metaboscape(file = ""), "File does not exist")
    
    file <- system.file("extdata", "metaboscape_test_file.xlsx", 
        package = "MatrixQCvisUtils")
    se <- metaboscape(file = file, sheet = 1)
    cols_samples <- c("sample_1", "sample_2", "sample_3", "sample_4", 
        "sample_5")
    expect_equal(colnames(se), cols_samples)
    expect_equal(rownames(se)[1:5], 
        c("feature_3", "feature_4", "feature_5", "feature_6", "feature_7"))
    expect_equal(dim(se), c(10, 5))
    
    ## check colData
    expect_equal(rownames(colData(se)), cols_samples)
    expect_equal(se$name, cols_samples)
    expect_equal(se$name_original, cols_samples)
    
    ## check rowData
    expect_equal(rownames(rowData(se)), rownames(se))
    expect_equal(rowData(se)$feature[1:5], 
        c("feature_3", "feature_4", "feature_5", "feature_6", "feature_7"))
    expect_equal(rowData(se)$rt_min[1:5], c(1, 2, 3, 4, 5))
    expect_equal(rowData(se)$ccs_a2[1:5], c(281, 198, 248.2, 263.8, 241.4))
    expect_equal(rowData(se)$deltaccs_percent, rep(NaN, 10))
    expect_equal(rowData(se)$mz[1:5], 
        c(684.5466, 608.8770, 968.7963, 1104.7707, 900.8089), tolerance = 1e-03)
    expect_equal(rowData(se)$molecular_mass[1:5], 
        c(683.5393, 607.8698, 967.7890, 1103.7634, 899.8016), tolerance = 1e-03)
    expect_equal(unique(rowData(se)$ions), "[M+H]+")
    expect_equal(rowData(se)$msms, c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        FALSE, TRUE, TRUE))
    expect_equal(unique(rowData(se)$qc_rsd_percent), NaN)
    
    ## check assay
    expect_equal(colnames(assay(se)), cols_samples)
    expect_equal(rownames(assay(se))[1:5], 
        c("feature_3", "feature_4", "feature_5", "feature_6", "feature_7"))
    expect_equal(sum(is.na(assay(se))), 6)
    expect_equal(sum(assay(se), na.rm = TRUE), 62593)
})

## function maxquant
test_that("maxquant", {
    expect_is(maxquant, "function")
    expect_error(maxquant(file = "", intensity = "foo"), "should be one of")
    suppressWarnings(expect_error(maxquant(file = "", intensity = "iBAQ", 
        type = "txt"), "no lines available"))
    expect_error(maxquant(file = "", intensity = "iBAQ", type = "xlsx"), 
        "File does not exist")
    suppressWarnings(expect_error(maxquant(file = "", intensity = "LFQ", 
        type = "txt"), "no lines available"))
    expect_error(maxquant(file = "", intensity = "LFQ", type = "xlsx"), 
        "File does not exist")
    expect_error(maxquant(file = "", intensity = "LFQ", type = "foo"),
        "should be one of")
    
    file <- system.file("extdata", "maxquant_test_file.xlsx", 
        package = "MatrixQCvisUtils")
    se <- maxquant(file = file, type = "xlsx", intensity = "LFQ", sheet = 1)
    cols_samples <- c("sample_1_LFQ", "sample_2_LFQ", "sample_3_LFQ", 
        "sample_4_LFQ", "sample_5_LFQ")
    expect_equal(colnames(se), cols_samples)
    expect_equal(rownames(se)[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(dim(se), c(10, 5))
    
    ## check colData
    expect_equal(rownames(colData(se)), cols_samples)
    expect_equal(se$name, cols_samples)
    expect_equal(se$name_original, cols_samples)
    expect_equal(se$name_cut, cols_samples)
    
    ## check rowData
    expect_equal(rownames(rowData(se)), rownames(se))
    expect_equal(rowData(se)$feature[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(unique(rowData(se)$count), "3")
    expect_equal(rowData(se)$gene_name[1:5], 
        c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5"))
    expect_equal(rowData(se)$majority_protein_ids[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(rowData(se)$number_of_proteins[1:5], 
        c("1", "1", "2", "2", "1"))
    expect_equal(rowData(se)$peptides[1:5], 
        c("1", "12", "6", "7", "5"))
    expect_equal(rowData(se)$protein_names[1:5], 
        c("protein_1_HUMAN", "protein_2_HUMAN", "protein_3_HUMAN", 
            "protein_4_HUMAN", "protein_5_HUMAN"))
    
    ## check assay
    expect_equal(colnames(assay(se)), cols_samples)
    expect_equal(rownames(assay(se))[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(sum(is.na(assay(se))), 6)
    expect_equal(sum(assay(se), na.rm = TRUE), 62593)
})

## function diann
test_that("diann", {
    expect_is(diann, "function")
    expect_error(suppressWarnings(diann(file = "")), "no lines available")
    
    file <- system.file("extdata", "diann_test_file.tsv", 
        package = "MatrixQCvisUtils")
    se <- diann(file = file)
    cols_samples <- c("sample_1", "sample_2", "sample_3", "sample_4", 
        "sample_5")
    expect_equal(colnames(se), cols_samples)
    expect_equal(rownames(se)[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(dim(se), c(10, 5))
    
    ## check colData
    expect_equal(rownames(colData(se)), cols_samples)
    expect_equal(se$name, cols_samples)
    expect_equal(se$name_original, cols_samples)
    expect_equal(se$name_cut, cols_samples)
    
    ## check rowData
    expect_equal(rownames(rowData(se)), rownames(se))
    expect_equal(rowData(se)$feature[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(rowData(se)$first_protein_description[1:5], 
        c("description_1", "description_2", "description_3", "description_4", 
            "description_5"))
    expect_equal(rowData(se)$genes[1:5], 
        c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5"))
    expect_equal(rowData(se)$protein_ids[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(rowData(se)$protein_names[1:5], 
        c("protein_1_HUMAN", "protein_2_HUMAN", "protein_3_HUMAN", 
            "protein_4_HUMAN", "protein_5_HUMAN"))
    
    ## check assay
    expect_equal(colnames(assay(se)), cols_samples)
    expect_equal(rownames(assay(se))[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(sum(is.na(assay(se))), 6)
    expect_equal(sum(assay(se), na.rm = TRUE), 62593)
})

## function spectronaut
test_that("spectronaut", {
    expect_is(spectronaut, "function")
    expect_error(spectronaut(file = ""), "File does not exist")
    
    file <- system.file("extdata", "spectronaut_test_file.xlsx", 
        package = "MatrixQCvisUtils")
    se <- suppressWarnings(spectronaut(file = file))
    cols_samples <- c("sample_1", "sample_2", "sample_3", "sample_4", 
        "sample_5")
    expect_equal(colnames(se), cols_samples)
    expect_equal(rownames(se)[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(dim(se), c(10, 5))
    
    ## check colData
    expect_equal(rownames(colData(se)), cols_samples)
    expect_equal(se$name, cols_samples)
    expect_equal(se$Condition, c("a", "a", "b", "b", "c"))
    expect_equal(se$Type, c("plasma", "plasma", "plasma", "plasma", "Plasma"))
    expect_equal(se$Time, c(0, 2, 0, 2, 0))
    expect_equal(se$Replicate, c(1, 1, 2, 2, 3))
    expect_equal(se$TechRep, c("T1", "T1", "T1", "T2", "T1"))
    expect_equal(se$Precursors, c(3307, 3112, 81, 3053, 3064))
    expect_equal(se$ModifiedSequences, c(2917, 2707, 80, 2673, 2682))
    expect_equal(se$Peptides, c(2694, 2476, 75, 2466, 2469))
    expect_equal(se$ProteinGroups, c(487, 424, 54, 461, 439))
    expect_equal(se$Proteins, c(729, 591, 71, 667, 631))
    expect_equal(se$Filter, c(NA, NA, "Exclude", NA, NA))
    expect_equal(se$name_original, c("sample_1", "sample_2", "sample_3", 
        "sample_4", "sample_5"))

    ## check rowData
    expect_equal(rownames(rowData(se)), rownames(se))
    expect_equal(rowData(se)$feature[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(rowData(se)$PG_Genes[1:5], 
        c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5"))
    
    ## check assay
    expect_equal(colnames(assay(se)), cols_samples)
    expect_equal(rownames(assay(se))[1:5], 
        c("protein_1", "protein_2", "protein_3", "protein_4", "protein_5"))
    expect_equal(sum(is.na(assay(se))), 6)
    expect_equal(sum(assay(se), na.rm = TRUE), 62593)
})

