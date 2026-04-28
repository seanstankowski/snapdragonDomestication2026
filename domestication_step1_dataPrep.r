################################################## data preprocessing scripts #######################################################

# Set working directory to the location of the data file
setwd("/Users/sean/projects/Snapdragons/pipeline_out")

# read in data
data <- read.csv("domestication_genotypes.csv", header = TRUE)

# identify valid line-style IDs once
has_line_id <- grepl("^L\\d{1,2}-\\d{1,2}$", data$PlantID)

# Add domesticated status based on PlantID
data$domesticated <- ifelse(has_line_id, "domesticated", "wild")

# Add a column indicating the line number
# If no valid line ID is present, set to NA
data$line <- NA_character_
data$line[has_line_id] <- sub("^L(\\d{1,2})-\\d{1,2}$", "L\\1", data$PlantID[has_line_id])

# Add a line replicate number
data$replicate_number <- NA_integer_
data$replicate_number[has_line_id] <- as.integer(
  sub("^L\\d{1,2}-(\\d{1,2})$", "\\1", data$PlantID[has_line_id])
)

# Some line-numbered samples are actually wild-derived comparison lines
wild_lines <- c("L8", "L9", "L15", "L18", "L24", "L38")
data$domesticated[data$line %in% wild_lines] <- "wild"

# add a species column indicating which species each sample belongs to
data$species <- "majus"
data$species[data$line == "L40"] <- "siculum"
data$species[data$line == "L27"] <- "molle"
data$species[data$line == "L5"]  <- "majus_ssp_tortuosum"

# Calculate proportion of missing data for each individual 
marker_cols <- 10:115
geno_mat <- as.matrix(data[, marker_cols])

data$prop_missing <- rowSums(geno_mat == -9 | geno_mat == -10) / ncol(geno_mat)

# Drop individuals with more than 20% missing data
data_filtered <- data[data$prop_missing <= 0.2 | is.na(data$prop_missing), ]

# remove individuals without valid GPS unless they have a valid line ID
data_filtered <- data_filtered[
  !is.na(data_filtered$line) |
    (!is.na(data_filtered$CorrectedLatitude) & !is.na(data_filtered$CorrectedLongitude)),
]

############### bringing in marker information

marker_info <- read.csv("SNPlociListUpdated_Dec_2019.csv", header = TRUE)

# keep only the rows in the marker info that correspond to the markers in the genotype data
marker_cols <- 10:115
geno_marker_names <- names(data_filtered)[marker_cols]

marker_info_filtered <- marker_info[marker_info$LocusName %in% geno_marker_names, ]

# reorder marker info to match the order of markers in the genotype data
marker_info_filtered <- marker_info_filtered[
  match(geno_marker_names, marker_info_filtered$LocusName),
]

# sanity check
stopifnot(all(marker_info_filtered$LocusName == geno_marker_names))

### designate marker effects

# start with NA
marker_info_filtered$effect <- NA

# marker sets
colour_markers <- c(
  "ros_assembly_473914",
  "ros_assembly_620992",
  "ros_assembly_653015",
  "ros_assembly_670530",
  "ros_assembly_737420",
  "ros_assembly_744403",
  "ros_assembly_748981",
  "ros_assembly_758578",
  "ros_assembly_715001",
  "s154_504353",
  "s2338_45429",
  "s829_8371",
  "s829_204463",
  "s316_93292",
  "s316_257789",
  "s444_38909",
  "s992_223854",
  "s148_425797",
  "s1140_224946",
  "s155_1201194",
  "ros_assembly_849332",
  "ros_assembly_543443",
  "ros_assembly_567004",
  "ros_assembly_575837",
  "ros_assembly_576271",
  "s91_78256",
  "s91_122561",
  "s91_181717",
  "s901_155439",
  "s1187_290152"
)

chloroplast_markers <- c(
  "scpDNA_seq_58467",
  "scpDNA_seq_21466"
)

neutral_markers <- c(
  "s816_1076784",
  "s200_39182",
  "s1152_21836",
  "s7_443011",
  "s1380_129448",
  "s1048_161750",
  "s320_60828",
  "s863_381450",
  "s1492_57885",
  "s787_264617",
  "s1135_562200",
  "s44_51095",
  "s1367_16317",
  "s899_432971",
  "s36_910340",
  "s282_483674",
  "s1668_18672",
  "s801_62707",
  "s334_46968",
  "s1099_632387",
  "s314_526456",
  "s720_870524",
  "s531_811478",
  "s260_949098",
  "s1180_371625",
  "s1056_246800",
  "s930_115738",
  "s382_741148",
  "s705_151285",
  "s235_82833",
  "s699_457945",
  "s549_1270640",
  "s751_467474",
  "s85_398904",
  "s743_443660",
  "s692_198224",
  "s16_391553",
  "s481_103209",
  "s1258_68093",
  "s682_573512",
  "s2786_119740",
  "s211_597773",
  "s822_183155",
  "s368_322607",
  "s89_95031",
  "s149_225968",
  "s494_499812",
  "s804_404293",
  "s88_362224",
  "s763_365058",
  "s442_292702",
  "s112_351672",
  "s758_80943",
  "s248_368997",
  "s427_390429",
  "s612_303555",
  "s1387_197749",
  "s1195_856423",
  "s1088_730874",
  "s50_706008",
  "s83_819752",
  "s432_102321",
  "s807_258952",
  "s121_351821",
  "s338_787544",
  "s445_464618",
  "s55_222383",
  "s24_510325",
  "s679_138960",
  "s866_450394",
  "s269_957088"
)

morphology_markers <- c(
  "s1346_227941",
  "s694_100451",
  "s699_530828"
)

# assign categories
marker_info_filtered$effect[marker_info_filtered$LocusName %in% colour_markers] <- "colour"
marker_info_filtered$effect[marker_info_filtered$LocusName %in% chloroplast_markers] <- "chloroplast"
marker_info_filtered$effect[marker_info_filtered$LocusName %in% neutral_markers] <- "neutral"
marker_info_filtered$effect[marker_info_filtered$LocusName %in% morphology_markers] <- "morphology"

# now add a column to designate gene names

# initialise as NA
marker_info_filtered$gene_id <- NA

# define groups
def_markers <- c("s1346_227941")
dich_markers <- c("s694_100451")
el_markers <- c("ros_assembly_715001")

fla_markers <- c(
  "s154_504353",
  "s2338_45429",
  "s829_8371",
  "s829_204463",
  "s316_93292",
  "s316_257789",
  "s444_38909",
  "s992_223854",
  "s148_425797",
  "s1140_224946",
  "s155_1201194"
)

mixta_markers <- c("s699_530828")

ros_markers <- c(
  "ros_assembly_473914",
  "ros_assembly_620992",
  "ros_assembly_653015",
  "ros_assembly_670530",
  "ros_assembly_737420",
  "ros_assembly_744403",
  "ros_assembly_748981",
  "ros_assembly_758578",
  "ros_assembly_849332",
  "ros_assembly_543443",
  "ros_assembly_567004",
  "ros_assembly_575837",
  "ros_assembly_576271"
)

sulf_markers <- c(
  "s91_78256",
  "s91_122561",
  "s91_181717"
)

ven_markers <- c("s901_155439")
cre_markers <- c("s1187_290152")

# assign gene IDs
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% def_markers] <- "def"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% dich_markers] <- "dich"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% el_markers] <- "el"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% fla_markers] <- "fla"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% mixta_markers] <- "mixta"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% ros_markers] <- "ros"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% sulf_markers] <- "sulf"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% ven_markers] <- "ven"
marker_info_filtered$gene_id[marker_info_filtered$LocusName %in% cre_markers] <- "cre"

########## now flower colour

colour_data <- read.csv("20250418_SnapPics_Apollonia.csv", header = TRUE)

## bring the relevant columns into the main data frame data_filtered

# make sure replicate_number is the same type in both data frames
data_filtered$replicate_number <- as.character(data_filtered$replicate_number)
colour_data$replicate_number   <- as.character(colour_data$replicate_number)

# keep only the columns needed from colour_data
colour_sub <- colour_data[, c("line", "replicate_number", "Red", "Yellow", "Peloric", "mean_hue_all_circles")]

# match rows in colour_data to data_filtered by line + replicate_number
match_idx <- match(
  paste(data_filtered$line, data_filtered$replicate_number),
  paste(colour_sub$line, colour_sub$replicate_number)
)

# fill Red and Yellow where matches exist; note that 5 of the individuals in colour data were not sequenced and one was dropped because it had > 0.2 missingness
data_filtered$Red[!is.na(match_idx)]    <- colour_sub$Red[match_idx[!is.na(match_idx)]]
data_filtered$Yellow[!is.na(match_idx)] <- colour_sub$Yellow[match_idx[!is.na(match_idx)]]

# add Peloric and mean_hue_all_circles columns
data_filtered$Peloric <- NA
data_filtered$mean_hue_all_circles <- NA

# fill new columns for matched rows
data_filtered$Peloric[!is.na(match_idx)] <- colour_sub$Peloric[match_idx[!is.na(match_idx)]]
data_filtered$mean_hue_all_circles[!is.na(match_idx)] <- colour_sub$mean_hue_all_circles[match_idx[!is.na(match_idx)]]

