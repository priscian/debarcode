CATALYST_debarcode <- function(
  x, # flowCore::flowFrame
  barcoding_key,
  assay = "exprs"
)
{
  if (is.null(barcoding_key))
    return (rep(basename(flowCore::description(x)$FILENAME), NROW(x)))

  ## TODO: Evaluate 'barcoding_key' if it's an expression...?
  bc_key <- barcoding_key

  row_data <- flowCore::pData(flowCore::parameters(x)) %>%
    dplyr::select(name, desc) %>%
    dplyr::rename(channel_name = "name", marker_name = "desc") %>%
    dplyr::mutate(
      dplyr::across(.fns = as.character),
      marker_name = dplyr::case_when(is.na(marker_name) ~ channel_name, TRUE ~ marker_name)
    ) %>%
    textshape::column_to_rownames(loc = "marker_name")

  e <- flowCore::exprs(x)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = t(e) %>% `rownames<-`(as.vector(rownames(.)))),
    rowData = row_data)

  sce <- CATALYST::assignPrelim(sce, bc_key = bc_key, assay = assay, verbose = TRUE)
  sce <- CATALYST::estCutoffs(sce)
  sce <- CATALYST::applyCutoffs(sce)
  #table(sce$bc_id)

  sce$bc_id
}


#' @export
debarcode <- function(
  input_path,
  key,
  output_dir = ".", create_output_dir = TRUE,
  b = 1/8, # asinh transformation parameter: FCM = 1/150, CyTOF = 1/8 (v. MetaCyto vignette)
  excluded_transform_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  output_transformed_events = FALSE,
  outfile_prefix = "", outfile_suffix = "", # also possibly 'NULL'
  filename_sample_sep = "-",
  ids_only = FALSE, # If TRUE, return only vector of sample IDs
  ...
)
{
  ff0 <- flowCore::read.FCS(input_path, transformation = FALSE, truncate_max_range = FALSE)

  temp <- flowCore::colnames(ff0); transformChannels <- temp[stringr::str_detect(temp, excluded_transform_channels_re, negate = TRUE)]
  asinhTrans <- flowCore::arcsinhTransform(transformationId = "flowpipe-transformation", a = 1, b = b, c = 0)
  transList <- flowCore::transformList(transformChannels, asinhTrans)
  tff <- flowCore::transform(ff0, transList)

  ## Debarcode
  sample_id <- CATALYST_debarcode(tff, key)
  table(sample_id) %>% print

  if (ids_only)
    return (sample_id)

  ## Split barcoded samples & remove unassigned "0" events
  ff <- ff0; if (output_transformed_events) ff <- tff
  ffs <- flowCore::split(ff, sample_id, flowSet = FALSE)# %>% `[[<-`("0", NULL)

  if (create_output_dir && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  fcsFilePaths <- sapply(names(ffs),
    function(i)
    {
      fcsFileName <-
        sprintf(paste0("%s", basename(input_path), "%s%s.fcs"), outfile_prefix, paste0(filename_sample_sep, i), outfile_suffix)
      fcsFilePath <- paste(output_dir, fcsFileName, sep = "/")
      ff <- ffs[[i]]
      flowCore::keyword(ff) <- list(`$FIL` = basename(fcsFilePath))
      cat(sprintf("Saving FCS file %s...", fcsFileName)); utils::flush.console()
      flowCore::write.FCS(ff, filename = fcsFilePath)
      cat(". Done.", fill = TRUE)

      fcsFilePath
    }, simplify = TRUE)

  structure(fcsFilePaths, transformed = output_transformed_events, sample_id = sample_id)
}


### Utility functions

## Simplify 'read.table()' from 'textConnection()'s; like SAS's 'CARDS' statement.
#' @export
cards <- function(
  x,
  header = TRUE,
  as.is = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  ...
)
{
  tab <- read.table(text = x, header = header, as.is = as.is, check.names = check.names, stringsAsFactors = stringsAsFactors, ...)

  return (tab)
}
