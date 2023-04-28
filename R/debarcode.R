CATALYST_debarcode <- function(
  x, # flowCore::flowFrame
  barcoding_key,
  assay = "exprs"
)
{
  if (is.null(barcoding_key))
    return (rep(basename(flowCore::description(x)$FILENAME), NROW(x)))

  bc_key <- barcoding_key

  keystone::poly_eval(barcoding_key$process_key)

  sce <- as_sce(x)

  sce <- CATALYST::assignPrelim(sce, bc_key = bc_key, assay = assay, verbose = TRUE)
  sce <- CATALYST::estCutoffs(sce)
  sce <- CATALYST::applyCutoffs(sce)
  #table(sce$bc_id)

  sce$bc_id
}


debarcode_single <- function(
  input_path,
  key,
  output_dir = ".", create_output_dir = TRUE,
  b = 1/8, # asinh transformation parameter: FCM = 1/150, CyTOF = 1/8 (v. MetaCyto vignette)
  excluded_transform_channels_re = stringr::regex("time|event_length", ignore_case = TRUE),
  output_transformed_events = FALSE,
  outfile_prefix = "", outfile_suffix = "", # also possibly 'NULL'
  filename_sample_sep = "-",
  ids_only = FALSE, # If TRUE, return only vector of sample IDs
  ignore_unassigned_files = TRUE,
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

      if (ignore_unassigned_files && i == "0")
        return (NULL)

      fcsFilePath
    }, simplify = FALSE) %>% unlist %>% as.vector

  ret_val <- structure(fcsFilePaths, transformed = output_transformed_events, sample_id = sample_id)

  ret_val
}


#' @export
## Good idea, but I need more control:
#debarcode <- Vectorize(FUN = debarcode_single, vectorize.args = c("input_path", "key", "output_dir"))
debarcode <- function(
  input_path, # Any vector of FCS files
  key, # Named list of barcoding keys w/ names exactly same as 'input_path' elements needing deconvolution
  output_dir = ".", # Vector of output directories, recycled to 'length(key)'
  ... # Arguments passed on to 'debarcode_single()'
)
{
  outputDirs <- structure(rep(output_dir, length.out = length(key)), .Names = names(key))

  pp <- keystone::psapply(input_path,
    function(a)
    {
      if (is.null(key[[a]]))
        return (a)

      debarcode_single(input_path = a, key = key[[a]], output_dir = outputDirs[a])
    }, simplify = FALSE)

  ret_val <- structure(pp %>% unlist %>% as.vector, details = pp)

  ret_val
}
