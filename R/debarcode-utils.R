#' @export
as_sce <- function(x, ...)
  UseMethod("as_sce")


#' @export
as_sce.flowFrame <- function(x, ...)
{
  row_data <- flowCore::pData(flowCore::parameters(x)) %>%
    dplyr::select(name, desc) %>%
    dplyr::rename(channel_name = "name", marker_name = "desc") %>%
    dplyr::mutate(
      dplyr::across(.cols = everything(), .fns = as.character),
      marker_name = dplyr::case_when(is.na(marker_name) ~ channel_name, TRUE ~ marker_name)
    ) %>%
    textshape::column_to_rownames(loc = "marker_name")

  e <- flowCore::exprs(x)
  sce <- SingleCellExperiment::SingleCellExperiment(
    #assays = list(exprs = t(e) %>% `rownames<-`(as.vector(rownames(.)))), # Unnecessary?
    assays = list(exprs = t(e) %>% `rownames<-`(NULL)),
    rowData = row_data)

  sce
}


## Drop-in update of 'premessa::concatenate_fcs_files()' to allow args to 'flowCore::read.FCS()'
#' @export
concatenate_fcs_files <- function(
  files.list,
  output.file = NULL,
  create_output_dir = TRUE,
  read.FCS... = list(), # Arguments to 'flowCore::read.FCS()'
  ... # Unused
)
{
  if (create_output_dir && !dir.exists(dirname(output.file)))
    dir.create(dirname(output.file), recursive = TRUE)

  read.FCSArgs <- list(
    filename = NULL,
    transformation = FALSE,
    truncate_max_range = FALSE
  )

  #m <- lapply(files.list, flowCore::read.FCS, ...)
  m <- lapply(files.list,
    function(a) {
      read_FCSArgs <- read.FCSArgs
      read_FCSArgs$filename = a
      read_FCSArgs <- utils::modifyList(read_FCSArgs, read.FCS..., keep.null = TRUE)

      do.call(flowCore::read.FCS, read_FCSArgs)
    })

  flow.frame <- m[[1]]
  m <- lapply(m, function(x) {
    flowCore::exprs(x)
  })
  m <- do.call(rbind, m)

  ret <- premessa::as_flowFrame(m, flow.frame)

  if (!is.null(output.file))
    premessa::write_flowFrame(ret, output.file)
  else return (ret)
}
