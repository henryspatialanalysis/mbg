#!/usr/bin/env Rscript
library(rmarkdown)
library(argparse)

# Function to rebuild a vignette
rebuild_vignette <- function(raw_rmd, out_file = NULL, out_dir = NULL){
  # The default output filename is the raw RMD basename without a leading underscore
  if(is.null(out_file)) out_file <- gsub('^_', '', basename(raw_rmd))
  # The default output directory is one directory above where the raw RMD file is stored
  if(is.null(out_dir)) out_dir <- dirname(dirname(raw_rmd))
  # If the rendered RMD already exists, remove it
  full_out_dir <- normalizePath(out_dir)
  full_out_path <- file.path(full_out_dir, out_file)
  if(file.exists(full_out_path)) file.remove(full_out_path)
  # Render the document
  rmarkdown::render(
    raw_rmd,
    output_file = out_file,
    output_dir = out_dir,
    output_format = rmarkdown::md_document(
      variant = 'gfm',
      preserve_yaml = TRUE,
      ext = '.Rmd'
    )
  )
  # Replace absolute paths with relative paths
  rendered_text <- readLines(full_out_path)
  if(!endsWith(full_out_dir, '/')) full_out_dir <- paste0(full_out_dir, '/')
  fixed_text <- gsub(full_out_dir, '', rendered_text)
  writeLines(fixed_text, full_out_path)
  invisible(TRUE)
}

# Parse the command line arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--raw_rmd", type = 'character')
parser$add_argument("--out_file", type = 'character', default = NULL)
parser$add_argument("--out_dir", type = 'character', default = NULL)
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Run rebuild_vignette() with the passed arguments
do.call(what = rebuild_vignette, args = args)
