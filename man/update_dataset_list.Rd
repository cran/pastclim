% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_dataset_list.R
\name{update_dataset_list}
\alias{update_dataset_list}
\title{Update the dataset list}
\usage{
update_dataset_list(on_cran = FALSE)
}
\arguments{
\item{on_cran}{boolean to make this function run on ci tests using tempdir}
}
\value{
TRUE if the dataset was updated
}
\description{
If a newer dataset list (which includes all the information about the files
storing the data for pastclim), download it and start using it as
'dataset_list_included.csv' in
\code{tools::R_user_dir("pastclim","config")}. If the latter is present, the last
column, named 'dataset_list_v', provides the version of this table, and the
most advanced table is used.
}
