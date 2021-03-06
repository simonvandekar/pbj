% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/niftiImage.R
\name{image.niftiImage}
\alias{image.niftiImage}
\title{Create images of a NiFTI object}
\usage{
\method{image}{niftiImage}(
  x,
  bgimg = NULL,
  thresh = 0,
  index = NULL,
  col = gray(0:64/64),
  colpos = pbj:::redyellow(64),
  colneg = pbj:::bluecyan(64),
  plane = c("axial", "coronal", "sagittal"),
  xlab = "",
  ylab = "",
  axes = FALSE,
  oma = rep(0, 4),
  mar = rep(0, 4),
  bg = "black",
  ...
)
}
\arguments{
\item{x}{the Nifti object to display images of}

\item{bgimg}{background image to use.}

\item{thresh}{A lower (an optionally upper length 2 vector) threshold to apply to the image, defaults to 0}

\item{index}{Any selected image planes. defaults to NULL}

\item{col}{a vector of colors to use for scaled intensities defaults to a grey scale.}

\item{colpos}{a vector of colors to use for positive values.}

\item{colneg}{a vector of colors to use for negative values.}

\item{plane}{the plane to display, can be axial, coronal or sagittal.}

\item{xlab}{a title for the x axis.}

\item{ylab}{a title for the y axis.}

\item{axes}{display axes, defaults to false}

\item{oma}{A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text. Default to 0.}

\item{mar}{A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. Defaults to 0.}

\item{bg}{background color, defaults to black.}

\item{...}{additional arguments passed to par}
}
\description{
modified from oro.nifti:::image.nifti
}
