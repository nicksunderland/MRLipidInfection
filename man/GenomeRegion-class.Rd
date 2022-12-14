% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome_region_class.R
\docType{class}
\name{GenomeRegion-class}
\alias{GenomeRegion-class}
\alias{GenomeRegion}
\title{GenomeRegion class}
\value{
a GenomeRegion object
}
\description{
This class defines an area of the genome and handles various parameters passed to the LD clumping
algorithms. This is useful when working with SNPs surrounding multiple gene targets of interest.
To perform normal genome wide SNP clumping simply create a GenomeRegion object with
'genome_wide=TRUE'.
}
\section{Slots}{

\describe{
\item{\code{name}}{character. The name of the gene or area of the genome referred to.}

\item{\code{assembly}}{character. The assembly information, e.g. GRCh37}

\item{\code{chromosome}}{numeric. The chromosome the gene is on.}

\item{\code{start}}{numeric. The start of the gene (number of bases from start of chromosome).}

\item{\code{end}}{numeric. The end of the gene (number of bases from start of chromosome).}

\item{\code{cis_tol}}{when defining SNPs in relation to a specific a gene, set a tolerance of how many
bases either side of the gene to include the SNP}

\item{\code{genome_wide}}{this is usually FALSE; if TRUE this is treated as a object placeholder and will
refer to the 'whole genome'}

\item{\code{clump_kb}}{numeric. Clumping window, default is 10000.}

\item{\code{clump_r2}}{numeric. Clumping r2 cutoff.}

\item{\code{clump_p1}}{numeric. Clumping sig level for index SNPs, default is 1.}

\item{\code{clump_p2}}{numeric. Clumping sig level for secondary SNPs, default is 1.}

\item{\code{pop}}{character. Super-population to use as reference panel. Default = "EUR". Options are
EUR, SAS, EAS, AFR, AMR. 'legacy' also available - which is a previously used version of the
EUR panel with a slightly different set of markers}
}}

