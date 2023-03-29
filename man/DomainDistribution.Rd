\name{DomainDistribution}
\alias{DomainDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DomainDistribution(hclust_result_df, datasets_label_combine, domain_num = ncol(hclust_result_df))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hclust_result_df}{
%%     ~~Describe \code{hclust_result_df} here~~
}
  \item{datasets_label_combine}{
%%     ~~Describe \code{datasets_label_combine} here~~
}
  \item{domain_num}{
%%     ~~Describe \code{domain_num} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (hclust_result_df, datasets_label_combine, domain_num = ncol(hclust_result_df)) 
{
    domain_result_each_list <- list()
    SO_tree <- hclust_result_df[, domain_num]
    names(SO_tree) <- row.names(hclust_result_df)
    domain_names <- names(table(SO_tree))
    datasets_cellNum <- table(datasets_label_combine)
    datasets_num <- length(datasets_cellNum)
    datasets_name <- names(datasets_cellNum)
    proportion_dataset_domain <- matrix(0, domain_num, datasets_num)
    Term_table <- table(SO_tree)
    row.names(proportion_dataset_domain) <- names(Term_table)
    colnames(proportion_dataset_domain) <- names(datasets_cellNum)
    for (j in 1:length(Term_table)) {
        Datasets_term_cellNum <- table(datasets_label_combine[names(SO_tree[SO_tree == 
            names(Term_table)[j]])])
        proportion_dataset_domain[j, ] <- Datasets_term_cellNum/datasets_cellNum
    }
    for (i in 1:datasets_num) {
        domain_result_each_list[[i]] <- hclust_result_df[names(datasets_label_combine[datasets_label_combine == 
            datasets_name[i]]), ]
    }
    names(domain_result_each_list) <- datasets_name
    return(list(domain_distribution = proportion_dataset_domain, 
        domain_result_each_list = domain_result_each_list))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
