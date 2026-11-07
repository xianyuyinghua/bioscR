% Documentation for R/group_order_set.R; keep in sync with its opening help comments.
\name{group_order_set}
\alias{group_order_set}
\title{
Set the order of group labels
}
\description{
Returns an order for group labels, either using a supplied order or
placing the first recognised control label before the remaining labels.
}
\usage{
group_order_set(vector, group_levels = NULL)
}
\arguments{
\item{vector}{A vector of group labels, typically a character vector or factor.
Required, with no default. The function applies \code{unique()} before
\code{as.character()}; unused factor levels are not included. Numeric labels
are therefore ordered as character strings. Missing and blank labels are
not explicitly validated or cleaned; see Details.}
\item{group_levels}{An optional vector specifying the desired order. Defaults
to \code{NULL}; any zero-length value selects automatic ordering. For intended
use, supply a character vector containing every distinct character label
in \code{vector}. Matching is case-sensitive. Extra labels and duplicates are
accepted, and a matching vector is returned unchanged, including its type
and attributes. There is no explicit type or uniqueness validation.}
}
\details{
In automatic mode, labels are deduplicated and converted to character.
The function searches for \code{"control"}, \code{"normal"}, or \code{"norm"} without regard
to case. If several labels match, the first one in input order is moved to
the front; the other labels retain their input order. There is no priority
among the three control names. Original spelling and case are retained,
and whitespace is not trimmed.

If no control label matches, \code{sort()} orders the character labels in
ascending order using the current locale and drops \code{NA} values by default.
When a control label matches, missing values instead remain among the
other labels. Blank strings are not removed. Empty input in automatic
mode returns \code{character(0)}.

With a non-empty \code{group_levels}, the function checks only whether every
processed input label occurs in \code{group_levels}, using \code{\%in\%}. It does not
filter extra levels, remove duplicates, or reorder the supplied vector.
If any label is absent, it prints \code{"group_levels Not Match the vector!"}
with ANSI colour escape sequences and returns \code{NULL} invisibly. This is
a console message produced by \code{cat()}, not a warning or an error.

The result specifies label order only. The function does not modify a
factor or metadata column, select observations, or perform comparisons;
it does not define a statistical reference group or comparison direction.

The implementation uses base R functions and the unqualified \code{\%>\%} pipe.
The pipe must be available in the function's evaluation environment, for
example by attaching dplyr before calling. dplyr is declared in the
package Imports, but the current NAMESPACE does not import \code{\%>\%}.
No external scripts, input files, or output files are used. The only
side effect is the console message on an explicit-order mismatch.
}
\value{
In automatic mode, a character vector of group labels in the computed
order. With a matching non-empty \code{group_levels}, that object is returned
unchanged, so its type, attributes, extra labels, and duplicates are
preserved. On an explicit-order mismatch, returns \code{NULL} invisibly after
printing a message. The result is a vector (or the supplied object), not
a list or an updated factor.
}
\examples{
library(dplyr) # Makes the unqualified pipe available.

# First recognised control label comes first; other labels keep input order.
group_order_set(c("treated", "Normal", "control", "treated"))
# c("Normal", "treated", "control")

# Without a control label, character labels are sorted.
group_order_set(c("B", "A", "B"))

# Apply an explicit order when constructing a factor.
groups <- c("treated", "control", "treated")
ordered_groups <- factor(
  groups,
  levels = group_order_set(groups, group_levels = c("control", "treated"))
)
levels(ordered_groups)

# Extra levels and duplicates are returned unchanged.
group_order_set(c("A", "B"), group_levels = c("B", "A", "unused", "A"))

# An incomplete explicit order prints a message and returns NULL.
is.null(group_order_set(c("A", "B"), group_levels = "A"))
}
