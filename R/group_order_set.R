#' Set the order of group labels
#'
#' Returns an order for group labels, either using a supplied order or
#' placing the first recognised control label before the remaining labels.
#'
#' @param vector A vector of group labels, typically a character vector or factor.
#' Required, with no default. The function applies `unique()` before
#' `as.character()`; unused factor levels are not included. Numeric labels
#' are therefore ordered as character strings. Missing and blank labels are
#' not explicitly validated or cleaned; see Details.
#' @param group_levels An optional vector specifying the desired order. Defaults
#' to `NULL`; any zero-length value selects automatic ordering. For intended
#' use, supply a character vector containing every distinct character label
#' in `vector`. Matching is case-sensitive. Extra labels and duplicates are
#' accepted, and a matching vector is returned unchanged, including its type
#' and attributes. There is no explicit type or uniqueness validation.
#'
#' @details
#' In automatic mode, labels are deduplicated and converted to character.
#' The function searches for `"control"`, `"normal"`, or `"norm"` without regard
#' to case. If several labels match, the first one in input order is moved to
#' the front; the other labels retain their input order. There is no priority
#' among the three control names. Original spelling and case are retained,
#' and whitespace is not trimmed.
#'
#' If no control label matches, `sort()` orders the character labels in
#' ascending order using the current locale and drops `NA` values by default.
#' When a control label matches, missing values instead remain among the
#' other labels. Blank strings are not removed. Empty input in automatic
#' mode returns `character(0)`.
#'
#' With a non-empty `group_levels`, the function checks only whether every
#' processed input label occurs in `group_levels`, using `%in%`. It does not
#' filter extra levels, remove duplicates, or reorder the supplied vector.
#' If any label is absent, it prints `"group_levels Not Match the vector!"`
#' with ANSI colour escape sequences and returns `NULL` invisibly. This is
#' a console message produced by `cat()`, not a warning or an error.
#'
#' The result specifies label order only. The function does not modify a
#' factor or metadata column, select observations, or perform comparisons;
#' it does not define a statistical reference group or comparison direction.
#'
#' The implementation uses base R functions and the unqualified `%>%` pipe.
#' The pipe must be available in the function's evaluation environment, for
#' example by attaching dplyr before calling. dplyr is declared in the
#' package Imports, but the current NAMESPACE does not import `%>%`.
#' No external scripts, input files, or output files are used. The only
#' side effect is the console message on an explicit-order mismatch.
#'
#' @return In automatic mode, a character vector of group labels in the computed
#' order. With a matching non-empty `group_levels`, that object is returned
#' unchanged, so its type, attributes, extra labels, and duplicates are
#' preserved. On an explicit-order mismatch, returns `NULL` invisibly after
#' printing a message. The result is a vector (or the supplied object), not
#' a list or an updated factor.
#'
#' @examples
#' library(dplyr) # Makes the unqualified pipe available.
#'
#' # First recognised control label comes first; other labels keep input order.
#' group_order_set(c("treated", "Normal", "control", "treated"))
#' # c("Normal", "treated", "control")
#'
#' # Without a control label, character labels are sorted.
#' group_order_set(c("B", "A", "B"))
#'
#' # Apply an explicit order when constructing a factor.
#' groups <- c("treated", "control", "treated")
#' ordered_groups <- factor(
#'   groups,
#'   levels = group_order_set(groups, group_levels = c("control", "treated"))
#' )
#' levels(ordered_groups)
#'
#' # Extra levels and duplicates are returned unchanged.
#' group_order_set(c("A", "B"), group_levels = c("B", "A", "unused", "A"))
#'
#' # An incomplete explicit order prints a message and returns NULL.
#' is.null(group_order_set(c("A", "B"), group_levels = "A"))
#' @export

group_order_set <- function(vector,group_levels = NULL){
    # 设置分组水平
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 获取列名
    vector <- unique(vector) %>% as.character()
    if(length(group_levels) == 0){
        if(tolower(vector) %in% c("control", "normal","norm") %>% any()){

            # 找出小写化后等于 "control" 或 "normal" 的列名
            target_col <- vector[tolower(vector) %in% c("control", "normal","norm")]
            
            # 如果存在符合条件的列
            if (length(target_col) > 0) {
              # 将目标列放在第一位，其余列保持原顺序但去掉目标列
              new_order <- c(target_col[1], setdiff(vector, target_col[1]))
              return(new_order)  
            }
        }else{
            return(sort(vector))  
        }
    }else{
        if(all(vector %in% group_levels)){
            return(group_levels)
        }else{
            cat("\n",yellow,"group_levels Not Match the vector!",reset,"\n")
        }
    }
}