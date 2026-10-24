% Documentation for R/Drug_DrugReflector.R; keep in sync with its opening help comments.
\name{Drug_DrugReflector}
\alias{Drug_DrugReflector}
\title{
Run DrugReflector from Seurat expression data and plot drug predictions
}
\description{
Exports two groups of cells from a Seurat object to a temporary AnnData file,
calls the installed Python bioscPy.drugreflector.run_drugreflector function,
reads its prediction tables, and saves drug-score and mechanism-of-action
(MOA) figures. Model inference and biological preprocessing are delegated
to the Python package.
}
\usage{
Drug_DrugReflector(seurat,
  assay = "RNA",
  exp_layer = "data",
  group_col,
  group_from,
  group_to,
  species,
  sample_col = NULL,
  min_landmark_coverage = 0.9,
  run_leave_one_sample_out = TRUE,
  topn_drug = 50,
  colors = basicR::get_colors(number = 10.1,package = "viridis",name = "D"),
  output_dir,
  checkpoint_dir = "/home/zuoanjian/software/drugreflector/checkpoints",
  human_mouse_ortholog_file = "/reference/genome/human_mouse/Ortholog/HOM_ProteinCoding.csv",
  drug_annotation_file = "/reference/functional/LINCS/LINCS2020/metadata/compoundinfo_beta.txt",
  python = "/home/zuoanjian/miniforge3/envs/scRNA_Base/bin/python",
  keep_h5ad = FALSE,
  formats = c("pdf","png"))
}
\arguments{
\item{seurat}{In-memory Seurat object, required with no default; file paths are not
accepted. Supply an expression assay with unique, non-missing, non-empty
gene-symbol row names and cell names matching the metadata. All cells in
the two requested groups are used; subset beforehand for a cell-type-specific
analysis. Normalised, log1p-transformed expression is expected; this wrapper
does not normalise or verify that transformation.}
\item{assay}{Single non-empty assay name. Defaults to "RNA". Must exist in seurat;
the selected assay must contain the requested expression layer.}
\item{exp_layer}{Single non-empty expression-layer name. Defaults to "data".
For Assay5 objects, split layers starting with this name followed by a dot
are joined locally. Having both the exact layer and split layers is an error.
There is no fallback to counts or scale.data. The selected matrix is written
to AnnData.X, and Python receives exp_layer = NULL.}
\item{group_col}{Single non-empty metadata column name, required with no default.
Values are converted to character. Cells with missing values or values
outside group_from and group_to are excluded.}
\item{group_from}{Single non-empty starting-group label, required with no default.
Must occur in group_col and differ from group_to. Passed unchanged after
character conversion to the Python backend as the starting group.}
\item{group_to}{Single non-empty target-group label, required with no default.
Must occur in group_col. The requested direction is group_from to group_to;
the R wrapper does not calculate or reverse a contrast.}
\item{species}{Single species value, required with no default. Converted to lowercase
and restricted to "human" or "mouse". Mouse input requires an existing
human_mouse_ortholog_file; conversion itself is delegated to Python.}
\item{sample_col}{Single non-empty sample metadata column name, or NULL (default).
With NULL, the first available column among "sample" and "orig.ident" is
selected, in that order, even when run_leave_one_sample_out = FALSE.
If no column is found, leave-one-sample-out analysis must be disabled.
Selected cells must have non-missing sample values that are not empty after
trimming whitespace; original values are passed without trimming. Use a
column distinct from group_col: identical names cause sample values to
overwrite group values in the exported metadata. The wrapper does not check
sample replication or the allocation of samples across groups.}
\item{min_landmark_coverage}{Single finite numeric value in [0, 1]. Defaults to 0.9.
Passed to Python as a landmark-coverage threshold. The R wrapper reports
failed coverage flags after prediction but does not stop or filter results
because of them; backend behaviour depends on the installed version.}
\item{run_leave_one_sample_out}{Single non-missing logical value. Defaults to TRUE.
Requests leave-one-sample-out analysis from Python and requires a sample
column. The R wrapper does not perform sample exclusion or summarise its
results itself.}
\item{topn_drug}{Number of leading prediction-table rows used for the MOA figure.
Defaults to 50. Supply a single non-negative integer suitable for
dplyr::slice_head(n = topn_drug); there is no wrapper validation. Selection
uses the original file order before filtering missing MOA values, without
sorting by rank. Does not limit predictions returned or the drug-score plot.}
\item{colors}{Vector of valid colours for continuous probability and MOA-count
scales. Defaults to basicR::get_colors(number = 10.1, package = "viridis",
name = "D"). Passed to ggplot2 gradient scales without wrapper validation.}
\item{output_dir}{Single non-empty writable output-directory path, required with no
default. Created recursively if necessary and normalised to an absolute
path. Existing outputs may be overwritten; use a separate directory per run.
No output cache is reused by the R wrapper.}
\item{checkpoint_dir}{Existing directory containing model_fold_0.pt,
model_fold_1.pt, and model_fold_2.pt. Defaults to
"/home/zuoanjian/software/drugreflector/checkpoints". Paths are normalised;
only file existence is checked in R, not checkpoint content or compatibility.}
\item{human_mouse_ortholog_file}{Ortholog file passed to Python. Defaults to
"/reference/genome/human_mouse/Ortholog/HOM_ProteinCoding.csv".
For mouse input it must exist and its path is normalised. R does not inspect
its columns or mapping rules; supply the format required by the installed
backend. For human input it is passed through without an existence check.}
\item{drug_annotation_file}{Path to a tab-delimited compound annotation file,
or NULL. Defaults to
"/reference/functional/LINCS/LINCS2020/metadata/compoundinfo_beta.txt".
A non-NULL file must exist and have pert_id and cmap_name columns.
The plotting stage additionally requires moa in the returned predictions.
Although NULL passes input validation, plotting unconditionally uses
cmap_name and moa, so unannotated output lacking these columns cannot
complete this wrapper. No annotation-column content checks are performed.}
\item{python}{Single non-empty interpreter path or executable name found on PATH.
Defaults to "/home/zuoanjian/miniforge3/envs/scRNA_Base/bin/python".
Must exist and run an environment with compatible bioscPy[drugreflector]
and its dependencies installed. The directory is normalised while preserving
the interpreter basename to avoid resolving a virtual-environment symlink.}
\item{keep_h5ad}{Single non-missing logical value. Defaults to FALSE. If TRUE,
copies the exported input to output_dir/input.h5ad after Python succeeds,
overwriting an existing file. Copy failure issues a warning. Temporary
conversion files are removed on function exit, including on errors.}
\item{formats}{Character vector of figure formats accepted by basicR::save_figure().
Defaults to c("pdf", "png"). Used for both figures; required graphics devices
must be available. There is no wrapper validation.}
}
\details{
Only cells in the two requested groups are retained, in metadata row order.
Both groups must contain cells, and the expression layer must cover all
selected cells. No minimum cell count, gene-expression filter, differential
expression test, or sample balancing is performed in R. Assay5 split layers
are joined on a local assay object; legacy assays use LayerData() under
SeuratObject >= 5.0.0 and GetAssayData(slot = exp_layer) otherwise.

The transposed expression matrix is exported as cells by genes in AnnData.X
with gzip compression. Only the grouping column and, when available, sample
column are exported in obs; gene symbols become var row names. A temporary
JSON configuration and two temporary Python scripts check and invoke the
installed backend. No user-supplied standalone Python script is required.
The initial API check tests only the presence of input_h5ad, group_col,
sample_col, checkpoint_dir, and exp_layer parameters, not compatibility with
all arguments subsequently passed or with the expected output schema.

Group labels are passed as group_from and group_to without reordering.
The internal comment describes a standardised mean(to) minus mean(from)
contrast, but that calculation is not implemented in this R function.
Contrast construction, ortholog mapping, landmark selection, model scoring,
rank generation, annotation matching, and leave-one-sample-out behaviour
are controlled by the installed Python backend. Their exact formulas and
filtering rules cannot be established from this wrapper alone. In particular,
prob is a backend score and is not validated here as a calibrated probability
of therapeutic efficacy.

After a successful Python process, the wrapper requires
04_landmark_coverage.tsv (including passes_threshold),
05_predictions_all.tsv, and 09_run_manifest.json (including output_files).
It reads 06_predictions_annotated.tsv instead of the raw prediction table
only when that filename occurs in manifest$output_files. Other columns and
backend-generated files are version-dependent. Empty strings and "NA" in
the prediction table are read as missing. Coverage flags are converted to
logical values by lowercasing their text and comparing with "true";
missing or malformed coverage output is not handled robustly.

For successful plotting, predictions must contain rank, logit, prob,
display_order, cmap_name, pert_id, and moa with suitable types. The score
plot includes all predictions, sorted for plotting by increasing rank,
decreasing logit, then decreasing prob. Colour represents prob. Labels use
cmap_name where display_order <= 5, excluding missing or empty labels;
this is independent of topn_drug and is not a test of rank <= 5. Text
repulsion uses seed = 123; no model-inference seed is set in R.

The MOA plot takes the first topn_drug rows of the original prediction table,
removes missing or empty moa, splits MOA strings on the literal pipe
character, trims whitespace, removes empty terms, deduplicates pert_id/moa
pairs, and counts distinct compounds per term. One compound can contribute
to several terms. Missing annotations are not replaced, and there is no
significance or score threshold in the R plotting stage. Plot objects, label
columns, and MOA counts are not added to the returned prediction table.

R writes drugreflector.log (Python stdout and stderr) and
10_r_input_manifest.json (original object class and dimensions, chosen assay
and layer, group/sample settings, interpreter, keep_h5ad, R version, and
versions of the four explicitly checked R packages). The recorded dimensions
refer to the original object, not the selected cells or necessarily the
selected assay. Optional input.h5ad contains the selected expression data.
Figure stems are 11_drugs_Sort (width = 6.5, height = 6) and 12_drugs_MOA
(width = 7.5, height = 6), with res = 600 passed to basicR::save_figure().
The files element is assembled before plotting and does not add these figure
paths. Additional backend files are listed in 09_run_manifest.json.

A writable R temporary directory, output directory, compatible Python
installation, model checkpoints, and appropriate reference files are needed.
The default interpreter and reference paths are machine-specific. R explicitly
checks SeuratObject, Matrix, anndataR, and jsonlite; plotting additionally uses
basicR, dplyr, ggplot2, ggrepel, tidyr, and forcats, with stats and utils for
file/result handling. DESCRIPTION currently omits Matrix, anndataR, jsonlite,
ggrepel, and tidyr. Several dplyr/ggplot2 helpers and the pipe operator are
called without namespace qualification and are not imported by the current
NAMESPACE. Attach dplyr and ggplot2 before calling this function, as in the
example; merely listing these packages in Imports does not import all helpers.

The function returns only after both figures have been saved. Python errors
stop execution with the last 30 log lines. Later file-reading or plotting
errors can leave completed predictions and partial outputs without a return
value, even after the completion message has been printed. Concurrent runs
have private temporary directories but must use different output directories
to avoid colliding result files.
}
\value{
A visible named list, returned only after successful plotting:
\describe{
\item{predictions}{Data frame read from the annotated prediction file when
listed in the Python manifest, otherwise from the raw file. All rows and
backend columns are retained in file order. rank, logit, prob, display_order,
cmap_name, pert_id, and moa are used for plotting; their values originate
from Python, not an R-side ranking or annotation calculation.}
\item{coverage}{Data frame read from 04_landmark_coverage.tsv, with
passes_threshold converted to logical; other columns are retained.}
\item{output_dir}{Normalised absolute output-directory path.}
\item{predictions_file}{Path to the table actually read into predictions.}
\item{raw_predictions_file}{Path to 05_predictions_all.tsv.}
\item{manifest_file}{Path to the Python manifest, 09_run_manifest.json.}
\item{r_manifest_file}{Path to 10_r_input_manifest.json.}
\item{log_file}{Path to drugreflector.log.}
\item{files}{Named list of paths for manifest$output_files, the log, the R
manifest, and input.h5ad if requested and present. Names are filenames.
Paths from the manifest are not all checked for existence. R-generated
figure paths are not appended to this list.}
\item{annotated_file}{Present only when the annotated table is listed in the
Python manifest; identical to predictions_file in that case.}
}
}
\examples{
\dontrun{
library(bioscR)
library(Seurat)
library(dplyr)
library(ggplot2)
seu <- readRDS("seurat_normalised.rds")
# Use an appropriate cell-type subset if needed. Replace paths and labels.
# RNA/data must already contain normalised, log1p-transformed expression.
# Annotation must yield cmap_name and moa; checkpoints must match bioscPy.
result <- Drug_DrugReflector(
  seurat = seu,
  group_col = "condition",
  group_from = "disease",
  group_to = "control",
  species = "human",
  sample_col = "sample",
  output_dir = "drugreflector_results",
  checkpoint_dir = "/path/to/drugreflector/checkpoints",
  drug_annotation_file = "/path/to/compoundinfo_beta.txt",
  python = "/path/to/environment/bin/python"
)
head(result$predictions)
result$coverage
# For mouse data also set species = "mouse" and provide a compatible
# human_mouse_ortholog_file. If sample metadata are unavailable, set
# run_leave_one_sample_out = FALSE; sample_col = NULL still auto-detects
# sample or orig.ident when either column exists.
}
}
