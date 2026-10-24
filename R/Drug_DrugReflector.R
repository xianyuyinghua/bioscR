#' Run DrugReflector from Seurat expression data and plot drug predictions
#'
#' Exports two groups of cells from a Seurat object to a temporary AnnData file,
#' calls the installed Python bioscPy.drugreflector.run_drugreflector function,
#' reads its prediction tables, and saves drug-score and mechanism-of-action
#' (MOA) figures. Model inference and biological preprocessing are delegated
#' to the Python package.
#'
#' @param seurat In-memory Seurat object, required with no default; file paths are not
#' accepted. Supply an expression assay with unique, non-missing, non-empty
#' gene-symbol row names and cell names matching the metadata. All cells in
#' the two requested groups are used; subset beforehand for a cell-type-specific
#' analysis. Normalised, log1p-transformed expression is expected; this wrapper
#' does not normalise or verify that transformation.
#' @param assay Single non-empty assay name. Defaults to "RNA". Must exist in seurat;
#' the selected assay must contain the requested expression layer.
#' @param exp_layer Single non-empty expression-layer name. Defaults to "data".
#' For Assay5 objects, split layers starting with this name followed by a dot
#' are joined locally. Having both the exact layer and split layers is an error.
#' There is no fallback to counts or scale.data. The selected matrix is written
#' to AnnData.X, and Python receives exp_layer = NULL.
#' @param group_col Single non-empty metadata column name, required with no default.
#' Values are converted to character. Cells with missing values or values
#' outside group_from and group_to are excluded.
#' @param group_from Single non-empty starting-group label, required with no default.
#' Must occur in group_col and differ from group_to. Passed unchanged after
#' character conversion to the Python backend as the starting group.
#' @param group_to Single non-empty target-group label, required with no default.
#' Must occur in group_col. The requested direction is group_from to group_to;
#' the R wrapper does not calculate or reverse a contrast.
#' @param species Single species value, required with no default. Converted to lowercase
#' and restricted to "human" or "mouse". Mouse input requires an existing
#' human_mouse_ortholog_file; conversion itself is delegated to Python.
#' @param sample_col Single non-empty sample metadata column name, or NULL (default).
#' With NULL, the first available column among "sample" and "orig.ident" is
#' selected, in that order, even when run_leave_one_sample_out = FALSE.
#' If no column is found, leave-one-sample-out analysis must be disabled.
#' Selected cells must have non-missing sample values that are not empty after
#' trimming whitespace; original values are passed without trimming. Use a
#' column distinct from group_col: identical names cause sample values to
#' overwrite group values in the exported metadata. The wrapper does not check
#' sample replication or the allocation of samples across groups.
#' @param min_landmark_coverage Single finite numeric value in [0, 1]. Defaults to 0.9.
#' Passed to Python as a landmark-coverage threshold. The R wrapper reports
#' failed coverage flags after prediction but does not stop or filter results
#' because of them; backend behaviour depends on the installed version.
#' @param run_leave_one_sample_out Single non-missing logical value. Defaults to TRUE.
#' Requests leave-one-sample-out analysis from Python and requires a sample
#' column. The R wrapper does not perform sample exclusion or summarise its
#' results itself.
#' @param topn_drug Number of leading prediction-table rows used for the MOA figure.
#' Defaults to 50. Supply a single non-negative integer suitable for
#' dplyr::slice_head(n = topn_drug); there is no wrapper validation. Selection
#' uses the original file order before filtering missing MOA values, without
#' sorting by rank. Does not limit predictions returned or the drug-score plot.
#' @param colors Vector of valid colours for continuous probability and MOA-count
#' scales. Defaults to basicR::get_colors(number = 10.1, package = "viridis",
#' name = "D"). Passed to ggplot2 gradient scales without wrapper validation.
#' @param output_dir Single non-empty writable output-directory path, required with no
#' default. Created recursively if necessary and normalised to an absolute
#' path. Existing outputs may be overwritten; use a separate directory per run.
#' No output cache is reused by the R wrapper.
#' @param checkpoint_dir Existing directory containing model_fold_0.pt,
#' model_fold_1.pt, and model_fold_2.pt. Defaults to
#' "/home/zuoanjian/software/drugreflector/checkpoints". Paths are normalised;
#' only file existence is checked in R, not checkpoint content or compatibility.
#' @param human_mouse_ortholog_file Ortholog file passed to Python. Defaults to
#' "/reference/genome/human_mouse/Ortholog/HOM_ProteinCoding.csv".
#' For mouse input it must exist and its path is normalised. R does not inspect
#' its columns or mapping rules; supply the format required by the installed
#' backend. For human input it is passed through without an existence check.
#' @param drug_annotation_file Path to a tab-delimited compound annotation file,
#' or NULL. Defaults to
#' "/reference/functional/LINCS/LINCS2020/metadata/compoundinfo_beta.txt".
#' A non-NULL file must exist and have pert_id and cmap_name columns.
#' The plotting stage additionally requires moa in the returned predictions.
#' Although NULL passes input validation, plotting unconditionally uses
#' cmap_name and moa, so unannotated output lacking these columns cannot
#' complete this wrapper. No annotation-column content checks are performed.
#' @param python Single non-empty interpreter path or executable name found on PATH.
#' Defaults to "/home/zuoanjian/miniforge3/envs/scRNA_Base/bin/python".
#' Must exist and run an environment with compatible bioscPy[drugreflector]
#' and its dependencies installed. The directory is normalised while preserving
#' the interpreter basename to avoid resolving a virtual-environment symlink.
#' @param keep_h5ad Single non-missing logical value. Defaults to FALSE. If TRUE,
#' copies the exported input to output_dir/input.h5ad after Python succeeds,
#' overwriting an existing file. Copy failure issues a warning. Temporary
#' conversion files are removed on function exit, including on errors.
#' @param formats Character vector of figure formats accepted by basicR::save_figure().
#' Defaults to c("pdf", "png"). Used for both figures; required graphics devices
#' must be available. There is no wrapper validation.
#'
#' @details
#' Only cells in the two requested groups are retained, in metadata row order.
#' Both groups must contain cells, and the expression layer must cover all
#' selected cells. No minimum cell count, gene-expression filter, differential
#' expression test, or sample balancing is performed in R. Assay5 split layers
#' are joined on a local assay object; legacy assays use LayerData() under
#' SeuratObject >= 5.0.0 and GetAssayData(slot = exp_layer) otherwise.
#'
#' The transposed expression matrix is exported as cells by genes in AnnData.X
#' with gzip compression. Only the grouping column and, when available, sample
#' column are exported in obs; gene symbols become var row names. A temporary
#' JSON configuration and two temporary Python scripts check and invoke the
#' installed backend. No user-supplied standalone Python script is required.
#' The initial API check tests only the presence of input_h5ad, group_col,
#' sample_col, checkpoint_dir, and exp_layer parameters, not compatibility with
#' all arguments subsequently passed or with the expected output schema.
#'
#' Group labels are passed as group_from and group_to without reordering.
#' The internal comment describes a standardised mean(to) minus mean(from)
#' contrast, but that calculation is not implemented in this R function.
#' Contrast construction, ortholog mapping, landmark selection, model scoring,
#' rank generation, annotation matching, and leave-one-sample-out behaviour
#' are controlled by the installed Python backend. Their exact formulas and
#' filtering rules cannot be established from this wrapper alone. In particular,
#' prob is a backend score and is not validated here as a calibrated probability
#' of therapeutic efficacy.
#'
#' After a successful Python process, the wrapper requires
#' 04_landmark_coverage.tsv (including passes_threshold),
#' 05_predictions_all.tsv, and 09_run_manifest.json (including output_files).
#' It reads 06_predictions_annotated.tsv instead of the raw prediction table
#' only when that filename occurs in manifest$output_files. Other columns and
#' backend-generated files are version-dependent. Empty strings and "NA" in
#' the prediction table are read as missing. Coverage flags are converted to
#' logical values by lowercasing their text and comparing with "true";
#' missing or malformed coverage output is not handled robustly.
#'
#' For successful plotting, predictions must contain rank, logit, prob,
#' display_order, cmap_name, pert_id, and moa with suitable types. The score
#' plot includes all predictions, sorted for plotting by increasing rank,
#' decreasing logit, then decreasing prob. Colour represents prob. Labels use
#' cmap_name where display_order <= 5, excluding missing or empty labels;
#' this is independent of topn_drug and is not a test of rank <= 5. Text
#' repulsion uses seed = 123; no model-inference seed is set in R.
#'
#' The MOA plot takes the first topn_drug rows of the original prediction table,
#' removes missing or empty moa, splits MOA strings on the literal pipe
#' character, trims whitespace, removes empty terms, deduplicates pert_id/moa
#' pairs, and counts distinct compounds per term. One compound can contribute
#' to several terms. Missing annotations are not replaced, and there is no
#' significance or score threshold in the R plotting stage. Plot objects, label
#' columns, and MOA counts are not added to the returned prediction table.
#'
#' R writes drugreflector.log (Python stdout and stderr) and
#' 10_r_input_manifest.json (original object class and dimensions, chosen assay
#' and layer, group/sample settings, interpreter, keep_h5ad, R version, and
#' versions of the four explicitly checked R packages). The recorded dimensions
#' refer to the original object, not the selected cells or necessarily the
#' selected assay. Optional input.h5ad contains the selected expression data.
#' Figure stems are 11_drugs_Sort (width = 6.5, height = 6) and 12_drugs_MOA
#' (width = 7.5, height = 6), with res = 600 passed to basicR::save_figure().
#' The files element is assembled before plotting and does not add these figure
#' paths. Additional backend files are listed in 09_run_manifest.json.
#'
#' A writable R temporary directory, output directory, compatible Python
#' installation, model checkpoints, and appropriate reference files are needed.
#' The default interpreter and reference paths are machine-specific. R explicitly
#' checks SeuratObject, Matrix, anndataR, and jsonlite; plotting additionally uses
#' basicR, dplyr, ggplot2, ggrepel, tidyr, and forcats, with stats and utils for
#' file/result handling. DESCRIPTION currently omits Matrix, anndataR, jsonlite,
#' ggrepel, and tidyr. Several dplyr/ggplot2 helpers and the pipe operator are
#' called without namespace qualification and are not imported by the current
#' NAMESPACE. Attach dplyr and ggplot2 before calling this function, as in the
#' example; merely listing these packages in Imports does not import all helpers.
#'
#' The function returns only after both figures have been saved. Python errors
#' stop execution with the last 30 log lines. Later file-reading or plotting
#' errors can leave completed predictions and partial outputs without a return
#' value, even after the completion message has been printed. Concurrent runs
#' have private temporary directories but must use different output directories
#' to avoid colliding result files.
#'
#' @return A visible named list, returned only after successful plotting:
#' \describe{
#' \item{predictions}{Data frame read from the annotated prediction file when
#' listed in the Python manifest, otherwise from the raw file. All rows and
#' backend columns are retained in file order. rank, logit, prob, display_order,
#' cmap_name, pert_id, and moa are used for plotting; their values originate
#' from Python, not an R-side ranking or annotation calculation.}
#' \item{coverage}{Data frame read from 04_landmark_coverage.tsv, with
#' passes_threshold converted to logical; other columns are retained.}
#' \item{output_dir}{Normalised absolute output-directory path.}
#' \item{predictions_file}{Path to the table actually read into predictions.}
#' \item{raw_predictions_file}{Path to 05_predictions_all.tsv.}
#' \item{manifest_file}{Path to the Python manifest, 09_run_manifest.json.}
#' \item{r_manifest_file}{Path to 10_r_input_manifest.json.}
#' \item{log_file}{Path to drugreflector.log.}
#' \item{files}{Named list of paths for manifest$output_files, the log, the R
#' manifest, and input.h5ad if requested and present. Names are filenames.
#' Paths from the manifest are not all checked for existence. R-generated
#' figure paths are not appended to this list.}
#' \item{annotated_file}{Present only when the annotated table is listed in the
#' Python manifest; identical to predictions_file in that case.}
#' }
#'
#' @examples
#' \dontrun{
#' library(bioscR)
#' library(Seurat)
#' library(dplyr)
#' library(ggplot2)
#' seu <- readRDS("seurat_normalised.rds")
#' # Use an appropriate cell-type subset if needed. Replace paths and labels.
#' # RNA/data must already contain normalised, log1p-transformed expression.
#' # Annotation must yield cmap_name and moa; checkpoints must match bioscPy.
#' result <- Drug_DrugReflector(
#'   seurat = seu,
#'   group_col = "condition",
#'   group_from = "disease",
#'   group_to = "control",
#'   species = "human",
#'   sample_col = "sample",
#'   output_dir = "drugreflector_results",
#'   checkpoint_dir = "/path/to/drugreflector/checkpoints",
#'   drug_annotation_file = "/path/to/compoundinfo_beta.txt",
#'   python = "/path/to/environment/bin/python"
#' )
#' head(result$predictions)
#' result$coverage
#' # For mouse data also set species = "mouse" and provide a compatible
#' # human_mouse_ortholog_file. If sample metadata are unavailable, set
#' # run_leave_one_sample_out = FALSE; sample_col = NULL still auto-detects
#' # sample or orig.ident when either column exists.
#' }
#' @export
Drug_DrugReflector <- function(seurat,
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
                              formats = c("pdf","png")
                             ){
    
    # group_from -> group_to: v = (mean(to) - mean(from)) / sqrt(var(to)+var(from)).
    # All cells in the supplied Seurat object belonging to these two groups are used.
    # Subset the Seurat object before calling if a specific cell type is needed.
    
    for (pkg in c("SeuratObject", "Matrix", "anndataR", "jsonlite")) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少 R 包：", pkg)
    }
    
    scalar_string <- function(x, label){
        if(length(x) != 1L || is.na(x) || !nzchar(as.character(x))){
            stop(label, " 必须是单个非空值")
        }
        as.character(x)
    }
    
    if (!inherits(seurat, "Seurat")) {
        stop("seurat 必须是内存中的 Seurat 对象（例如 subseurat），不接受文件路径")
    }
    
    group_col <- scalar_string(group_col, "group_col")
    group_from <- scalar_string(group_from, "group_from")
    group_to <- scalar_string(group_to, "group_to")
    species <- tolower(scalar_string(species, "species"))
    
    if (!species %in% c("human", "mouse")) stop("species 仅支持 human / mouse")
    if (identical(group_from, group_to)) stop("group_from 和 group_to 必须不同")
    if (!is.numeric(min_landmark_coverage) || length(min_landmark_coverage) != 1L ||
        !is.finite(min_landmark_coverage) || min_landmark_coverage < 0 || min_landmark_coverage > 1) {
        stop("min_landmark_coverage 必须为 0–1 的数值；它仅用于提示，不阻止预测")
    }
    
    for (x in list(run_leave_one_sample_out, keep_h5ad)) {
        if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("开关参数必须为 TRUE / FALSE")
    }
    
    assay <- scalar_string(assay, "assay")
    exp_layer <- scalar_string(exp_layer, "exp_layer")
    python <- path.expand(scalar_string(python, "python（请指定安装 bioscPy 的 Python）"))
    if (!file.exists(python)) python <- unname(Sys.which(python))
    if (!nzchar(python) || !file.exists(python)) stop("找不到 Python，请传入解释器的完整路径")
    
    # Preserve the interpreter basename: resolving a venv symlink can select base Python.
    python <- file.path(normalizePath(dirname(python), mustWork = TRUE), basename(python))
    checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
    for (fold in 0:2) {
        p <- file.path(checkpoint_dir, paste0("model_fold_", fold, ".pt"))
        if (!file.exists(p)) stop("缺少模型：", p)
    }
    
    if (species == "mouse") {
        human_mouse_ortholog_file <- normalizePath(human_mouse_ortholog_file, mustWork = TRUE)
    }
    if (!is.null(drug_annotation_file)) {
        drug_annotation_file <- normalizePath(
            scalar_string(drug_annotation_file, "drug_annotation_file"), mustWork = TRUE
        )
        annotation_columns <- names(utils::read.delim(drug_annotation_file, nrows = 0L, check.names = FALSE))
        if (!all(c("pert_id", "cmap_name") %in% annotation_columns)) {
            stop("药物注释文件必须包含 pert_id 和 cmap_name（药物名称）列")
        }
    }
    
    output_dir <- path.expand(scalar_string(output_dir, "output_dir"))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    # A private temporary directory prevents concurrent calls sharing conversion files.
    work_dir <- tempfile("drugreflector_")
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

    check_script <- file.path(work_dir, "check_bioscpy.py")
    writeLines(c(
        "import inspect",
        "from bioscPy.drugreflector import run_drugreflector",
        "required = {'input_h5ad', 'group_col', 'sample_col', 'checkpoint_dir', 'exp_layer'}",
        "if not required.issubset(inspect.signature(run_drugreflector).parameters):",
        "    raise RuntimeError('Installed bioscPy has an incompatible run_drugreflector API')"
    ), check_script)
    
    check_log <- file.path(work_dir, "check.log")
    check_status <- system2(python, shQuote(check_script), stdout = check_log, stderr = check_log)
    
    if (check_status != 0L) {
        stop("指定的 Python 无法导入 bioscPy.drugreflector.run_drugreflector：\n",
             paste(readLines(check_log, warn = FALSE), collapse = "\n"),
             "\n请在该 Python 环境安装 bioscPy[drugreflector]。")
    }

    message("使用传入的 Seurat 对象")
    obj <- seurat
    input_summary <- list(class = class(obj), n_cells = ncol(obj), n_genes = nrow(obj))
    if (!assay %in% names(obj@assays)) stop("Seurat 中没有 assay：", assay)
    metadata <- obj[[]]
    if (!group_col %in% colnames(metadata)) stop("缺少分组列：", group_col)
    groups <- as.character(metadata[[group_col]])
    selected <- !is.na(groups) & groups %in% c(group_from, group_to)
    if (!all(c(group_from, group_to) %in% groups[selected])) stop("起始组或目标组没有细胞")
    cells <- rownames(metadata)[selected]
    metadata <- metadata[cells, , drop = FALSE]
    if (is.null(sample_col)) {
        candidates <- intersect(c("sample", "orig.ident"), colnames(metadata))
        if (length(candidates)) sample_col <- candidates[[1]]
    }
    if (!is.null(sample_col)) {
        sample_col <- scalar_string(sample_col, "sample_col")
        if (!sample_col %in% colnames(metadata)) stop("缺少样本列：", sample_col)
    } else if (run_leave_one_sample_out) {
        stop("留一样本分析需要样本列；请指定 sample_col，或设置 run_leave_one_sample_out = FALSE")
    }
    message("分组方向：", group_from, " -> ", group_to,"；样本列：", if (is.null(sample_col)) "未指定" else sample_col)
        
    # Join split Seurat v5 layers before extraction. Never substitute counts/scale.data.
    a <- obj[[assay]]
    if (inherits(a, "Assay5")) {
        layers <- SeuratObject::Layers(a, search = NA)
        split_layers <- layers[startsWith(layers, paste0(exp_layer, "."))]
        if (length(split_layers)) {
        if (exp_layer %in% layers) stop("同时存在完整和分拆的表达层，请先整理：", exp_layer)
        a <- SeuratObject::JoinLayers(a, layers = exp_layer, new = exp_layer)
        }
        if (!exp_layer %in% SeuratObject::Layers(a, search = NA)) {
        stop("缺少表达层：", exp_layer, "；需要已归一化且 log1p 转换的 RNA data")
        }
        x <- SeuratObject::LayerData(a, layer = exp_layer)
    }else if (utils::packageVersion("SeuratObject") >= "5.0.0") {
        # Legacy Assay objects loaded under SeuratObject v5 also require the layer API.
        x <- SeuratObject::LayerData(a, layer = exp_layer)
    } else {
        x <- SeuratObject::GetAssayData(obj, assay = assay, slot = exp_layer)
    }
            
    if (!all(cells %in% colnames(x))) stop("表达层未覆盖所选细胞")
    x <- x[, cells, drop = FALSE]
    if (!nrow(x) || !ncol(x) || anyDuplicated(rownames(x))) stop("表达为空或基因名重复")
    if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
        stop("表达矩阵必须使用非空基因符号作为行名")
    }
    # Preserve actual metadata names; pass these names to the Python package.
    obs <- data.frame(row.names = cells)
    obs[[group_col]] <- as.character(metadata[[group_col]])
    if (!is.null(sample_col)) {
        obs[[sample_col]] <- as.character(metadata[[sample_col]])
        if (anyNA(obs[[sample_col]]) || any(!nzchar(trimws(obs[[sample_col]])))) {
            stop("所选细胞的样本列包含空值")
        }
    }
    h5ad <- file.path(work_dir, "input.h5ad")
    adata <- anndataR::AnnData(X = Matrix::t(x), obs = obs,var = data.frame(row.names = rownames(x)))
    anndataR::write_h5ad(adata, h5ad, compression = "gzip")
    rm(obj, a, x, adata)
    invisible(gc())
    # Only Python function arguments go into this JSON; no analysis is reimplemented in R.
    parameters <- list(input_h5ad = h5ad,
                       output_dir = output_dir,
                       group_col = group_col, 
                       group_from = group_from, 
                       group_to = group_to,
                       species = species, 
                       sample_col = sample_col,
                       human_mouse_ortholog_file = human_mouse_ortholog_file,
                       min_landmark_coverage = min_landmark_coverage,
                       run_leave_one_sample_out = run_leave_one_sample_out,
                       checkpoint_dir = checkpoint_dir,
                       drug_annotation_file = drug_annotation_file,
                       exp_layer = NULL  # R exported the chosen Seurat layer into AnnData.X.
                      )
    config_path <- file.path(work_dir, "parameters.json")
    jsonlite::write_json(parameters, config_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    # Temporary bridge calls the installed package directly. No project Python script needed.
    bridge <- file.path(work_dir, "call_bioscpy.py")
    writeLines(c("import json, sys",
                 "from bioscPy.drugreflector import run_drugreflector",
                 "with open(sys.argv[1], encoding='utf-8') as handle:",
                 "    parameters = json.load(handle)",
                 "run_drugreflector(**parameters)"
                ), bridge)
    log_file <- file.path(output_dir, "drugreflector.log")
    message("调用已安装的 bioscPy 包；日志：", log_file)
            
    status <- system2(python, args = c(shQuote(bridge), shQuote(config_path)),stdout = log_file, stderr = log_file)
    if (status != 0L) {
        tail_log <- tail(readLines(log_file, warn = FALSE), 30)
        stop("DrugReflector 运行失败（exit ", status, "）：\n", paste(tail_log, collapse = "\n"))
    }
    if (keep_h5ad && !file.copy(h5ad, file.path(output_dir, "input.h5ad"), overwrite = TRUE)) {
        warning("预测已完成，但未能保留 input.h5ad")
    }
            
    r_manifest_file <- file.path(output_dir, "10_r_input_manifest.json")
    r_manifest <- list(
        input_type = "Seurat", input_summary = input_summary,
        assay = assay, expression_layer = exp_layer,
        group_col = group_col, group_from = group_from, group_to = group_to,
        sample_col = sample_col, python = python, keep_h5ad = keep_h5ad,
        r_version = R.version.string,
        r_packages = lapply(c("SeuratObject", "Matrix", "anndataR", "jsonlite"),
                           function(p) list(package = p, version = as.character(utils::packageVersion(p))))
    )
        
    jsonlite::write_json(r_manifest, r_manifest_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
    coverage <- read.delim(file.path(output_dir, "04_landmark_coverage.tsv"), check.names = FALSE)
    coverage$passes_threshold <- tolower(as.character(coverage$passes_threshold)) == "true"
    if (any(!coverage$passes_threshold)) {
        message("提示：landmark 覆盖率低于设定阈值，但已继续完成全部预测。")
    }
        
    manifest_file <- file.path(output_dir, "09_run_manifest.json")
    manifest <- jsonlite::read_json(manifest_file, simplifyVector = TRUE)
    raw_predictions_file <- file.path(output_dir, "05_predictions_all.tsv")
    has_annotation <- "06_predictions_annotated.tsv" %in% manifest$output_files
    predictions_file <- if (has_annotation) {
        file.path(output_dir, "06_predictions_annotated.tsv")
    } else raw_predictions_file
    output_files <- c(manifest$output_files, "drugreflector.log", "10_r_input_manifest.json")
    if (keep_h5ad && file.exists(file.path(output_dir, "input.h5ad"))) output_files <- c(output_files, "input.h5ad")
        
    result <- list(predictions = read.delim(predictions_file, check.names = FALSE, na.strings = c("", "NA")),
                 coverage = coverage,
                 output_dir = output_dir,
                 predictions_file = predictions_file, 
                 raw_predictions_file = raw_predictions_file,
                 manifest_file = manifest_file, 
                 r_manifest_file = r_manifest_file,
                 log_file = log_file,
                 files = stats::setNames(as.list(file.path(output_dir, output_files)), output_files)
                )
        
    if (has_annotation) {
        result$annotated_file <- predictions_file
        message("药物名称列：cmap_name；未匹配到注释的化合物仍保留。")
    }
        
    message("完成：", nrow(result$predictions), " 个化合物；结果目录：", output_dir)
        
    
        
    ##################################################    Plot   #########################################################
    df_result <- result$predictions
    df_result$label <- ifelse(df_result$display_order <= 5,df_result$cmap_name,NA)    

    df_plot <- df_result %>% arrange(rank, desc(logit), desc(prob))
    df_label <- df_plot %>% filter(!is.na(label), label != "")
    
    p <- ggplot(df_plot, aes(x = rank, y = logit)) +
          geom_point(aes(color = prob),alpha = 0.7,shape = 16,size = 2) +
          ggrepel::geom_text_repel(data = df_label,aes(label = label),size = 4,fontface = "bold",
                                   box.padding = 0.8,point.padding = 0.4,force = 5,force_pull = 0.1,
                                   max.iter = 50000,max.time = 5,max.overlaps = Inf,min.segment.length = 0,seed = 123) +
          scale_color_gradientn(colours = colors)+
          scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
          labs(x = "Drug rank",y = "DrugReflector logit score",color = "Probability") +
          theme_classic(base_size = 16) +
          theme(
            axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            plot.margin = margin(t = 25, r = 20, b = 20, l = 30)
          )
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, "11_drugs_Sort"),
                        width = 6.5,
                        height = 6,
                        res = 600,
                        formats = formats
                       ) 

        
    moa_count <- df_result %>% slice_head(n = topn_drug) %>%
          dplyr::filter(!is.na(moa), nzchar(moa)) %>%
          tidyr::separate_rows(moa, sep = "\\|") %>%
          dplyr::mutate(moa = trimws(moa)) %>%
          dplyr::filter(nzchar(moa)) %>%
          dplyr::distinct(pert_id, moa) %>%
          dplyr::count(moa, name = "n_drugs", sort = TRUE)

    df_plot <- moa_count %>% mutate(moa = forcats::fct_reorder(moa, n_drugs))
    
    p <- ggplot(df_plot,aes(x = n_drugs, y = moa)) +
      geom_col(aes(fill = n_drugs),width = 0.7) +
      scale_fill_gradientn(colours = colors) +
      geom_text(aes(label = n_drugs),hjust = -0.3,size = 4,fontface = "bold") +
      labs(x = "Number of candidate drugs",y = NULL) +
      theme_classic(base_size = 16) +
      theme(axis.title.x = element_text(face = "bold"),
            axis.ticks.y = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
          )  
    basicR::save_figure(obj = p,
                        filename = file.path(output_dir, "./12_drugs_MOA"),
                        width = 7.5,
                        height = 6,
                        res = 600,
                        formats = formats
                       )    
        
    return(result)
  
}
