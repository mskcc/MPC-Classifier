# MPC Classifier

A command-line tool for classifying patients from a tumor registry as having a **Single Primary Cancer (SPC)** or **Multiple Primary Cancers (MPC)**, based on ICD-O-3 morphology and ICD-10 topology codes.

---

## Dependencies

| Package   | Version |
|-----------|---------|
| Python    | ≥ 3.8   |
| pandas    | 3.0.3   |
| numpy     | 2.4.6   |
| PyYAML    | 6.0.2   |

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Input

A tab-separated file (`.tsv` or `.txt`) with one row per tumor record. The following columns are required:


| Column | Type | Description |
|--------|------|-------------|
| `PATIENT_ID` | string | Unique patient identifier |
| `START_DATE` | integer | Numeric date offset — must be consistent within patients to allow ordering and time difference calculations (e.g. days from a reference date) |
| `DX_DESCRIPTION` | string | Tumor registry diagnosis string containing ICD-O-3 morphology and ICD-10 topology codes in the format: `DESCRIPTION (M{morphology}/{behavior} \| C{topology})` |
| `STAGE_CDM_DERIVED_GRANULAR` | numeric | Tumour stage — used for colon/rectum and prostate/bladder deduplication logic |
| `LATERALITY` | integer | Side of paired organ (0 = not applicable/unknown, 1 = right, 2 = left) 


Optional columns (used if present, ignored if absent):

| Column | Type | Description |
|--------|------|-------------|
| `SUMMARY` | string | Free-text summary field passed through to output


See `toy_dataset.txt` for a minimal working example.

---

## Configuration

All classification logic is driven by `mpc_config.yaml`, which must be present in the **working directory** when the script is run. The config defines:


## How to Run

```bash
cd mpc_classifer

python mpc_classifier.py <input_file> <output_label_name> [output_dir]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `input_file` | Yes | Path to input `.tsv` file |
| `run_name` | Yes | Label used in output filenames |
| `output_dir` | No | Output directory (default: `outputs/`) |

**Example:**
```bash
python mpc_classifier.py my_cohort.tsv cohort_v1 results/
```

---

## Output

Four files are written to the output directory, all date-stamped:

| File | Description |
|------|-------------|
| `MPC_{run_name}_{date}.tsv` | One row per patient with ≥2 distinct primary cancers |
| `SP_{run_name}_{date}.tsv` | One row per patient with exactly 1 primary cancer |
| `NO_PC_{run_name}_{date}.txt` | Patient IDs excluded due to no valid tumor records |
| `start_date_exceptions.txt` | Patient IDs where a time-based exception was applied |


See the `outputs` directory for minimal outputs.


# MPC-PREDICT Models

Code for running a mock example of second primary cancer risk prediction (breast-ovary),
contained in the `mpc_predict_models/` directory. It applies a pre-fitted
penalized Fine–Gray model (saved coefficients) to toy data
to demonstrate a simple application of the models. For the full methodology, refer to the Materials and Methods and Supplementary Methods of the manuscript.
Models for the remaining cancer pairs are available from the authors upon reasonable request.

## Dependencies

Requires **R** and the **tidyverse** package.

Verify R is installed:

​```bash
R --version
​```

Install tidyverse if missing (run the top of the `.Rmd` block)

## Feature Set

The model uses the following features, derived from each patient's first primary
cancer record and linked germline/registry data. All predictors are measured at
or before the model landmark time (within one year of first cancer diagnosis).

| Feature | Type | Description |
|---------|------|-------------|
| `Age` | numeric | Age (years) at diagnosis of the first primary cancer |
| `BRCA1` | binary | BRCA1 pathogenic variant status (1 = carrier, 0 = non-carrier) |
| `Hormone Therapy` | binary | Received hormone therapy prior to t₀ (1 = yes, 0 = no) |
| `Yost Index` | numeric | Yost index (neighborhood socioeconomic status) |
| `Breast Adenocarcinoma` | binary | First primary is a breast adenocarcinoma (1 = yes, 0 = no) |
| `Ethnicity` | categorical | Patient ethnicity Latino/White/Unknown |
| `Stage` | binary | First primary diagnosed at a late stage (1 = III/IV, 0 = I–II) |
| `Breast PRS` | numeric | Breast cancer polygenic risk score (z-standardized) |
| `Ovarian PRS` | numeric | Ovarian cancer polygenic risk score (z-standardized) |
| `Prostate PRS` | numeric | Prostate cancer polygenic risk score (z-standardized) |
| `Pancreas PRS` | numeric | Pancreatic cancer polygenic risk score (z-standardized) |
| `Endometrial PRS` | numeric | Endometrial cancer polygenic risk score (z-standardized) |
| `Testicular PRS` | numeric | Testicular cancer polygenic risk score (z-standardized) |
| `Colorectal PRS` | numeric | Colorectal cancer polygenic risk score (z-standardized) |
| `Cervical PRS` | numeric | Cervical cancer polygenic risk score (z-standardized) |
| `Melanoma PRS` | numeric | Melanoma polygenic risk score (z-standardized) |
