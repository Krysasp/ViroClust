# ViroClust

**Auto-align sequences and generate consensus based on cluster files**

ViroClust is a comprehensive bioinformatics pipeline for viral sequence analysis. It automatically aligns sequences based on cluster definitions (e.g., from CD-HIT) and generates high-quality consensus sequences accepting IUPAC ambiguity codes, gap filtering, and quality scoring.

## Features

- **Modular Architecture**: Clean separation of concerns with reusable components
- **Smart Alignment**: Reference sequence prioritization and extended region handling
- **Quality Scoring**: Position-level quality metrics with configurable thresholds
- **Checkpoint/Resume**: Automatic progress saving for long-running analyses
- **Progress Tracking**: Real-time progress bars and status updates
- **Comprehensive Output**: Multiple consensus types, detailed reports, and visualization
- **Flexible Configuration**: JSON-based configuration with command-line overrides

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/yourusername/ViroClust.git
cd ViroClust

# Install in development mode
pip install -e .

# Or install with development dependencies
pip install -e .[dev]
```

### Direct Usage

```bash
python3 bin/viroclust.py -f input.fasta -c clusters.clstr -o output_dir
```

### With Executable Script

```bash
# Make executable (if not already)
chmod +x bin/viroclust

# Run using the executable
./bin/viroclust -f input.fasta -c clusters.clstr -o output_dir
```

## Quick Start

```bash
# Basic usage with defaults
python3 bin/viroclust.py -f sequences.fasta -c clusters.clstr -o output_dir

# With custom configuration
python3 bin/viroclust.py -f sequences.fasta -c clusters.clstr -o output_dir --config config.json

# Test run with first 3 clusters
python3 bin/viroclust.py -f sequences.fasta -c clusters.clstr -o test_output --test-clusters 3

# With custom thresholds and parallel processing
python3 bin/viroclust.py -f sequences.fasta -c clusters.clstr -o output_dir \
  --threads 16 \
  --conservation 0.95 \
  --inclusion 0.80 \
  --quality-threshold 0.7
```

## Project Structure

```
ViroClust/
├── bin/
│   ├── __init__.py          # Package initialization
│   ├── viroclust            # Executable CLI wrapper
│   ├── viroclust.py         # Main CLI entry point
│   └── example_config.json  # Example configuration
├── src/
│   ├── __init__.py          # Package initialization
│   ├── alignment.py         # MAFFT alignment module
│   ├── checkpoint.py        # Checkpoint/restart management
│   ├── config.py            # Configuration management
│   ├── consensus.py         # Consensus generation
│   ├── disparity.py         # Length disparity detection
│   ├── output.py            # Output file generation
│   ├── parsers.py           # FASTA and cluster parsing
│   ├── progress.py          # Progress tracking
│   ├── quality.py           # Quality scoring
│   └── visualizer.py        # Visualization module
├── tests/
│   ├── test_config.py       # Unit tests
│   ├── test_metadata_10clusters/  # Test data
│   └── test_metadata_missing/     # Test data
├── test_data/
│   ├── sequences.fasta      # Example FASTA file
│   └── clusters.clstr       # Example cluster file
├── setup.py                 # Package installation
├── requirements.txt         # Dependencies
├── README.md                # This file
└── MODULE_IMPROVEMENTS.md   # Detailed improvement documentation
```

## Command-Line Options

### Required Arguments

- `-f, --fasta` - Input FASTA file containing sequences
- `-c, --cluster` - Cluster file (.clstr format, e.g., from CD-HIT)
- `-o, --output` - Output directory for results

### Optional Arguments

| Option | Default | Description |
|--------|---------|-------------|
| `--config` | None | Path to configuration JSON file |
| `--conservation` | 0.90 | Conservation threshold for consensus |
| `--inclusion` | 0.70 | Inclusion threshold for consensus |
| `--threads` | CPU count | Number of parallel threads |
| `--test-clusters` | None | Process only first N clusters (testing) |
| `--test-run` | False | Limit to first 10 sequences (testing) |
| `--filter-gaps` | True | Enable gap filtering |
| `--no-filter-gaps` | False | Disable gap filtering |
| `--gap-threshold` | 0.10 | Gap threshold for position removal |
| `--length-disparity-threshold` | 0.20 | Length disparity threshold |
| `--quality-threshold` | 0.5 | Minimum quality score |

## Configuration

ViroClust uses JSON-based configuration. Create a config file or use command-line arguments.

### Example Configuration (config.json)

```json
{
  "conservation": 0.90,
  "inclusion": 0.70,
  "threads": 16,
  "filter_gaps": true,
  "gap_threshold": 0.10,
  "length_disparity_threshold": 0.20,
  "timeout_base": 120,
  "timeout_per_base": 5000,
  "max_timeout": 1800,
  "quality_threshold": 0.5,
  "enable_progress_bars": true,
  "enable_checkpointing": true,
  "checkpoint_interval": 10
}
```

### Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `conservation` | 0.90 | Minimum frequency of dominant base (0.0-1.0) |
| `inclusion` | 0.70 | Minimum fraction of sequences for consensus (0.0-1.0) |
| `threads` | CPU count | Number of parallel processing threads |
| `filter_gaps` | true | Enable gap position filtering |
| `gap_threshold` | 0.10 | Fraction of sequences with gaps to trigger removal |
| `length_disparity_threshold` | 0.20 | Ratio threshold for extended region handling |
| `quality_threshold` | 0.5 | Minimum acceptable quality score (0.0-1.0) |
| `enable_checkpointing` | true | Enable checkpoint/restart capability |
| `checkpoint_interval` | 10 | Number of clusters between checkpoints |

## Output Files

### FASTA Files

| File | Description |
|------|-------------|
| `aligned_sequences.fasta` | All aligned sequences with cluster prefixes |
| `consensus_sequences.fasta` | Gap-filtered consensus sequences |
| `iupac_consensus_sequences.fasta` | IUPAC consensus sequences |
| `degap_consensus_sequences.fasta` | Degap-filtered consensus sequences |
| `renamed_sequences.fasta` | Sequences with standardized names |

### CSV Reports

| File | Description |
|------|-------------|
| `gap_filter_analysis.csv` | Gap filtering statistics |
| `iupac_analysis.csv` | IUPAC base analysis by cluster |
| `oligo_analysis.csv` | Oligo design impact analysis |
| `degap_analysis.csv` | Degap filtering analysis |
| `length_disparity_analysis.csv` | Length disparity statistics |
| `quality_summary.csv` | Quality scores per cluster |
| `sequence_metadata.csv` | Sequence metadata information |

### Text Reports

| File | Description |
|------|-------------|
| `gap_analysis_summary.txt` | Gap analysis summary |
| `missing_sequence_report.txt` | Sequences not found in FASTA |
| `quality_reports.txt` | Detailed quality reports |
| `visualization_report.txt` | Alignment visualizations |
| `checkpoint_summary.txt` | Checkpoint status summary |
| `sequence_rename_map.txt` | Original to renamed sequence mapping |
| `viroclust_config.json` | Final configuration used |

## Key Features Explained

### 1. Quality Scoring

Position-level quality scores (0.0-1.0) based on:
- Agreement percentage with consensus
- IUPAC ambiguity level
- Gap frequency

**Quality Categories**:
- **High (≥0.9)**: Strong agreement, single base or low ambiguity
- **Medium (≥0.7)**: Moderate agreement, 2-way ambiguity
- **Low (≥0.5)**: Weak agreement, 3-way ambiguity
- **Very Low (<0.5)**: Poor agreement, high ambiguity or gaps

### 2. Extended Region Handling

For clusters with significant length disparity (≥20% ratio):
- Reference sequences (marked with '*') are identified and placed first
- Extended regions (leading/trailing) use per-position conservation
- Bypass standard inclusion threshold for better coverage
- Separate quality scoring for extended regions

### 3. Checkpoint/Resume

- Automatic checkpoint every N clusters (configurable)
- Resume interrupted runs from last checkpoint
- Atomic writes to prevent corruption
- Detailed checkpoint summary reports

### 4. Progress Tracking

- Real-time progress bars using tqdm
- Cluster-by-cluster status updates
- Graceful degradation if tqdm unavailable
- Progress persistence across restarts

## Input Format

### FASTA File

Standard FASTA format with descriptive headers:

### Cluster File (.clstr)

CD-HIT cluster format:
```
>Cluster 0
0	36511aa, >PQ164815|Human... *
>Cluster 1
0	35927aa, >PX778758|Human... at 99.16%
1	35947aa, >PX778759|Human... at 99.04%
```

The `*` marker indicates the reference sequence.

## Dependencies

### Core Dependencies

- **tqdm** (≥4.62.0) - Progress bars

### Optional Dependencies

- **biopython** (≥1.80) - Advanced sequence analysis
- **numpy** (≥1.20) - Numerical operations

### Development Dependencies

- **pytest** (≥6.0) - Testing framework
- **pytest-cov** (≥2.0) - Test coverage
- **black** (≥21.0) - Code formatting
- **flake8** (≥3.8) - Linting

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=src

# Run specific test file
pytest tests/test_config.py -v
```

## GitHub Deployment

### Setup Steps

1. **Initialize Git Repository** (if not already done):
   ```bash
   cd ViroClust
   git init
   ```

2. **Add Files**:
   ```bash
   git add .
   git commit -m "Initial ViroClust release v1.0.0"
   ```

3. **Create .gitignore**:
   ```bash
   # Create .gitignore file
   cat > .gitignore << EOF
   # Python
   __pycache__/
   *.py[cod]
   *.egg-info/
   dist/
   build/

   # Output directories
   output_*/
   tests/output_*/

   # IDE
   .idea/
   .vscode/
   *.swp
   *.swo

   # Test data (optional - uncomment if needed)
   # test_data/*.fasta
   # test_data/*.clstr
   EOF
   ```

4. **Create GitHub Repository**:
   ```bash
   # Create repository on GitHub at https://github.com/yourusername/ViroClust
   git remote add origin https://github.com/yourusername/ViroClust.git
   git branch -M main
   git push -u origin main
   ```

5. **Optional: Set up GitHub Actions CI/CD**:
   Create `.github/workflows/ci.yml`:
   ```yaml
   name: CI/CD

   on: [push, pull_request]

   jobs:
     test:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v2
         - name: Set up Python
           uses: actions/setup-python@v2
           with:
             python-version: '3.8'
         - name: Install dependencies
           run: |
             pip install -e .[dev]
         - name: Run tests
           run: pytest tests/ -v --cov=src
   ```

## Citation

```bibtex
@software{viroclust2026,
  author = {IHCM/NGS},
  title = {ViroClust: Auto-align sequences and generate consensus based on cluster files},
  year = {2026},
  url = {https://github.com/Krysasp/ViroClust},
  version = {1.0.0}
}
```


## Changelog

### v1.0.0 (2026-06-13)

- **Initial release** with full modular architecture
- **Core features**: Alignment, consensus generation, quality scoring
- **Enhanced modules**: Checkpoint/restart, progress tracking, visualization
- **Comprehensive output**: Multiple consensus types, detailed reports
- **Flexible configuration**: JSON-based with CLI overrides
- **Test suite**: Unit tests for core modules

---

**Built for viral sequence analysis**
