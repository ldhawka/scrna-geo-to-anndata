import pandas as pd
import os
from os import listdir
import glob as glob_module
import anndata
import scanpy as sc
import numpy as np
import scipy.sparse as sp
from tqdm import tqdm
import gc
import random
import warnings


class GEOAnndataCompiler:
    """
    Compile AnnData objects from GEO raw data.

    Handles three data formats:
    - 'simple': One expression file per sample (CSV/TSV)
    - 'with_cell_metadata': Paired count + cell metadata files per sample
    - '10x_mtx': Cell Ranger MTX format (matrix.mtx + barcodes.tsv + features.tsv)

    The pipeline stops after QC filtering, normalization, and HVG detection.
    Downstream analysis (PCA, UMAP, Leiden) should be performed interactively.
    See docs/downstream_analysis_guide.md.
    """

    def __init__(self, config):
        """
        Initialize with configuration dictionary.

        Required keys:
            raw_data_dir: path to GEO RAW directory
            metadata_file: path to metadata file (CSV or Excel)
            output_file: path for output h5ad file
            sample_id_column: column name for sample IDs in metadata

        Optional keys:
            max_cells_per_sample: max cells per sample (None = no subsampling, default: None)
            target_sum: normalization target (default: 1e4)
            n_top_genes: HVGs to detect (default: 3000)
            delimiter: file delimiter for CSV/TSV formats (default: 'whitespace')
            random_state: seed (default: 42)
            data_format: 'auto', 'simple', 'with_cell_metadata', '10x_mtx' (default: 'auto')
            min_genes: min genes per cell for QC filtering (default: 200)
            min_cells: min cells per gene for QC filtering (default: 3)
            max_mito_pct: max mitochondrial % for QC filtering (default: 20.0)
            mito_prefix: mitochondrial gene prefix (default: 'MT-' for human, 'mt-' for mouse)
            counts_pattern: glob for count files (default: None, auto-detect)
            cell_metadata_pattern: glob for cell metadata files (default: None)
            metadata_columns: list of metadata columns to include (default: None = all)
        """
        self.config = config
        self.adata_list = []
        self.metadata_df = None

        # Set defaults
        self.config.setdefault('max_cells_per_sample', None)
        self.config.setdefault('target_sum', 1e4)
        self.config.setdefault('n_top_genes', 3000)
        self.config.setdefault('delimiter', 'whitespace')
        self.config.setdefault('random_state', 42)

        # Data format handling
        self.config.setdefault('data_format', 'auto')
        self.config.setdefault('counts_pattern', None)
        self.config.setdefault('cell_metadata_pattern', None)
        self.config.setdefault('metadata_columns', None)

        # QC filtering thresholds
        self.config.setdefault('min_genes', 200)
        self.config.setdefault('min_cells', 3)
        self.config.setdefault('max_mito_pct', 20.0)
        self.config.setdefault('mito_prefix', 'MT-')

        # Set random seeds
        random.seed(self.config['random_state'])
        np.random.seed(self.config['random_state'])

    def load_metadata(self):
        """Load metadata from CSV or Excel file."""
        metadata_file = self.config['metadata_file']

        if metadata_file.endswith('.xlsx') or metadata_file.endswith('.xls'):
            self.metadata_df = pd.read_excel(metadata_file)
        elif metadata_file.endswith('.csv'):
            self.metadata_df = pd.read_csv(metadata_file)
        else:
            raise ValueError("Metadata file must be CSV or Excel format")

        print(f"Loaded metadata with {len(self.metadata_df)} rows")
        return self.metadata_df

    def detect_data_format(self):
        """Auto-detect data format: simple, with_cell_metadata, or 10x_mtx."""
        raw_dir = self.config['raw_data_dir']
        files = listdir(raw_dir)

        # Check for 10x MTX format — subdirectory layout
        subdirs_with_mtx = []
        for item in files:
            item_path = os.path.join(raw_dir, item)
            if os.path.isdir(item_path):
                subfiles = listdir(item_path)
                if any('matrix.mtx' in f for f in subfiles):
                    subdirs_with_mtx.append(item)

        if subdirs_with_mtx:
            print(f"Detected data format: 10x_mtx (subdirectory layout, {len(subdirs_with_mtx)} samples)")
            return '10x_mtx'

        # Check for 10x MTX format — prefixed flat file layout
        has_mtx = any('matrix.mtx' in f for f in files)
        has_barcodes = any('barcodes.tsv' in f for f in files)
        has_features = any('features.tsv' in f or 'genes.tsv' in f for f in files)

        if has_mtx and has_barcodes and has_features:
            print("Detected data format: 10x_mtx (prefixed file layout)")
            return '10x_mtx'

        # Check for with_cell_metadata format
        has_counts = any('_counts' in f or '_expression' in f for f in files)
        has_cell_metadata = any('_metadata' in f or '_barcodes' in f for f in files)

        if has_counts and has_cell_metadata:
            print("Detected data format: with_cell_metadata")
            if not self.config['counts_pattern']:
                if any('_counts.csv' in f for f in files):
                    self.config['counts_pattern'] = '*_counts.csv*'
                elif any('_expression' in f for f in files):
                    self.config['counts_pattern'] = '*_expression*'

            if not self.config['cell_metadata_pattern']:
                if any('_metadata.csv' in f for f in files):
                    self.config['cell_metadata_pattern'] = '*_metadata.csv*'
                elif any('_barcodes' in f for f in files):
                    self.config['cell_metadata_pattern'] = '*_barcodes*'

            return 'with_cell_metadata'

        print("Detected data format: simple")
        return 'simple'

    def extract_sample_id(self, filename):
        """Extract sample ID from filename. Override this method for custom extraction."""
        return filename.split('_')[0]

    def get_sample_metadata(self, sample_id):
        """Get metadata for a specific sample ID."""
        sample_id_col = self.config['sample_id_column']
        matching_rows = self.metadata_df.loc[self.metadata_df[sample_id_col] == sample_id]

        if matching_rows.empty:
            return None

        metadata_dict = matching_rows.iloc[0].to_dict()

        if self.config['metadata_columns'] is not None:
            columns_to_keep = set(self.config['metadata_columns']) | {sample_id_col}
            metadata_dict = {k: v for k, v in metadata_dict.items() if k in columns_to_keep}

        return metadata_dict

    # -------------------------------------------------------------------------
    # Data readers for each format
    # -------------------------------------------------------------------------

    def process_sample_file(self, filepath):
        """Process a single CSV/TSV sample file into a DataFrame."""
        if self.config['delimiter'] == 'whitespace':
            data = pd.read_csv(filepath, sep=r'\s+').T
        else:
            data = pd.read_csv(filepath, delimiter=self.config['delimiter']).T

        # Set gene names as columns
        data.rename(columns=data.iloc[0], inplace=True)
        data = data.iloc[1:]

        # Subsample if configured
        max_cells = self.config['max_cells_per_sample']
        if max_cells is not None and len(data) > max_cells:
            data = data.sample(n=max_cells, random_state=self.config['random_state'])

        return data

    def create_anndata(self, data, metadata):
        """Create AnnData object with metadata."""
        adata = sc.AnnData(data)
        for key, value in metadata.items():
            adata.obs[key] = value
        return adata

    def load_cell_metadata_file(self, filepath):
        """Load cell metadata from file."""
        if self.config['delimiter'] == 'whitespace':
            cell_meta = pd.read_csv(filepath, sep=r'\s+', index_col=0)
        else:
            cell_meta = pd.read_csv(filepath, delimiter=self.config['delimiter'], index_col=0)
        return cell_meta

    def load_counts_file(self, filepath):
        """Load counts matrix from file."""
        if self.config['delimiter'] == 'whitespace':
            counts = pd.read_csv(filepath, sep=r'\s+', index_col=0)
        else:
            counts = pd.read_csv(filepath, delimiter=self.config['delimiter'], index_col=0)
        return counts

    def process_sample_with_cell_metadata(self, sample_id, counts_file, metadata_file):
        """Process a sample that has both counts and cell metadata files."""
        counts_path = os.path.join(self.config['raw_data_dir'], counts_file)
        counts = self.load_counts_file(counts_path)

        metadata_path = os.path.join(self.config['raw_data_dir'], metadata_file)
        cell_metadata = self.load_cell_metadata_file(metadata_path)

        adata = anndata.AnnData(X=counts.T)

        common_cells = adata.obs.index.intersection(cell_metadata.index)
        if len(common_cells) < len(adata.obs.index):
            print(f"Warning: Only {len(common_cells)}/{len(adata.obs.index)} cells have metadata")

        for col in cell_metadata.columns:
            adata.obs[col] = cell_metadata.loc[adata.obs.index, col] if col in cell_metadata.columns else np.nan

        adata.obs['sample_id'] = sample_id

        sample_metadata = self.get_sample_metadata(sample_id)
        if sample_metadata:
            for key, value in sample_metadata.items():
                adata.obs[key] = value

        # Subsample if configured
        max_cells = self.config['max_cells_per_sample']
        if max_cells is not None and adata.n_obs > max_cells:
            indices = np.random.choice(adata.n_obs, max_cells, replace=False)
            adata = adata[sorted(indices)].copy()

        return adata

    def read_10x_mtx_sample(self, sample_path, sample_id=None):
        """
        Read a single 10x Cell Ranger MTX sample using scanpy.

        Args:
            sample_path: Directory containing matrix.mtx, barcodes.tsv, features.tsv.
            sample_id: Sample identifier. If None, extracted from directory name.

        Returns:
            AnnData object, or None if reading fails.
        """
        if sample_id is None:
            sample_id = os.path.basename(sample_path)

        try:
            adata = sc.read_10x_mtx(
                sample_path,
                var_names='gene_symbols',
                make_unique=True,
                gex_only=True,
            )
        except Exception as e:
            print(f"Error reading 10x MTX for {sample_id}: {e}")
            return None

        # Make cell barcodes unique across samples
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        adata.obs['sample_id'] = sample_id

        # Subsample if configured
        max_cells = self.config['max_cells_per_sample']
        if max_cells is not None and adata.n_obs > max_cells:
            indices = np.random.choice(adata.n_obs, max_cells, replace=False)
            adata = adata[sorted(indices)].copy()

        # Attach sample-level metadata
        sample_metadata = self.get_sample_metadata(sample_id)
        if sample_metadata:
            for key, value in sample_metadata.items():
                adata.obs[key] = value
        else:
            print(f"Warning: No metadata found for sample {sample_id}")

        return adata

    def _group_10x_prefixed_files(self):
        """
        Group 10x MTX files by sample prefix for flat-file GEO layouts.

        Returns:
            dict mapping sample_id -> {'prefix': str, 'matrix': str, 'barcodes': str, 'features': str}
        """
        raw_dir = self.config['raw_data_dir']
        files = listdir(raw_dir)

        mtx_files = [f for f in files if 'matrix.mtx' in f]

        samples = {}
        for mtx_file in mtx_files:
            prefix = mtx_file.split('matrix.mtx')[0]
            sample_id = self.extract_sample_id(mtx_file)

            barcodes = [f for f in files if f.startswith(prefix) and 'barcodes' in f]
            features = [f for f in files if f.startswith(prefix) and ('features' in f or 'genes' in f)]

            if barcodes and features:
                samples[sample_id] = {
                    'prefix': prefix,
                    'matrix': mtx_file,
                    'barcodes': barcodes[0],
                    'features': features[0],
                }
            else:
                print(f"Warning: Incomplete 10x file set for prefix '{prefix}', skipping")

        return samples

    # -------------------------------------------------------------------------
    # Main processing pipeline
    # -------------------------------------------------------------------------

    def process_all_samples(self):
        """Process all samples in the raw data directory."""
        raw_dir = self.config['raw_data_dir']

        if not os.path.exists(raw_dir):
            raise FileNotFoundError(f"Raw data directory not found: {raw_dir}")

        self.load_metadata()

        if self.config['data_format'] == 'auto':
            self.config['data_format'] = self.detect_data_format()

        file_names = listdir(raw_dir)
        print(f"Found {len(file_names)} files/directories in raw data directory")

        processed_count = 0

        if self.config['data_format'] == '10x_mtx':
            # Check for subdirectory layout first
            subdirs = [
                f for f in file_names
                if os.path.isdir(os.path.join(raw_dir, f))
                and any('matrix.mtx' in sf for sf in listdir(os.path.join(raw_dir, f)))
            ]

            if subdirs:
                for subdir in tqdm(subdirs, desc="Processing 10x samples"):
                    try:
                        sample_path = os.path.join(raw_dir, subdir)
                        sample_id = self.extract_sample_id(subdir)
                        adata = self.read_10x_mtx_sample(sample_path, sample_id)
                        if adata is not None:
                            self.adata_list.append(adata)
                            processed_count += 1
                    except Exception as e:
                        print(f"Error processing {subdir}: {e}")
                        continue
            else:
                # Prefixed flat file layout
                sample_groups = self._group_10x_prefixed_files()
                print(f"Found {len(sample_groups)} 10x samples from prefixed files")

                for sample_id, file_info in tqdm(sample_groups.items(), desc="Processing 10x samples"):
                    try:
                        adata = sc.read_10x_mtx(
                            raw_dir,
                            var_names='gene_symbols',
                            make_unique=True,
                            gex_only=True,
                            prefix=file_info['prefix'],
                        )
                        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
                        adata.obs['sample_id'] = sample_id

                        max_cells = self.config['max_cells_per_sample']
                        if max_cells is not None and adata.n_obs > max_cells:
                            indices = np.random.choice(adata.n_obs, max_cells, replace=False)
                            adata = adata[sorted(indices)].copy()

                        sample_metadata = self.get_sample_metadata(sample_id)
                        if sample_metadata:
                            for key, value in sample_metadata.items():
                                adata.obs[key] = value
                        else:
                            print(f"Warning: No metadata found for sample {sample_id}")

                        self.adata_list.append(adata)
                        processed_count += 1
                    except Exception as e:
                        print(f"Error processing {sample_id}: {e}")
                        continue

        elif self.config['data_format'] == 'with_cell_metadata':
            counts_pattern = self.config['counts_pattern'] or '*_counts*'
            counts_files = glob_module.glob(os.path.join(raw_dir, counts_pattern))
            print(f"Found {len(counts_files)} count files matching pattern: {counts_pattern}")

            for counts_file in tqdm(counts_files, desc="Processing samples"):
                try:
                    basename = os.path.basename(counts_file)
                    sample_id = self.extract_sample_id(basename)

                    expected_meta = basename.replace('_counts', '_metadata')
                    metadata_file = expected_meta if expected_meta in file_names else None

                    if not metadata_file:
                        meta_files = glob_module.glob(os.path.join(raw_dir, f"*{sample_id}*metadata*"))
                        if meta_files:
                            metadata_file = os.path.basename(meta_files[0])

                    if metadata_file and os.path.exists(os.path.join(raw_dir, metadata_file)):
                        adata = self.process_sample_with_cell_metadata(sample_id, basename, metadata_file)
                        self.adata_list.append(adata)
                        processed_count += 1
                    else:
                        print(f"Warning: No cell metadata file found for {sample_id}, processing as simple format")
                        metadata = self.get_sample_metadata(sample_id)
                        if metadata:
                            data = self.process_sample_file(counts_file)
                            adata = self.create_anndata(data, metadata)
                            self.adata_list.append(adata)
                            processed_count += 1

                except Exception as e:
                    print(f"Error processing {counts_file}: {str(e)}")
                    continue

        else:
            # Simple format
            for filename in tqdm(file_names, desc="Processing samples"):
                try:
                    sample_id = self.extract_sample_id(filename)

                    metadata = self.get_sample_metadata(sample_id)
                    if metadata is None:
                        print(f"No metadata found for sample: {sample_id}")
                        continue

                    filepath = os.path.join(raw_dir, filename)
                    data = self.process_sample_file(filepath)

                    adata = self.create_anndata(data, metadata)
                    self.adata_list.append(adata)
                    processed_count += 1

                except Exception as e:
                    print(f"Error processing {filename}: {str(e)}")
                    continue

        print(f"Successfully processed {processed_count} samples")
        return self.adata_list

    def merge_datasets(self):
        """Merge all AnnData objects into one."""
        if len(self.adata_list) == 0:
            raise ValueError("No samples were processed successfully")

        print("Merging datasets...")
        adata_merged = anndata.concat(self.adata_list)

        del self.adata_list
        gc.collect()

        # Convert to float32 for memory efficiency (handles both sparse and dense)
        if sp.issparse(adata_merged.X):
            adata_merged.X = adata_merged.X.astype(np.float32)
        else:
            adata_merged.X = adata_merged.X.astype(np.float32)

        # Convert object columns to string
        for col in adata_merged.obs.columns:
            if adata_merged.obs[col].dtype == 'object':
                adata_merged.obs[col] = adata_merged.obs[col].astype(str)

        print(f"Merged data shape: {adata_merged.shape}")
        return adata_merged

    def preprocess_data(self, adata):
        """
        Standard preprocessing: QC filtering, normalization, HVG detection.

        Steps:
            1. Calculate QC metrics (including mitochondrial %)
            2. Filter cells and genes by QC thresholds
            3. Store raw counts in adata.layers['counts']
            4. Normalize (CPM to target_sum + log1p)
            5. Detect highly variable genes (flagged, not subsetted)

        Returns:
            Preprocessed AnnData with raw counts in layers['counts'].
        """
        print("Preprocessing data...")

        # Step 1: QC metrics
        mito_prefix = self.config['mito_prefix']
        adata.var['mt'] = adata.var_names.str.startswith(mito_prefix)

        n_mito = adata.var['mt'].sum()
        if n_mito == 0:
            print(f"  WARNING: No mitochondrial genes found with prefix '{mito_prefix}'. "
                  f"If using mouse data, set mito_prefix='mt-' in config.")

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
        )

        n_cells_before = adata.n_obs
        n_genes_before = adata.n_vars

        # Step 2: QC filtering
        min_genes = self.config['min_genes']
        min_cells = self.config['min_cells']
        max_mito_pct = self.config['max_mito_pct']

        sc.pp.filter_cells(adata, min_genes=min_genes)
        print(f"  Filtered cells by min_genes={min_genes}: {n_cells_before:,} -> {adata.n_obs:,}")

        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(f"  Filtered genes by min_cells={min_cells}: {n_genes_before:,} -> {adata.n_vars:,}")

        n_before_mito = adata.n_obs
        adata = adata[adata.obs['pct_counts_mt'] < max_mito_pct].copy()
        print(f"  Filtered cells by max_mito_pct={max_mito_pct}%: {n_before_mito:,} -> {adata.n_obs:,}")

        # Step 3: Preserve raw counts
        adata.layers['counts'] = adata.X.copy()

        # Step 4: Normalize
        sc.pp.normalize_total(adata, target_sum=self.config['target_sum'])
        sc.pp.log1p(adata)

        # Step 5: HVG detection (flag only, keep all genes)
        sc.pp.highly_variable_genes(
            adata, flavor='seurat',
            n_top_genes=self.config['n_top_genes'],
            subset=False,
        )

        n_hvgs = np.sum(adata.var['highly_variable'])
        print(f"  Flagged {n_hvgs} highly variable genes")
        print(f"  Raw counts preserved in adata.layers['counts']")
        print(f"  Final shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        return adata

    def save_data(self, adata):
        """Save the processed AnnData object."""
        output_file = self.config['output_file']
        print(f"Saving data to {output_file}")
        adata.write(output_file)

    def run_full_pipeline(self, **kwargs):
        """
        Run the data compilation pipeline:
            1. Process all sample files
            2. Merge into single AnnData
            3. QC filtering and preprocessing
            4. Save h5ad

        Returns:
            Preprocessed AnnData, ready for downstream analysis (PCA, UMAP, clustering).
        """
        if 'plot_colors' in kwargs:
            warnings.warn(
                "plot_colors parameter has been removed in v0.2.0. "
                "Plotting is now a downstream step. See docs/downstream_analysis_guide.md",
                DeprecationWarning, stacklevel=2,
            )

        try:
            self.process_all_samples()

            adata_merged = self.merge_datasets()

            adata_merged = self.preprocess_data(adata_merged)

            self.save_data(adata_merged)

            print(f"\nPipeline completed successfully!")
            print(f"Output saved to: {self.config['output_file']}")
            print(f"Final dataset: {adata_merged.n_obs:,} cells x {adata_merged.n_vars:,} genes")
            print(f"\nNext steps: See docs/downstream_analysis_guide.md for PCA, UMAP, and clustering.")

            return adata_merged

        except Exception as e:
            print(f"Pipeline failed: {str(e)}")
            raise


def create_config_template():
    """Create a template configuration dictionary with documented defaults."""
    return {
        # Required
        'raw_data_dir': './GSE123456_RAW',
        'metadata_file': './metadata.csv',
        'output_file': './processed_data.h5ad',
        'sample_id_column': 'sample_id',

        # Data format
        'data_format': 'auto',           # 'auto', 'simple', 'with_cell_metadata', '10x_mtx'

        # Subsampling (None = keep all cells)
        'max_cells_per_sample': None,

        # Preprocessing
        'target_sum': 1e4,
        'n_top_genes': 3000,

        # QC filtering
        'min_genes': 200,                # Min genes per cell
        'min_cells': 3,                  # Min cells per gene
        'max_mito_pct': 20.0,            # Max mitochondrial %
        'mito_prefix': 'MT-',            # 'MT-' for human, 'mt-' for mouse

        # File parsing
        'delimiter': 'whitespace',
        'random_state': 42,

        # Optional (for with_cell_metadata format)
        'counts_pattern': None,
        'cell_metadata_pattern': None,
        'metadata_columns': None,
    }


if __name__ == "__main__":
    config = {
        'raw_data_dir': './GSE244263_RAW',
        'metadata_file': './GSE244263_metadata.csv',
        'output_file': './compiled_data.h5ad',
        'sample_id_column': 'Sample_ID',
        'min_genes': 200,
        'min_cells': 3,
        'max_mito_pct': 20.0,
    }

    compiler = GEOAnndataCompiler(config)
    adata = compiler.run_full_pipeline()
