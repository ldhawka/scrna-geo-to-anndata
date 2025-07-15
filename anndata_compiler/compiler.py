import pandas as pd
import os
from os import listdir
import anndata
import scanpy as sc
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from tqdm import tqdm
import gc
import random


class GEOAnndataCompiler:
    """
    Reusable template for compiling AnnData objects from GEO raw data.
    """
    
    def __init__(self, config):
        """
        Initialize with:
        - raw_data_dir: path to GEO RAW directory
        - metadata_file: path to metadata file (CSV or Excel)
        - output_file: path for output h5ad file
        - sample_id_column: column name for sample IDs in metadata
        - max_cells_per_sample: maximum cells to sample per file
        - target_sum: normalization target sum
        - n_top_genes: number of highly variable genes
        - delimiter: file delimiter (default: whitespace)
        - random_state: random seed for reproducibility
        """
        self.config = config
        self.adata_list = []
        self.metadata_df = None
        
        # Set defaults
        self.config.setdefault('max_cells_per_sample', 750)
        self.config.setdefault('target_sum', 1e4)
        self.config.setdefault('n_top_genes', 3000)
        self.config.setdefault('delimiter', 'whitespace')
        self.config.setdefault('random_state', 42)
        self.config.setdefault('optimize_params', False)
        
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
    
    def extract_sample_id(self, filename):
        """Extract sample ID from filename. Override this method for custom extraction."""
        return filename.split('_')[0]
    
    def get_sample_metadata(self, sample_id):
        """Get metadata for a specific sample ID."""
        sample_id_col = self.config['sample_id_column']
        matching_rows = self.metadata_df.loc[self.metadata_df[sample_id_col] == sample_id]
        
        if matching_rows.empty:
            return None
        
        return matching_rows.iloc[0].to_dict()
    
    def process_sample_file(self, filepath):
        """Process a single sample file into AnnData format."""
        if self.config['delimiter'] == 'whitespace':
            data = pd.read_csv(filepath, delim_whitespace=True).T
        else:
            data = pd.read_csv(filepath, delimiter=self.config['delimiter']).T
        
        # Set gene names as columns
        data.rename(columns=data.iloc[0], inplace=True)
        data = data.iloc[1:]
        
        # Sample subset if necessary
        max_cells = self.config['max_cells_per_sample']
        if len(data) > max_cells:
            data = data.sample(n=max_cells, random_state=self.config['random_state'])
        
        return data
    
    def create_anndata(self, data, metadata):
        """Create AnnData object with metadata."""
        adata = sc.AnnData(data)
        
        # Add metadata to observations
        for key, value in metadata.items():
            adata.obs[key] = value
        
        return adata
    
    def process_all_samples(self):
        """Process all samples in the raw data directory."""
        raw_dir = self.config['raw_data_dir']
        
        if not os.path.exists(raw_dir):
            raise FileNotFoundError(f"Raw data directory not found: {raw_dir}")
        
        # Load metadata
        self.load_metadata()
        
        # Get all files
        file_names = listdir(raw_dir)
        print(f"Found {len(file_names)} files to process")
        
        processed_count = 0
        
        for filename in tqdm(file_names, desc="Processing samples"):
            try:
                # Extract sample ID
                sample_id = self.extract_sample_id(filename)
                
                # Get metadata
                metadata = self.get_sample_metadata(sample_id)
                if metadata is None:
                    print(f"No metadata found for sample: {sample_id}")
                    continue
                
                # Process file
                filepath = os.path.join(raw_dir, filename)
                data = self.process_sample_file(filepath)
                
                # Create AnnData
                adata = self.create_anndata(data, metadata)
                self.adata_list.append(adata)
                processed_count += 1
                
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue
        
        print(f"Successfully processed {processed_count} samples")
        return self.adata_list
    
    def merge_datasets(self):
        """Merge all AnnData objects."""
        if len(self.adata_list) == 0:
            raise ValueError("No samples were processed successfully")
        
        print("Merging datasets...")
        adata_merged = anndata.concat(self.adata_list)
        
        # Clear memory
        del self.adata_list
        gc.collect()
        
        # Convert data types
        adata_merged.X = adata_merged.X.astype(float)
        
        # Convert categorical columns to string
        for col in adata_merged.obs.columns:
            if adata_merged.obs[col].dtype == 'object':
                adata_merged.obs[col] = adata_merged.obs[col].astype(str)
        
        print(f"Merged data shape: {adata_merged.shape}")
        return adata_merged
    
    def optimize_parameters(self, adata):
        """Optimize PCA and neighbor parameters."""
        print("Optimizing parameters...")
        
        # PCA analysis
        max_comps = min(100, min(adata.shape) - 1)
        sc.pp.pca(adata, svd_solver='arpack', n_comps=max_comps)
        variance_ratio = adata.uns['pca']['variance_ratio']
        cumulative_variance_ratio = np.cumsum(variance_ratio)
        
        # Find optimal PCs (90% variance, fallback to 80% or elbow point)
        pcs_90 = np.where(cumulative_variance_ratio >= 0.9)[0]
        pcs_80 = np.where(cumulative_variance_ratio >= 0.8)[0]
        
        if len(pcs_90) > 0:
            optimal_pcs = pcs_90[0] + 1
        elif len(pcs_80) > 0:
            optimal_pcs = pcs_80[0] + 1
        else:
            # Use elbow method - find where slope changes significantly
            diffs = np.diff(variance_ratio)
            optimal_pcs = np.where(diffs > np.mean(diffs) + np.std(diffs))[0]
            optimal_pcs = optimal_pcs[0] + 1 if len(optimal_pcs) > 0 else max_comps // 2
        
        optimal_pcs = min(optimal_pcs, 50)  # Cap at 50
        optimal_pcs = max(optimal_pcs, 10)  # Minimum 10
        
        print(f"Variance explained with {optimal_pcs} PCs: {cumulative_variance_ratio[optimal_pcs-1]:.3f}")
        
        # Optimize neighbors
        pca_matrix = adata.obsm['X_pca'][:, :optimal_pcs]
        scaler = StandardScaler()
        pca_scaled = scaler.fit_transform(pca_matrix)
        
        # Adjust k_range based on dataset size
        max_k = min(50, adata.shape[0] // 10)
        k_range = range(5, max_k + 1, 5)
        avg_distances = []
        
        for k in k_range:
            if k >= adata.shape[0]:
                break
            nbrs = NearestNeighbors(n_neighbors=k).fit(pca_scaled)
            distances, _ = nbrs.kneighbors(pca_scaled)
            avg_distances.append(np.mean(distances[:, -1]))
        
        # Find elbow point
        if len(avg_distances) > 1:
            diffs = np.diff(avg_distances)
            optimal_neighbors_idx = np.where(diffs > np.mean(diffs) + np.std(diffs))[0]
            optimal_neighbors = list(k_range)[optimal_neighbors_idx[0]] if len(optimal_neighbors_idx) > 0 else 15
        else:
            optimal_neighbors = 15
        
        optimal_neighbors = max(optimal_neighbors, 5)  # Minimum 5
        
        print(f"Optimal parameters: n_pcs={optimal_pcs}, n_neighbors={optimal_neighbors}")
        return optimal_pcs, optimal_neighbors
    
    def preprocess_data(self, adata):
        """Standard preprocessing pipeline."""
        print("Preprocessing data...")
        
        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
        # Normalization
        sc.pp.normalize_total(adata, target_sum=self.config['target_sum'])
        sc.pp.log1p(adata)
        
        # Highly variable genes
        sc.pp.highly_variable_genes(adata, flavor='seurat', 
                                  n_top_genes=self.config['n_top_genes'], 
                                  subset=False)
        
        print(f"Selected {np.sum(adata.var['highly_variable'])} highly variable genes")
        return adata
    
    def run_analysis(self, adata):
        """Run clustering and dimensionality reduction."""
        print("Running analysis...")
        
        if self.config['optimize_params']:
            optimal_pcs, optimal_neighbors = self.optimize_parameters(adata)
        else:
            optimal_pcs = 40
            optimal_neighbors = 15
        
        # Use only highly variable genes for downstream analysis
        if 'highly_variable' in adata.var.columns:
            adata_hvg = adata[:, adata.var['highly_variable']].copy()
        else:
            adata_hvg = adata.copy()
        
        # PCA
        sc.tl.pca(adata_hvg, svd_solver='arpack')
        
        # Neighbors and clustering
        sc.pp.neighbors(adata_hvg, n_neighbors=optimal_neighbors, n_pcs=optimal_pcs)
        sc.tl.leiden(adata_hvg)
        sc.tl.umap(adata_hvg)
        
        # Transfer results back to original object
        adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
        adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
        adata.obs['leiden'] = adata_hvg.obs['leiden']
        adata.uns.update(adata_hvg.uns)
        adata.obsp.update(adata_hvg.obsp)
        
        return adata
    
    def generate_plots(self, adata, color_by=None):
        """Generate standard plots."""
        if color_by is None:
            color_by = ['leiden']
        
        print("Generating plots...")
        for color in color_by:
            if color in adata.obs.columns:
                sc.pl.umap(adata, color=color, save=f'_{color}.pdf')
    
    def save_data(self, adata):
        """Save the processed AnnData object."""
        output_file = self.config['output_file']
        print(f"Saving data to {output_file}")
        adata.write(output_file)
    
    def run_full_pipeline(self, plot_colors=None):
        """Run the complete processing pipeline."""
        try:
            # Process all samples
            self.process_all_samples()
            
            # Merge datasets
            adata_merged = self.merge_datasets()
            
            # Preprocess
            adata_merged = self.preprocess_data(adata_merged)
            
            # Run analysis
            adata_merged = self.run_analysis(adata_merged)
            
            # Generate plots
            if plot_colors:
                self.generate_plots(adata_merged, plot_colors)
            
            # Save data
            self.save_data(adata_merged)
            
            print("Pipeline completed successfully!")
            return adata_merged
            
        except Exception as e:
            print(f"Pipeline failed: {str(e)}")
            raise


def create_config_template():
    """Create a template configuration dictionary."""
    return {
        'raw_data_dir': './GSE123456_RAW',
        'metadata_file': './metadata.csv',
        'output_file': './processed_data.h5ad',
        'sample_id_column': 'sample_id',
        'max_cells_per_sample': 750,
        'target_sum': 1e4,
        'n_top_genes': 3000,
        'delimiter': 'whitespace',
        'random_state': 42,
        'optimize_params': False
    }


# Example usage
if __name__ == "__main__":
    # Create configuration
    config = {
        'raw_data_dir': './GSE244263_RAW',
        'metadata_file': './GSE244263_metadata.csv',
        'output_file': './compiled_data.h5ad',
        'sample_id_column': 'Sample_ID',
        'max_cells_per_sample': 750,
        'target_sum': 1e4,
        'n_top_genes': 3000,
        'optimize_params': True
    }
    
    # Initialize compiler
    compiler = GEOAnndataCompiler(config)
    
    # Run pipeline
    adata = compiler.run_full_pipeline(plot_colors=['leiden', 'Sample_ID', 'Disease_state'])