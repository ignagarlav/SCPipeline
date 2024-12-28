import numpy as np
import pandas as pd
import anndata
import logging
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
import warnings
warnings.filterwarnings('ignore')


@dataclass
class  HVGConfig:
    hvg_flavor: str = 'seurat'
    n_top_genes: int = 3000

@dataclass
class CountsConfig:
    low_percentile: int = 1
    high_percentile: int = 99
    min_threshold: int = 500
    max_threshold: int = 40000

@dataclass
class GenesConfig:
    low_percentile: int = 1
    high_percentile: int = 99
    min_threshold: int = 500
    max_threshold: int = 10000
    min_cells_percent: float = 0.01
    modelos_anotacion: list = None
    mt_pattern: str = '^MT-'

@dataclass 
class CellConfig:
    max_pct_mt: int = 30
    min_complexity: float = 0.1
    expected_doublet_rate: float = 0.1

@dataclass
class NormalizationConfig:
    target_sum: int = 1e4
    max_value: int = 10

@dataclass
class IntegrationConfig:
    batch_key: str = "batch"
    n_pcs: int = 30

@dataclass
class DimReductionConfig:
    resolutions: List[float] = field(default_factory=list)  
    min_distances: List[float] = field(default_factory=list)
    n_neighbors: int = 15
    n_pcs: int = 30
    spread: float = 1.0
    n_components: int = 2
    representation: str = 'X_pca_harmony'
    embedding_keys: str = None
    

@dataclass
class SCDataConfig:
    
    data_dir: str = None
    
    samples: List[str] = field(default_factory=list)
    
    all_percentiles: List[int] = field(default_factory=lambda: [1, 5, 95, 99])
    hvg: Any = field(default_factory=lambda: HVGConfig())
    counts: Any = field(default_factory=lambda: CountsConfig())
    genes: Any = field(default_factory=lambda: GenesConfig())
    integration: Any = field(default_factory=lambda: IntegrationConfig())
    dim_reduction: Any = field(default_factory=lambda: DimReductionConfig())
    cells: Any = field(default_factory=lambda: CellConfig())
    normalization: Any = field(default_factory=lambda: NormalizationConfig())
    
    logger: Any = None



  
class SCData:

    def __init__(self,
                 Adata_config: SCDataConfig, 
                 data_dict: Optional[Dict] = None, 
                 adata: Optional[Any] = None, ) -> None:
        
        self.config = Adata_config
        self._data_dict: Dict = data_dict if data_dict is not None else {}
        self.adata: Optional[Any] = adata
        self.all_percentiles: Dict[str, Dict[str, float]] = {}  
    # Properties  
    
    @property
    def obsm_names(self):
        """Dynamically get obsm names from adata."""
        if self.adata is not None:
            return list(self.adata.obsm.keys())
        return []

    @obsm_names.setter
    def obsm_names(self, value):
        """Prevent direct setting of obsm_names."""
        raise AttributeError("obsm_names is derived from adata and cannot be set directly.")

    @property 
    def layer_names(self):
        """Dynamically get layer names from adata."""
        if self.adata is not None:
            return list(self.adata.layers.keys())
        return []
    @layer_names.setter
    def layer_names(self, value):
        """Prevent direct setting of layer_names."""
        raise AttributeError("layer_names is derived from adata and cannot be set directly.")
    
    @property
    def graph_keys(self):
        """Dynamically get graph keys from adata."""
        if self.adata is not None:
            return list(self.adata.obsp.keys())
        return []
    @graph_keys.setter
    def graph_keys(self, value):
        """Prevent direct setting of graph_keys."""
        raise AttributeError("graph_keys is derived from adata and cannot be set directly.")

    @property
    def adata_list(self) -> List:
        # Generate the adata_list dynamically from data_dict
        return list(self._data_dict.values())

    @adata_list.setter
    def adata_list(self, value):
        raise AttributeError("Cannot set adata_list directly. Update data_dict instead.")

    @property
    def data_dict(self) -> Dict:
        return self._data_dict

    @data_dict.setter
    def data_dict(self, value: Dict):
        self._data_dict = value


    #Preprocessing 
    def create_adata_list(self):
        import scanpy as sc
        import os 
                   
        for sample in self.config.samples:
            sample_dir = os.path.join(self.config.data_dir, sample)
            if self.config.logger:
                self.config.logger.info(f"Loading data for sample: {sample} from {sample_dir}")
            try:     
                adata = sc.read_10x_mtx(
                    sample_dir,
                    var_names='gene_symbols',
                    cache=True
                )
                # Add sample metadata (the name of the sample)
                adata.obs['sample'] = sample
                
                self.data_dict[sample] = adata
            
            except Exception as e:
                if self.config.logger:
                    self.config.logger.error(f"Error loading data for sample '{sample}': {e}")
                raise e
        if self.config.logger:   
            self.config.logger.info(f"Successfully created AnnData objects for {len(self.config.samples)} samples.")

    def compute_qc_metrics(self, **kwargs):
        
        import scanpy as sc

        if self.config.logger: 
            self.config.logger.info("Starting QC metrics computation.")
        
         # Validate that data_dict is populated
        if not self.data_dict:
            raise ValueError("No AnnData objects found in data_dict. Run 'create_adata_list' first.")

        for sample, adata in self.data_dict.items():
            if 'mt' in adata.var.columns:
                if self.config.logger:
                    self.config.logger.warning(f"'mt' column already exists in AnnData for sample '{sample}'. Overwriting.")
        
            # Identify mitochondrial genes
            adata.var['mt'] = adata.var_names.str.match(self.config.genes.      mt_pattern)
            if self.config.logger:
                self.config.logger.info(f"Identified {adata.var['mt'].sum()} mitochondrial genes in {adata.n_vars} total genes.")

            # Compute QC metrics
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], log1p=False, inplace=True)
            adata.obs['complexity'] = adata.obs['n_genes_by_counts'] / adata.obs['total_counts']

            if self.config.logger:
                self.config.logger.info("Computed QC metrics: total counts and % mitochondrial genes.")

    def calculate_percentiles(self, **kwargs):

        '''
        Calculates the percentiles to select tresholds for QC metrics
        Checks whether the all_percentiles attribute exists in the self object (typically an instance of a class). If it does not exist, it initializes self.all_percentiles as an empty list.
        If a percentile_dir is provided, the percentiles are saved as a csv file
        '''
        import numpy as np
        import pandas as pd

        percentiles = self.config.all_percentiles  # [1, 5, 95, 99]
        metrics = ['total_counts', 'n_genes_by_counts', 'total_counts_mt', 'complexity', 'doublet_scores']
        
        if not isinstance(self.all_percentiles, dict):
            self.all_percentiles = {}

        for sample, adata in self.data_dict.items():
            percentiles_muestra: Dict[str, float] = {"Muestra": sample}
            
            for metric in metrics:
                # Check if the metric exists in adata.obs
                if metric not in adata.obs.columns:
                    if self.config.logger:
                        self.config.logger.warning(f"Metric '{metric}' not found in adata.obs for sample '{sample}'. Skipping.")
                    continue
            
                quantile_values = adata.obs[metric].quantile([p / 100.0 for p in percentiles]).to_dict()
                
                for p in percentiles:
                        key = f"{metric} P{p}"
                        percentiles_muestra[key] = quantile_values.get(p / 100.0, np.nan)
                        
            self.all_percentiles[sample] = percentiles_muestra

                    
        percentile_dir = kwargs.get('percentile_dir',None)
        if percentile_dir:
            import os
            os.makedirs(percentile_dir, exist_ok=True)
            percentile_file = os.path.join(percentile_dir, "percentiles.csv")

            # Convert the dictionary of percentiles into a DataFrame for saving
            if self.all_percentiles:
                all_percentiles_df = pd.DataFrame.from_dict(self.all_percentiles, orient='index')
                all_percentiles_df.to_csv(percentile_file, index=False)
                if self.config.logger:
                    self.config.logger.info(f"Saved percentiles to '{percentile_file}'.")
            else:
                if self.config.logger:
                    self.config.logger.warning("No percentiles were calculated. Skipping saving to CSV.")
        else:  
            if self.config.logger:
                self.config.logger.info("Percentiles computed but not saved to CSV.")
            
    def calculate_thresholds(self, **kwargs):
        import pandas as pd
   
        counts_config = self.config.counts
        genes_config = self.config.genes

        for sample, adata in self.data_dict.items():
            
            percentiles_muestra = self.all_percentiles.get(sample, {})
            
            if not percentiles_muestra:
                
                if self.config.logger:
                    self.config.logger.warning(f"No percentiles computed for sample '{sample}'. Skipping threshold computation.")
                continue
        
            # Compute counts thresholds
            min_counts = max(
                percentiles_muestra.get(f"total_counts P{counts_config.low_percentile}", 0),
                counts_config.min_threshold
            )
            max_counts = min(
                percentiles_muestra.get(f"total_counts P{counts_config.high_percentile}", np.inf),
                counts_config.max_threshold
            )

            # Compute genes thresholds
            min_genes = max(
                percentiles_muestra.get(f"n_genes_by_counts P{genes_config.low_percentile}", 0),
                genes_config.min_threshold
            )
            max_genes = min(
                percentiles_muestra.get(f"n_genes_by_counts P{genes_config.high_percentile}", np.inf),
                genes_config.max_threshold
            )

            # Assign thresholds to adata or use them as needed
            adata.uns['thresholds'] = {
            'min_counts': min_counts,
            'max_counts': max_counts,
            'min_genes': min_genes,
            'max_genes': max_genes
            }

            if self.config.logger:
                self.config.logger.info(
                "Computed thresholds: "
                f"min_counts={min_counts:.2f}" 
                f"max_counts={max_counts:.2f}" 
                f"min_genes={min_genes:.2f}" 
                f"max_genes={max_genes:.2f}" 
            )        
        if self.config.logger:        
            self.config.logger.info("Finished calculating dynamic thresholds.")

    def filter_low_quality_cells(self, **kwargs):   
        
        import numpy as np
        import scanpy as sc
        
        for i, (sample, adata) in enumerate(self.data_dict.items()):
            if self.config.logger:
                self.config.logger.info(f"Starting cell filtering for sample '{sample}'.")
            
            # Store the initial number of cells
            initial_cells = adata.n_obs

            # Get dynamic thresholds

            thresholds = adata.uns['thresholds']
            min_counts = int(thresholds['min_counts'])
            max_counts = int(thresholds['max_counts'])
            min_genes = int(thresholds['min_genes'])
            max_genes = int(thresholds['max_genes'])


            # Filter cells based on QC metrics
            if self.config.logger:
               self.config.logger.info(f"Initial number of cells: {adata.n_obs}")
            
            sc.pp.filter_cells(adata, min_genes=min_genes)
            if self.config.logger:
                self.config.logger.info(f"After min_genes filtering: {adata.n_obs} cells remaining.")
            
            sc.pp.filter_cells(adata, max_genes=max_genes)
            if self.config.logger:
                self.config.logger.info(f"After max_genes filtering: {adata.n_obs} cells remaining.")
            
            sc.pp.filter_cells(adata, min_counts=min_counts)
            if self.config.logger:
                self.config.logger.info(f"After min_counts filtering: {adata.n_obs} cells remaining.")
            
            sc.pp.filter_cells(adata, max_counts=max_counts)
            if self.config.logger:
                self.config.logger.info(f"After max_counts filtering: {adata.n_obs} cells remaining.")
            
            adata = adata[adata.obs['pct_counts_mt'] <= self.config.cells.max_pct_mt, :].copy()
            if self.config.logger:
                self.config.logger.info(f"After mitochondrial filtering: {adata.n_obs} cells remaining.")
            
            adata = adata[adata.obs['complexity'] > self.config.cells.min_complexity , :].copy()  
            if self.config.logger:
                self.config.logger.info(f"After complexity filtering: {adata.n_obs} cells remaining.")

            # Log results
            if self.config.logger:
                self.config.logger.info(f"{initial_cells - adata.n_obs} have been removed.")
                self.config.logger.info(f"There are {adata.n_obs} cells remaining.")

            self.data_dict[sample] = adata

    def filter_low_quality_genes(self, **kwargs):
   
        import scanpy as sc
        import pandas as pd
        import numpy as np
        
        for sample, adata in self.data_dict.items():
            if self.config.logger:
                self.config.logger.info("Starting low quality gene filtering...")

            # Store the initial number of genes
            initial_genes = adata.n_vars
            if self.config.logger:
                self.config.logger.info(f"Initial number of genes: {initial_genes}")

            # Calculate the minimum number of cells a gene must be expressed in
            min_cells = int(np.ceil(adata.n_obs * self.config.genes.min_cells_percent))
            if self.config.logger:
                self.config.logger.info(f"Filtering genes expressed in fewer than {min_cells} cells")

            # Apply gene filtering  
            sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)
            sc.pp.filter_genes(adata, min_counts=3, inplace=True)


            # Store the number of genes after filtering
            filtered_genes = adata.n_vars
            if self.config.logger:
                self.config.logger.info(f"Number of genes after filtering: {filtered_genes}")

            # Calculate the number of genes removed
            removed_genes = initial_genes - filtered_genes
            if self.config.logger:
                self.config.logger.info(f"Number of genes removed: {removed_genes}")

          # Save original counts in adata.layers['counts']
            adata.layers["counts"] = adata.X.copy()
            if self.config.logger:
                self.config.logger.info("Saved original counts in adata.layers['counts'].")

    def detect_and_remove_doublets(self, **kwargs):  


        cell_config = self.config.cells
        import scrublet as scr
        if self.config.logger:
            self.config.logger.info("Starting doublet detection with Scrublet...")
        
        for sample, adata in self.data_dict.items():
            
            # Get the number of samples (cells) and features (genes)
            n_samples, n_features = adata.shape

            # Determine the maximum number of PCA components
            n_components = min(30, n_samples, n_features) 
            
            if self.config.logger:
                self.config.logger.info(f"Using n_components={n_components} for PCA.")


            # Initialize Scrublet
            scrub = scr.Scrublet(adata.raw.X if adata.raw is not None else adata.X,                               expected_doublet_rate=cell_config.expected_doublet_rate)
            
            # Perform doublet detection    
            doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=n_components)

            # Assign doublet metrics to adata.obs
            adata.obs['doublet_scores'] = doublet_scores
            adata.obs['predicted_doublets'] = predicted_doublets
            
            # Count doublets
            n_doublets = predicted_doublets.sum()
            if self.config.logger:
                self.config.logger.info(f"Detected {n_doublets} doublets out of {adata.n_obs} cells.")
                self.config.logger.info(f"Doublet rate: {n_doublets / adata.n_obs:.2%}")

            remove = kwargs.get('remove', False)
            if remove:
                # Remove doublets with copy to avoid SettingWithCopyWarning
                adata_filtered = adata[~adata.obs['predicted_doublets'], :].copy()
                if self.config.logger:
                    self.config.logger.info(f"Removed {n_doublets} doublets. {adata_filtered.n_obs} cells remaining.")
                    self.data_dict[sample] = adata_filtered
                
    def normalize(self, **kwargs): 
        
        import scanpy as sc

        adata = kwargs.get('adata', self.adata)
        data_dict = kwargs.get('data_dict', self.data_dict)
        target_sum = kwargs.get('target_sum', self.config.normalization.target_sum)
        extra_name = kwargs.get('extra_name', '')


        if data_dict and not adata:
            if not isinstance(data_dict, dict):
                raise ValueError("data_dict must be a dictionary of AnnData objects.")
            
            for sample, adata in data_dict.items():
                
                if self.config.logger:
                    self.config.logger.info(f"Starting normalization for sample '{sample}'.")
                
                # Normalize total counts per cell
                if self.config.logger:
                    self.config.logger.info("Normalizing total counts.")
                sc.pp.normalize_total(adata, target_sum=target_sum)

                if extra_name:
                    adata.layers[f"norm_{extra_name}"] = adata.X.copy()
                else:
                    if "norm" in adata.layers: 
                        if self.config.logger:
                            self.config.logger.warning("Overwriting existing 'norm' layer.")
                            adata.layers["norm"] = adata.X.copy()
                    else:
                        if self.config.logger:
                            self.config.logger.info("Creating 'norm' layer.")
                            adata.layers["norm"] = adata.X.copy()
                        
                
                # Log1p transformation
                if self.config.logger:
                    self.config.logger.info("Applying log1p transformation.")
                sc.pp.log1p(adata)

                if extra_name:
                    adata.layers[f"log1p_norm_{extra_name}"] = adata.X.copy()
                else: 
                    if "log1p_norm" in adata.layers:
                        if self.config.logger:
                            self.config.logger.warning("Overwriting existing 'log1p_norm' layer.")
                            adata.layers["log1p_norm"] = adata.X.copy()
                    else:
                        if self.config.logger:
                            self.config.logger.info("Creating 'log1p_norm' layer.")
                            adata.layers["log1p_norm"] = adata.X.copy()
                    if self.config.logger:
                        self.config.logger.info(f"Normalization complete. Data shape: {adata.shape} (Cells x Genes).")
        
        elif adata and not data_dict:
            if self.config.logger:
                self.config.logger.info("Starting normalization.")
            
            # Normalize total counts per cell
            if self.config.logger:
                self.config.logger.info("Normalizing total counts.")
            sc.pp.normalize_total(adata, target_sum=target_sum)

            if extra_name:
                    adata.layers[f"norm_{extra_name}"] =  adata.X.copy()
            else:
                if "norm" in adata.layers: 
                    if self.config.logger:
                        self.config.logger.warning("Overwriting existing 'norm' layer.")
                        adata.layers["norm"] = adata.X.copy()
                else:
                    if self.config.logger:
                        self.config.logger.info("Creating 'norm' layer.")
                        adata.layers["norm"] = adata.X.copy()
            
            # Log1p transformation
            if self.config.logger:
                self.config.logger.info("Applying log1p transformation.")
            sc.pp.log1p(adata)

            if extra_name:
                adata.layers[f"log1p_norm_{extra_name}"] = adata.X.copy()
            else:
                if "log1p_norm" in adata.layers:
                    if self.config.logger:
                        self.config.logger.warning("Overwriting existing 'log1p_norm' layer.")
                        adata.layers["log1p_norm"] = adata.X.copy()
                else:
                    if self.config.logger:
                        self.config.logger.info("Creating 'log1p_norm' layer.")

            if self.config.logger:
                self.config.logger.info(f"Normalization complete. Data shape: {adata.shape} (Cells x Genes).")      
        
        elif adata and data_dict:
            raise ValueError("Both AnnData object and data_dict provided. Please provide only one.")
        else:
            raise ValueError("No AnnData object or data_dict provided.")
        
    def identify_hvg(self, **kwargs):
        import scanpy as sc

        adata = kwargs.get('adata', self.adata)
        data_dict = kwargs.get('data_dict', self.data_dict)
        hvg_flavor = kwargs.get('hvg_flavor', self.config.hvg.hvg_flavor)
        n_top_genes = kwargs.get('n_top_genes', self.config.hvg.n_top_genes)
        batch_key = kwargs.get('batch_key', self.config.integration.batch_key)
        what_to_process = kwargs.get('process', 'both') 

        def hvgs_multiple(data_dict, hvg_flavor, n_top_genes, batch_key):  
                for sample, adata in data_dict.items():
                    if self.config.logger:
                        self.config.logger.info(f"Identifying top {n_top_genes} highly variable genes using flavor='{hvg_flavor}' for sample '{sample}'.")
                        
                    sc.pp.highly_variable_genes(
                        adata,
                        n_top_genes=n_top_genes,
                        subset=False,
                        flavor=hvg_flavor,
                        batch_key=batch_key
                    )
        def hvgs_single(adata, hvg_flavor, n_top_genes, batch_key):
                if self.config.logger:
                    self.config.logger.info(f"Identifying top {n_top_genes} highly variable genes using flavor='{hvg_flavor}'.")
            
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes,
                    subset=False,
                    flavor=hvg_flavor,
                    batch_key=batch_key
                )

        if data_dict and not adata:
            if not isinstance(data_dict, dict):
                raise ValueError("data_dict must be a dictionary of AnnData objects.")
            
            
            hvgs_multiple(data_dict, hvg_flavor, n_top_genes, batch_key)
            
        elif adata and not data_dict:
            
            hvgs_single(adata, hvg_flavor, n_top_genes, batch_key)
        
        elif adata and data_dict:
            if what_to_process == 'both':
                if self.config.logger:
                    self.config.logger.info("Identifying highly variable genes for both AnnData object and data_dict.")
                hvgs_single(adata, hvg_flavor, n_top_genes, batch_key)
                hvgs_multiple(data_dict, hvg_flavor, n_top_genes, batch_key)
            elif what_to_process == 'adata':
                if self.config.logger:
                    self.config.logger.info("Identifying highly variable genes for AnnData object.")
                hvgs_single(adata, hvg_flavor, n_top_genes, batch_key)
            elif what_to_process == 'data_dict':
                if self.config.logger:
                    self.config.logger.info("Identifying highly variable genes for data_dict.")
                hvgs_multiple(data_dict, hvg_flavor, n_top_genes, batch_key)
            
        else:
            raise ValueError("No AnnData object or data_dict provided.")

    def scale_data(self, **kwargs):
        import scanpy as sc

        adataprocessed = kwargs.get('adata', self.adata)
        data_dict = kwargs.get('data_dict', self.data_dict)
        max_value = kwargs.get('max_value', self.config.normalization.max_value)
        extra_name = kwargs.get('extra_name', '')

        if data_dict and not adataprocessed:
            if not isinstance(data_dict, dict):
                raise ValueError("data_dict must be a dictionary of AnnData objects.")
            
            for sample, adata in data_dict.items():
                if self.config.logger:
                    self.config.logger.info(f"Scaling data for sample '{sample}'.")
                sc.pp.scale(adata, max_value=max_value)

                if extra_name:
                    adata.layers[f"scaled_{extra_name}"] = adata.X.copy()
                
                else:
                    if "scaled" in adata.layers:
                        if self.config.logger:
                            self.config.logger.warning("Overwriting existing 'scaled' layer.")
                            adata.layers["scaled"] = adata.X.copy()
                    else:
                        if self.config.logger:
                            self.config.logger.info("Creating 'scaled' layer.")
                            adata.layers["scaled"] = adata.X.copy()
                if self.config.logger:
                    self.config.logger.info(f"Preprocessing complete. Data shape: {adata.shape} (Cells x Genes).")
        
        elif adata and not data_dict:
            if self.config.logger:
                self.config.logger.info("Scaling data.")
            sc.pp.scale(adataprocessed, max_value=max_value)

            if extra_name:
                adataprocessed.layers[f"scaled_{extra_name}"] = adataprocessed.X.copy()
            else:
                if "scaled" in adataprocessed.layers:
                    if self.config.logger:
                        self.config.logger.warning("Overwriting existing 'scaled' layer.")
                        adataprocessed.layers["scaled"] = adataprocessed.X.copy()
                else:
                    if self.config.logger:
                        self.config.logger.info("Creating 'scaled' layer.")
                        adataprocessed.layers["scaled"] = adataprocessed.X.copy()
            if self.config.logger:
                self.config.logger.info(f"Scalation complete. Data shape: {adataprocessed.shape} (Cells x Genes).")
            
        
        elif adata and data_dict:
            raise ValueError("Both AnnData object and data_dict provided. Please provide only one.")
        else:
            raise ValueError("No AnnData object or data_dict provided.")

    def annotate_celltypist(self, **kwargs):
        '''
        The expression matrix is pre-processed (and required) as log1p normalised expression to 10,000 counts per cell. 
        A raw count matrix (reads or UMIs) is required. Non-expressed genes (if you are sure of their expression absence in your data) are suggested to be included in the input table as well, as they point to the negative transcriptomic signatures when compared with the model used.
        '''

        import celltypist
        import anndata as ad
        import scanpy as sc
        import copy
        import os
        from celltypist import models

        modelos_deseados = kwargs.get('modelos_deseados', self.config.genes.modelos_anotacion)
        models.download_models(model = modelos_deseados)
        
        for sample, adata in self.data_dict.items():
            adata_celltypist = copy.deepcopy(adata)
            sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
            sc.pp.log1p(adata_celltypist)

            # Convertir a dense instead of sparse matrix
            adata_celltypist.X = adata_celltypist.X.toarray()


            # Cargar el/los modelo(s) de interés
            modelos = {os.path.splitext(modelo)[0]:models.Model.load(model = modelo) for modelo in modelos_deseados}

            #  Anotar las células con los modelos de interés
            for nombre_modelo, modelo in modelos.items(): 
                # Hacer predicción con datos crudos
                prediccion = celltypist.annotate(adata_celltypist, model=modelo, majority_voting=True)
                prediccion_adata = prediccion.to_adata()
                
                etiqueta_labels = f'celltypist_labels_{nombre_modelo}'
                etiqueta_scores = f'celltypist_scores_{nombre_modelo}'
                

                # Meter las predicciones en el objecto que se seguirá trabajando
                adata.obs[etiqueta_labels] = prediccion_adata.obs.loc[adata.obs.index, "majority_voting"]
                adata.obs[etiqueta_scores] = prediccion_adata.obs.loc[adata.obs.index, "conf_score"]

            del adata_celltypist
               
    def concatenate_samples(self, **kwargs):
        
        import anndata as ad

        if len(self.adata_list) != len(self.config.samples):
            if self.config.logger:
                self.config.logger.error("The number of samples does not match the number of AnnData objects.")
            raise ValueError("The number of samples must match the number of AnnData objects.")

        self.adata = ad.concat(
            self.adata_list, 
            join='outer', 
            label='batch', 
            merge='unique',
            keys = self.config.samples, 
            fill_value=0, 
            index_unique='-')
        
        self.adata.obs_names_make_unique()

        #all_var = [x.var for x in self.data_dict.values()]
        #all_var = pd.concat(all_var, join="outer")
        #all_var = all_var[~all_var.index.duplicated()]
        #self.adata.var = all_var
        if self.config.logger:
            self.config.logger.info(f"Concatenated data has {self.adata.n_obs} cells and {self.adata.n_vars} genes.")
    
    # Integration 
    def integrate_with_harmony(self, **kwargs):
        
        import harmonypy as hm
        import scanpy as sc

        adata = kwargs.get('adata', self.adata)
        batch_key = kwargs.get('batch_key', self.config.integration.batch_key)
        n_pcs = kwargs.get('n_pcs', self.config.integration.n_pcs)
        
        if self.config.logger:
            self.config.logger.info("Starting Harmony integration...")

        # Perform PCA if not already done
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata, svd_solver='arpack', n_comps=n_pcs)
            if self.config.logger:
                self.config.logger.info(f"Performed PCA with {n_pcs} components.")

        # Extract PCA embeddings
        data_mat = adata.obsm["X_pca"]
        metadata = adata.obs
        varsused = [batch_key]

        # Run Harmony
        ho = hm.run_harmony(
            data_mat=data_mat,
            meta_data=metadata, 
            vars_use=varsused)
        
        adata.obsm['X_pca_harmony'] = ho.Z_corr.T
        if self.config.logger:
            self.config.logger.info("Completed Harmony integration.")
    
    # Dimensionality Reduction 
    def pca(self): 
        import scanpy as sc
        
        if self.config.logger:
            self.config.logger.info("Performing PCA with {} components...".format(self.config.integration.n_pcs))

        sc.tl.pca(self.adata, n_comps=self.config.integration.n_pcs, svd_solver='arpack')  
        
        if self.config.logger:
            self.config.logger.info("Performed PCA.")

        
         
    def neighbors(self, **kwargs):
        
        import scanpy as sc

        representation = kwargs.get('representation', self.config.dim_reduction.representation)

        if self.config.logger:
            self.config.logger.info(f"Computing {self.config.dim_reduction.n_neighbors} neighborhood graph using {representation} representation...")

        sc.pp.neighbors(
            self.adata, 
            n_neighbors=self.config.dim_reduction.n_neighbors,
            n_pcs=self.config.dim_reduction.n_pcs, 
            use_rep= representation)  
        
        if self.config.logger:
            self.config.logger.info("Computed neighborhood graph.")

    def clustering(self, **kwargs):
        
        import scanpy as sc
        resolutions = kwargs.get('listaResoluciones', self.config.dim_reduction.resolutions)
        
        # Perform Leiden clustering using the constructed neighbor graph
        for res in resolutions:
            if self.config.logger:
                self.config.logger.info("Performing Leiden clustering with resolution of {}...".format(res))
            
            sc.tl.leiden(self.adata, resolution=res, key_added=f'leiden_res{res}')

        if self.config.logger:
            self.config.logger.info("Performed Leiden clustering.")

    def umap(self, **kwargs):
        import scanpy as sc
               
        min_distances = kwargs.get('lista_distancias',self.config.dim_reduction.min_distances)
        
        for min_dist in min_distances:
            umap_key = f'X_umap_minDist{min_dist}'
            
            if self.config.logger:
                self.config.logger.info(f"Computing UMAP with min_dist={min_dist}, spread={self.config.dim_reduction.spread}, n_components={self.config.dim_reduction.n_components}...")
            
            sc.tl.umap(self.adata, min_dist=min_dist, spread=self.config.dim_reduction.spread, n_components=self.config.dim_reduction.n_components)
            
            self.adata.obsm[umap_key] = self.adata.obsm['X_umap']  # Save UMAP embedding with unique name
            if self.config.logger:
                self.config.logger.info(f"Stored UMAP embedding in obsm with key '{umap_key}'.")
