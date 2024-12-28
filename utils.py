import os
import logging
import sys
def setup_logging(logger_name, log_file, level=logging.INFO):

    logger = logging.getLogger(logger_name)
    
    # Avoid adding multiple handlers if the logger already exists
    if logger.hasHandlers():
        return logger
    
    # Set the logging level
    logger.setLevel(level)

    # Ensure the log file's directory exists
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Create file handler for logging to a file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(level)

    # Create console handler for logging to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    # Create a logging format
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)



    return logger



def exportar_datosV1(adata, seurat_dir): 

    # Importar los paquetes necesarios
    import os 
    import datetime
    from anndata import AnnData
    import pandas as pd
    from scipy.io import mmwrite
    import scipy
    # Crear el directorio base si no existe
    current_date = datetime.datetime.now().strftime("%Y-%m-%d")
    os.makedirs(seurat_dir, exist_ok=True)
    fullPath = os.path.join(seurat_dir, current_date)
    os.makedirs(fullPath, exist_ok=True)

    # Save .var (gene metadata)
    var_path = os.path.join(fullPath, "var.csv")
    adata.var.to_csv(var_path)
    print(f"Saved .var to {var_path}")


    # Export metadata (e.g., cell type annotations, clustering)
    file_name = f"{current_date}_metadata.csv"
    file_path = os.path.join(fullPath, file_name)
    adata.obs.to_csv(file_path)
 
    # Save each layer in .layers
    for layer_name in adata.layers.keys():
        layer_path = os.path.join(fullPath, f"layer_{layer_name}_sparse.mtx")

        if isinstance(adata.layers.get(layer_name), scipy.sparse.spmatrix):
            mmwrite(layer_path, adata.layers[layer_name])
        else:
            layer_sparse = scipy.sparse.csr_matrix(adata.layers[layer_name])
            mmwrite(layer_path, layer_sparse)

        print(f"Saved {layer_name} in MTX format to {layer_path}")

   
    # Export UMAP coordinates
    if any(key.startswith('X_umap') for key in adata.obsm.keys()):

        for key in adata.obsm.keys():
                # Check if the key corresponds to a UMAP embedding
                if key.startswith('X_umap'):
                    # Extract UMAP coordinates
                    umap_coords = pd.DataFrame(
                        adata.obsm[key],
                        index=adata.obs_names,
                        columns=[f'{key}_1', f'{key}_2']
                    )
                    # Construct file name and path
                    file_name = f"{current_date}_{key}_coords.csv"
                    file_path = os.path.join(fullPath, file_name)
                
                    # Save to CSV
                    umap_coords.to_csv(file_path)







