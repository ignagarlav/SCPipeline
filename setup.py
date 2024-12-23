from setuptools import setup, find_packages

setup(
    name="SCPipeline",                 # Nombre del proyecto
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=find_packages(),           # Detecta automáticamente los módulos y paquetes
    install_requires=[],                # Lista de dependencias externas si aplica
    description="Herramientas para el análisis de single cell RNA-Seq",  # Descripción del proyecto
    author="Jose Ignacio Garzón",                 # Autor del proyecto
    author_email="igarzonalva@hpc-2.unav.es",  # Email de contacto
    url="https://github.com/ignagarlalv/SCPipeline",  # URL del repositorio si aplica
)
