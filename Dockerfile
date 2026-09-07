FROM condaforge/miniforge3:latest

WORKDIR /app

# Create a dedicated conda environment with MAGeCK from Bioconda
RUN conda create -n mageck-env -c bioconda -c conda-forge mageck=0.5.9.5 python=3.11 -y \
    && conda clean -afy

# Make the conda environment's binaries (including mageck) available on PATH
ENV PATH="/opt/conda/envs/mageck-env/bin:$PATH"

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app/ ./app/

EXPOSE 8501

HEALTHCHECK CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

ENTRYPOINT ["streamlit", "run", "app/main.py", "--server.port=8501", "--server.address=0.0.0.0"]