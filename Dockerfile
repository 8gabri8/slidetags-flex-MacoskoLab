# Image for the slidetags_flex Terra workflow:
# Python 3.11 + numpy/pandas/scipy/matplotlib, scanpy + h5py (section summaries), gcloud (gs:// paths),
# pigz (fast .gz reading), and the three scripts.
#
# Built by .github/workflows/docker.yml (GitHub -> Docker Hub). Or with Google Cloud Build:
#   gcloud builds submit --tag us-central1-docker.pkg.dev/<PROJECT>/<REPOSITORY>/slidetags-flex:v3.0.0 .
FROM python:3.11-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl gnupg ca-certificates pigz \
 && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
 && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends google-cloud-cli \
 && rm -rf /var/lib/apt/lists/*

# scanpy and everything it uses, pinned to versions that work together (scanpy 1.10.4 needs
# scipy < 1.14: newer scipy removed a function scanpy calls for sections > 4,096 nuclei)
RUN pip install --no-cache-dir numpy==2.0.2 pandas==2.2.3 scipy==1.13.1 matplotlib==3.9.2 \
      h5py==3.12.1 scanpy==1.10.4 anndata==0.11.1 scikit-learn==1.5.2 umap-learn==0.5.7 \
      pynndescent==0.5.13 numba==0.60.0 igraph==0.11.8 leidenalg==0.10.2

COPY slidetags_flex.py merge_sections.py report_section.py /opt/
ENV PYTHONUNBUFFERED=1 PYTHONPATH=/opt NUMBA_CACHE_DIR=/tmp/numba MPLCONFIGDIR=/tmp/matplotlib

# Checks at build time: the command line works, and scanpy (UMAP, Leiden) + the .h5 reader work
# on a small made-up dataset. If one fails, the image is not built.
RUN python3 /opt/slidetags_flex.py --help > /dev/null \
 && python3 /opt/merge_sections.py --self-test
