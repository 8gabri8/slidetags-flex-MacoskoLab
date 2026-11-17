# Image for the slidetags_flex Terra workflow:
# Python 3.11 + numpy/pandas/scipy/matplotlib, gcloud (gs:// paths), pigz (fast .gz reading), slidetags_flex.py
#
# Built automatically by .github/workflows/docker.yml (GitHub -> Docker Hub). Or with Google Cloud Build:
#   gcloud builds submit --tag us-central1-docker.pkg.dev/<PROJECT>/<REPOSITORY>/slidetags-flex:v2.0.0 .
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

RUN pip install --no-cache-dir numpy==2.1.3 pandas==2.2.3 scipy==1.14.1 matplotlib==3.9.2

COPY slidetags_flex.py /opt/slidetags_flex.py
ENV PYTHONUNBUFFERED=1

# Smoke test at build time
RUN python3 /opt/slidetags_flex.py --help > /dev/null
