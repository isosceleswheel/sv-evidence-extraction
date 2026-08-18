# Starting point for the sv_evidence_extraction WDL task's runtime image.
# Build from this repo's root, e.g.:
#   docker build -t us-central1-docker.pkg.dev/PROJECT_ID/sv-evidence-extraction/sv-evidence-extraction:latest .
#   docker push us-central1-docker.pkg.dev/PROJECT_ID/sv-evidence-extraction/sv-evidence-extraction:latest
# then point the WDL's `docker` input at the pushed image.
#
# pysam's PyPI wheel already bundles an htslib built with GCS/libcurl
# support -- confirmed working against a public gs:// tabix file with no
# extra setup -- so no gcloud SDK or manual GCS_OAUTH_TOKEN plumbing is
# needed for the primary access path. gcloud is still installed as the
# fallback path's dependency, in case a given runtime's htslib build
# ever lacks GCS support.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl gnupg ca-certificates \
    && curl -sSL https://sdk.cloud.google.com | bash -s -- --disable-prompts --install-dir=/opt \
    && apt-get purge -y curl gnupg \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/google-cloud-sdk/bin:${PATH}"

RUN pip install --no-cache-dir pysam pandas pyarrow

COPY sv_evidence_extraction/ /opt/sv_evidence_extraction/
ENV PYTHONPATH="/opt:${PYTHONPATH}"

WORKDIR /work
