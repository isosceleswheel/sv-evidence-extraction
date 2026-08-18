# Starting point for the sv_evidence_extraction WDL task's runtime image.
# Build and push from this repo's root via Cloud Build (no local Docker
# daemon required) -- this is the actual command used for the image the
# WDL currently defaults to:
#   gcloud builds submit \
#     --project=talkowski-sv-gnomad \
#     --tag us.gcr.io/talkowski-sv-gnomad/dam/sv-evidence-extraction:latest \
#     .
# Swap the project/path for your own if you're maintaining a fork --
# "dam" here is just the pushing user's initials, used as a namespacing
# folder within the shared talkowski-sv-gnomad registry.
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
