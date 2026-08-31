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
# pysam's PyPI wheel bundles an htslib built with GCS/libcurl support,
# but htslib's GCS backend does NOT discover Application Default
# Credentials on its own -- confirmed by hand against a real private
# bucket, where pysam.TabixFile() failed until GCS_OAUTH_TOKEN was set
# explicitly. sv_evidence_extraction.core refreshes that env var via
# `gcloud auth print-access-token` before opening each gs:// file, so
# both gcloud AND the `tabix` CLI (the fallback if pysam still can't
# open a file) need to actually be present in this image.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl gnupg ca-certificates tabix \
    && curl -sSL https://sdk.cloud.google.com | bash -s -- --disable-prompts --install-dir=/opt \
    && apt-get purge -y curl gnupg \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/google-cloud-sdk/bin:${PATH}"

RUN pip install --no-cache-dir pysam pandas pyarrow

COPY sv_evidence_extraction/ /opt/sv_evidence_extraction/
ENV PYTHONPATH="/opt:${PYTHONPATH}"

WORKDIR /work
