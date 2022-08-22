###########
# BUILDER #
###########
FROM clinicalgenomics/python3.8-venv:1.0 AS python-builder

ENV PATH="/venv/bin:$PATH"

WORKDIR /app

RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y install -y --no-install-recommends gcc && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/*

# reqs
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

#########
# FINAL #
#########
FROM python:3.8-slim

LABEL about.home="https://github.com/NBISweden/vcf2cytosure"
LABEL about.documentation="https://github.com/NBISweden/vcf2cytosure/blob/master/README.md"
LABEL about.tags="WGS,WES,Rare diseases,VCF,SV,Coverage,Structural variation,CNV"
LABEL about.license="MIT License (MIT)"

# Do not upgrade to the latest pip version to ensure more reproducible builds
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"
RUN echo export PATH="/venv/bin:\$PATH" > /etc/profile.d/venv.sh

# Create a non-root user to run commands
RUN groupadd --gid 1000 worker && useradd -g worker --uid 1000 --shell /usr/sbin/nologin --create-home worker

# Copy virtual environment from builder
COPY --chown=worker:worker --from=python-builder /venv /venv

WORKDIR /home/worker/app
COPY --chown=worker:worker . /home/worker/app

# Install only vcf2cytosure
RUN pip install --no-cache-dir --editable .

# Run the app as non-root user
USER worker

WORKDIR /data/