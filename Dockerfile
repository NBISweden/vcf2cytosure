###########
# BUILDER #
###########
FROM clinicalgenomics/python3.8-venv:1.0

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

LABEL about.home="https://github.com/NBISweden/vcf2cytosure"
LABEL about.documentation="https://github.com/NBISweden/vcf2cytosure/blob/master/README.md"
LABEL about.tags="WGS,WES,Rare diseases,VCF,SV,Coverage,Structural variation,CNV"
LABEL about.license="MIT License (MIT)"

# Do not upgrade to the latest pip version to ensure more reproducible builds
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"
RUN echo export PATH="/venv/bin:\$PATH" > /etc/profile.d/venv.sh

COPY . /app

# Install only vcf2cytosure
RUN pip install --no-cache-dir --editable .

WORKDIR /data/