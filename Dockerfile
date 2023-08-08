FROM us.gcr.io/landerlab-atacseq-200218/gcloud-samtools:0.1
FROM python:3.8
RUN apt install zlib1g

WORKDIR /app
ENV PATH=$PATH:/app

COPY requirements.txt .
RUN python3 -m pip install -r requirements.txt

COPY donor_assignment ./donor_assignment
