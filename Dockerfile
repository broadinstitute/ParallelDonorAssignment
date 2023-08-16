FROM us.gcr.io/landerlab-atacseq-200218/gcloud-samtools:0.1

RUN apt install zlib1g

ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt
RUN rm requirements.txt  # cleanup

WORKDIR /app
