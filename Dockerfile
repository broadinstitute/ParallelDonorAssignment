FROM us.gcr.io/landerlab-atacseq-200218/gcloud-samtools:0.1
RUN apt install zlib1g

WORKDIR /app
ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt

COPY donor_assignment ./donor_assignment
COPY cisvar ./cisvar
COPY monitor_script.sh .
RUN chmod a+rx /app/monitor_script.sh
