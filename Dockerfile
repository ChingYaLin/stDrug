FROM tensorflow/tensorflow:2.19.0
RUN apt-get update && apt-get install -y git unzip wget curl
RUN pip3 install --upgrade pip cmake
WORKDIR /stDrug
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
COPY . .

## scMatch
RUN git clone https://github.com/asrhou/scMatch.git /opt/scMatch

RUN unzip /opt/scMatch/refDB/FANTOM5/10090_HID.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/10090_HID.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/10090_symbol.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/10090_symbol.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/9606_HID.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/9606_HID.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/9606_symbol.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/9606_symbol.csv.zip 

RUN sed -i 's/\.ix/.loc/g' /opt/scMatch/scMatch.py
RUN sed -i 's/loc\[commonRows, ].fillna(0\.0)/reindex(commonRows, axis="index", fill_value=0.0)/g' /opt/scMatch/scMatch.py

# survival analysis
RUN mkdir -p data \
 && curl -fL --retry 5 --retry-all-errors --max-time 600 \
      -H "User-Agent: Mozilla/5.0" \
      "https://figshare.com/ndownloader/files/35612942?download=1" \
      -o data/TCGA.zip \
 && unzip -q data/TCGA.zip -d data \
 && rm data/TCGA.zip

## CaDRReS-Sc
RUN git clone https://github.com/CSB5/CaDRReS-Sc.git /opt/CaDRReS-Sc
RUN pip install gdown && \
    gdown --id 19d5kIP2YlChqLFZZo4aWoU2maspa9oe1 -O /opt/CaDRReS-Sc/data/GDSC/GDSC_exp.tsv && \
    mkdir -p /opt/CaDRReS-Sc/data/CCLE && \
    gdown --id 19lU6RnCjx57Oj0UZlpIwieHYTYdxc7MN -O //opt/CaDRReS-Sc/data/CCLE/CCLE_expression.csv
RUN mkdir -p /opt/CaDRReS-Sc/preprocessed_data/PRISM
RUN wget -q --no-check-certificate 'https://docs.google.com/uc?export=download&id=1AhjNbT88--PmG8qZSOcF_C2tYm-FoC2c' -O /opt/CaDRReS-Sc/preprocessed_data/PRISM/PRISM_drug_info.csv
RUN wget -q --no-check-certificate 'https://docs.google.com/uc?export=download&id=1_TCD-OO-l1dsnwLoLlb8eD92grE6_ZPU' -O /opt/CaDRReS-Sc/preprocessed_data/PRISM/feature_genes.txt

RUN sed -i 's/import tensorflow as tf/import tensorflow.compat.v1 as tf\ntf.disable_v2_behavior()/g' /opt/CaDRReS-Sc/cadrres_sc/model.py
RUN sed -i 's/import tensorflow\.python\.util\.deprecation as deprecation/from tensorflow.python.util import deprecation/g' /opt/CaDRReS-Sc/cadrres_sc/model.py

CMD [ "/bin/bash" ]


