FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit
RUN pip install tqdm==4.64.1
RUN pip install ratelimit 
RUN pip install standardiser
RUN pip install pandas 

WORKDIR /repo
COPY . /repo
