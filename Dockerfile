FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit==2023.3.3
RUN pip install standardiser==0.1.9
RUN pip install requests==2.31.0
RUN pip install six==1.16
RUN pip install numpy==1.24.4

WORKDIR /repo
COPY . /repo
