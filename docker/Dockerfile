FROM ubuntu:20.04

LABEL autor "Romuald Marin"
# Romuald Marin

ENV DEBIAN_FRONTEND noninteractive
ENV http_proxy http://ftprox.ibcp.fr:3128/
ENV https_proxy http://ftprox.ibcp.fr:3128/


RUN apt-get update && apt-get install -y \
    git \
    python \
    python3.9 \
    python3-pip

RUN pip3 install pip --upgrade
RUN pip3 install rdflib --upgrade \
    && pip3 install cython --upgrade \
    && pip3 install numpy --upgrade \
    && pip3 install Pillow --upgrade \
    && pip3 install matplotlib --upgrade \
    && pip3 install pandas --upgrade \
    && pip3 install --upgrade setuptools 

  
RUN pip install biopython

# LIGHTDOCK PART
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'
RUN pip3 install lightdock
RUN export LIGHTDOCK_HOME="$(which lightdock)"
RUN export PATH=$PATH:$LIGHTDOCK_HOME/bin
RUN export PYTHONPATH=$PYTHONPATH:$LIGHTDOCK_HOME 

#MDAnalysys
RUN pip install --upgrade MDAnalysis

# DOCKQ
#need to modify the index
#grep if len(model_chains) > 2 or len(native_chains)> 2:
#and change 2 to 1 
# " if len(model_chains) > 1 or len(native_chains)> 1:"

RUN git clone https://github.com/bjornwallner/DockQ/
WORKDIR "/DockQ"
RUN make
RUN echo 'alias DockQ="./DockQ/DockQ.py"' >> ~/.bashrc


WORKDIR /
#Vorona score
#https://github.com/kliment-olechnovic/voronota
RUN apt install voronota


#JWALK score
RUN git clone https://github.com/Topf-Lab/XLM-Tools
WORKDIR "/XLM-Tools"
RUN python3 setup.py install
RUN echo 'alias jwalk="python XLM-Tools/Jwalk.v2.1.py"' >> ~/.bashrc

WORKDIR /

RUN rm -r /usr/bin/python &&\
    ln -s /usr/bin/python3 /usr/bin/python 

ENV PATH="${PATH}:/DockQ"

RUN apt-get install -y pymol

RUN echo 'alias jwalk="python XLM-Tools/Jwalk.v2.1.py"' >> /root/.bashrc
RUN alias jwalk="python XLM-Tools/Jwalk.v2.1.py"

RUN apt-get install bc