FROM ubuntu:18.04

RUN apt update -y && \
    apt install -y python3 python3-pip python3-venv && \
    apt install -y g++ && \
    apt install -y zlib1g-dev && \
    apt install -y wget

RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar xzf cd-hit-v4.8.1-2019-0228.tar.gz && \
    rm -rf cd-hit-v4.8.1-2019-0228.tar.gz && \
    cd cd-hit-v4.8.1-2019-0228 && \
    make && \
    make install && \
    cd .. && \
    rm -rf cd-hit-v4.8.1-2019-0228

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.450-1_amd64.deb && \
    dpkg -i mafft_7.450-1_amd64.deb && \
    rm -rf mafft_7.450-1_amd64.deb

RUN mkdir /app
WORKDIR /app

COPY requirements.txt /app
RUN python3 -m venv .venv && \
    . .venv/bin/activate && \
    python -m pip install -U pip && \
    pip install -U setuptools && \
    pip install -U wheel && \
    pip install -r requirements.txt

COPY . /app
RUN cd UnitePSSMs && \
    g++ *.cpp -std=c++11 -O3 -o UnitePSSMs && \
    cd ../PSSM_score_Peptide && \
    g++ *.cpp -std=c++11 -O3 -o PSSM_score_Peptide && \
    cd ../hits_cpp && \
    g++ *.cpp -std=c++11 -O3 -o hits && \
    cd ../tfidf && \
    g++ *.cpp -std=c++11 -O3 -o tfidf

ENV APP_FILE IgOmeProfiling_pipeline.py
RUN chmod +x entrypoint.sh
ENTRYPOINT ["./entrypoint.sh"]
CMD ["-h"]