# IgOme Profiling pipeline

## Overview
IgOme Profiling pipeline is ...

## Getting started

### OS
The pipeline runs on Linux environment only.  
For Windows the pipeline can run using [WSL and Visual Studio Code](https://code.visualstudio.com/remote-tutorials/wsl/getting-started).

Tested on Ubuntu 18.04 via Windows.

### Prerequisites
The following are required in order to run the code:  
* git, e.g.:
  ```bash
  sudo apt install git
  ```
* Python 3.6+, .e.g:
  ```bash
  sudo apt install python3 python3-pip python3-venv
  ```
* g++, e.g.:
  ```
  sudo apt install g++
  ```
* CD-HIT:
  * Install using instructions: http://weizhongli-lab.org/cd-hit/
  * Make sure cd-hit is installed globally, e.g. /usr/local/bin
  * Test by running:
    ```bash
    # From project directory
    cd-hit -h
    ```
* MAFFT:
  * Install using instructions: https://mafft.cbrc.jp/alignment/software/
  * Make sure mafft is installed globally, e.g. /usr/bin
  * Test by running:
    ```bash
    # From project directory
    mafft -h
    ```

It is advised to update the packages manager before installing, e.g.:
```bash
sudo apt update
```

## Repository
Configure git user and email (if not done yet):
```bash
git config --global user.name "FIRST_NAME LAST_NAME"
git config --global user.email "you@mail.com"
```

Clone the repository:
```bash
# From projects/workspace directory
git clone https://github.com/Webiks/IgomeProfiling.git
cd IgomeProfiling
```

### Compile
Most of the code is in Python but some of the code is in C++ and requires compilation:
* UnitePSSMs:
  ```bash
  # From project directory
  cd UnitePSSMs
  g++ *.cpp -std=c++11 -O3 -o UnitePSSMs
  ```
* PSSM_score_Peptide:
  ```bash
  # From project directory
  cd PSSM_score_Peptide
  g++ *.cpp -std=c++11 -O3 -o PSSM_score_Peptide
  ```
* Hits
  ```bash
  # From project directory
  cd hits_cpp
  g++ *.cpp -std=c++11 -O3 -o hits
  ```
* TF-IDF/Merge
  ```bash
  # From project directory
  cd tfidf
  g++ *.cpp -std=c++11 -O3 -o tfidf
  ```

### Python 
Python3 should be installed in previous steps.  
Use the following to install the required packages:
* Create a virtual environment for the project:
  ```bash
  # From project directory
  python3 -m venv .venv
  source .venv/bin/activate
  ```
* Install the dependencies:
  ```bash
  pip3 install -r requirements.txt
  ```

## Running
The pipeline is modular, each step is a python file which can be run by itself by passing command line arguments.  
To run the entire pipeline the entry point is ```IgOmeProfiling_pipeline.py```.  
Run the pipeline using the following:  
```bash
python3 IgOmeProfiling_pipeline.py fastq_file_path barcodes_to_samples_path samples_to_biological_conditions_path analysis_results_directory run_logs_directory
```

For more information and options run the following:
```bash
python3 IgOmeProfiling_pipeline.py -h
```

The code assumes existence of ```bash``` in order to run.  
This can be changed via ```local_command_prefix``` variable in ```global_params.py```.

## Important flags:
The following is are important flags, see help for details:
* left_construct: Depends on experiment, default is for EXP12
* right_construct: Depends on experiment, default is for EXP12
* allowed_gap_frequency: Threshold for gappy columns
* concurrent_cutoffs: If set then cutoffs are calculated concurrently, if not each BC run as a whole.
* meme_split_size: Batch size of memes for cutoffs (if concurrent_cutoffs is set) and pValues calculation.
* number_of_random_pssms: Number of random PSSMs for calculating pValues.

## Testing
The code contains mock data for testing.

Set configuration in ```global_params.py```.
To run local set:
* run_using_celery to False
* run_local_in_parallel_mode to True

Run the following:
```bash
# From project directory
python3 IgOmeProfiling_pipeline.py mock_data/exp12_10M_rows.fastq.gz mock_data/barcode2samplename.txt mock_data/samplename2biologicalcondition.txt output/analysis output/logs
```

The entire pipeline might run a few hours on mock data on 96 cores.  
The output of the pipeline are heatmaps of most important motifs per biological condition va all samples.

## Docker
The code can be containerized using Docker:

Set configuration in ```global_params.py```.
To run local set:
* run_using_celery to True if using multiple machines
* run_local_in_parallel_mode to True

```bash
# From project directory
./build_docker.sh
```

The image is built without data and command argument of "-h".
To run the mock-data using docker:
```bash
# From project directory
docker run --name igome --rm -v ./some_data:/data -v ./test:/output webiks/igome-profile /data/exp12_10M_rows.fastq.gz /data/barcode2samplename.txt /data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

The mock data is included in the docker, to run it:
```bash
docker run --name igome --rm -v ./test:/output webiks/igome-profile ./mock_data/exp12_10M_rows.fastq.gz ./mock_data/barcode2samplename.txt ./mock_data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

Upload to AWS (using aws-cli with credentials set):
```bash
./aws_upload.sh
```

In AWS machine (with aws-cli credentials set):
```bash
$(aws ecr get-login --no-include-email --region us-west-2)
docker pull 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest
docker run --name igome --rm -v ./test:/output 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest ./mock_data/exp12_10M_rows.fastq.gz ./mock_data/barcode2samplename.txt ./mock_data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

## Running on multiple machines
Running on multiple machines is support via [Celery](http://www.celeryproject.org/) and [RabbitMQ](https://www.rabbitmq.com/).  
Pipeline steps synchronize via files. Therefore, the same paths must be mounted on all machines. For example if using AWS mount an EFS to /data on all machines.  

Run RabbitMQ + Flower via docker-compose:
```yaml
version: '3'

services:
  rabbit:
    image: rabbitmq:management
    restart: always
    ports:
     - 5672:5672
     - 15672:15672
  flower:
    image: mher/flower
    restart: always
    command: flower --broker=pyamqp://guest@rabbit//
    ports:
     - 5555:5555
```
Open [RabbitMQ Management](http://localhost:15672) and/or [Flower](http://localhost:5555) to monitor the pipeline tasks.  

On a worker instance:
* Mount shared drive to directory, e.g. /data
* Pull worker image from repository
* Run the worker
  ```bash
  docker run -d --restart always -e CELERY_BROKER_URL="pyamqp://guest@[broker-ip]//" -v /data:/data --name igomeworker 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile-worker:latest
  ```
  * Replace [broker-ip] and guest with RabbitMQ host and user/password (if defined).
  * Change volume mount to shared drive

Run the pipeline of some machine:
* Mount shared drive, e.g. /data
* Pull pipeline image from repository
* Run the pipeline, e.g.
  ```bash
  docker run --name igome --rm -e CELERY_BROKER_URL="pyamqp://guest@[broker-ip]//" -d -v /data:/data 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile /data/mock_data/exp12_10M_rows.fastq.gz /data/mock_data/barcode2samplename.txt /data/mock_data/samplename2biologicalcondition.txt /data/output/analysis /data/output/logs
  ```
* Use ```docker logs --tail 10 -f igome``` to see over all progress


## License
TBD