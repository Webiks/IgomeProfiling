#!/bin/bash
docker build . -t webiks/igome-profile:latest
docker build . -f DockerfileWorker -t webiks/igome-profile-worker:latest
