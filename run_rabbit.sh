#!/bin/bash
docker-compose -f rabbit-docker-compose.yml down
docker-compose -f rabbit-docker-compose.yml up -d
