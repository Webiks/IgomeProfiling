#/bin/bash
$(aws ecr get-login --no-include-email --region us-west-2)
docker tag webiks/igome-profile:latest 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest
docker push 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest

docker tag webiks/igome-profile-worker:latest 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile-worker:latest
docker push 686447933053.dkr.ecr.us-west-2.amazonaws.com/igome-profile-worker:latest
