#build:
docker build -t gemli-r:v0 -f Dockerfile_v0 .

#push to registry
docker login #user and pass from dockerhub here
docker tag gemli-r:v0 fenix07/gemli-r:v0
docker push fenix07/gemli-r:v0

#pull with docker
docker pull fenix07/gemli-r:v0
docker run -it --rm --privileged fenix07/gemli-r:v0
