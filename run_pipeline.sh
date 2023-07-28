# REPLACE WITH CONTAINER ID
CONTAINER="3a456ba31724a7d6dd4dc65be4e12815d0" # e.g 3a456ba31724a7d6dd4dc65be4e12815d0

docker run \
    -v $PWD:$PWD \
    -w $PWD $CONTAINER \
    bash nextflow_command.sh