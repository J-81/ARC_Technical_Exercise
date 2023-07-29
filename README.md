# ARC_Technical_Exercise

## Building Dockerfile

This repository contains required files to build a docker image for the ARC Technical Exercise.

The following build script should take about no more than 30 minutes to complete, depending on your internet connection.

```bash
bash build_docker.sh
```

This will output a docker image sha256 hash that will be used for running the workflow in the next step.

Example
> => writing image sha256:cc8d7ea8b1363b109f7c769513f17ad54cf8e3e0fbae44eaee57607069c237b

## Running The Workflow Via Container

Once built the docker image can be run with the following command:

```bash
bash run_pipeline.sh <DOCKER BUILD SHA>
# e.g.
# bash run_pipeline.sh cc8d7ea8b1363b109f7c769513f17ad54cf8e3e0fbae44eaee57607069c237b
```

This will launch the container and run the workflow. The output will be written to the `final_test` directory.

The read files associated with the exercise are located in the `test_data` directory and used as input for the workflow.

## Output Files

The output files are located in the `final_test` directory include organized relevant results for the workflow.

Most importantly, the requested text report is located at `full_report/report.txt` and the requested plot is located at `full_report/Quality_scores.png`.

## Known Issues

A permissions issue within the container currently prevents the container from being run internally as the host user.

The planned fix is simply updating the permissions within the container for the mOTU databases; however, storage limitations in the development environment prevented this from being completed in time.
