#! /usr/bin/env bash

kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
kraken2-build --download-taxonomy --db kraken2-human-db/
kraken2-build --build --db kraken2-human-db/ --threads 30
kraken2-build --clean --db kraken2-human-db/