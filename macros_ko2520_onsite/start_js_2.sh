#!/usr/bin/env bash
cd "$(dirname "$0")" || exit
source ../../lilak.sh
lilak js -I 192.168.1.102:9092 -D data_reco2
